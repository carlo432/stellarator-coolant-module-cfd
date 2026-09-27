#!/usr/bin/env python3
"""Build first CFD-ready OpenFOAM case from anchor-exact Ferrero Fig. 2.1 contours."""

from __future__ import annotations

import argparse
import json
import math
import stat
from pathlib import Path

from build_arc_blanket_2d_case import MU, NU, RHO, THICKNESS, foam_header, write


TANK_SPEED = 0.313
CHANNEL_SPEED = 1.98
L_REF = 0.02
DEFAULT_PATCH_TARGETS = Path("references/digitization/ferrero_arc/arc_fig21_cfd_patch_targets_default.json")


def chmod_x(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def signed_area_mm(points: list[list[float]]) -> float:
    return 0.5 * sum(
        x0 * y1 - x1 * y0
        for (x0, y0), (x1, y1) in zip(points, points[1:] + points[:1])
    )


def split_axis_segment(points: list[list[float]], split_points: list[list[float]]) -> list[list[float]]:
    """Insert split points that lie on existing straight axis-aligned segments."""
    out: list[list[float]] = []
    n = len(points)
    for i, a in enumerate(points):
        b = points[(i + 1) % n]
        out.append(a)
        inserts = []
        ax, ay = a
        bx, by = b
        for p in split_points:
            px, py = p
            on_x = abs(ax - bx) < 1e-9 and abs(px - ax) < 1e-9 and min(ay, by) < py < max(ay, by)
            on_y = abs(ay - by) < 1e-9 and abs(py - ay) < 1e-9 and min(ax, bx) < px < max(ax, bx)
            if on_x or on_y:
                if abs(bx - ax) >= abs(by - ay):
                    t = (px - ax) / (bx - ax)
                else:
                    t = (py - ay) / (by - ay)
                inserts.append((t, p))
        for _, p in sorted(inserts):
            out.append(p)
    return out


def mm_to_local_m(points: list[list[float]], mesh_size: float, inner: bool = False) -> list[tuple[float, float, float]]:
    h = mesh_size * 0.55 if inner else mesh_size
    return [((r - 3300.0) / 1000.0, z / 1000.0, h) for r, z in points]


def dist2(a: tuple[float, float], b: tuple[float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def line_mid(points: list[tuple[float, float, float]], i: int) -> tuple[float, float]:
    a = points[i]
    b = points[(i + 1) % len(points)]
    return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5)


def find_line_by_endpoints(
    points: list[tuple[float, float, float]],
    target_a: tuple[float, float],
    target_b: tuple[float, float],
) -> int:
    best_i = -1
    best = 1e99
    for i in range(len(points)):
        a = points[i]
        b = points[(i + 1) % len(points)]
        score = dist2((a[0], a[1]), target_a) + dist2((b[0], b[1]), target_b)
        score_rev = dist2((a[0], a[1]), target_b) + dist2((b[0], b[1]), target_a)
        score = min(score, score_rev)
        if score < best:
            best = score
            best_i = i
    if best > 1e-10:
        raise ValueError(f"could not find line endpoints {target_a} {target_b}; best={best:g}")
    return best_i


def closest_inner_line(inner: list[tuple[float, float, float]], target_rz_mm: tuple[float, float]) -> int:
    target = ((target_rz_mm[0] - 3300.0) / 1000.0, target_rz_mm[1] / 1000.0)
    return min(range(len(inner)), key=lambda i: dist2(line_mid(inner, i), target))


def mm_endpoint(point: list[float]) -> tuple[float, float]:
    return ((point[0] - 3300.0) / 1000.0, point[1] / 1000.0)


def load_patch_targets(path: Path) -> dict:
    data = json.loads(path.read_text())
    required = ("outer_split_points_mm", "inner_targets_mm", "basis")
    missing = [k for k in required if k not in data]
    if missing:
        raise ValueError(f"patch target file {path} missing keys: {missing}")
    for name in ("inlet_ch_a", "inlet_ch_b", "inlet_main40"):
        if name not in data["inner_targets_mm"]:
            raise ValueError(f"patch target file {path} missing inner target {name}")
    return data


def outward_normal(points: list[tuple[float, float, float]], i: int, boundary: str) -> tuple[float, float]:
    a = points[i]
    b = points[(i + 1) % len(points)]
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    L = math.hypot(dx, dy)
    if L == 0:
        return (0.0, 0.0)
    # Outer is CCW: outward is right normal. Inner contour is CW in this case:
    # outward from fluid into the exclusion is also right normal for the line order.
    return (dy / L, -dx / L)


def inlet_vec(points: list[tuple[float, float, float]], idx: int, speed: float) -> tuple[float, float, float]:
    nx, ny = outward_normal(points, idx, "boundary")
    return (-speed * nx, -speed * ny, 0.0)


def write_geo(
    case: Path,
    contours: Path,
    mesh_size: float,
    patch_targets_path: Path,
) -> dict[str, list[int] | dict[str, tuple[float, float, float]]]:
    data = json.loads(contours.read_text())
    targets = load_patch_targets(patch_targets_path)
    outer_mm = data["tank_inner_wall"]
    inner_mm = data["vv_exclusion_union_contour"]
    if signed_area_mm(outer_mm) < 0:
        outer_mm = list(reversed(outer_mm))
    if signed_area_mm(inner_mm) > 0:
        inner_mm = list(reversed(inner_mm))

    outer_mm = split_axis_segment(
        outer_mm,
        targets["outer_split_points_mm"],
    )
    outer = mm_to_local_m(outer_mm, mesh_size)
    inner = mm_to_local_m(inner_mm, mesh_size, inner=True)

    outer_splits = targets["outer_split_points_mm"]
    outlet_i = find_line_by_endpoints(outer, mm_endpoint(outer_splits[0]), mm_endpoint(outer_splits[1]))
    if len(outer_splits) >= 4:
        aux_i = find_line_by_endpoints(outer, mm_endpoint(outer_splits[2]), mm_endpoint(outer_splits[3]))
    else:
        aux_i = find_line_by_endpoints(outer, ((4605 - 3300) / 1000, 3.560), ((4605 - 3300) / 1000, 3.660))
    inner_targets = targets["inner_targets_mm"]
    ch_a_i = closest_inner_line(inner, tuple(inner_targets["inlet_ch_a"]))
    ch_b_i = closest_inner_line(inner, tuple(inner_targets["inlet_ch_b"]))
    main_i = closest_inner_line(inner, tuple(inner_targets["inlet_main40"]))

    geom = case / "geometry"
    geom.mkdir(parents=True, exist_ok=True)
    all_points = outer + inner
    lines = []
    for i, (x, y, h) in enumerate(all_points, 1):
        lines.append(f"Point({i}) = {{{x:.8g}, {y:.8g}, 0, {h:.8g}}};")

    outer_ids = []
    for i in range(len(outer)):
        lid = i + 1
        lines.append(f"Line({lid}) = {{{i + 1}, {(i + 1) % len(outer) + 1}}};")
        outer_ids.append(lid)
    inner_ids = []
    offset = len(outer)
    for i in range(len(inner)):
        lid = len(outer) + i + 1
        a = offset + i + 1
        b = offset + ((i + 1) % len(inner)) + 1
        lines.append(f"Line({lid}) = {{{a}, {b}}};")
        inner_ids.append(lid)

    patch_ids = {
        "outlet180": [outer_ids[outlet_i]],
        "inlet_aux100": [outer_ids[aux_i]],
        "inlet_ch_a": [inner_ids[ch_a_i]],
        "inlet_ch_b": [inner_ids[ch_b_i]],
        "inlet_main40": [inner_ids[main_i]],
    }
    used = {idx for ids in patch_ids.values() for idx in ids}
    all_ids = outer_ids + inner_ids
    wall_ids = [idx for idx in all_ids if idx not in used]

    def surf_refs(curve_indices: list[int]) -> str:
        return ", ".join(f"ext[{2 + (idx - 1)}]" for idx in curve_indices)

    geo = "\n".join(
        [
            'SetFactory("Built-in");',
            f"Mesh.CharacteristicLengthMin = {mesh_size * 0.25:.8g};",
            f"Mesh.CharacteristicLengthMax = {mesh_size:.8g};",
            "Mesh.Algorithm = 6;",
            "Mesh.RecombineAll = 0;",
            *lines,
            f"Curve Loop(1) = {{{', '.join(str(i) for i in outer_ids)}}};",
            f"Curve Loop(2) = {{{', '.join(str(-i) for i in reversed(inner_ids))}}};",
            "Plane Surface(1) = {1, 2};",
            f"ext[] = Extrude {{0, 0, {THICKNESS:.8g}}} {{ Surface{{1}}; Layers{{1}}; Recombine; }};",
            'Physical Surface("frontAndBack") = {1, ext[0]};',
            f'Physical Surface("outlet180") = {{{surf_refs(patch_ids["outlet180"])}}};',
            f'Physical Surface("inlet_aux100") = {{{surf_refs(patch_ids["inlet_aux100"])}}};',
            f'Physical Surface("inlet_ch_a") = {{{surf_refs(patch_ids["inlet_ch_a"])}}};',
            f'Physical Surface("inlet_main40") = {{{surf_refs(patch_ids["inlet_main40"])}}};',
            f'Physical Surface("inlet_ch_b") = {{{surf_refs(patch_ids["inlet_ch_b"])}}};',
            f'Physical Surface("walls") = {{{surf_refs(wall_ids)}}};',
            'Physical Volume("fluid") = {ext[1]};',
            "",
        ]
    )
    write(geom / "arc_fig21_anchor_exact_cfd.geo", geo)

    u = {
        "inlet_aux100": inlet_vec(outer, aux_i, TANK_SPEED),
        "inlet_ch_a": inlet_vec(inner, ch_a_i, CHANNEL_SPEED),
        "inlet_ch_b": inlet_vec(inner, ch_b_i, CHANNEL_SPEED),
        "inlet_main40": inlet_vec(inner, main_i, CHANNEL_SPEED),
    }
    rows = [
        "# ARC Fig. 2.1 Anchor-Exact CFD Patch Manifest",
        "",
        f"Source contours: `{contours}`",
        f"Patch target file: `{patch_targets_path}`",
        "",
        "| Patch | Curve IDs | Basis | Initial U [m/s] |",
        "|---|---:|---|---:|",
    ]
    bases = targets["basis"]
    for name in ("outlet180", "inlet_aux100", "inlet_ch_a", "inlet_main40", "inlet_ch_b"):
        vec = u.get(name, (0.0, 0.0, 0.0))
        rows.append(f"| `{name}` | `{patch_ids[name]}` | {bases[name]} | `({vec[0]:.6g} {vec[1]:.6g} {vec[2]:.6g})` |")
    rows += [
        "",
        "Patch placement follows the supplied target JSON. Change the patch target JSON, regenerate, and rerun `scripts/audit_arc2d_inlet_bc.py` after meshing.",
        "",
    ]
    write(geom / "arc_fig21_anchor_exact_cfd_patch_manifest.md", "\n".join(rows))
    return {"patch_ids": patch_ids, "u": u}


def write_case(case: Path, u: dict[str, tuple[float, float, float]], end_time: float, delta_t: float) -> None:
    write_interval = min(0.2, max(end_time / 4.0, delta_t))
    for sub in ("0", "constant", "system"):
        (case / sub).mkdir(parents=True, exist_ok=True)
    write(
        case / "constant/transportProperties",
        foam_header("dictionary", "transportProperties")
        + f"""transportModel Newtonian;
nu [0 2 -1 0 0 0 0] {NU:.12g};
// rho={RHO:.6g} kg/m3, mu={MU:.6g} Pa s.
""",
    )
    write(
        case / "constant/turbulenceProperties",
        foam_header("dictionary", "turbulenceProperties")
        + """simulationType LES;
LES
{
    LESModel WALE;
    turbulence on;
    printCoeffs on;
    delta cubeRootVol;
    cubeRootVolCoeffs { deltaCoeff 1; }
}
""",
    )

    def vec(name: str) -> str:
        v = u[name]
        return f"({v[0]:.12g} {v[1]:.12g} {v[2]:.12g})"

    write(
        case / "0/U",
        foam_header("volVectorField", "U")
        + f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inlet_aux100 {{ type fixedValue; value uniform {vec("inlet_aux100")}; }}
    inlet_ch_a   {{ type fixedValue; value uniform {vec("inlet_ch_a")}; }}
    inlet_main40 {{ type fixedValue; value uniform {vec("inlet_main40")}; }}
    inlet_ch_b   {{ type fixedValue; value uniform {vec("inlet_ch_b")}; }}
    outlet180    {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    walls        {{ type noSlip; }}
    frontAndBack {{ type empty; }}
}}
""",
    )
    write(
        case / "0/p",
        foam_header("volScalarField", "p")
        + """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    "(inlet_.*)" { type zeroGradient; }
    outlet180    { type fixedValue; value uniform 0; }
    walls        { type zeroGradient; }
    frontAndBack { type empty; }
}
""",
    )
    write(
        case / "0/nut",
        foam_header("volScalarField", "nut")
        + """dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{
    "(inlet_.*|outlet180)" { type calculated; value uniform 0; }
    walls        { type nutUSpaldingWallFunction; value uniform 0; }
    frontAndBack { type empty; }
}
""",
    )
    write(
        case / "system/controlDict",
        foam_header("dictionary", "controlDict")
        + f"""application pimpleFoam;
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {end_time:.12g};
deltaT {delta_t:.12g};
writeControl runTime;
writeInterval {write_interval:.12g};
purgeWrite 0;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
adjustTimeStep yes;
maxCo 0.7;
maxDeltaT {delta_t:.12g};
functions
{{
    fieldAverage1
    {{
        type fieldAverage;
        libs ("libfieldFunctionObjects.so");
        timeStart {max(end_time / 2.0, delta_t):.12g};
        writeControl writeTime;
        fields
        (
            U {{ mean on; prime2Mean on; base time; }}
            p {{ mean on; prime2Mean on; base time; }}
        );
    }}
}}
""",
    )
    write(
        case / "system/fvSchemes",
        foam_header("dictionary", "fvSchemes")
        + """ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss LUST grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""",
    )
    write(
        case / "system/fvSolution",
        foam_header("dictionary", "fvSolution")
        + """solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.05; }
    pFinal { $p; relTol 0; }
    "(U|k)(Final)?" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; }
}
PIMPLE
{
    nCorrectors 2;
    nNonOrthogonalCorrectors 1;
    nOuterCorrectors 1;
    pRefCell 0;
    pRefValue 0;
}
""",
    )
    write(
        case / "system/changeDictionaryDict",
        foam_header("dictionary", "changeDictionaryDict")
        + """boundary
{
    frontAndBack { type empty; }
    walls { type wall; }
}
""",
    )
    allpre = case / "Allpre.arc_fig21_cfd"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
gmsh -3 geometry/arc_fig21_anchor_exact_cfd.geo -format msh2 -o arc_fig21_anchor_exact_cfd.msh > log.gmsh 2>&1
gmshToFoam arc_fig21_anchor_exact_cfd.msh > log.gmshToFoam 2>&1
changeDictionary > log.changeDictionary 2>&1
checkMesh > log.checkMesh 2>&1
""",
    )
    chmod_x(allpre)
    allrun = case / "Allrun.arc_fig21_cfd"
    write(
        allrun,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d constant/polyMesh ]; then
    ./Allpre.arc_fig21_cfd
fi
pimpleFoam > log.pimpleFoam 2>&1
postProcess -func writeCellCentres -latestTime > log.cellCentres 2>&1 || true
touch case.foam
""",
    )
    chmod_x(allrun)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_fig21_anchor_exact_cfd_candidate")
    ap.add_argument("--contours", default="references/digitization/ferrero_arc/arc_fig21_anchor_exact_mesh_contours.json")
    ap.add_argument("--patch-targets", default=str(DEFAULT_PATCH_TARGETS))
    ap.add_argument("--mesh-size", type=float, default=0.045)
    ap.add_argument("--end-time", type=float, default=0.005)
    ap.add_argument("--delta-t", type=float, default=1e-4)
    args = ap.parse_args()
    case = Path(args.case)
    case.mkdir(parents=True, exist_ok=True)
    meta = write_geo(case, Path(args.contours), args.mesh_size, Path(args.patch_targets))
    write_case(case, meta["u"], args.end_time, args.delta_t)  # type: ignore[arg-type]
    print(f"wrote {case}")


if __name__ == "__main__":
    main()

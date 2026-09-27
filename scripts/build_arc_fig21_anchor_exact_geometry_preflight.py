#!/usr/bin/env python3
"""Build a geometry-only OpenFOAM mesh preflight from anchor-exact ARC contours."""

from __future__ import annotations

import argparse
import json
import stat
from pathlib import Path

from build_arc_blanket_2d_case import THICKNESS, foam_header, write


def chmod_x(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def mm_to_local_m(points: list[list[float]]) -> list[tuple[float, float, float]]:
    out = []
    for r_mm, z_mm in points:
        out.append(((r_mm - 3300.0) / 1000.0, z_mm / 1000.0, 0.0))
    return out


def poly_area(points: list[tuple[float, float, float]]) -> float:
    xy = [(p[0], p[1]) for p in points]
    return abs(
        0.5
        * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(xy, xy[1:] + xy[:1]))
    )


def write_geo(case: Path, contours: Path, mesh_size: float) -> None:
    data = json.loads(contours.read_text())
    outer = mm_to_local_m(data["tank_inner_wall"])
    inner = mm_to_local_m(data["vv_exclusion_union_contour"])
    geom = case / "geometry"
    geom.mkdir(parents=True, exist_ok=True)

    all_points = outer + inner
    lines = []
    for i, (x, y, z) in enumerate(all_points, 1):
        h = mesh_size * 0.55 if i > len(outer) else mesh_size
        lines.append(f"Point({i}) = {{{x:.8g}, {y:.8g}, {z:.8g}, {h:.8g}}};")

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

    all_boundary_ids = list(range(1, len(outer_ids) + len(inner_ids) + 1))

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
            f'Physical Surface("walls") = {{{surf_refs(all_boundary_ids)}}};',
            'Physical Volume("fluid") = {ext[1]};',
            "",
        ]
    )
    write(geom / "arc_fig21_anchor_exact.geo", geo)

    outer_area = poly_area(outer)
    inner_area = poly_area(inner)
    fluid_area = outer_area - inner_area
    write(
        geom / "arc_fig21_anchor_exact_preflight_manifest.md",
        "\n".join(
            [
                "# ARC Fig. 2.1 Anchor-Exact Geometry Preflight Manifest",
                "",
                f"Source contours: `{contours}`",
                "",
                "This is a geometry-only mesh preflight. It has no inlet/outlet patch split yet.",
                "",
                f"- outer points: `{len(outer)}`",
                f"- inner contour points: `{len(inner)}`",
                f"- outer area: `{outer_area:.4f} m2`",
                f"- inner exclusion area: `{inner_area:.4f} m2`",
                f"- fluid area: `{fluid_area:.4f} m2`",
                "",
            ]
        ),
    )


def write_preflight_case(case: Path) -> None:
    for sub in ("0", "constant", "system"):
        (case / sub).mkdir(parents=True, exist_ok=True)
    write(
        case / "system/controlDict",
        foam_header("dictionary", "controlDict")
        + """application simpleFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1;
deltaT 1;
writeControl timeStep;
writeInterval 1;
purgeWrite 0;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
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
    write(
        case / "system/fvSchemes",
        foam_header("dictionary", "fvSchemes")
        + """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
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
    p { solver GAMG; tolerance 1e-7; relTol 0.1; smoother GaussSeidel; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 0;
}
""",
    )
    allpre = case / "Allpre.arc_fig21"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
gmsh -3 geometry/arc_fig21_anchor_exact.geo -format msh2 -o arc_fig21_anchor_exact.msh > log.gmsh 2>&1
gmshToFoam arc_fig21_anchor_exact.msh > log.gmshToFoam 2>&1
changeDictionary > log.changeDictionary 2>&1
checkMesh > log.checkMesh 2>&1
touch case.foam
""",
    )
    chmod_x(allpre)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_fig21_anchor_exact_geometry_preflight")
    ap.add_argument(
        "--contours",
        default="references/digitization/ferrero_arc/arc_fig21_anchor_exact_mesh_contours.json",
    )
    ap.add_argument("--mesh-size", type=float, default=0.045)
    args = ap.parse_args()

    case = Path(args.case)
    case.mkdir(parents=True, exist_ok=True)
    write_geo(case, Path(args.contours), args.mesh_size)
    write_preflight_case(case)
    print(f"wrote {case}")


if __name__ == "__main__":
    main()

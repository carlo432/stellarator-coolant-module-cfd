#!/usr/bin/env python3
"""Build a local OpenFOAM approximation of Leffler's ARC blanket 2D case.

The geometry is a slide-traced approximation of the ARC blanket tank in
Leffler_Christopher_10705.pptx, not Ferrero's original CAD.  The numerical
values that are explicit in the deck are preserved: 20/40/100/180 mm inlet and
outlet widths, FLiBe properties at 834 K, and U_ref = 2.7 m/s.
"""
from __future__ import annotations

import argparse
import math
import os
import stat
from pathlib import Path


RHO = 2005.781
MU = 0.01047
NU = MU / RHO
U_REF = 2.7
L_REF = 0.02
THICKNESS = 0.01


def foam_header(cls: str, obj: str) -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*/\n"
        f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
    )


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def chmod_x(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def dist(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def add_polyline(points: list[tuple[float, float, float]], start_line: int) -> tuple[list[str], list[int]]:
    lines = []
    ids = []
    n = len(points)
    for i in range(n):
        lid = start_line + i
        lines.append(f"Line({lid}) = {{{i + 1}, {(i + 1) % n + 1}}};")
        ids.append(lid)
    return lines, ids


def build_geo(case: Path, mesh_size: float) -> None:
    geom = case / "geometry"
    geom.mkdir(parents=True, exist_ok=True)

    # Outer tank outline, traced from Leffler slide 8/10.  Coordinates are in m
    # in the 2D cross-section plane.  The outline preserves the broad Ferrero
    # tank silhouette, not every CAD fillet.
    outer = [
        (-1.30, 1.94, mesh_size),
        (0.10, 1.94, mesh_size),
        (0.28, 1.94, mesh_size),  # 180 mm outlet segment before this point
        (1.05, 1.94, mesh_size),
        (1.42, 1.74, mesh_size),
        (1.42, 1.64, mesh_size * 0.5),  # 100 mm auxiliary inlet segment
        (1.55, 1.55, mesh_size),
        (1.55, 0.62, mesh_size),
        (1.88, 0.62, mesh_size),
        (2.15, 0.35, mesh_size),
        (2.35, -0.05, mesh_size),
        (2.35, -0.25, mesh_size),
        (2.15, -0.55, mesh_size),
        (1.86, -0.75, mesh_size),
        (1.55, -0.75, mesh_size),
        (1.55, -1.72, mesh_size),
        (1.20, -1.94, mesh_size),
        (-0.75, -1.94, mesh_size),
        (-1.08, -1.72, mesh_size),
        (-1.55, -1.72, mesh_size),
        (-1.55, -0.62, mesh_size),
        (-1.88, -0.62, mesh_size),
        (-1.88, 0.62, mesh_size),
        (-1.55, 0.62, mesh_size),
        (-1.55, 1.55, mesh_size),
    ]

    # Inner vessel/divertor exclusion.  This is intentionally polyline-based so
    # the trace can be edited later from slide measurements.  Small line
    # segments with inlet labels are included where Leffler marks the channel
    # inlets.
    inner = [
        (0.30, 1.72, mesh_size),
        (0.52, 1.60, mesh_size),
        (0.56, 1.38, mesh_size),
        (0.43, 1.18, mesh_size),
        (0.22, 1.12, mesh_size * 0.35),
        (0.22, 1.10, mesh_size * 0.35),  # 20 mm inlet_ch_a
        (0.05, 0.80, mesh_size),
        (0.40, 0.55, mesh_size),
        (0.68, 0.18, mesh_size),
        (0.78, -0.25, mesh_size),
        (0.68, -0.65, mesh_size),
        (0.43, -0.95, mesh_size),
        (0.25, -1.16, mesh_size * 0.35),
        (0.27, -1.17, mesh_size * 0.35),  # 20 mm inlet_ch_b
        (0.45, -1.20, mesh_size),
        (0.65, -1.36, mesh_size),
        (0.66, -1.60, mesh_size),
        (0.50, -1.78, mesh_size),
        (0.25, -1.82, mesh_size),
        (0.06, -1.66, mesh_size),
        (0.02, -1.42, mesh_size),
        (0.12, -1.25, mesh_size),
        (-0.39, -1.05, mesh_size * 0.5),
        (-0.35, -1.05, mesh_size * 0.5),  # 40 mm inlet_main40 group
        (-0.65, -0.65, mesh_size),
        (-0.78, -0.20, mesh_size),
        (-0.75, 0.35, mesh_size),
        (-0.55, 0.75, mesh_size),
        (-0.25, 0.98, mesh_size),
        (0.05, 1.05, mesh_size),
    ]

    all_points = outer + inner
    lines = []
    for i, (x, y, h) in enumerate(all_points, 1):
        lines.append(f"Point({i}) = {{{x:.8g}, {y:.8g}, 0, {h:.8g}}};")

    outer_lines = []
    for i in range(len(outer)):
        lid = i + 1
        lines.append(f"Line({lid}) = {{{i + 1}, {(i + 1) % len(outer) + 1}}};")
        outer_lines.append(lid)

    offset = len(outer)
    inner_lines = []
    for i in range(len(inner)):
        lid = len(outer) + i + 1
        a = offset + i + 1
        b = offset + ((i + 1) % len(inner)) + 1
        # The hole loop is reversed below in the curve loop.
        lines.append(f"Line({lid}) = {{{a}, {b}}};")
        inner_lines.append(lid)

    n_outer = len(outer_lines)

    # Labels by curve index in the combined boundary order.  Line ids are
    # one-based; ext[2 + zero_based_curve_index] is the lateral extruded surface.
    outlet180 = [2]          # outer point 2 -> 3, length 0.18 m
    inlet_aux100 = [5]       # outer point 5 -> 6, length 0.10 m
    inlet_ch_a = [n_outer + 5]      # inner point 5 -> 6, length ~0.02 m
    inlet_ch_b = [n_outer + 13]     # inner point 13 -> 14, length ~0.02 m
    inlet_main40 = [n_outer + 23]  # one 40 mm segment

    inlet_or_outlet = set(outlet180 + inlet_aux100 + inlet_ch_a + inlet_main40 + inlet_ch_b)
    all_curve_indices = list(range(1, n_outer + len(inner_lines) + 1))
    wall_indices = [idx for idx in all_curve_indices if idx not in inlet_or_outlet]

    def surf_refs(curve_indices: list[int]) -> str:
        return ", ".join(f"ext[{2 + (idx - 1)}]" for idx in curve_indices)

    geo = "\n".join(
        [
            'SetFactory("Built-in");',
            f"Mesh.CharacteristicLengthMin = {mesh_size * 0.3:.8g};",
            f"Mesh.CharacteristicLengthMax = {mesh_size:.8g};",
            "Mesh.Algorithm = 6;",
            "Mesh.RecombineAll = 0;",
            *lines,
            f"Curve Loop(1) = {{{', '.join(str(i) for i in outer_lines)}}};",
            f"Curve Loop(2) = {{{', '.join(str(-i) for i in reversed(inner_lines))}}};",
            "Plane Surface(1) = {1, 2};",
            f"ext[] = Extrude {{0, 0, {THICKNESS:.8g}}} {{ Surface{{1}}; Layers{{1}}; Recombine; }};",
            'Physical Surface("frontAndBack") = {1, ext[0]};',
            f'Physical Surface("outlet180") = {{{surf_refs(outlet180)}}};',
            f'Physical Surface("inlet_aux100") = {{{surf_refs(inlet_aux100)}}};',
            f'Physical Surface("inlet_ch_a") = {{{surf_refs(inlet_ch_a)}}};',
            f'Physical Surface("inlet_main40") = {{{surf_refs(inlet_main40)}}};',
            f'Physical Surface("inlet_ch_b") = {{{surf_refs(inlet_ch_b)}}};',
            f'Physical Surface("walls") = {{{surf_refs(wall_indices)}}};',
            'Physical Volume("fluid") = {ext[1]};',
            "",
        ]
    )
    write(geom / "arc_blanket_2d.geo", geo)

    manifest = [
        "# ARC 2D Geometry Trace Manifest",
        "",
        "Source: Leffler slide deck pages 8-12.",
        "",
        "This is a slide-traced approximation, not Ferrero CAD.",
        "",
        "| Patch | Intended deck width | Trace length [m] |",
        "|---|---:|---:|",
        f"| outlet180 | 0.180 | {dist(outer[1], outer[2]):.6f} |",
        f"| inlet_aux100 | 0.100 | {dist(outer[4], outer[5]):.6f} |",
        f"| inlet_ch_a | 0.020 | {dist(inner[4], inner[5]):.6f} |",
        f"| inlet_main40 | 0.040 | {dist(inner[22], inner[23]):.6f} |",
        f"| inlet_ch_b | 0.020 | {dist(inner[12], inner[13]):.6f} |",
        "",
    ]
    write(geom / "arc_blanket_2d_trace_manifest.md", "\n".join(manifest))


def inlet_vectors(inlet_mode: str) -> dict[str, tuple[float, float, float]]:
    if inlet_mode == "legacy-global":
        return {
            "inlet_aux100": (0.0, -U_REF, 0.0),
            "inlet_ch_a": (U_REF, 0.0, 0.0),
            "inlet_main40": (U_REF, 0.0, 0.0),
            "inlet_ch_b": (U_REF, 0.0, 0.0),
        }
    if inlet_mode == "normal":
        # Inward patch-normal vectors for the current slide-traced geometry.
        # These are derived from the OpenFOAM patch face-area normals and are
        # documented by scripts/audit_arc2d_inlet_bc.py.
        return {
            "inlet_aux100": (-U_REF, 0.0, 0.0),
            "inlet_ch_a": (U_REF, 0.0, 0.0),
            "inlet_main40": (0.0, U_REF, 0.0),
            "inlet_ch_b": (U_REF / math.sqrt(5.0), 2.0 * U_REF / math.sqrt(5.0), 0.0),
        }
    if inlet_mode == "ferrero-average-normal":
        # Source-informed diagnostic only.  Ferrero's thesis summary reports
        # average tank velocity 0.313 m/s and average channel velocity
        # 1.98 m/s.  These are not exact Leffler/Nek5000 inlet data, but they
        # are a useful sensitivity against the equal-2.7 m/s normal-inlet mode.
        tank_u = 0.313
        channel_u = 1.98
        return {
            "inlet_aux100": (-tank_u, 0.0, 0.0),
            "inlet_ch_a": (channel_u, 0.0, 0.0),
            "inlet_main40": (0.0, channel_u, 0.0),
            "inlet_ch_b": (channel_u / math.sqrt(5.0), 2.0 * channel_u / math.sqrt(5.0), 0.0),
        }
    raise ValueError(f"unknown inlet mode {inlet_mode}")


def write_openfoam_case(case: Path, end_time: float, delta_t: float, inlet_mode: str) -> None:
    for sub in ("0", "constant", "system"):
        (case / sub).mkdir(parents=True, exist_ok=True)

    write(
        case / "constant/transportProperties",
        foam_header("dictionary", "transportProperties")
        + f"""transportModel Newtonian;
nu [0 2 -1 0 0 0 0] {NU:.12g};
// Leffler deck: rho={RHO:.6g} kg/m3, mu={MU:.6g} Pa s at 834 K, Uref={U_REF} m/s, Lref={L_REF} m.
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
    u_in = inlet_vectors(inlet_mode)

    def vec(name: str) -> str:
        u = u_in[name]
        return f"({u[0]:.12g} {u[1]:.12g} {u[2]:.12g})"

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
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {end_time:.12g};
deltaT {delta_t:.12g};
writeControl runTime;
writeInterval {max(end_time / 4.0, delta_t):.12g};
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
        restartOnOutput false;
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
        case / "system/decomposeParDict",
        foam_header("dictionary", "decomposeParDict")
        + """numberOfSubdomains 4;
method scotch;
""",
    )
    write(
        case / "system/changeDictionaryDict",
        foam_header("dictionary", "changeDictionaryDict")
        + """boundary
{
    frontAndBack
    {
        type empty;
    }
    walls
    {
        type wall;
    }
}
""",
    )

    allpre = case / "Allpre.arc2d"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
gmsh -3 geometry/arc_blanket_2d.geo -format msh2 -o arc_blanket_2d.msh > log.gmsh 2>&1
gmshToFoam arc_blanket_2d.msh > log.gmshToFoam 2>&1
changeDictionary > log.changeDictionary 2>&1
checkMesh > log.checkMesh 2>&1
""",
    )
    chmod_x(allpre)

    allrun = case / "Allrun.arc2d"
    write(
        allrun,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d constant/polyMesh ]; then
    ./Allpre.arc2d
fi
pimpleFoam > log.pimpleFoam 2>&1
postProcess -func writeCellCentres -latestTime > log.cellCentres 2>&1 || true
""",
    )
    chmod_x(allrun)

    inlet_lines = "\n".join(
        f"- `{name}`: `({u[0]:.8g} {u[1]:.8g} {u[2]:.8g}) m/s`"
        for name, u in u_in.items()
    )

    write(
        case / "README.arc_blanket_2d.md",
        f"""# ARC Blanket 2D OpenFOAM Case

Purpose: local OpenFOAM reproduction-style case for Leffler's ARC blanket 2D
velocity-only LES workflow.

This is not a Nek5000 reproduction and not Ferrero CAD.  It is a slide-traced
OpenFOAM approximation using the explicit deck values:

- `L_ref = {L_REF} m`
- `U_ref = {U_REF} m/s`
- `rho = {RHO} kg/m3`
- `mu = {MU} Pa s`
- `nu = {NU:.8g} m2/s`
- `Re = rho*U_ref*L_ref/mu = {RHO * U_REF * L_REF / MU:.2f}`
- reference velocity for nondimensionalization: `{U_REF} m/s`
- inlet vector mode: `{inlet_mode}`
- outlet pressure: `0` kinematic gauge pressure in OpenFOAM's incompressible form

Inlet velocity vectors:

{inlet_lines}

The full Leffler run reports fields at `t* = 7566.2`.  With
`t* = U_ref*t/L_ref`, that corresponds to about
`{7566.2 * L_REF / U_REF:.2f} s`.  This local case starts as a short smoke run;
extend `system/controlDict` for production-style averaging.
""",
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d")
    ap.add_argument("--mesh-size", type=float, default=0.022)
    ap.add_argument("--end-time", type=float, default=0.01)
    ap.add_argument("--delta-t", type=float, default=1.0e-4)
    ap.add_argument(
        "--inlet-mode",
        choices=("legacy-global", "normal", "ferrero-average-normal"),
        default="legacy-global",
    )
    args = ap.parse_args()

    case = Path(args.case)
    case.mkdir(parents=True, exist_ok=True)
    build_geo(case, args.mesh_size)
    write_openfoam_case(case, args.end_time, args.delta_t, args.inlet_mode)
    print(f"wrote {case}")
    print(f"nu={NU:.8g} Re={RHO * U_REF * L_REF / MU:.2f}")
    print(f"inlet_mode={args.inlet_mode}")


if __name__ == "__main__":
    main()

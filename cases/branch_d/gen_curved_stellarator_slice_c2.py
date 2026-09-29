#!/usr/bin/env python3
"""C2: compact-FLiBe-stellarator curved-slice OpenFOAM scaffold.

This creates a multi-block curved duct with a gently helical/wavy centerline and
stellarator-relevant patch names.  It is intentionally still a local module, not
real stellarator CAD.  The first proof is blockMesh/checkMesh plus a render.
"""
from __future__ import annotations

import argparse
import json
import math
import stat
from pathlib import Path


def foam_header(cls: str, obj: str) -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*/\n"
        f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
    )


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def vadd(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vscale(a: tuple[float, float, float], s: float) -> tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def fmt(v: tuple[float, float, float]) -> str:
    return f"({v[0]:.10g} {v[1]:.10g} {v[2]:.10g})"


def centerline(t: float, radius: float, arc: float, y_amp: float) -> tuple[float, float, float]:
    theta = arc * t
    return (
        radius * math.sin(theta),
        y_amp * math.sin(2.0 * math.pi * t),
        radius * (1.0 - math.cos(theta)),
    )


def radial_normal(t: float, arc: float) -> tuple[float, float, float]:
    theta = arc * t
    # Points approximately from plasma-facing first wall toward blanket back wall.
    return (math.sin(theta), 0.0, -math.cos(theta))


def station_vertices(
    t: float,
    radius: float,
    arc: float,
    y_amp: float,
    radial_width: float,
    height: float,
) -> list[tuple[float, float, float]]:
    c = centerline(t, radius, arc, y_amp)
    n = radial_normal(t, arc)
    ey = (0.0, 1.0, 0.0)
    # A/D are heated first-wall side; B/C are blanket-back side.
    a = vadd(vadd(c, vscale(n, -0.5 * radial_width)), vscale(ey, -0.5 * height))
    b = vadd(vadd(c, vscale(n, 0.5 * radial_width)), vscale(ey, -0.5 * height))
    cpt = vadd(vadd(c, vscale(n, 0.5 * radial_width)), vscale(ey, 0.5 * height))
    d = vadd(vadd(c, vscale(n, -0.5 * radial_width)), vscale(ey, 0.5 * height))
    return [a, b, cpt, d]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="c2_curved_stellarator_slice_preflight")
    ap.add_argument("--segments", type=int, default=24)
    ap.add_argument("--stream-cells", type=int, default=4)
    ap.add_argument("--radial-cells", type=int, default=10)
    ap.add_argument("--height-cells", type=int, default=8)
    ap.add_argument("--radius", type=float, default=0.18)
    ap.add_argument("--arc-deg", type=float, default=95.0)
    ap.add_argument("--y-amp", type=float, default=0.010)
    ap.add_argument("--radial-width", type=float, default=0.030)
    ap.add_argument("--height", type=float, default=0.016)
    ap.add_argument("--u-in", type=float, default=1.1598)
    ap.add_argument("--nu", type=float, default=3.09278e-6)
    ap.add_argument("--iterations", type=int, default=120)
    args = ap.parse_args()

    case = Path(args.case)
    arc = math.radians(args.arc_deg)
    stations = args.segments + 1
    verts: list[tuple[float, float, float]] = []
    for i in range(stations):
        verts.extend(
            station_vertices(
                i / args.segments,
                args.radius,
                arc,
                args.y_amp,
                args.radial_width,
                args.height,
            )
        )

    def idx(station: int, corner: int) -> int:
        return 4 * station + corner

    blocks = []
    heated_faces = []
    back_faces = []
    side_faces = []
    for i in range(args.segments):
        # Block order maps local x=stream, y=height, z=back->heated.
        # Using back-side vertices as z=0 keeps the hex right-handed at t=0.
        # 0 B_i, 1 B_j, 2 C_j, 3 C_i, 4 A_i, 5 A_j, 6 D_j, 7 D_i.
        v0, v1 = idx(i, 1), idx(i + 1, 1)
        v2, v3 = idx(i + 1, 2), idx(i, 2)
        v4, v5 = idx(i, 0), idx(i + 1, 0)
        v6, v7 = idx(i + 1, 3), idx(i, 3)
        blocks.append((v0, v1, v2, v3, v4, v5, v6, v7))
        back_faces.append((v0, v1, v2, v3))
        heated_faces.append((v4, v7, v6, v5))
        side_faces.append((v0, v4, v5, v1))
        side_faces.append((v3, v2, v6, v7))

    inlet_face = (idx(0, 1), idx(0, 2), idx(0, 3), idx(0, 0))
    outlet_face = (
        idx(args.segments, 1),
        idx(args.segments, 0),
        idx(args.segments, 3),
        idx(args.segments, 2),
    )
    cells = args.segments * args.stream_cells * args.radial_cells * args.height_cells

    vertices_txt = "\n".join(f"    {fmt(v)}" for v in verts)
    blocks_txt = "\n".join(
        "    hex (%s) (%d %d %d) simpleGrading (1 1 1)"
        % (" ".join(str(x) for x in b), args.stream_cells, args.height_cells, args.radial_cells)
        for b in blocks
    )

    def faces_txt(faces: list[tuple[int, ...]]) -> str:
        return "\n".join("        (%s)" % " ".join(str(x) for x in face) for face in faces)

    write(
        case / "system/blockMeshDict",
        foam_header("dictionary", "blockMeshDict")
        + f"""scale 1;
vertices
(
{vertices_txt}
);
blocks
(
{blocks_txt}
);
edges ();
boundary
(
    inlet
    {{
        type patch;
        faces
        (
            ({' '.join(str(x) for x in inlet_face)})
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            ({' '.join(str(x) for x in outlet_face)})
        );
    }}
    heated_first_wall
    {{
        type wall;
        faces
        (
{faces_txt(heated_faces)}
        );
    }}
    blanket_back_wall
    {{
        type wall;
        faces
        (
{faces_txt(back_faces)}
        );
    }}
    side_walls
    {{
        type wall;
        faces
        (
{faces_txt(side_faces)}
        );
    }}
);
mergePatchPairs ();
""",
    )

    for sub in ("0", "constant", "system"):
        (case / sub).mkdir(parents=True, exist_ok=True)
    write(
        case / "constant/transportProperties",
        foam_header("dictionary", "transportProperties")
        + f"""transportModel Newtonian;
nu [0 2 -1 0 0 0 0] {args.nu:.12g};
DT [0 2 -1 0 0 0 0] {args.nu / 14.4:.12g};
""",
    )
    write(
        case / "constant/turbulenceProperties",
        foam_header("dictionary", "turbulenceProperties")
        + """simulationType RAS;
RAS
{
    RASModel kOmegaSST;
    turbulence on;
    printCoeffs on;
}
""",
    )
    write(
        case / "system/controlDict",
        foam_header("dictionary", "controlDict")
        + f"""application simpleFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {args.iterations};
deltaT 1;
writeControl timeStep;
writeInterval {max(args.iterations, 1)};
purgeWrite 0;
writeFormat ascii;
writePrecision 6;
writeCompression off;
timeFormat general;
timePrecision 6;
runTimeModifiable true;
functions
{{
    yPlus1
    {{
        type yPlus;
        libs ("libfieldFunctionObjects.so");
        writeControl writeTime;
    }}
}}
""",
    )
    write(
        case / "system/fvSchemes",
        foam_header("dictionary", "fvSchemes")
        + """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    div(phi,U) bounded Gauss upwind;
    div(phi,k) bounded Gauss upwind;
    div(phi,omega) bounded Gauss upwind;
    div(phi,nut) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
""",
    )
    write(
        case / "system/fvSolution",
        foam_header("dictionary", "fvSolution")
        + """solvers
{
    p { solver GAMG; tolerance 1e-7; relTol 0.01; smoother GaussSeidel; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
    "(k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;
}
relaxationFactors
{
    fields { p 0.3; }
    equations { U 0.7; k 0.7; omega 0.7; }
}
""",
    )
    write(
        case / "0/U",
        foam_header("volVectorField", "U")
        + f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform ({args.u_in:.12g} 0 0);
boundaryField
{{
    inlet {{ type fixedValue; value uniform ({args.u_in:.12g} 0 0); }}
    outlet {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    "(heated_first_wall|blanket_back_wall|side_walls)" {{ type noSlip; }}
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
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    "(heated_first_wall|blanket_back_wall|side_walls)" { type zeroGradient; }
}
""",
    )
    write(
        case / "0/k",
        foam_header("volScalarField", "k")
        + """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0.00505;
boundaryField
{
    inlet { type fixedValue; value uniform 0.00505; }
    outlet { type inletOutlet; inletValue uniform 0.00505; value uniform 0.00505; }
    "(heated_first_wall|blanket_back_wall|side_walls)" { type kqRWallFunction; value uniform 0.00505; }
}
""",
    )
    write(
        case / "0/omega",
        foam_header("volScalarField", "omega")
        + """dimensions [0 0 -1 0 0 0 0];
internalField uniform 24.0;
boundaryField
{
    inlet { type fixedValue; value uniform 24.0; }
    outlet { type inletOutlet; inletValue uniform 24.0; value uniform 24.0; }
    "(heated_first_wall|blanket_back_wall|side_walls)" { type omegaWallFunction; value uniform 24.0; }
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
    inlet { type calculated; value uniform 0; }
    outlet { type calculated; value uniform 0; }
    "(heated_first_wall|blanket_back_wall|side_walls)" { type nutkWallFunction; value uniform 0; }
}
""",
    )
    write(
        case / "0/T",
        foam_header("volScalarField", "T")
        + """dimensions [0 0 0 1 0 0 0];
internalField uniform 800;
boundaryField
{
    inlet { type fixedValue; value uniform 800; }
    outlet { type inletOutlet; inletValue uniform 800; value uniform 800; }
    heated_first_wall
    {
        type fixedGradient;
        gradient uniform 100000;
        value uniform 800;
    }
    "(blanket_back_wall|side_walls)" { type zeroGradient; }
}
""",
    )

    allpre = case / "Allpre.c2"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1
""",
    )
    allpre.chmod(allpre.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    write(
        case / "geometry_c2.json",
        json.dumps(
            {
                "vertices": verts,
                "blocks": blocks,
                "patches": {
                    "inlet": [inlet_face],
                    "outlet": [outlet_face],
                    "heated_first_wall": heated_faces,
                    "blanket_back_wall": back_faces,
                    "side_walls": side_faces,
                },
                "centerline": [
                    centerline(i / args.segments, args.radius, arc, args.y_amp)
                    for i in range(stations)
                ],
                "cells": cells,
                "parameters": vars(args),
            },
            indent=2,
        ),
    )
    write(
        case / "README.c2.md",
        f"""# C2 Curved Stellarator Slice Preflight

Generated by `gen_curved_stellarator_slice_c2.py`.

- cells: `{cells}`
- segments: `{args.segments}`
- arc: `{args.arc_deg} deg`
- radial coolant-pocket width: `{args.radial_width} m`
- vertical helical/wavy amplitude: `{args.y_amp} m`

This is a geometry/RANS scaffold, not a full stellarator blanket and not CHT.
Use `./Allpre.c2` for local `blockMesh` / `checkMesh`, then run
`simpleFoam > log.simpleFoam 2>&1` for the isothermal RANS smoke.  From the
parent folder, run `python3 postprocess_curved_stellarator_slice_c2.py {case.name}`
to extract owner-cell pressure-drop, heated-wall-adjacent speed, and y+ hooks.
""",
    )
    print(f"wrote {case}/ cells={cells} segments={args.segments} arc={args.arc_deg:g}deg")


if __name__ == "__main__":
    main()

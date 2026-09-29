#!/usr/bin/env python3
"""C2 straight-control OpenFOAM scaffold.

This matches the curved C2 slice in cross-section, cell count, inlet speed, and
heat-flux-style scalar boundary condition.  It is the first control case for the
curved-vs-straight flushing-on-curvature diagnostic.
"""
from __future__ import annotations

import argparse
import json
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


def fmt(v: tuple[float, float, float]) -> str:
    return f"({v[0]:.10g} {v[1]:.10g} {v[2]:.10g})"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="c2_straight_control_preflight")
    ap.add_argument("--length", type=float, default=0.3016540039235351)
    ap.add_argument("--radial-width", type=float, default=0.030)
    ap.add_argument("--height", type=float, default=0.016)
    ap.add_argument("--stream-cells", type=int, default=96)
    ap.add_argument("--radial-cells", type=int, default=10)
    ap.add_argument("--height-cells", type=int, default=8)
    ap.add_argument("--u-in", type=float, default=1.1598)
    ap.add_argument("--nu", type=float, default=3.09278e-6)
    ap.add_argument("--iterations", type=int, default=120)
    args = ap.parse_args()

    case = Path(args.case)
    l = args.length
    w = args.radial_width
    h = args.height

    # Same block vertex ordering as the C2 curved generator:
    # B_in, B_out, C_out, C_in, A_in, A_out, D_out, D_in.
    verts = [
        (0.0, -0.5 * h, -0.5 * w),
        (l, -0.5 * h, -0.5 * w),
        (l, 0.5 * h, -0.5 * w),
        (0.0, 0.5 * h, -0.5 * w),
        (0.0, -0.5 * h, 0.5 * w),
        (l, -0.5 * h, 0.5 * w),
        (l, 0.5 * h, 0.5 * w),
        (0.0, 0.5 * h, 0.5 * w),
    ]
    block = (0, 1, 2, 3, 4, 5, 6, 7)
    patches = {
        "inlet": [(0, 3, 7, 4)],
        "outlet": [(1, 5, 6, 2)],
        "heated_first_wall": [(4, 7, 6, 5)],
        "blanket_back_wall": [(0, 1, 2, 3)],
        "side_walls": [(0, 4, 5, 1), (3, 2, 6, 7)],
    }
    cells = args.stream_cells * args.height_cells * args.radial_cells

    vertices_txt = "\n".join(f"    {fmt(v)}" for v in verts)

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
    hex ({' '.join(str(x) for x in block)}) ({args.stream_cells} {args.height_cells} {args.radial_cells}) simpleGrading (1 1 1)
);
edges ();
boundary
(
    inlet
    {{
        type patch;
        faces
        (
{faces_txt(patches['inlet'])}
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
{faces_txt(patches['outlet'])}
        );
    }}
    heated_first_wall
    {{
        type wall;
        faces
        (
{faces_txt(patches['heated_first_wall'])}
        );
    }}
    blanket_back_wall
    {{
        type wall;
        faces
        (
{faces_txt(patches['blanket_back_wall'])}
        );
    }}
    side_walls
    {{
        type wall;
        faces
        (
{faces_txt(patches['side_walls'])}
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
                "blocks": [block],
                "patches": patches,
                "centerline": [(0.0, 0.0, 0.0), (l, 0.0, 0.0)],
                "cells": cells,
                "parameters": vars(args),
            },
            indent=2,
        ),
    )
    write(
        case / "README.c2.md",
        f"""# C2 Straight Control Preflight

Generated by `gen_straight_control_c2.py`.

- cells: `{cells}`
- length: `{args.length} m`
- radial coolant-pocket width: `{args.radial_width} m`
- height: `{args.height} m`

This is a same-scale control for the C2 curved slice, not a reactor CAD model.
Use `./Allpre.c2`, then run `simpleFoam > log.simpleFoam 2>&1`, then run the
C2 postprocessor from the parent folder.
""",
    )
    print(f"wrote {case}/ cells={cells} length={args.length:g}m")


if __name__ == "__main__":
    main()

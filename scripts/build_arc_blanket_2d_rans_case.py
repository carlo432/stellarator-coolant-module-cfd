#!/usr/bin/env python3
"""Build a steady RANS baseline for the local Leffler ARC 2D case."""
from __future__ import annotations

import argparse
import stat
from pathlib import Path

from build_arc_blanket_2d_case import L_REF, MU, NU, RHO, U_REF, build_geo, foam_header, write


TURB_INTENSITY = 0.05
K0 = 1.5 * (TURB_INTENSITY * U_REF) ** 2
CMU_025 = 0.09**0.25
L_TURB = 0.07 * L_REF
OMEGA0 = (K0**0.5) / (CMU_025 * L_TURB)


def chmod_x(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_rans_case(case: Path, end_time: int) -> None:
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
        case / "0/U",
        foam_header("volVectorField", "U")
        + f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inlet_aux100 {{ type fixedValue; value uniform (0 {-U_REF:.12g} 0); }}
    inlet_ch_a   {{ type fixedValue; value uniform ({U_REF:.12g} 0 0); }}
    inlet_main40 {{ type fixedValue; value uniform ({U_REF:.12g} 0 0); }}
    inlet_ch_b   {{ type fixedValue; value uniform ({U_REF:.12g} 0 0); }}
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
        case / "0/k",
        foam_header("volScalarField", "k")
        + f"""dimensions [0 2 -2 0 0 0 0];
internalField uniform {K0:.12g};
boundaryField
{{
    "(inlet_.*)" {{ type fixedValue; value uniform {K0:.12g}; }}
    outlet180    {{ type inletOutlet; inletValue uniform {K0:.12g}; value uniform {K0:.12g}; }}
    walls        {{ type kqRWallFunction; value uniform {K0:.12g}; }}
    frontAndBack {{ type empty; }}
}}
""",
    )
    write(
        case / "0/omega",
        foam_header("volScalarField", "omega")
        + f"""dimensions [0 0 -1 0 0 0 0];
internalField uniform {OMEGA0:.12g};
boundaryField
{{
    "(inlet_.*)" {{ type fixedValue; value uniform {OMEGA0:.12g}; }}
    outlet180    {{ type inletOutlet; inletValue uniform {OMEGA0:.12g}; value uniform {OMEGA0:.12g}; }}
    walls        {{ type omegaWallFunction; value uniform {OMEGA0:.12g}; }}
    frontAndBack {{ type empty; }}
}}
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
    walls        { type nutkWallFunction; value uniform 0; }
    frontAndBack { type empty; }
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
endTime {end_time};
deltaT 1;
writeControl timeStep;
writeInterval {max(50, end_time // 4)};
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
        case / "system/fvSchemes",
        foam_header("dictionary", "fvSchemes")
        + """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss upwind;
    div(phi,k) bounded Gauss upwind;
    div(phi,omega) bounded Gauss upwind;
    div(phi,nut) bounded Gauss upwind;
    div(nuEff*dev2(T(grad(U)))) Gauss linear;
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
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.01; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
    "(k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;
    pRefCell 0;
    pRefValue 0;
}
relaxationFactors
{
    fields { p 0.3; }
    equations { U 0.7; k 0.7; omega 0.7; }
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

    allpre = case / "Allpre.arc2d_rans"
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

    allrun = case / "Allrun.arc2d_rans"
    write(
        allrun,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d constant/polyMesh ]; then
    ./Allpre.arc2d_rans
fi
simpleFoam > log.simpleFoam 2>&1
postProcess -func writeCellCentres -latestTime > log.cellCentres 2>&1 || true
touch case.foam
""",
    )
    chmod_x(allrun)

    write(
        case / "README.arc_blanket_2d_rans.md",
        f"""# ARC Blanket 2D RANS Baseline

Purpose: steady `simpleFoam`/`kOmegaSST` baseline on the same slide-traced
Leffler ARC 2D geometry as `cases/arc_blanket_2d`.

This is for RANS-vs-unsteady visual comparison only. It is not Leffler's
Nek5000 method and not Ferrero CAD.

- `U_ref = {U_REF} m/s`
- `nu = {NU:.8g} m2/s`
- `Re = {RHO * U_REF * L_REF / MU:.2f}`
- turbulence intensity seed: `{TURB_INTENSITY:.2%}`
- `k0 = {K0:.6g} m2/s2`
- `omega0 = {OMEGA0:.6g} 1/s`
""",
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d_rans_baseline")
    ap.add_argument("--mesh-size", type=float, default=0.022)
    ap.add_argument("--end-time", type=int, default=300)
    args = ap.parse_args()

    case = Path(args.case)
    case.mkdir(parents=True, exist_ok=True)
    build_geo(case, args.mesh_size)
    write_rans_case(case, args.end_time)
    print(f"wrote {case}")
    print(f"k0={K0:.6g} omega0={OMEGA0:.6g} Re={RHO * U_REF * L_REF / MU:.2f}")


if __name__ == "__main__":
    main()

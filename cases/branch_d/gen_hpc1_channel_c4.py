#!/usr/bin/env python3
"""C4: parametric HPC-1 straight-channel LES case generator.

This is the production-readiness sibling of the O/C1 generators.  It writes a
periodic straight FLiBe-like channel, optional high-Pr passive scalar, and a
simple-decomposition `decomposeParDict` that avoids the local OF2512 gotchas:

* do not rely on `decomposePar -numberOfSubdomains N`
* do not assume Scotch exists

Use tiny/smaller cases for local preflight and production-scale counts only on
Grace or another machine with enough wall time.
"""
from __future__ import annotations

import argparse
import math
import stat
from pathlib import Path


def header(cls: str, obj: str) -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*/\n"
        f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
    )


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def choose_decomp(ncores: int) -> tuple[int, int, int]:
    """Pick a simple geometric decomposition biased toward x, then z, then y."""
    best = (ncores, 1, 1)
    best_score = float("inf")
    for nx in range(1, ncores + 1):
        if ncores % nx:
            continue
        rem = ncores // nx
        for ny in range(1, rem + 1):
            if rem % ny:
                continue
            nz = rem // ny
            # Prefer x >= z >= y and avoid extreme imbalance.
            penalty = 0
            if nx < nz:
                penalty += 10
            if nz < ny:
                penalty += 10
            score = penalty + abs(nx / max(nz, 1) - 2.0) + abs(nz / max(ny, 1) - 2.0)
            if score < best_score:
                best = (nx, ny, nz)
                best_score = score
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="hpc1_channel_c4_preflight")
    ap.add_argument("--nx", type=int, default=128)
    ap.add_argument("--ny", type=int, default=96)
    ap.add_argument("--nz", type=int, default=96)
    ap.add_argument("--cores", type=int, default=32)
    ap.add_argument("--decomp", nargs=3, type=int, metavar=("DX", "DY", "DZ"))
    ap.add_argument("--delta", type=float, default=0.01)
    ap.add_argument("--u-bar", type=float, default=1.1598)
    ap.add_argument("--nu", type=float, default=3.09278e-6)
    ap.add_argument("--pr", type=float, default=14.4)
    ap.add_argument("--sgs-prt", type=float, default=0.9, help="LES SGS turbulent Prandtl number for scalarTransport alphaDt")
    ap.add_argument("--grad-y", type=float, default=16.0)
    ap.add_argument("--end-time", type=float, default=4.0)
    ap.add_argument("--delta-t", type=float, default=5.0e-5)
    ap.add_argument("--max-delta-t", type=float, default=2.5e-4)
    ap.add_argument("--write-interval", type=float, default=0.25)
    ap.add_argument("--avg-start", type=float, default=1.5)
    ap.add_argument("--no-thermal", action="store_true")
    args = ap.parse_args()

    case = Path(args.case)
    cells = args.nx * args.ny * args.nz
    decomp = tuple(args.decomp) if args.decomp else choose_decomp(args.cores)
    if math.prod(decomp) != args.cores:
        raise SystemExit(f"decomp product {math.prod(decomp)} != cores {args.cores}")

    d = args.delta
    lx, ly, lz = 2 * math.pi * d, 2 * d, math.pi * d
    alpha = args.nu / args.pr
    alpha_d = 1.0 / args.pr
    alpha_dt = 1.0 / args.sgs_prt
    thermal = not args.no_thermal
    cells_per_core = cells / args.cores

    for sub in ("system", "constant", "0"):
        (case / sub).mkdir(parents=True, exist_ok=True)

    write(
        case / "system/blockMeshDict",
        header("dictionary", "blockMeshDict")
        + f"""scale 1;
vertices
(
    (0 0 0) ({lx:.12g} 0 0) ({lx:.12g} {ly:.12g} 0) (0 {ly:.12g} 0)
    (0 0 {lz:.12g}) ({lx:.12g} 0 {lz:.12g}) ({lx:.12g} {ly:.12g} {lz:.12g}) (0 {ly:.12g} {lz:.12g})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({args.nx} {args.ny} {args.nz})
    simpleGrading (1 ( (0.5 0.5 {args.grad_y:.12g}) (0.5 0.5 {1.0 / args.grad_y:.12g}) ) 1)
);
edges ();
boundary
(
    bottomWall {{ type wall;   faces ( (0 1 5 4) ); }}
    topWall    {{ type wall;   faces ( (3 7 6 2) ); }}
    inlet      {{ type cyclic; neighbourPatch outlet; faces ( (0 4 7 3) ); }}
    outlet     {{ type cyclic; neighbourPatch inlet;  faces ( (1 2 6 5) ); }}
    front      {{ type cyclic; neighbourPatch back;   faces ( (4 5 6 7) ); }}
    back       {{ type cyclic; neighbourPatch front;  faces ( (0 3 2 1) ); }}
);
mergePatchPairs ();
""",
    )

    write(
        case / "constant/transportProperties",
        header("dictionary", "transportProperties")
        + f"""transportModel Newtonian;
nu [0 2 -1 0 0 0 0] {args.nu:.12g};
DT [0 2 -1 0 0 0 0] {alpha:.12g};
""",
    )
    write(
        case / "constant/turbulenceProperties",
        header("dictionary", "turbulenceProperties")
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
    write(
        case / "constant/fvOptions",
        header("dictionary", "fvOptions")
        + f"""momentumSource
{{
    type meanVelocityForce;
    selectionMode all;
    fields (U);
    Ubar ({args.u_bar:.12g} 0 0);
}}
""",
    )

    write(
        case / "0/U",
        header("volVectorField", "U")
        + f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform ({args.u_bar:.12g} 0 0);
boundaryField
{{
    "(bottomWall|topWall)" {{ type noSlip; }}
    "(inlet|outlet)" {{ type cyclic; }}
    "(front|back)" {{ type cyclic; }}
}}
""",
    )
    write(
        case / "0/p",
        header("volScalarField", "p")
        + """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    "(bottomWall|topWall)" { type zeroGradient; }
    "(inlet|outlet)" { type cyclic; }
    "(front|back)" { type cyclic; }
}
""",
    )
    write(
        case / "0/nut",
        header("volScalarField", "nut")
        + """dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{
    "(bottomWall|topWall)" { type nutUSpaldingWallFunction; value uniform 0; }
    "(inlet|outlet)" { type cyclic; }
    "(front|back)" { type cyclic; }
}
""",
    )
    if thermal:
        write(
            case / "0/T",
            header("volScalarField", "T")
            + """dimensions [0 0 0 1 0 0 0];
internalField uniform 0;
boundaryField
{
    bottomWall { type fixedGradient; gradient uniform 1; value uniform 0; }
    topWall    { type fixedGradient; gradient uniform -1; value uniform 0; }
    "(inlet|outlet)" { type cyclic; }
    "(front|back)" { type cyclic; }
}
""",
        )

    functions = f"""fieldAverage1
    {{
        type fieldAverage;
        libs ("libfieldFunctionObjects.so");
        timeStart {args.avg_start:.12g};
        writeControl writeTime;
        restartOnOutput false;
        fields
        (
            U {{ mean on; prime2Mean on; base time; }}
            p {{ mean on; prime2Mean on; base time; }}
{('            T { mean on; prime2Mean on; base time; }' if thermal else '')}
        );
    }}

    yPlus
    {{
        type yPlus;
        libs ("libfieldFunctionObjects.so");
        writeControl writeTime;
    }}"""
    if thermal:
        functions = f"""Ttransport
    {{
        type scalarTransport;
        libs ("libsolverFunctionObjects.so");
        field T;
        schemesField T;
        // Use turbulence-model scalar diffusivity. Setting D directly forces
        // molecular-only transport and drops the LES SGS thermal flux.
        // D_eff = alphaD*nu + alphaDt*nut = nu/Pr + nut/Prt_sgs.
        alphaD {alpha_d:.12g};
        alphaDt {alpha_dt:.12g};
        nCorr 1;
        resetOnStartUp false;
        writeControl writeTime;
    }}

    {functions}"""

    write(
        case / "system/controlDict",
        header("dictionary", "controlDict")
        + f"""application pimpleFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {args.end_time:.12g};
deltaT {args.delta_t:.12g};
writeControl runTime;
writeInterval {args.write_interval:.12g};
purgeWrite 3;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 9;
runTimeModifiable true;
adjustTimeStep yes;
maxCo 0.5;
maxDeltaT {args.max_delta_t:.12g};

functions
{{
    {functions}
}}
""",
    )
    write(
        case / "system/fvSchemes",
        header("dictionary", "fvSchemes")
        + """ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss LUST grad(U);
    div(phi,T) Gauss limitedLinear 1;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""",
    )
    write(
        case / "system/fvSolution",
        header("dictionary", "fvSolution")
        + """solvers
{
    p { solver GAMG; smoother DICGaussSeidel; tolerance 1e-6; relTol 0.05; }
    pFinal { $p; relTol 0; }
    "(U|T|k)(Final)?" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }
}
PIMPLE
{
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    nOuterCorrectors 1;
    pRefCell 0;
    pRefValue 0;
}
""",
    )
    write(
        case / "system/decomposeParDict",
        header("dictionary", "decomposeParDict")
        + f"""numberOfSubdomains {args.cores};
method simple;
coeffs
{{
    n ( {decomp[0]} {decomp[1]} {decomp[2]} );
    delta 0.001;
}}
""",
    )
    t_expr = """

    T
    {
        field T;
        dimensions [0 0 0 1 0 0 0];
        expression #{ 0.02*sin(35.0*pos().x())*sin(30.0*pos().z()) #};
    }""" if thermal else ""
    write(
        case / "system/setExprFieldsDict",
        header("dictionary", "setExprFieldsDict")
        + f"""expressions
(
    U
    {{
        field U;
        dimensions [0 1 -1 0 0 0 0];
        expression #{{ vector
        (
            1.74*(1.0 - sqr(pos().y()/0.01 - 1.0)) + 0.35*sin(40.0*pos().z())*cos(3000.0*pos().y()),
            0.35*sin(35.0*pos().x())*sin(28.0*pos().z()),
            0.35*cos(32.0*pos().x())*sin(2500.0*pos().y())
        ) #}};
    }}{t_expr}
);
""",
    )

    allpre = case / "Allpre.hpc1"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1
setExprFields > log.setExprFields 2>&1
decomposePar -force > log.decomposePar 2>&1
""",
    )
    allpre.chmod(allpre.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    allrun = case / "Allrun.hpc1"
    write(
        allrun,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
: "${MPI_NP:?Set MPI_NP to match system/decomposeParDict numberOfSubdomains}"
MPI_ROOT_ARGS=()
if [ "$(id -u)" -eq 0 ]; then
    export OMPI_ALLOW_RUN_AS_ROOT=1
    export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
    MPI_ROOT_ARGS=(--allow-run-as-root)
fi
mpirun "${MPI_ROOT_ARGS[@]}" -np "$MPI_NP" pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
postProcess -func writeCellCentres -latestTime > log.cc 2>&1 || true
""",
    )
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    write(
        case / "README.hpc1.md",
        f"""# C4 HPC-1 Channel Case

Generated by `gen_hpc1_channel_c4.py`.

- cells: `{cells}`
- decomposition: `{args.cores}` subdomains, simple `n ({decomp[0]} {decomp[1]} {decomp[2]})`
- cells/core: `{cells_per_core:.0f}`
- thermal scalar: `{'on' if thermal else 'off'}`
- Pr: `{args.pr}`
- molecular alpha: `{alpha:.8g} m^2/s`
- scalarTransport coefficients: `alphaD = 1/Pr = {alpha_d:.8g}`, `alphaDt = 1/Prt_sgs = {alpha_dt:.8g}` with `Prt_sgs = {args.sgs_prt:.8g}`

This is a readiness artifact.  Use `./Allpre.hpc1` locally for mesh/decompose
checks.  Use the Grace SLURM templates for real queued runs.
""",
    )

    print(
        f"wrote {case}/ cells={cells} cores={args.cores} "
        f"decomp={decomp} cells/core={cells_per_core:.0f}"
    )
    if not 25_000 <= cells_per_core <= 50_000:
        print("warning: cells/core outside 25k-50k planning band")


if __name__ == "__main__":
    main()

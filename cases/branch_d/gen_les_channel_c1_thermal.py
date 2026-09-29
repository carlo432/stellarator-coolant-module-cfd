#!/usr/bin/env python3
"""C1: periodic-channel LES with a high-Pr passive thermal scalar.

This extends the refined O channel setup with OpenFOAM's scalarTransport
function object.  The default mesh is deliberately a smoke-test mesh; use the
O-like 64x64x64 settings only for a real run.

The default C1 thermal condition is an opposed fixed-gradient scalar on the two
walls.  That gives a statistically stationary periodic scalar problem for
plumbing and profile work.  An opt-in Kawamura transformed-temperature mode is
available for source-term smoke testing, but it is not certification until it
is run long enough and compared against digitized DNS profiles.
"""
from __future__ import annotations

import argparse
import json
import math
import os
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


def kader_beta(pr: float) -> float:
    return (3.85 * pr ** (1.0 / 3.0) - 1.3) ** 2 + 2.12 * math.log(pr)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="les_channel_c1_thermal_smoke")
    parser.add_argument("--nx", type=int, default=24)
    parser.add_argument("--ny", type=int, default=32)
    parser.add_argument("--nz", type=int, default=24)
    parser.add_argument("--delta", type=float, default=0.01, help="channel half-height [m]")
    parser.add_argument("--u-bar", type=float, default=1.1598)
    parser.add_argument("--nu", type=float, default=3.09278e-6)
    parser.add_argument("--pr", type=float, default=14.4)
    parser.add_argument("--sgs-prt", type=float, default=0.9, help="LES SGS turbulent Prandtl number for scalarTransport alphaDt")
    parser.add_argument("--grad-y", type=float, default=12.0)
    parser.add_argument("--t-grad", type=float, default=1.0, help="opposed wall scalar gradient magnitude")
    parser.add_argument(
        "--thermal-wall-bc",
        choices=("opposed_gradient", "kader_exprmixed"),
        default="opposed_gradient",
        help="ordinary-T wall BC; kader_exprmixed is an OPT-1 expression-BC plumbing prototype",
    )
    parser.add_argument("--kader-yplus-ref", type=float, default=30.0)
    parser.add_argument("--kader-theta-scale", type=float, default=0.1)
    parser.add_argument(
        "--thermal-mode",
        choices=("opposed_gradient", "kawamura_transformed"),
        default="opposed_gradient",
        help="thermal scalar setup; default preserves the existing C1 scaffold",
    )
    parser.add_argument(
        "--kawamura-source-re-tau",
        type=float,
        default=180.0,
        help="Re_tau used to scale S_H=(u_tau/delta)*(U_x/U_bulk) in transformed mode",
    )
    parser.add_argument("--end-time", type=float, default=0.01)
    parser.add_argument("--delta-t", type=float, default=1.0e-4)
    parser.add_argument("--max-delta-t", type=float, default=2.5e-4)
    parser.add_argument("--write-interval", type=float, default=0.002)
    parser.add_argument("--avg-start", type=float, default=0.004)
    args = parser.parse_args()

    case = Path(args.case)
    d = args.delta
    lx, ly, lz = 2 * math.pi * d, 2 * d, math.pi * d
    alpha = args.nu / args.pr
    alpha_d = 1.0 / args.pr
    alpha_dt = 1.0 / args.sgs_prt
    cells = args.nx * args.ny * args.nz
    transformed_mode = args.thermal_mode == "kawamura_transformed"
    scalar_field = "H" if transformed_mode else "T"
    scalar_dims = "[0 0 0 0 0 0 0]" if transformed_mode else "[0 0 0 1 0 0 0]"
    kawamura_source_scale = args.kawamura_source_re_tau * args.nu / (d * d)

    for sub in ("system", "constant", "0"):
        (case / sub).mkdir(parents=True, exist_ok=True)

    write(
        case / "system/blockMeshDict",
        foam_header("dictionary", "blockMeshDict")
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
        foam_header("dictionary", "transportProperties")
        + f"""transportModel Newtonian;
nu [0 2 -1 0 0 0 0] {args.nu:.12g};
// Passive scalar diffusivity used by scalarTransport.
DT [0 2 -1 0 0 0 0] {alpha:.12g};
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

    write(
        case / "constant/fvOptions",
        foam_header("dictionary", "fvOptions")
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
        foam_header("volVectorField", "U")
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
        foam_header("volScalarField", "p")
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
        foam_header("volScalarField", "nut")
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

    if transformed_mode:
        scalar_initial = f"""dimensions {scalar_dims};
internalField uniform 0;
boundaryField
{{
    "(bottomWall|topWall)" {{ type fixedValue; value uniform 0; }}
    "(inlet|outlet)" {{ type cyclic; }}
    "(front|back)" {{ type cyclic; }}
}}
"""
    else:
        if args.thermal_wall_bc == "opposed_gradient":
            wall_block = f"""    bottomWall
    {{
        type fixedGradient;
        gradient uniform {args.t_grad:.12g};
        value uniform 0;
    }}
    topWall
    {{
        type fixedGradient;
        gradient uniform {-args.t_grad:.12g};
        value uniform 0;
    }}
"""
            wall_treatment = f"opposed fixedGradient +/- {args.t_grad:.8g}"
        else:
            theta_plus = 2.12 * math.log(args.kader_yplus_ref) + kader_beta(args.pr)
            wall_value = args.kader_theta_scale * theta_plus
            wall_block = f"""    "(bottomWall|topWall)"
    {{
        type exprMixed;
        // OPT-1 Kader expression-BC plumbing prototype.
        // Pr={args.pr:.8g}, yPlusRef={args.kader_yplus_ref:.8g}, thetaPlus={theta_plus:.8g}.
        // valueFraction=1 makes this a fixed-value parser/runnability proof, not final physics.
        valueExpr #{{ {wall_value:.12g} #}};
        gradientExpr #{{ 0 #}};
        fractionExpr #{{ 1 #}};
        value uniform {wall_value:.12g};
    }}
"""
            wall_treatment = (
                f"exprMixed Kader prototype: yPlusRef={args.kader_yplus_ref:.8g}, "
                f"thetaPlus={theta_plus:.8g}, value={wall_value:.8g}"
            )
        scalar_initial = f"""dimensions {scalar_dims};
internalField uniform 0;
boundaryField
{{
{wall_block}
    "(inlet|outlet)" {{ type cyclic; }}
    "(front|back)" {{ type cyclic; }}
}}
"""
    if transformed_mode:
        wall_treatment = "fixedValue H=0 with inline Kawamura source"

    write(
        case / f"0/{scalar_field}",
        foam_header("volScalarField", scalar_field) + scalar_initial,
    )

    if transformed_mode:
        source_block = f"""
        fvOptions
        {{
            kawamuraSource
            {{
                type scalarSemiImplicitSource;
                volumeMode specific;
                selectionMode all;

                sources
                {{
                    H
                    {{
                        explicit
                        {{
                            type exprField;
                            expression #{{ {kawamura_source_scale:.12g}*(U.x()/{args.u_bar:.12g}) #}};
                        }}
                    }}
                }}
            }}
        }}
"""
        statistics_window_start = args.avg_start + args.write_interval
        gate_contract = {
            "version": 7,
            "version_note": (
                "Generated at the final C9 v7 standard: late-window mean "
                "balance from fieldAverage totalTime accumulators, "
                "volume-weighted flatness, and stale-endpoint enforcement. "
                "A short smoke is machinery-only and is expected not to pass."
            ),
            "balance_window_start": statistics_window_start,
            "flatness_mean": "volume",
            "volume_field": "V",
            "balance_model": "closed_periodic_channel_with_two_fixed-value_wall_sinks",
            "coordinate_field": "Cy",
            "walls_m": [0.0, 2 * d],
            "wall_values": [0.0, 0.0],
            "source_to_sink_area_length_m": d,
            "near_wall_fit_levels": 3,
            "coordinate_round_decimals": 10,
            "balance_tolerance_fraction": 0.05,
            "flatness_time_start": statistics_window_start,
            "flatness_tolerance_fraction": 0.01,
            "minimum_flatness_snapshots": 3,
            "scalars": [
                {
                    "field": "H",
                    "balance_field": "HMean",
                    "molecular_diffusivity_m2_s": alpha,
                    "source": {
                        "type": "constant_mean",
                        "mean_value_per_s": kawamura_source_scale,
                        "basis": "meanVelocityForce fixes volume-mean U.x/Ubar to one",
                    },
                }
            ],
        }
        write(
            case / "constant/sourceDrivenScalarGate.json",
            json.dumps(gate_contract, indent=2) + "\n",
        )
        transfer_contract = {
            "version": 1,
            "registered_before_data": True,
            "scope": "straight_periodic_channel_high_pr_transfer",
            "prandtl_number": args.pr,
            "benchmark_parent": "qualified C9 Pr=5 channel-method certification",
            "smoke_interpretation": "machinery_only",
            "smoke_expected_physical_gate_status": "FAIL_OR_NOT_READY",
            "smoke_acceptance": {
                "mesh_and_solver_complete": True,
                "alphaD": alpha_d,
                "alphaDt": alpha_dt,
                "required_fields": [
                    "H", "HMean", "HPrime2Mean", "U", "UMean",
                    "UPrime2Mean", "yPlus", "Cy", "V",
                ],
                "all_required_fields_finite": True,
            },
            "production_requirements": {
                "first_cell_thermal_resolution": (
                    "demonstrate yPlus_1*Pr^(1/3) or an equivalent theta+ "
                    "conductive-sublayer check; do not inherit Pr=5 grading"
                ),
                "relative_sublayer_thinning_vs_pr5": (args.pr / 5.0) ** (1.0 / 3.0),
                "default_initialization": "Kader-profile seeded near equilibrium",
                "equilibration_budget_unseeded_s": [25.0, 30.0],
                "physical_gate": {
                    "balance_tolerance_fraction": 0.05,
                    "flatness_tolerance_fraction": 0.01,
                    "balance_mode": "fieldAverage totalTime accumulator difference",
                    "flatness_mean": "volume",
                    "stale_endpoint_invariant": True,
                    "two_wall_reporting": True,
                    "independent_window_stability": True,
                },
            },
            "claim_limits": {
                "kawamura_p1_p4_certification_allowed": math.isclose(args.pr, 5.0),
                "pr14_status": "transfer/readiness only without a registered matching reference",
                "stellarator_geometry_validation": False,
            },
        }
        write(
            case / "constant/highPrTransferContract.json",
            json.dumps(transfer_contract, indent=2) + "\n",
        )
    else:
        source_block = ""

    write(
        case / "system/controlDict",
        foam_header("dictionary", "controlDict")
        + f"""application pimpleFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {args.end_time:.12g};
deltaT {args.delta_t:.12g};
writeControl runTime;
writeInterval {args.write_interval:.12g};
purgeWrite 0;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
adjustTimeStep yes;
maxCo 0.5;
maxDeltaT {args.max_delta_t:.12g};

functions
{{
    {scalar_field}transport
    {{
        type scalarTransport;
        libs ("libsolverFunctionObjects.so");
        field {scalar_field};
        schemesField {scalar_field};
        // Use turbulence-model scalar diffusivity. Setting D directly forces
        // molecular-only transport and drops the LES SGS thermal flux.
        // D_eff = alphaD*nu + alphaDt*nut = nu/Pr + nut/Prt_sgs.
        alphaD {alpha_d:.12g};
        alphaDt {alpha_dt:.12g};
        nCorr 1;
        resetOnStartUp false;
        writeControl writeTime;
{source_block}
    }}

    fieldAverage1
    {{
        type fieldAverage;
        libs ("libfieldFunctionObjects.so");
        timeStart {args.avg_start:.12g};
        writeControl writeTime;
        restartOnOutput false;
        fields
        (
            U {{ mean on; prime2Mean on; base time; }}
            {scalar_field} {{ mean on; prime2Mean on; base time; }}
            p {{ mean on; prime2Mean on; base time; }}
        );
    }}

    yPlus
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
        + """ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss LUST grad(U);
    div(phi,T) Gauss limitedLinear 1;
    div(phi,H) Gauss limitedLinear 1;
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
    p { solver GAMG; smoother DICGaussSeidel; tolerance 1e-6; relTol 0.05; }
    pFinal { $p; relTol 0; }
    "(U|T|H|k)(Final)?" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }
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
        case / "system/setExprFieldsDict",
        foam_header("dictionary", "setExprFieldsDict")
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
    }}

    {scalar_field}
    {{
        field {scalar_field};
        dimensions {scalar_dims};
        expression #{{ 0.02*sin(35.0*pos().x())*sin(30.0*pos().z()) #}};
    }}
);
""",
    )

    allrun = case / "Allrun.c1"
    write(
        allrun,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1
setExprFields > log.setExprFields 2>&1
pimpleFoam > log.pimpleFoam 2>&1
postProcess -func writeCellCentres -latestTime > log.cc 2>&1 || true
postProcess -func writeCellVolumes -latestTime > log.volumes 2>&1 || true
""",
    )
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    allpre = case / "Allpre.c1"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1
setExprFields > log.setExprFields 2>&1
postProcess -func writeCellCentres -latestTime > log.cc 2>&1 || true
postProcess -func writeCellVolumes -latestTime > log.volumes 2>&1 || true
""",
    )
    allpre.chmod(allpre.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    write(
        case / "README.c1_thermal.md",
        f"""# C1 Thermal LES Channel Case

Generated by `gen_les_channel_c1_thermal.py`.

Purpose: OpenFOAM plumbing case for a high-Pr passive scalar thermal LES in the
periodic channel lineage of O.

Key values:

- cells: `{cells}`
- thermal mode: `{args.thermal_mode}`
- transported scalar field: `{scalar_field}`
- `nu = {args.nu:.8g} m^2/s`
- `Pr = {args.pr:.8g}`
- molecular scalar diffusivity `nu/Pr = {alpha:.8g} m^2/s`
- scalarTransport coefficients: `alphaD = 1/Pr = {alpha_d:.8g}`, `alphaDt = 1/Prt_sgs = {alpha_dt:.8g}` with `Prt_sgs = {args.sgs_prt:.8g}`
- wall scalar treatment: `{wall_treatment}`
- Kawamura source scale: `{kawamura_source_scale:.8g} 1/s`

Important limitation: the physical balance/flatness gate is necessary but not
sufficient for a benchmark label. Only `Pr=5` may use the project's registered
Kawamura P1-P4 comparison. Other Prandtl numbers, including `Pr=14.4`, receive
transfer/readiness status unless a matching reference and acceptance criteria
are separately preregistered. A short smoke is machinery-only and is expected
not to pass the physical gate.

For `Pr=14.4` production, `highPrTransferContract.json` registers the thermal
sublayer-resolution check, Kader-profile seeded default, 25-30 s unseeded
equilibration budget, v7 accumulator/volume gate, two-wall reporting, and
independent-window stability requirement.
""",
    )

    print(
        f"wrote {case}/ cells={cells} mode={args.thermal_mode} "
        f"field={scalar_field} Pr={args.pr:g} alphaD={alpha_d:.6g} alphaDt={alpha_dt:.6g} "
        f"endTime={args.end_time:g}"
    )


if __name__ == "__main__":
    main()

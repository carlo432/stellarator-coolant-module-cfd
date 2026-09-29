#!/usr/bin/env python3
"""Decomposition variant: CURVED duct WITHOUT the sinusoidal span waviness (y_amp = 0).

The straight control (second-geometry pass) removed curvature AND waviness at once, so the
+35%/+16% penalty attribution to Dean secondary flow carries a confound. This case isolates
it: same arc/radius/section/cells as geom_pipeline (fuse-recorded parameters, y_amp -> 0),
same DFSEM inlet (inlet plane at t=0 is identical: sin(0)=0), same annular OpenMC zones
(independent of y), same wall BCs including the unchanged atan2 plasma map.

  straight            = no curvature, no waviness
  nowave (this)       = curvature only
  geom_pipeline family = curvature + waviness

superheat(nowave) - superheat(straight) = curvature; pipeline - nowave = waviness.
"""
from __future__ import annotations
import argparse, json, shutil, stat, subprocess
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--case", default="geom_nowave_pipeline")
ap.add_argument("--ref", default="geom_pipeline")
ap.add_argument("--manifest", default="figs/openmc_depth_source.json")
ap.add_argument("--qwall", type=float, default=0.5e6)
ap.add_argument("--u-bulk", type=float, default=1.99)
ap.add_argument("--end-time", type=float, default=1.3)
ap.add_argument("--avg-start", type=float, default=0.6)
args = ap.parse_args()

NU = 3.09278e-6
ALPHAD = 0.069444444
ALPHADT = 1.1111111
m = json.loads(Path(args.manifest).read_text())
zones = m["zones"]
rhocp = m["rhocp_J_m3K"]
grad_molec = args.qwall / (rhocp * ALPHAD * NU)

case, ref = Path(args.case), Path(args.ref)
scratch = Path("scratch_nowave_gen")
for p in (case, scratch):
    if p.exists():
        shutil.rmtree(p)

# 1) mesh dict from the C2 generator with the fuse-recorded parameters, y_amp = 0
subprocess.run(["python3", "gen_curved_stellarator_slice_c2.py", "--case", str(scratch),
                "--segments", "24", "--stream-cells", "4", "--radial-cells", "24",
                "--height-cells", "24", "--radius", "0.18", "--arc-deg", "95.0",
                "--y-amp", "0.0", "--radial-width", "0.030", "--height", "0.016"],
               check=True, capture_output=True)

# 2) assemble the pipeline case (same skeleton as make_straight_pipeline_case)
(case / "system").mkdir(parents=True)
(case / "0").mkdir()
shutil.copy(scratch / "system/blockMeshDict", case / "system/blockMeshDict")
shutil.copytree(ref / "constant", case / "constant")
if (case / "constant/polyMesh").exists():
    shutil.rmtree(case / "constant/polyMesh")
for d in ("fvSchemes", "fvSolution", "decomposeParDict"):
    shutil.copy(ref / "system" / d, case / "system" / d)
shutil.copy(ref / "system/topoSetDict.qzones", case / "system/topoSetDict.qzones")
shutil.rmtree(scratch)

def H_(cls, obj):
    return ("/*--------------------------------*- C++ -*----------------------------------*/\n"
            f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n")

(case / "0/U").write_text(
    H_("volVectorField","U") + "dimensions [0 1 -1 0 0 0 0];\n"
    f"internalField uniform ({args.u_bulk:g} 0 0);\n"
    "boundaryField\n{\n"
    "    inlet\n    {\n        type turbulentDFSEMInlet;\n        delta 0.0099919;\n"
    "        U mappedFile; UCoeffs {}\n        R mappedFile; RCoeffs {}\n        L mappedFile; LCoeffs {}\n"
    f"        value uniform ({args.u_bulk:g} 0 0);\n    }}\n"
    "    outlet { type pressureInletOutletVelocity; value uniform (0 0 0); }\n"
    '    "(heated_first_wall|blanket_back_wall|side_walls)" { type noSlip; }\n}\n')
(case / "0/p").write_text(
    H_("volScalarField","p") + "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
    "boundaryField\n{\n    inlet { type zeroGradient; }\n    outlet { type fixedValue; value uniform 0; }\n"
    '    "(heated_first_wall|blanket_back_wall|side_walls)" { type zeroGradient; }\n}\n')
(case / "0/nut").write_text(
    H_("volScalarField","nut") + "dimensions [0 2 -1 0 0 0 0];\ninternalField uniform 0;\n"
    "boundaryField\n{\n    inlet { type calculated; value uniform 0; }\n    outlet { type calculated; value uniform 0; }\n"
    '    "(heated_first_wall|blanket_back_wall|side_walls)" { type nutkWallFunction; blending stepwise; value uniform 0; }\n}\n')
gexpr = f"{args.qwall:g} / ({rhocp:g} * ({ALPHAD:g} * {NU:g} + {ALPHADT:g} * nut))"
(case / "0/T").write_text(
    H_("volScalarField","T") + "dimensions [0 0 0 1 0 0 0];\ninternalField uniform 900;\n"
    "boundaryField\n{\n"
    "    inlet { type fixedValue; value uniform 900; }\n"
    "    outlet { type inletOutlet; inletValue uniform 900; value uniform 900; }\n"
    "    heated_first_wall\n    {\n        type exprMixed;\n"
    "        refValue uniform 900;\n"
    f"        refGradient uniform {grad_molec:.6g};\n"
    "        valueFraction uniform 0;\n        source uniform 0;\n        value uniform 900;\n"
    "        valueExpr #{ 900 #};\n"
    f"        gradientExpr #{{ {gexpr} #}};\n"
    "        fractionExpr #{ 0 #};\n        variables ();\n    }\n"
    '    "(blanket_back_wall|side_walls)" { type zeroGradient; }\n}\n')
(case / "0/Ctrit").write_text(
    H_("volScalarField","Ctrit") + "dimensions [0 -3 0 0 1 0 0];\ninternalField uniform 0;\n"
    "boundaryField\n{\n    inlet { type fixedValue; value uniform 0; }\n"
    "    outlet { type inletOutlet; inletValue uniform 0; value uniform 0; }\n"
    '    "(heated_first_wall|blanket_back_wall|side_walls)" { type zeroGradient; }\n}\n')

def src_blocks(field, key):
    return "\n".join(
        f"            qzone{z['zone']}_{field} {{ type scalarSemiImplicitSource; selectionMode cellZone; "
        f"cellZone qzone{z['zone']}; volumeMode specific; sources {{ {field} ({z[key]:.6g} 0); }} }}"
        for z in zones)
def transport(field, key):
    return (f"    transport{field}\n    {{\n        type scalarTransport; libs (solverFunctionObjects);\n"
            f"        field {field}; schemesField U; resetOnStartUp false;\n"
            f"        alphaD {ALPHAD}; alphaDt {ALPHADT};\n        fvOptions\n        {{\n{src_blocks(field, key)}\n        }}\n    }}")
def patch_fo(name, patch, op):
    return (f"    {name}\n    {{\n        type surfaceFieldValue; libs (fieldFunctionObjects);\n"
            f"        log false; writeControl timeStep; writeInterval 20; writeFields false;\n"
            f"        regionType patch; name {patch};\n        {op}\n    }}")
zone_fos = "\n".join(
    f"    qzone{z['zone']}_stats\n    {{\n        type volFieldValue; libs (fieldFunctionObjects);\n"
    f"        log false; writeControl timeStep; writeInterval 100; writeFields false;\n"
    f"        regionType cellZone; name qzone{z['zone']}; operation volAverage; writeVolume true;\n"
    f"        fields (T Ctrit);\n    }}" for z in zones)
(case / "system/controlDict").write_text(
    H_("dictionary","controlDict") +
    f"""application pimpleFoam;
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {args.end_time};
deltaT 0.0002;
writeControl adjustableRunTime;
writeInterval 0.1;
purgeWrite 3;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
adjustTimeStep yes;
maxCo 0.6;
maxDeltaT 0.001;
functions
{{
    fieldAverage1
    {{
        type fieldAverage; libs (fieldFunctionObjects);
        writeControl writeTime;
        timeStart {args.avg_start};
        fields ( U {{ mean on; prime2Mean on; base time; }} T {{ mean on; prime2Mean on; base time; }}
                 Ctrit {{ mean on; prime2Mean off; base time; }} );
    }}
{transport("T", "S_T_K_per_s")}
{transport("Ctrit", "strit_mol_m3_s")}
    yPlus1 {{ type yPlus; libs (fieldFunctionObjects); writeControl writeTime; }}
{patch_fo("inletPhi", "inlet", "operation sum; fields (phi);")}
{patch_fo("outletPhi", "outlet", "operation sum; fields (phi);")}
{patch_fo("inletScal", "inlet", "operation weightedAverage; weightField phi; fields (T Ctrit);")}
{patch_fo("outletScal", "outlet", "operation weightedAverage; weightField phi; fields (T Ctrit);")}
{patch_fo("fwT", "heated_first_wall", "operation areaAverage; writeArea true; fields (T);")}
{zone_fos}
}}
""")

rs = case / "run_pipeline.sh"
rs.write_text("""#!/usr/bin/env bash
set -e
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
cd "$(dirname "$0")"
blockMesh > log.blockMesh 2>&1
topoSet -dict system/topoSetDict.qzones > log.topoSet 2>&1
decomposePar -force > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)

subprocess.run(["bash", "-c",
    f"source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1; cd {case}; "
    "blockMesh > log.blockMesh 2>&1 && topoSet -dict system/topoSetDict.qzones > log.topoSet 2>&1"],
    check=True)
print(f"built {case}: curved (r=0.18, 95deg) WITHOUT waviness; endTime {args.end_time}, avg {args.avg_start}")

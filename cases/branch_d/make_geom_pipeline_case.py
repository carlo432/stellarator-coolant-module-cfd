#!/usr/bin/env python3
"""Task 1.1/1.2: build geom_pipeline from the geom_heat recipe + OpenMC depth sources.

Reads figs/openmc_depth_source.json (from openmc_depth_source.py) and emits a warm-started
pimpleFoam LES case:
  - 6 annular depth cellZones (topoSetDict.qzones, cylinderAnnulusToCell about the arc axis)
  - scalarTransport T with per-zone scalarSemiImplicitSource S_i = q'''_i/(rho*cp)
  - scalarTransport Ctrit with per-zone tritium sources S_T,i [mol/m3/s]
  - surface 0.5 MW/m2 on heated_first_wall via exprMixed: snGrad = q''/(rhoCp*(alphaD*nu+alphaDt*nut))
    -- divides by the LOCAL D_eff incl. wall-function nut. The old fixedGradient 500650 (molecular-only
    conversion) over-delivers by 1+16*(nut_w/nu) wherever y+ > ~11.5; at this mesh's y+ ~ 15 that is
    a multi-x error. This case corrects it; energy balance FOs verify delivery.
  - balance/verification functionObjects (patch phi sums, phi-weighted T/Ctrit, zone volumes+averages)
  - fieldAverage timeStart AFTER warm-start re-equilibration (Task 1.4 hygiene)

Usage: make_geom_pipeline_case.py [--case geom_pipeline] [--smoke] [--src geom_heat]
       --smoke: endTime 0.82, no averaging -- BC/source shakedown only.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import stat
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--case", default="geom_pipeline")
ap.add_argument("--src", default="geom_heat")
ap.add_argument("--manifest", default="figs/openmc_depth_source.json")
ap.add_argument("--smoke", action="store_true")
ap.add_argument("--qwall", type=float, default=0.5e6, help="plasma surface flux [W/m2]")
ap.add_argument("--end-time", type=float, default=2.2)
ap.add_argument("--avg-start", type=float, default=1.5)
args = ap.parse_args()

m = json.loads(Path(args.manifest).read_text())
zones = m["zones"]
g = m["geometry"]
rhocp = m["rhocp_J_m3K"]
NU = 3.09278e-6
ALPHAD = 0.069444444
ALPHADT = 1.1111111
grad_molec = args.qwall / (rhocp * ALPHAD * NU)

src = Path(args.src)
case = Path(args.case)
if case.exists():
    shutil.rmtree(case)
case.mkdir()
shutil.copytree(src / "constant", case / "constant")
shutil.copytree(src / "system", case / "system")
start = "0.8"
shutil.copytree(src / start, case / start)
# reset function-object state (old averaging window) but keep restart time/deltaT
fop = case / start / "uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((case / start).glob("*Mean")) | set((case / start).glob("*Prime2Mean")):
    f.unlink()

# --- topoSetDict: 6 annular depth zones about the arc axis (0,y,0.18) ---
ax, ad = g["axis_point"], g["axis_dir"]
p1 = f"({ax[0]:g} {ax[1] - 1:g} {ax[2]:g})"
p2 = f"({ax[0]:g} {ax[1] + 1:g} {ax[2]:g})"
acts = []
for z in zones:
    i = z["zone"]
    r_in = z["r_inner_m"] - (0.005 if i == 0 else 0.0)              # slack past the walls
    r_out = z["r_outer_m"] + (0.005 if i == len(zones) - 1 else 0.0)
    acts.append(
        f"    {{ name qzone{i}c; type cellSet; action new; source cylinderAnnulusToCell;\n"
        f"      p1 {p1}; p2 {p2}; innerRadius {r_in:g}; outerRadius {r_out:g}; }}\n"
        f"    {{ name qzone{i}; type cellZoneSet; action new; source setToCellZone; set qzone{i}c; }}"
    )
(case / "system/topoSetDict.qzones").write_text(
    "FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }\n"
    "actions\n(\n" + "\n".join(acts) + "\n);\n"
)

# --- 0.8/T: heated_first_wall -> exprMixed flux BC (local-D_eff-correct) ---
tfile = case / start / "T"
txt = tfile.read_text()
old = re.search(r"    heated_first_wall\n    \{[^}]*\}\n", txt).group(0)
expr = (
    "    heated_first_wall\n"
    "    {\n"
    "        type            exprMixed;\n"
    "        refValue        uniform 900;\n"
    f"        refGradient     uniform {grad_molec:.6g};\n"
    "        valueFraction   uniform 0;\n"
    "        source          uniform 0;\n"
    "        value           uniform 900;\n"
    "        valueExpr       #{ 900 #};\n"
    f"        gradientExpr    #{{ {args.qwall:g} / ({rhocp:g} * ({ALPHAD:g} * {NU:g} + {ALPHADT:g} * nut)) #}};\n"
    "        fractionExpr    #{ 0 #};\n"
    "        variables       ();\n"
    "    }\n"
)
tfile.write_text(txt.replace(old, expr))

# --- 0.8/Ctrit: tritium concentration, starts empty ---
(case / start / "Ctrit").write_text(
    "FoamFile { version 2.0; format ascii; class volScalarField; object Ctrit; }\n"
    "dimensions [0 -3 0 0 1 0 0];\n"
    "internalField uniform 0;\n"
    "boundaryField\n{\n"
    "    inlet { type fixedValue; value uniform 0; }\n"
    "    outlet { type inletOutlet; inletValue uniform 0; value uniform 0; }\n"
    '    "(heated_first_wall|blanket_back_wall|side_walls)" { type zeroGradient; }\n'
    "}\n"
)

# --- fvSolution: give Ctrit the same solver as U/T ---
fvs = case / "system/fvSolution"
fvs.write_text(fvs.read_text().replace('"(U|T)"', '"(U|T|Ctrit)"').replace('"(U|T)Final"', '"(U|T|Ctrit)Final"'))

# --- controlDict ---
def src_blocks(field: str, key: str) -> str:
    return "\n".join(
        f"            qzone{z['zone']}_{field} {{ type scalarSemiImplicitSource; selectionMode cellZone; "
        f"cellZone qzone{z['zone']}; volumeMode specific; sources {{ {field} ({z[key]:.6g} 0); }} }}"
        for z in zones
    )

def transport(field: str, key: str) -> str:
    return f"""    transport{field}
    {{
        type scalarTransport; libs (solverFunctionObjects);
        field {field}; schemesField U; resetOnStartUp false;
        alphaD {ALPHAD}; alphaDt {ALPHADT};
        fvOptions
        {{
{src_blocks(field, key)}
        }}
    }}"""

def patch_fo(name: str, patch: str, op_lines: str) -> str:
    return f"""    {name}
    {{
        type surfaceFieldValue; libs (fieldFunctionObjects);
        log false; writeControl timeStep; writeInterval 20; writeFields false;
        regionType patch; name {patch};
        {op_lines}
    }}"""

zone_fos = "\n".join(
    f"""    qzone{z['zone']}_stats
    {{
        type volFieldValue; libs (fieldFunctionObjects);
        log false; writeControl timeStep; writeInterval 100; writeFields false;
        regionType cellZone; name qzone{z['zone']}; operation volAverage; writeVolume true;
        fields (T Ctrit);
    }}""" for z in zones
)

avg = "" if args.smoke else f"""    fieldAverage1
    {{
        type fieldAverage; libs (fieldFunctionObjects);
        writeControl writeTime;
        timeStart {args.avg_start};
        fields ( U {{ mean on; prime2Mean on; base time; }} T {{ mean on; prime2Mean on; base time; }}
                 Ctrit {{ mean on; prime2Mean off; base time; }} );
    }}
"""

end_time = 0.82 if args.smoke else args.end_time
(case / "system/controlDict").write_text(f"""FoamFile
{{
    version 2.0; format ascii; class dictionary; object controlDict;
}}
application pimpleFoam;
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {end_time};
deltaT 0.0002;
writeControl adjustableRunTime;
writeInterval 0.05;
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
{avg}{transport("T", "S_T_K_per_s")}
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

# --- run script ---
rs = case / "run_pipeline.sh"
rs.write_text(f"""#!/usr/bin/env bash
set -e
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
cd "$(dirname "$0")"
topoSet -dict system/topoSetDict.qzones > log.topoSet 2>&1
grep -E "qzone[0-9]  " log.topoSet || true
decomposePar -force > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)

exp_vol = sum(z["power_W"] for z in zones)
print(f"built {case} (smoke={args.smoke})  endTime={end_time}  avgStart={args.avg_start if not args.smoke else '-'}")
print(f"  expected volumetric power {exp_vol:.1f} W (analytic zone volumes)")
print(f"  wall flux {args.qwall/1e6:g} MW/m2 via exprMixed (molecular-only grad would be {grad_molec:.6g} K/m)")
print(f"  expected tritium generation {sum(z['trit_gen_mol_s'] for z in zones):.4g} mol/s")

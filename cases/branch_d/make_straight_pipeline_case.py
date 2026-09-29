#!/usr/bin/env python3
"""Second-geometry pipeline pass: a matched STRAIGHT control through the full pipeline.

Builds a straight duct dimensionally locked to the curved geom_pipeline family (arc length
L = 0.2984513 m = 95 deg x R_mid 0.18; radial gap 30 mm along z; span 16 mm along y;
96 x 24 x 24 = 55,296 cells) and drops the SAME pipeline onto it:
  - DFSEM inlet: the curved case's boundaryData transplants VERBATIM (identical inlet plane
    x=0, y +/-8 mm, z +/-15 mm) -> identical delivered bulk / Re.
  - 6 OpenMC depth zones as z-slabs from the heated wall (z=+0.015) inward (the 1D depth
    column is native to a slab), same specific W/m3 and mol/m3/s as the annulus.
  - heated_first_wall = z=+0.015 plane; exact-flux exprMixed BC (uniform 0.5 MW/m2, or the
    plasma map with s = pos().x()).
The curved heated wall is the INNER (r=0.165) wall, where Dean secondary flow steers fast
fluid AWAY -> curvature disadvantages heated-wall flushing. The straight duct removes exactly
that mechanism, isolating load-shape x mixing from geometry. Cold start (uniform internal +
DFSEM inflow), not warm-start, so establish 4 flow-throughs before averaging.

Usage: make_straight_pipeline_case.py [--case geom_straight_pipeline] [--smoke]
"""
from __future__ import annotations
import argparse, json, shutil, stat, subprocess
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--case", default="geom_straight_pipeline")
ap.add_argument("--ref", default="geom_pipeline", help="curved case to lift constant/boundaryData from")
ap.add_argument("--manifest", default="figs/openmc_depth_source.json")
ap.add_argument("--qwall", type=float, default=0.5e6)
ap.add_argument("--u-bulk", type=float, default=1.99)
ap.add_argument("--end-time", type=float, default=1.3)
ap.add_argument("--avg-start", type=float, default=0.6)
ap.add_argument("--smoke", action="store_true")
args = ap.parse_args()

L = 0.2984513          # arc length of the curved duct = 95deg * 0.18, matched exactly
W = 0.030              # radial gap (z)
H = 0.016              # span (y)
NX, NY, NZ = 96, 24, 24
Z_WALL = 0.5 * W       # heated wall plane z = +0.015
NU = 3.09278e-6
ALPHAD = 0.069444444
ALPHADT = 1.1111111

m = json.loads(Path(args.manifest).read_text())
zones = m["zones"]
rhocp = m["rhocp_J_m3K"]
grad_molec = args.qwall / (rhocp * ALPHAD * NU)

case = Path(args.case)
ref = Path(args.ref)
if case.exists():
    shutil.rmtree(case)
(case / "system").mkdir(parents=True)
(case / "0").mkdir()
# constant/: lift transportProperties, turbulenceProperties (WALE), and the DFSEM boundaryData verbatim
shutil.copytree(ref / "constant", case / "constant")
for junk in ("polyMesh",):
    if (case / "constant" / junk).exists():
        shutil.rmtree(case / "constant" / junk)

def H_(cls, obj):
    return ("/*--------------------------------*- C++ -*----------------------------------*/\n"
            f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n")

# --- blockMeshDict: straight duct, heated wall at z=+0.015 (matches gen_straight_control_c2) ---
v = [(0,-.5*H,-.5*W),(L,-.5*H,-.5*W),(L,.5*H,-.5*W),(0,.5*H,-.5*W),
     (0,-.5*H,.5*W),(L,-.5*H,.5*W),(L,.5*H,.5*W),(0,.5*H,.5*W)]
verts = "\n".join(f"    ({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})" for p in v)
pat = {"inlet":[(0,3,7,4)], "outlet":[(1,5,6,2)],
       "heated_first_wall":[(4,7,6,5)], "blanket_back_wall":[(0,1,2,3)],
       "side_walls":[(0,4,5,1),(3,2,6,7)]}
def faces(fs): return "\n".join("            (%s)" % " ".join(map(str,f)) for f in fs)
btxt = ""
for name, fs in pat.items():
    typ = "patch" if name in ("inlet","outlet") else "wall"
    btxt += f"    {name}\n    {{\n        type {typ};\n        faces\n        (\n{faces(fs)}\n        );\n    }}\n"
(case / "system/blockMeshDict").write_text(
    H_("dictionary","blockMeshDict") +
    f"scale 1;\nvertices\n(\n{verts}\n);\nblocks\n(\n"
    f"    hex (0 1 2 3 4 5 6 7) ({NX} {NY} {NZ}) simpleGrading (1 1 1)\n);\nedges ();\n"
    f"boundary\n(\n{btxt});\nmergePatchPairs ();\n")

# --- 0/ fields (cold start: uniform internal + DFSEM inflow) ---
(case / "0/U").write_text(
    H_("volVectorField","U") +
    "dimensions [0 1 -1 0 0 0 0];\n"
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
# heated wall exact-flux exprMixed; walls else adiabatic; inlet 900, outlet inletOutlet
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

# --- topoSetDict: 6 z-slab depth zones from the heated wall (z=+0.015) inward ---
acts = []
for z in zones:
    i = z["zone"]
    z_hi = Z_WALL - z["depth0_m"] + (0.001 if i == 0 else 0.0)          # slack past the wall
    z_lo = Z_WALL - z["depth1_m"] - (0.001 if i == len(zones) - 1 else 0.0)
    acts.append(
        f"    {{ name qzone{i}c; type cellSet; action new; source boxToCell;\n"
        f"      box (-0.01 -0.01 {z_lo:.6g}) ({L+0.01:.6g} 0.01 {z_hi:.6g}); }}\n"
        f"    {{ name qzone{i}; type cellZoneSet; action new; source setToCellZone; set qzone{i}c; }}")
(case / "system/topoSetDict.qzones").write_text(
    "FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }\n"
    "actions\n(\n" + "\n".join(acts) + "\n);\n")

# --- system dicts lifted from ref (schemes/solution identical), controlDict rebuilt ---
for d in ("fvSchemes", "fvSolution", "decomposeParDict"):
    shutil.copy(ref / "system" / d, case / "system" / d)

def src_blocks(field, key):
    return "\n".join(
        f"            qzone{z['zone']}_{field} {{ type scalarSemiImplicitSource; selectionMode cellZone; "
        f"cellZone qzone{z['zone']}; volumeMode specific; sources {{ {field} ({z[key]:.6g} 0); }} }}"
        for z in zones)
def transport(field, key):
    return (f"    transport{field}\n    {{\n"
            f"        type scalarTransport; libs (solverFunctionObjects);\n"
            f"        field {field}; schemesField U; resetOnStartUp false;\n"
            f"        alphaD {ALPHAD}; alphaDt {ALPHADT};\n        fvOptions\n        {{\n"
            f"{src_blocks(field, key)}\n        }}\n    }}")
def patch_fo(name, patch, op):
    return (f"    {name}\n    {{\n        type surfaceFieldValue; libs (fieldFunctionObjects);\n"
            f"        log false; writeControl timeStep; writeInterval 20; writeFields false;\n"
            f"        regionType patch; name {patch};\n        {op}\n    }}")
zone_fos = "\n".join(
    f"    qzone{z['zone']}_stats\n    {{\n        type volFieldValue; libs (fieldFunctionObjects);\n"
    f"        log false; writeControl timeStep; writeInterval 100; writeFields false;\n"
    f"        regionType cellZone; name qzone{z['zone']}; operation volAverage; writeVolume true;\n"
    f"        fields (T Ctrit);\n    }}" for z in zones)
avg = "" if args.smoke else (
    "    fieldAverage1\n    {\n        type fieldAverage; libs (fieldFunctionObjects);\n"
    "        writeControl writeTime;\n"
    f"        timeStart {args.avg_start};\n"
    "        fields ( U { mean on; prime2Mean on; base time; } T { mean on; prime2Mean on; base time; }\n"
    "                 Ctrit { mean on; prime2Mean off; base time; } );\n    }\n")
end_time = 0.4 if args.smoke else args.end_time
(case / "system/controlDict").write_text(
    H_("dictionary","controlDict") +
    f"""application pimpleFoam;
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {end_time};
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

rs = case / "run_pipeline.sh"
rs.write_text(f"""#!/usr/bin/env bash
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

# build the mesh now so cell/zone counts are verifiable before the long run
subprocess.run(["bash", "-c",
    f"source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1; cd {case}; "
    "blockMesh > log.blockMesh 2>&1 && topoSet -dict system/topoSetDict.qzones > log.topoSet 2>&1"],
    check=True)
zc = subprocess.run(["bash","-c", f"grep -c . /dev/null; grep 'Added' {case}/log.topoSet | grep zone"],
    capture_output=True, text=True)
print(f"built {case} (smoke={args.smoke}) endTime={end_time} avgStart={args.avg_start if not args.smoke else '-'}")
print(f"  L={L} m, {NX}x{NY}x{NZ}={NX*NY*NZ} cells; heated wall z=+{Z_WALL} plane")
print(f"  wall flux {args.qwall/1e6:g} MW/m2 via exprMixed; molecular-only grad would be {grad_molec:.6g} K/m")
print(f"  expected volumetric power {sum(z['power_W'] for z in zones):.1f} W (annular analytic; straight zone vols differ, totals reported at run)")

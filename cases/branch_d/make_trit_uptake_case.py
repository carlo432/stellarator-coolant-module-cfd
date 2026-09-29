#!/usr/bin/env python3
"""Tritium wall-uptake case: geom_pipeline_plasma endpoint + 6 extra Ctrit scalars, each
with a first-order Robin uptake BC on heated_first_wall (see TRITIUM_UPTAKE_PLAN.md).

One run carries the whole k_w sweep because the scalars are passive and independent.
Robin -D_eff dC/dn = k_w C via exprMixed: refValue 0, refGradient 0,
fractionExpr 1/(1 + D_eff*deltaCoeffs/k_w), deltaCoeffs = 1600 (ungraded 1.25 mm cells).
New scalars start from the SOLVED washout Ctrit field (fast equilibration).
"""
from __future__ import annotations
import argparse, json, re, shutil, stat
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--src", default="geom_pipeline_plasma")
ap.add_argument("--case", default="geom_trit_uptake")
ap.add_argument("--manifest", default="figs/openmc_depth_source.json")
ap.add_argument("--end-time", type=float, default=5.0)
ap.add_argument("--avg-start", type=float, default=4.3)
args = ap.parse_args()

NU = 3.09278e-6
ALPHAD = 0.069444444
ALPHADT = 1.1111111
DELTAC = 1600.0                     # 1/m, wall-to-cell-center 0.625 mm, ungraded mesh
KWS = {f"Ctrit_w{n}": 10.0 ** (-n) for n in (5, 4, 3, 2, 1, 0)}   # m/s

m = json.loads(Path(args.manifest).read_text())
zones = m["zones"]

src, case = Path(args.src), Path(args.case)
t0s = "3.6"
assert (src / t0s).is_dir(), f"expected endpoint {src}/{t0s}"
if case.exists():
    shutil.rmtree(case)
case.mkdir()
shutil.copytree(src / "constant", case / "constant")
shutil.copytree(src / "system", case / "system")
shutil.copytree(src / t0s, case / t0s)
fop = case / t0s / "uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((case / t0s).glob("*Mean")) | set((case / t0s).glob("*Prime2Mean")):
    f.unlink()

# --- new scalar fields: washout Ctrit internal field + Robin wall BC ---
ct = (case / t0s / "Ctrit").read_text()
internal = ct[:ct.index("boundaryField")]
for fld, kw in KWS.items():
    frac = f"1 / (1 + ({ALPHAD:g} * {NU:g} + {ALPHADT:g} * nut) * {DELTAC:g} / {kw:g})"
    (case / t0s / fld).write_text(
        re.sub(r"object\s+Ctrit;", f"object {fld};", internal) +
        "boundaryField\n{\n"
        "    inlet { type fixedValue; value uniform 0; }\n"
        "    outlet { type inletOutlet; inletValue uniform 0; value uniform 0; }\n"
        "    heated_first_wall\n    {\n"
        "        type exprMixed;\n"
        "        refValue uniform 0;\n"
        "        refGradient uniform 0;\n"
        "        valueFraction uniform 0;\n"
        "        source uniform 0;\n"
        "        value uniform 0;\n"
        "        valueExpr #{ 0 #};\n"
        "        gradientExpr #{ 0 #};\n"
        f"        fractionExpr #{{ {frac} #}};\n"
        "        variables ();\n    }\n"
        '    "(blanket_back_wall|side_walls)" { type zeroGradient; }\n}\n')

# --- fvSolution: solver selector covers the new fields ---
fvs = case / "system/fvSolution"
fvs.write_text(fvs.read_text().replace('"(U|T|Ctrit)"', '"(U|T|Ctrit.*)"')
                              .replace('"(U|T|Ctrit)Final"', '"(U|T|Ctrit.*)Final"'))

# --- controlDict regenerated (same generator as make_geom_pipeline_case, extended) ---
scalars = ["Ctrit"] + list(KWS)
def src_blocks(field):
    return "\n".join(
        f"            qzone{z['zone']}_{field} {{ type scalarSemiImplicitSource; selectionMode cellZone; "
        f"cellZone qzone{z['zone']}; volumeMode specific; sources {{ {field} ({z['strit_mol_m3_s']:.6g} 0); }} }}"
        for z in zones)
def transport_T():
    blocks = "\n".join(
        f"            qzone{z['zone']}_T {{ type scalarSemiImplicitSource; selectionMode cellZone; "
        f"cellZone qzone{z['zone']}; volumeMode specific; sources {{ T ({z['S_T_K_per_s']:.6g} 0); }} }}"
        for z in zones)
    return ("    transportT\n    {\n        type scalarTransport; libs (solverFunctionObjects);\n"
            "        field T; schemesField U; resetOnStartUp false;\n"
            f"        alphaD {ALPHAD}; alphaDt {ALPHADT};\n        fvOptions\n        {{\n{blocks}\n        }}\n    }}")
def transport_C(field):
    return (f"    transport{field}\n    {{\n        type scalarTransport; libs (solverFunctionObjects);\n"
            f"        field {field}; schemesField U; resetOnStartUp false;\n"
            f"        alphaD {ALPHAD}; alphaDt {ALPHADT};\n        fvOptions\n        {{\n{src_blocks(field)}\n        }}\n    }}")
def patch_fo(name, patch, op):
    return (f"    {name}\n    {{\n        type surfaceFieldValue; libs (fieldFunctionObjects);\n"
            f"        log false; writeControl timeStep; writeInterval 20; writeFields false;\n"
            f"        regionType patch; name {patch};\n        {op}\n    }}")
scal_list = " ".join(scalars)
zone_fos = "\n".join(
    f"    qzone{z['zone']}_stats\n    {{\n        type volFieldValue; libs (fieldFunctionObjects);\n"
    f"        log false; writeControl timeStep; writeInterval 100; writeFields false;\n"
    f"        regionType cellZone; name qzone{z['zone']}; operation volAverage; writeVolume true;\n"
    f"        fields (T {scal_list});\n    }}" for z in zones)
avg_fields = ("fields ( U { mean on; prime2Mean on; base time; } T { mean on; prime2Mean on; base time; }\n"
              + "                 "
              + " ".join(f"{s} {{ mean on; prime2Mean off; base time; }}" for s in scalars) + " );")
transports = "\n".join([transport_T()] + [transport_C(s) for s in scalars])
(case / "system/controlDict").write_text(f"""FoamFile
{{
    version 2.0; format ascii; class dictionary; object controlDict;
}}
application pimpleFoam;
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
        {avg_fields}
    }}
{transports}
    yPlus1 {{ type yPlus; libs (fieldFunctionObjects); writeControl writeTime; }}
{patch_fo("inletPhi", "inlet", "operation sum; fields (phi);")}
{patch_fo("outletPhi", "outlet", "operation sum; fields (phi);")}
{patch_fo("inletScal", "inlet", f"operation weightedAverage; weightField phi; fields (T {scal_list});")}
{patch_fo("outletScal", "outlet", f"operation weightedAverage; weightField phi; fields (T {scal_list});")}
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
decomposePar -force -time {t0s} > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)

gen = sum(z["trit_gen_mol_s"] for z in zones)
deff_mol = ALPHAD * NU * DELTAC
print(f"built {case}: 7 tritium scalars (baseline + kw {sorted(KWS.values())}), 3.6 -> {args.end_time}, avg {args.avg_start}")
print(f"  generation {gen:.4g} mol/s; molecular D*deltaCoeffs = {deff_mol:.3g} m/s (sweep straddles the turbulent value)")

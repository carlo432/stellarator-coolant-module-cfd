#!/usr/bin/env python3
"""THERMAL-RESOLUTION TEST: is C9's thermal failure a GRID problem or an SGS-MODEL problem?

THE QUESTION (this locks the Grace production mesh spec -- C4's missing "thermal-y+ spec")
  C9 (64x160x64, Re_tau=180, Pr=5) passes every MOMENTUM check but over-predicts the thermal
  profile by ~25%, with the error injected in the BUFFER layer (5<y+<30). Two candidate causes:
    (A) GRID: at Pr>1 the scalar cascades BELOW the momentum scales (Batchelor scale = eta/sqrt(Pr)),
        so resolving the scalar needs a wall-parallel grid ~sqrt(Pr) finer than momentum needs.
        C9's grid is dx+ = 17.7, dz+ = 8.8 -- comfortably momentum-resolved, but for the SCALAR at
        Pr=5 it would need dz+ <~ 10/sqrt(5) = 4.5. It is ~2x under-resolved FOR THE SCALAR.
    (B) SGS MODEL: the constant-Prt_sgs scalar flux model simply under-delivers in the buffer.
  These have OPPOSITE consequences for Grace. If (A), the production mesh must be sqrt(14.4)=3.8x
  finer in x-z than a momentum-resolved mesh (~14x the cells, ~35x the cost). If (B), the grid is
  fine and the fix is a better scalar model -- vastly cheaper.

THE TEST (exact, and cheap)
  The scalar is PASSIVE (one-way coupled), so THREE scalars at THREE Prandtl numbers can ride the
  SAME velocity field in ONE run. Identical flow, identical grid, identical turbulence, identical
  source -- the ONLY thing that differs is Pr, i.e. the scalar's scale RELATIVE TO THE GRID.
    H  : alphaD 0.2  -> Pr = 5   (C9's case; grid ~2.0x under-resolved for the scalar)
    H1 : alphaD 1.0  -> Pr = 1   (grid ADEQUATE for the scalar)
    H2 : alphaD 0.5  -> Pr = 2   (grid ~1.2x under-resolved)
  Same alphaDt (Prt_sgs = 0.9) and the same Pr-independent source 5.567004*(U.x()/0.86) for all.

*** PRE-REGISTERED PREDICTION (written BEFORE the run, as the O certification was) ***
  Scalar-adequate needs dz+ <~ 10/sqrt(Pr). C9 has dz+ = 8.8, so the under-resolution ratio is:
      Pr=1 -> needs 10.0, have 8.8  -> ADEQUATE   => thermal error should be SMALL
      Pr=2 -> needs  7.1, have 8.8  -> 1.24x under => error MODERATE
      Pr=5 -> needs  4.5, have 8.8  -> 1.96x under => error LARGE (the ~+25% already measured)
  IF the error grows monotonically with Pr  => cause (A) GRID. Production mesh must be refined in
     x-z by sqrt(Pr); the Grace cost estimate is ~35x too low and MUST be revised before the ask.
  IF the error is FLAT across Pr (large even at Pr=1) => cause (B) SGS MODEL. The grid is fine, the
     Grace cost stands, and the fix is a better scalar flux model.
  Either answer is decisive. A null (no Pr dependence) would REFUTE my own Batchelor argument.

Usage: make_c9_pr_sweep.py   then  bash c9_pr_sweep/run.sh
"""
from __future__ import annotations
import re, shutil, stat
from pathlib import Path

SRC = Path("c9_thermal_cert_run")
CASE = Path("c9_pr_sweep")
END, AVG = 2.0, 1.6          # run 1.2 -> 2.0 (~11 flow-throughs); average the last ~5.5

times = sorted((p.name for p in SRC.iterdir() if re.fullmatch(r"[0-9.]+", p.name)), key=float)
t0 = times[-1]
print(f"cloning {CASE} from {SRC}/{t0}")
if CASE.exists():
    shutil.rmtree(CASE)
CASE.mkdir()
shutil.copytree(SRC / "constant", CASE / "constant")
shutil.copytree(SRC / "system", CASE / "system")
shutil.copytree(SRC / t0, CASE / t0)

# --- reset averaging state; drop stale Mean/Prime2Mean ---
fop = CASE / t0 / "uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((CASE / t0).glob("*Mean")) | set((CASE / t0).glob("*Prime2Mean")):
    f.unlink()

# --- seed H1 (Pr=1) and H2 (Pr=2) from the developed H field ---
Hf = (CASE / t0 / "H").read_text()
for new in ("H1", "H2"):
    (CASE / t0 / new).write_text(re.sub(r"object\s+H\s*;", f"object {new};", Hf))
    print(f"  seeded {new} from the developed H field")

# --- fvSchemes: div schemes for the new scalars ---
fs = CASE / "system/fvSchemes"
t = fs.read_text()
t = t.replace("    div(phi,H) Gauss limitedLinear 1;",
              "    div(phi,H) Gauss limitedLinear 1;\n"
              "    div(phi,H1) Gauss limitedLinear 1;\n"
              "    div(phi,H2) Gauss limitedLinear 1;")
fs.write_text(t)

# --- fvSolution: register H1/H2 with the existing smoothSolver ---
fv = CASE / "system/fvSolution"
t = fv.read_text()
t = t.replace('"(U|T|H|k)(Final)?"', '"(U|T|H|H1|H2|k)(Final)?"')
fv.write_text(t)

# --- controlDict: two more scalarTransport FOs, extend fieldAverage, new window ---
cd = CASE / "system/controlDict"
t = cd.read_text()

def block(field, alphaD, Pr):
    return f"""
    {field}transport
    {{
        type scalarTransport;
        libs ("libsolverFunctionObjects.so");
        field {field};
        schemesField {field};
        // Pr = {Pr}:  D_eff = alphaD*nu + alphaDt*nut = nu/Pr + nut/Prt_sgs
        alphaD {alphaD};
        alphaDt 1.11111111111;
        nCorr 1;
        resetOnStartUp false;
        writeControl writeTime;
        fvOptions
        {{
            kawamuraSource
            {{
                type scalarSemiImplicitSource;
                volumeMode specific;
                selectionMode all;
                sources
                {{
                    {field}
                    {{
                        explicit
                        {{
                            type exprField;
                            expression #{{ 5.567004*(U.x()/0.86) #}};
                        }}
                    }}
                }}
            }}
        }}
    }}
"""

# insert the two new transports right after the existing Htransport block
m = re.search(r"(    Htransport\s*\{.*?\n    \}\n)", t, re.S)
t = t[:m.end(1)] + block("H1", "1.0", 1) + block("H2", "0.5", 2) + t[m.end(1):]
# fieldAverage: add H1/H2, restart the window
t = t.replace("            H { mean on; prime2Mean on; base time; }",
              "            H { mean on; prime2Mean on; base time; }\n"
              "            H1 { mean on; prime2Mean on; base time; }\n"
              "            H2 { mean on; prime2Mean on; base time; }")
t = re.sub(r"timeStart\s+[0-9.]+;", f"timeStart {AVG};", t)
t = re.sub(r"endTime\s+[0-9.]+;", f"endTime {END};", t)
cd.write_text(t)

rs = CASE / "run.sh"
rs.write_text(f"""#!/usr/bin/env bash
set -e
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
cd "$(dirname "$0")"
decomposePar -force -time {t0} > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)
print(f"  built {CASE}: 3 scalars (Pr=5,1,2) on ONE velocity field, {t0} -> {END}, avg from {AVG}")
print("  PRE-REGISTERED: error should GROW with Pr if the cause is GRID; FLAT if the cause is SGS MODEL.")

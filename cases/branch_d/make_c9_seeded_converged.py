#!/usr/bin/env python3
"""FIX THE THERMAL CERTIFICATION BUG: seed the scalar at its EXACT equilibrium amplitude.

THE BUG (found 2026-07-13)
  The Kawamura transformed temperature H is driven by a volumetric source and removed at the walls
  (H=0). At statistical steady state the source balance is EXACT and reference-free:
      integral(S)dV = 2*A*q_w ,  S = 5.567004*(U.x()/0.86) ,  V = 2*delta*A
      => q_w = S_mean * delta = 0.05567   (kinematic; INDEPENDENT of Pr, grid and model)
  MEASURED: C9 ended at q_w = 0.00450 -- EIGHT PERCENT of equilibrium. The scalar had barely begun
  to fill. Its time constant is tau = delta*H/q_w ~ 8 s; C9 ran to t=1.2 s and averaged from 0.45 s.
  And theta+ is NOT invariant during the fill (it drifts +3..+4% per 0.8 s, and the error vs Kader
  SHRINKS as the scalar fills: core +30.8% at 8% fill -> +25.8% at 17% fill).
  => C9's "+25% error / SGS-NECESSARY-NOT-SUFFICIENT", C1's "+28%", the buffer-layer
     under-resolution diagnosis, AND my own Batchelor/Pr test were ALL measured on an
     unconverged field. They are artifacts of under-development, not physics.

THE FIX (kills the slow mode instead of waiting ~25 s / ~250 h for it)
  Seed H(y) = theta_tau * theta+(y+) with theta_tau = q_w_exact/u_tau taken from the EXACT source
  balance (1.0771), and theta+ from Kader. In the conductive sublayer Kader gives theta+ = Pr*y+,
  so dH/dy|wall = theta_tau*Pr*u_tau/nu and hence
      q_w = (nu/Pr)*dH/dy = theta_tau*u_tau = 0.05567   EXACTLY the target, BY CONSTRUCTION.
  The source is therefore balanced from step one: the slow bulk-filling mode is GONE, and only the
  fast turbulent SHAPE relaxation remains (a few flow-throughs).

WHY SEEDING WITH KADER DOES NOT BIAS THE ANSWER
  The convergence criteria are REFERENCE-FREE: (1) q_w -> 0.05567 (exact source balance) and
  (2) <H> flat. The LES then relaxes the profile SHAPE to ITS OWN equilibrium under its own grid and
  SGS model; the seed only fixes the amplitude. Any departure of the converged shape from Kader is
  the LES's own answer. (Seed-independence should still be confirmed with a perturbed seed.)

Usage: make_c9_seeded_converged.py  then  bash c9_thermal_converged/run.sh
"""
from __future__ import annotations
import json
import math, re, shutil, stat
from pathlib import Path
import numpy as np

SRC = Path("c9_pr_sweep"); T0 = "2.000015945742716"
CASE = Path("c9_thermal_converged")
NU, DELTA = 3.09278e-6, 0.01
U_TAU = 0.05169                     # measured from the momentum field
QW_EXACT = 0.05567                  # = S_mean * delta, exact source balance
THETA_TAU = QW_EXACT / U_TAU        # = 1.0771
END, AVG = 2.8, 2.4                 # relax 2.0->2.4 (~5 flow-throughs), average 2.4->2.8

def kader(yp, pr):
    b = (3.85*pr**(1/3) - 1.3)**2 + 2.12*math.log(pr)
    o = np.zeros_like(yp)
    for i, y in enumerate(yp):
        G = 0.01*(pr*y)**4/(1 + 5*pr**3*y)
        o[i] = pr*y*math.exp(-G) + (2.12*math.log(1+y) + b)*math.exp(-1/(G+1e-30))
    return o

def rd(p):
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;', open(p).read(), re.S)
    return np.array([float(x) for x in m.group(1).split()])

print(f"theta_tau (EXACT, from source balance) = {THETA_TAU:.4f}   q_w target = {QW_EXACT}")
if CASE.exists():
    shutil.rmtree(CASE)
CASE.mkdir()
shutil.copytree(SRC/"constant", CASE/"constant")
shutil.copytree(SRC/"system", CASE/"system")
shutil.copytree(SRC/T0, CASE/T0)
fop = CASE/T0/"uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((CASE/T0).glob("*Mean")) | set((CASE/T0).glob("*Prime2Mean")):
    f.unlink()

# --- seed each scalar at its exact equilibrium profile ---
Cy = rd(CASE/T0/"Cy")
yw = np.minimum(Cy, 2*DELTA - Cy)          # wall distance
yp = yw*U_TAU/NU
for fld, pr in (("H", 5.0), ("H1", 1.0), ("H2", 2.0)):
    Heq = THETA_TAU * kader(yp, pr)
    p = CASE/T0/fld
    txt = p.read_text()
    body = "\n".join(f"{v:.8g}" for v in Heq)
    txt = re.sub(r'(internalField\s+nonuniform\s+List<scalar>\s*)\d+\s*\(.*?\)\s*;',
                 lambda m: f"{m.group(1)}{len(Heq)}\n(\n{body}\n)\n;", txt, flags=re.S)
    p.write_text(txt)
    # verify the seeded wall flux reproduces the exact balance
    o = np.argsort(yw); n0 = 6
    dHdy = np.polyfit(np.r_[0, yw[o][:n0]], np.r_[0, Heq[o][:n0]], 1)[0]
    qw = (NU/pr)*abs(dHdy)
    print(f"  seeded {fld:2s} (Pr={pr:.0f}): q_w = {qw:.5f}  ({100*qw/QW_EXACT:.1f}% of exact target)")

cd = CASE/"system/controlDict"
t = cd.read_text()
t = re.sub(r"endTime\s+[0-9.]+;", f"endTime {END};", t)
t = re.sub(r"timeStart\s+[0-9.]+;", f"timeStart {AVG};", t)
t = re.sub(r"startFrom\s+\w+;", "startFrom       latestTime;", t)
cd.write_text(t)

gate = {
    "version": 1,
    "balance_model": "closed_periodic_channel_with_two_fixed-value_wall_sinks",
    "coordinate_field": "Cy",
    "walls_m": [0.0, 2 * DELTA],
    "wall_values": [0.0, 0.0],
    "source_to_sink_area_length_m": DELTA,
    "near_wall_fit_levels": 3,
    "coordinate_round_decimals": 9,
    "balance_tolerance_fraction": 0.05,
    "flatness_time_start": AVG,
    "flatness_tolerance_fraction": 0.01,
    "minimum_flatness_snapshots": 3,
    "scalars": [
        {
            "field": field,
            "molecular_diffusivity_m2_s": NU / pr,
            "source": {"type": "constant_mean", "mean_value_per_s": QW_EXACT / DELTA},
        }
        for field, pr in (("H", 5.0), ("H1", 1.0), ("H2", 2.0))
    ],
}
(CASE / "constant/sourceDrivenScalarGate.json").write_text(
    json.dumps(gate, indent=2) + "\n"
)

rs = CASE/"run.sh"
rs.write_text(f"""#!/usr/bin/env bash
set -e
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
cd "$(dirname "$0")"
decomposePar -force -time {T0} > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)
print(f"  built {CASE}: {T0} -> {END}, relax to {AVG}, average {AVG}-{END}")
print(f"  CONVERGENCE GATE (reference-free): q_w -> {QW_EXACT} AND <H> flat. Only then compare to Kawamura.")

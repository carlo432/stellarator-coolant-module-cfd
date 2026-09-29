#!/usr/bin/env python3
"""PHASE SWEEP: does it matter WHERE along the coolant channel the flux peak lands?

WHY THIS EXISTS (found by auditing C1, 2026-07-12)
  The load map compresses one whole toroidal field period onto a 0.3 m duct. WHERE the peak falls
  along the channel is therefore a MODELLING CHOICE, not a physical fact:
     v1 (used by EVERY CFD case we ever ran) centred the peak at s = 146 mm (mid-duct).
     A2/C1 used the natural zeta=0 -> s=0 mapping, putting the peak at the INLET.
  Nobody checked whether that choice matters. It plausibly does: wall superheat is set by the
  ACCUMULATED UPSTREAM thermal boundary layer, not local q/h -- a peak dumped at the thermal
  entrance meets cold fluid and a thin BL; the same peak mid-duct meets pre-warmed fluid and a
  thick one. Same pk/mean, different physics.

WHAT THIS TESTS
  A real first wall has MANY channels at MANY phases, so the design-relevant number is the
  WORST-PHASE hotspot. This sweep finds it -- and audits whether the project's 2.6x headline
  (measured at ONE arbitrary phase) is robust or lucky.

CONSTRUCTION (exact single-variable test)
  The standoff d(s) is PERIODIC over the duct (d(0)=0.27156, d(L)=0.27178, 0.23 mm apart), so a
  phase shift is a pure roll that preserves the load SHAPE, MEAN (0.5 MW/m2) and PEAK/MEAN (3.1585)
  EXACTLY -- only the peak's position moves. We fit d(s) with a FOURIER series (K=16, the natural
  periodic basis: 0.185 mm error -> 0.25% of peak load; no Runge, and the shift s -> s+delta is
  exact with wrap handled automatically by cos/sin -- v1's degree-10 polynomial could not do this).
  Everything else is byte-identical to the v1 leg: same A, B, lambda, d_min, same denominator.

Usage: make_phase_sweep_case.py --peak-s 0     (peak at the inlet)
       make_phase_sweep_case.py --peak-s 75    (quarter)
       make_phase_sweep_case.py --peak-s 222   (three-quarter)   [146 = the existing baseline]
"""
from __future__ import annotations
import argparse, json, re, shutil, stat
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--peak-s", type=float, required=True, help="target peak location along the duct [mm]")
ap.add_argument("--src", default="geom_straight_pipeline_plasma")
ap.add_argument("--qmap", default="figs/plasma_qmap.json")
ap.add_argument("--K", type=int, default=16)
ap.add_argument("--end-time", type=float, default=4.0)
ap.add_argument("--avg-start", type=float, default=3.3)
args = ap.parse_args()

L = 0.2984513
A_CONST, B_CONST, D_MIN, LAM = 250000.0, 1330253.9, 0.029960662, 0.03
DENOM = "(4.65e+06 * (0.0694444 * 3.09278e-06 + 1.11111 * nut))"

q = json.loads(Path(args.qmap).read_text())
s = np.array(q["s_m"]); d = np.array(q["d_m"])
S_PEAK0 = s[d.argmin()]                      # 0.1461 m -- where v1 put it
target = args.peak_s * 1e-3
delta = S_PEAK0 - target                     # d_shift(s) = d(s + delta) moves the min to `target`

# --- Fourier fit of the smooth periodic standoff ---
M = [np.ones_like(s)]
for k in range(1, args.K + 1):
    M.append(np.cos(2*np.pi*k*s/L)); M.append(np.sin(2*np.pi*k*s/L))
c, *_ = np.linalg.lstsq(np.array(M).T, d, rcond=None)

def d_fit(x):                                 # numpy evaluation, for verification
    v = np.full_like(np.asarray(x, float), c[0])
    for k in range(1, args.K + 1):
        v += c[2*k-1]*np.cos(2*np.pi*k*np.asarray(x)/L) + c[2*k]*np.sin(2*np.pi*k*np.asarray(x)/L)
    return v

# --- emit the OpenFOAM expression for d(s + delta), s = pos().x() (straight duct) ---
X = f"(pos().x() + {delta:.9g})"
terms = [f"{c[0]:.9g}"]
for k in range(1, args.K + 1):
    w = 2*np.pi*k/L
    terms.append(f"{c[2*k-1]:+.9g}*cos({w:.9g}*{X})")
    terms.append(f"{c[2*k]:+.9g}*sin({w:.9g}*{X})")
d_expr = "(" + " ".join(terms) + ")"
gexpr = f"({A_CONST:g} + {B_CONST:g}*exp(-({d_expr} - {D_MIN:.9g})/{LAM:g})) / {DENOM}"

# --- verify the shifted load: shape/mean/peaking must be PRESERVED, only position moves ---
ds = d_fit(s + delta)
qs = A_CONST + B_CONST*np.exp(-(ds - D_MIN)/LAM)
q0 = A_CONST + B_CONST*np.exp(-(d - D_MIN)/LAM)
print(f"phase shift delta = {delta*1e3:+.1f} mm  -> peak target s = {args.peak_s:.0f} mm")
print(f"  emitted load: peak at s = {s[qs.argmax()]*1e3:6.1f} mm | mean {qs.mean()/1e6:.4f} MW/m2 "
      f"| pk/mean {qs.max()/qs.mean():.4f}")
print(f"  v1 reference: peak at s = {s[q0.argmax()]*1e3:6.1f} mm | mean {q0.mean()/1e6:.4f} MW/m2 "
      f"| pk/mean {q0.max()/q0.mean():.4f}")
print(f"  -> shape/power preserved: mean drift {100*(qs.mean()/q0.mean()-1):+.2f}%, "
      f"pk/mean drift {100*(qs.max()/qs.mean())/(q0.max()/q0.mean())-100:+.2f}%")

# --- clone the straight plasma baseline at its endpoint and swap the wall BC ---
src = Path(args.src)
case = Path(f"geom_phase_{int(args.peak_s):03d}")
times = sorted((p.name for p in src.iterdir() if re.fullmatch(r"[0-9.]+", p.name)), key=float)
t0 = times[-1]
if case.exists():
    shutil.rmtree(case)
case.mkdir()
shutil.copytree(src / "constant", case / "constant")
shutil.copytree(src / "system", case / "system")
shutil.copytree(src / t0, case / t0)
fop = case / t0 / "uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((case / t0).glob("*Mean")) | set((case / t0).glob("*Prime2Mean")):
    f.unlink()

tf = case / t0 / "T"
txt = tf.read_text()
txt = re.sub(r"gradientExpr\s+#\{.*?#\};", f"gradientExpr #{{ {gexpr} #}};", txt, flags=re.S)
tf.write_text(txt)

cd = case / "system/controlDict"
t = cd.read_text()
t = re.sub(r"endTime\s+[0-9.]+;", f"endTime {args.end_time};", t)
t = re.sub(r"timeStart\s+[0-9.]+;", f"timeStart {args.avg_start};", t)
cd.write_text(t)

rs = case / "run.sh"
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
print(f"  built {case} (warm-start {t0} -> {args.end_time}, avg from {args.avg_start})")

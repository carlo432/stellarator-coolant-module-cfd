#!/usr/bin/env python3
"""GRID or SGS MODEL? The test that locks the Grace thermal mesh spec.

Three passive scalars ride the SAME velocity field, SAME grid, SAME turbulence, SAME source.
The ONLY variable is Pr -- i.e. the scalar's scale relative to the mesh:
    H1  alphaD 1.0  -> Pr = 1     C9 grid (dz+ = 8.8) is ADEQUATE for the scalar
    H2  alphaD 0.5  -> Pr = 2     ~1.24x under-resolved
    H   alphaD 0.2  -> Pr = 5     ~1.96x under-resolved   (C9's certification case)

PRE-REGISTERED (before the run):
  GRID cause  => thermal error GROWS with Pr (small at 1, large at 5).
  MODEL cause => thermal error is FLAT in Pr (large even at Pr=1). This would REFUTE the
                 Batchelor/sqrt(Pr) argument and mean the Grace cost estimate stands.
Reference: Kader (1981), which was already cross-checked against the digitized Kawamura Pr=5 DNS.
"""
import math, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE = Path("c9_pr_sweep"); T = "2.000015945742716"
NU, DELTA = 3.09278e-6, 0.01
SCALARS = [("H1Mean", 1.0, "#2a5d8f"), ("H2Mean", 2.0, "#c07830"), ("HMean", 5.0, "#7a1f1f")]

def rd_scalar(p):
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;', open(p).read(), re.S)
    return np.array([float(x) for x in m.group(1).split()])

def rd_vecx(p):
    m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;', open(p).read(), re.S)
    return np.array([float(n.split()[0]) for n in re.findall(r'\(([^()]*)\)', m.group(1))])

def kader(yp, pr):
    beta = (3.85*pr**(1/3) - 1.3)**2 + 2.12*math.log(pr)
    o = np.zeros_like(yp)
    for i, y in enumerate(yp):
        G = 0.01*(pr*y)**4/(1 + 5*pr**3*y)
        o[i] = pr*y*math.exp(-G) + (2.12*math.log(1+y) + beta)*math.exp(-1/(G+1e-30))
    return o

Cy = rd_scalar(CASE/T/"Cy"); Ux = rd_vecx(CASE/T/"UMean")
lv = np.unique(np.round(Cy, 9))
Up = np.array([Ux[np.round(Cy, 9) == y].mean() for y in lv])
yw = np.minimum(lv, 2*DELTA - lv)
o = np.argsort(yw); ywo, Uwo = yw[o], Up[o]

# u_tau from near-wall shear (U(0)=0)
n0 = 6
dUdy = np.polyfit(np.r_[0, ywo[:n0]], np.r_[0, Uwo[:n0]], 1)[0]
u_tau = math.sqrt(NU*abs(dUdy)); Re_tau = u_tau*DELTA/NU
yp = ywo*u_tau/NU
dzp = (math.pi*DELTA/64)/(NU/u_tau); dxp = (2*math.pi*DELTA/64)/(NU/u_tau)
print(f"u_tau={u_tau:.5g}  Re_tau={Re_tau:.1f}   GRID: dx+={dxp:.1f}  dz+={dzp:.1f}\n")

def region(err, lo, hi):
    m = (yp >= lo) & (yp < hi)
    return err[m].mean() if m.any() else float("nan")

print(f"{'Pr':>4s} {'dz+ needed':>11s} {'under-res':>10s} | {'sublayer':>9s} {'BUFFER':>8s} {'log':>8s} {'core':>8s}   (LES vs Kader, %)")
res = {}
for fld, pr, col in SCALARS:
    H = rd_scalar(CASE/T/fld)
    Hp = np.array([H[np.round(Cy, 9) == y].mean() for y in lv])[o]
    alpha = NU/pr
    dHdy = np.polyfit(np.r_[0, ywo[:n0]], np.r_[0, Hp[:n0]], 1)[0]
    th_tau = alpha*abs(dHdy)/u_tau
    thp = Hp/th_tau
    kd = kader(yp, pr)
    err = 100*(thp - kd)/np.maximum(kd, 1e-9)
    need = 10/math.sqrt(pr); under = dzp/need
    r = dict(sub=region(err, 0, 5), buf=region(err, 5, 30), log=region(err, 30, 100),
             core=region(err, 100, 1e9), thp=thp, kd=kd, under=under, col=col)
    res[pr] = r
    print(f"{pr:4.0f} {need:11.1f} {under:9.2f}x | {r['sub']:+8.1f} {r['buf']:+7.1f} {r['log']:+7.1f} {r['core']:+7.1f}")

print("\n-- VERDICT --")
b1, b2, b5 = res[1]["buf"], res[2]["buf"], res[5]["buf"]
grew = abs(b5) > abs(b1) + 5
print(f"  BUFFER error:  Pr=1 {b1:+.1f}%   Pr=2 {b2:+.1f}%   Pr=5 {b5:+.1f}%")
if grew:
    print("  => GROWS with Pr. CAUSE = GRID (scalar under-resolution). PRE-REGISTERED PREDICTION HELD.")
    print(f"     Production at Pr=14.4 needs dz+ <= {10/math.sqrt(14.4):.1f} (vs {dzp:.1f} now)"
          f" -> refine x-z by {math.sqrt(14.4):.2f}x -> ~{14.4:.0f}x cells. GRACE COST MUST BE REVISED.")
else:
    print("  => FLAT in Pr. CAUSE = SGS MODEL, not the grid. MY BATCHELOR ARGUMENT IS REFUTED.")
    print("     The grid is adequate; Grace cost estimate STANDS; fix the scalar-flux model instead.")

fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.8), facecolor="white")
a = ax[0]
for fld, pr, col in SCALARS:
    r = res[pr]
    a.semilogx(yp[yp > 0], r["thp"][yp > 0], color=col, lw=2, label=f"LES Pr={pr:.0f}")
    a.semilogx(yp[yp > 0], r["kd"][yp > 0], color=col, lw=1.2, ls="--", alpha=.7)
a.set_xlabel("$y^+$"); a.set_ylabel(r"$\theta^+$")
a.set_title("LES (solid) vs Kader (dashed) — same grid, same flow, only Pr differs", fontsize=10, loc="left")
a.legend(fontsize=8.5); a.grid(alpha=.3)
a = ax[1]
prs = [1, 2, 5]; bufs = [res[p]["buf"] for p in prs]; unders = [res[p]["under"] for p in prs]
a.plot(prs, bufs, "o-", color="#111", ms=9, lw=2)
for p, b, u in zip(prs, bufs, unders):
    a.annotate(f"{u:.2f}x under", (p, b), textcoords="offset points", xytext=(6, -12), fontsize=8, color="0.4")
a.axhline(0, color="0.6", lw=1, ls=":")
a.set_xlabel("Prandtl number"); a.set_ylabel("buffer-layer error vs Kader [%]")
a.set_title("If GRID: error grows with Pr.  If MODEL: flat.", fontsize=10, loc="left")
a.grid(alpha=.3)
fig.tight_layout(); fig.savefig("figs/c9_pr_sweep_verdict.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/c9_pr_sweep_verdict.png")

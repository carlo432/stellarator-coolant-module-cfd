#!/usr/bin/env python3
"""Bent (v23) vs straight (v14) heated-wall temperature profile comparison.

Parses the ASCII VTP heated_wall export for per-face x (centroid) and CellData T,
bins T(x) along the heated wall, and overlays the straight v14 profile (existing CSV).
"""
import os, re, csv, glob
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import numpy as np
import matplotlib.pyplot as plt

TIN = 800.0
ROOT = str(Path(__file__).resolve().parents[2])
V23_VTP = sorted(glob.glob(f"{ROOT}/cases/branch_d/v23_*/VTK/**/boundary/heated_wall.vtp", recursive=True))[0]
V14_CSV = f"{ROOT}/results/wall_profiles/v14_flibe_thermal_refined_100kw_heated_wall_surface_profile.csv"
OUT = f"{ROOT}/cases/branch_d/figs/wall_dT_bent_vs_straight.png"


def _ascii_block(txt, tag_re):
    m = re.search(tag_re, txt, re.DOTALL)
    if not m:
        raise ValueError(f"missing {tag_re}")
    return np.fromstring(m.group(1), sep=" ")


def parse_vtp(path):
    txt = open(path).read()
    pts = _ascii_block(txt, r"Name=['\"]Points['\"][^>]*>(.*?)</DataArray>").reshape(-1, 3)
    conn = _ascii_block(txt, r"Name=['\"]connectivity['\"][^>]*>(.*?)</DataArray>").astype(int)
    offs = _ascii_block(txt, r"Name=['\"]offsets['\"][^>]*>(.*?)</DataArray>").astype(int)
    T = _ascii_block(txt, r"<CellData>.*?Name=['\"]T['\"][^>]*>(.*?)</DataArray>")
    xc = np.empty(len(offs))
    start = 0
    for i, end in enumerate(offs):
        xc[i] = pts[conn[start:end], 0].mean()
        start = end
    return xc, T


def binned(x, T, lo=0.10, hi=0.30, nb=40):
    edges = np.linspace(lo, hi, nb + 1)
    xs, mean, mx = [], [], []
    for i in range(nb):
        m = (x >= edges[i]) & (x < edges[i + 1])
        if m.sum() == 0:
            continue
        xs.append(0.5 * (edges[i] + edges[i + 1]))
        mean.append(T[m].mean() - TIN)
        mx.append(T[m].max() - TIN)
    return np.array(xs), np.array(mean), np.array(mx)


def load_v14(path):
    x, mean, mx = [], [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            x.append(float(r["x"])); mean.append(float(r["dT_mean"])); mx.append(float(r["dT_max"]))
    o = np.argsort(x)
    return np.array(x)[o], np.array(mean)[o], np.array(mx)[o]


xc, T = parse_vtp(V23_VTP)
bx, bmean, bmax = binned(xc, T)
vx, vmean, vmax = load_v14(V14_CSV)

fig, ax = plt.subplots(figsize=(8.0, 4.6), dpi=180)
ax.plot(vx, vmean, "-",  color="tab:blue", lw=1.8, label="straight v14  mean ΔT")
ax.plot(vx, vmax,  "--", color="tab:blue", lw=1.2, label="straight v14  max ΔT")
ax.plot(bx, bmean, "-",  color="tab:red",  lw=1.8, label="bent v23  mean ΔT")
ax.plot(bx, bmax,  "--", color="tab:red",  lw=1.2, label="bent v23  max ΔT")
ax.set_xlabel("x along heated wall [m]  (flow →)")
ax.set_ylabel("wall ΔT = T − 800 K  [K]")
ax.set_title("Heated-wall temperature rise: bent (recirculation) vs straight\n"
             "100 kW/m², FLiBe-like, Re=10,000 — same heated area & heat load")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(OUT)
print("wrote", OUT)
print(f"v23 bent  : wall ΔT mean={bmean.mean():.2f} K, peak(max-bin)={bmax.max():.2f} K")
print(f"v14 strt  : wall ΔT mean={vmean.mean():.2f} K, peak(max-bin)={vmax.max():.2f} K")
print(f"peak shift: bent peak at x={bx[np.argmax(bmax)]:.3f} m vs straight peak at x={vx[np.argmax(vmax)]:.3f} m")

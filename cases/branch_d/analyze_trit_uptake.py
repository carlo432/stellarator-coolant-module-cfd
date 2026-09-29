#!/usr/bin/env python3
"""Tritium wall-uptake verdict (see TRITIUM_UPTAKE_PLAN.md).

Double-entry bookkeeping per scalar over the averaging window:
  FO route   : uptake = (outflow_baseline - outflow_k), PAIRED against the k_w=0 scalar in
               the SAME flow so turbulent outflow scatter cancels common-mode (uptake is a
               sub-0.1% effect; gen - outflow alone would drown in window noise)
  wall route : uptake = sum over heated-wall faces of k_w * CMean_face * A_face
Then the S-curve (uptake fraction vs k_w) and the spatial map (does uptake co-locate with
the thermal hotspot?).
"""
import json
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE = "geom_trit_uptake"
T_END = 5.0
AVG0 = 4.3
KWS = {"Ctrit": 0.0, **{f"Ctrit_w{n}": 10.0 ** (-n) for n in (5, 4, 3, 2, 1, 0)}}
GEN = sum(z["trit_gen_mol_s"] for z in json.load(open("figs/openmc_depth_source.json"))["zones"])
ARC, R_MID = np.radians(95.0), 0.18

def window_mean(path, cols):
    """time-windowed mean of selected 0-based data columns of a FO .dat file."""
    rows = [list(map(float, ln.split())) for ln in open(path) if not ln.startswith("#")]
    a = np.array([r for r in rows if r[0] >= AVG0])
    return a[:, cols].mean(axis=0)

import glob
outQ = window_mean(glob.glob(f"{CASE}/postProcessing/outletPhi/*/surfaceFieldValue.dat")[0], [1])[0]
# outletScal columns: time, T, Ctrit, Ctrit_w5, w4, w3, w2, w1, w0  (order of the FO fields list)
names = ["Ctrit", "Ctrit_w5", "Ctrit_w4", "Ctrit_w3", "Ctrit_w2", "Ctrit_w1", "Ctrit_w0"]
cvals = window_mean(glob.glob(f"{CASE}/postProcessing/outletScal/*/surfaceFieldValue.dat")[0],
                    list(range(2, 2 + len(names))))
outflow = {n: float(cvals[i] * outQ) for i, n in enumerate(names)}

# wall route from the endpoint mean fields
open(f"{CASE}/case.foam", "w").close()
rd = pv.OpenFOAMReader(f"{CASE}/case.foam")
rd.set_active_time_value(T_END)
hw = rd.read()["boundary"]["heated_first_wall"]
hw = hw.compute_cell_sizes(length=False, volume=False)
A = np.asarray(hw.cell_data["Area"])
Tm = np.asarray(hw.cell_data["TMean"])

print(f"generation {GEN:.4g} mol/s | outlet Q {outQ:.6g} m3/s | window {AVG0}-{T_END}")
print(f"absolute closure of the baseline: (gen - outflow_0)/gen = {(GEN - outflow['Ctrit'])/GEN*100:+.2f}% (window-limited)")
print(f"{'field':10s} {'k_w':>8s} {'uptake(FO,paired) %':>20s} {'uptake(wall) %':>15s}")
rows = []
for n in names:
    kw = KWS[n]
    up_fo = (outflow["Ctrit"] - outflow[n]) / GEN * 100
    cf = np.asarray(hw.cell_data[f"{n}Mean"])
    up_wall = float((kw * cf * A).sum()) / GEN * 100
    rows.append(dict(name=n, kw=kw, up_fo=up_fo, up_wall=up_wall))
    print(f"{n:10s} {kw:8.0e} {up_fo:19.4f} {up_wall:14.4f}")

# --- figure ---
fig = plt.figure(figsize=(13.2, 5.6), facecolor="white")
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.35], hspace=0.5, wspace=0.22)

ax = fig.add_subplot(gs[:, 0])
ks = [r["kw"] for r in rows[1:]]
fo = [r["up_fo"] for r in rows[1:]]
wl = [r["up_wall"] for r in rows[1:]]
ax.semilogx(ks, wl, "o-", color="#1a6b52", lw=2, ms=8, label="wall integral  $\\sum k_w C_f A$")
ax.semilogx(ks, fo, "s--", color="#4a6d8c", lw=1.4, ms=7, label="outflow deficit vs $k_w$=0 (paired)")
ceiling = wl[-1]
ax.axhline(ceiling, color="0.6", lw=1, ls=":")
ax.text(1.5e-5, ceiling * 1.04, f"transport-limited ceiling ~{ceiling:.3f}%", fontsize=9, color="0.35")
ax.set_xlabel("uptake rate constant $k_w$ [m/s]")
ax.set_ylabel("first-wall tritium uptake  [% of generation]")
ax.set_title("Uptake S-curve: kinetics-limited $\\to$ transport-limited", fontsize=11)
ax.legend(fontsize=9, loc="upper left")
ax.grid(alpha=0.3, which="both")
for sp in ax.spines.values(): sp.set_color("0.75")

# spatial strips: uptake flux at kw=1e-3 and 1e-1, unrolled
pts = hw.cell_centers().points
th = np.arctan2(pts[:, 0], 0.18 - pts[:, 2])
s = th * R_MID
t = th / ARC
eta = pts[:, 1] - 0.01 * np.sin(2 * np.pi * t)
for i, (n, kw) in enumerate((("Ctrit_w3", 1e-3), ("Ctrit_w1", 1e-1))):
    axs = fig.add_subplot(gs[i, 1])
    J = kw * np.asarray(hw.cell_data[f"{n}Mean"])          # mol/m2/s
    M, xe, ye = np.histogram2d(s, eta, bins=[96, 24], weights=J)
    N, _, _ = np.histogram2d(s, eta, bins=[xe, ye])
    M = (M / np.maximum(N, 1)).T
    pm = axs.pcolormesh(xe * 1e3, ye * 1e3, M, cmap="Greens", shading="flat")
    cb = fig.colorbar(pm, ax=axs, pad=0.01)
    cb.set_label("$J_T$ [mol/m$^2$/s]", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    # thermal hotspot marker (TMean p99 faces), core span only
    hot = Tm >= np.percentile(Tm, 99)
    axs.plot(s[hot] * 1e3, eta[hot] * 1e3, ".", color="#7a1f1f", ms=2, alpha=0.6)
    axs.set_title(f"uptake flux, $k_w$={kw:g} m/s  (red dots = thermal hotspot faces)", fontsize=9.5, loc="left")
    axs.set_ylabel("span [mm]", fontsize=9)
    for sp in axs.spines.values(): sp.set_color("0.75")
axs.set_xlabel("streamwise arc length s [mm]", fontsize=9)

fig.suptitle("Tritium first-wall uptake: how much can the wall take, and where?", fontsize=13, fontweight="bold")
fig.text(0.5, 0.005,
         "linear first-order surrogate for Sievert's-law uptake, k_w swept over decades | heated first wall only | "
         "plasma-shaped load design point | wall-modeled LES tier, indicative",
         ha="center", fontsize=8, color="0.35")
fig.savefig("figs/trit_uptake_verdict.png", dpi=140, bbox_inches="tight")
print("wrote figs/trit_uptake_verdict.png")

# hotspot co-location stat
for n, kw in (("Ctrit_w3", 1e-3), ("Ctrit_w1", 1e-1)):
    J = kw * np.asarray(hw.cell_data[f"{n}Mean"])
    hot = Tm >= np.percentile(Tm, 99)
    print(f"{n}: uptake flux at thermal hotspot / wall mean = {J[hot].mean() / J.mean():.2f}")

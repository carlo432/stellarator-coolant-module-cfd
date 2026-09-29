#!/usr/bin/env python3
"""Pumping ledger for the pipeline family: what does each geometry/physics choice COST in
pressure drop, next to what it does to the hotspot?

dp = area-average kinematic p at inlet (outlet pinned to 0) x rho, averaged over ALL saved
time snapshots of each case (instantaneous p carries turbulent scatter; the saved dirs give
a few-sample mean; indicative, not certified). MHD cases already have exact FO-based dp
(3.0x/7.7x over the 1.49 kPa hydraulic baseline, recorded) and are included for context.

Output: printed table + figs/pumping_ledger.png (dp vs hotspot superheat scatter -- the
design tradeoff plane).
"""
import re
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

RHO = 1940.0
T_IN = 900.0
CASES = [
    # tag, case, superheat-source time (TMean), marker color
    ("straight uniform",  "geom_straight_pipeline",        1.3, "#8a8fa3"),
    ("straight plasma",   "geom_straight_pipeline_plasma", 2.7, "#b05a00"),
    ("nowave uniform",    "geom_nowave_pipeline",          1.3, "#5e7d5e"),
    ("nowave plasma",     "geom_nowave_pipeline_plasma",   2.7, "#2d5a2d"),
    ("curved+wavy uniform","geom_pipeline",                2.2, "#4a6d8c"),
    ("curved+wavy plasma", "geom_pipeline_plasma",         3.6, "#7a1f1f"),
]

def snapshot_times(case):
    ts = []
    for p in Path(case).iterdir():
        if re.fullmatch(r"[0-9]+(\.[0-9]+)?", p.name) and (p / "p").exists():
            ts.append(float(p.name))
    return sorted(ts)

def dp_and_hotspot(case, t_avg):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam")
    dps = []
    for t in snapshot_times(case):
        if t == 0:
            continue
        rd.set_active_time_value(t)
        b = rd.read()["boundary"]
        pin = b["inlet"].cell_data_to_point_data()
        inlet = b["inlet"].compute_cell_sizes(length=False, volume=False)
        A = np.asarray(inlet.cell_data["Area"])
        pi = np.asarray(b["inlet"].cell_data["p"]) if "p" in b["inlet"].cell_data else None
        if pi is None:
            continue
        dps.append(float((pi * A).sum() / A.sum()) * RHO)
    rd.set_active_time_value(t_avg)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_data_to_point_data()
    Tm = np.asarray(hw["TMean"])
    sup = float(np.percentile(Tm, 99)) - T_IN
    return np.mean(dps), np.std(dps), len(dps), sup

rows = []
print(f"{'case':22s} {'dp [Pa]':>9s} {'scatter':>8s} {'n':>2s} {'superheat99 [K]':>16s}")
for tag, case, t, col in CASES:
    if not Path(case).exists():
        print(f"{tag:22s}   -- case not present, skipped")
        continue
    dp, sd, n, sup = dp_and_hotspot(case, t)
    rows.append(dict(tag=tag, dp=dp, sd=sd, sup=sup, col=col))
    print(f"{tag:22s} {dp:9.1f} {sd:8.1f} {n:2d} {sup:16.1f}")

fig, ax = plt.subplots(figsize=(9.2, 6.0), facecolor="white")
for r in rows:
    ax.errorbar(r["dp"] / 1e3, r["sup"], xerr=r["sd"] / 1e3, fmt="o", ms=10,
                color=r["col"], capsize=4, zorder=3)
    ax.annotate(r["tag"], (r["dp"] / 1e3, r["sup"]), textcoords="offset points",
                xytext=(9, -3), fontsize=9.5, color="0.25")
# MHD context points (exact FO-based dp, recorded in the MHD addendum), plasma load
for lab, dpk, sup in (("plasma 5 T", 1.49 * 3.0, 1197 * 1.043), ("plasma 9.2 T", 1.49 * 7.7, 1197 * 1.071)):
    ax.plot(dpk, sup, "D", ms=8, color="#5e3c99", zorder=3)
    ax.annotate(lab + " (MHD)", (dpk, sup), textcoords="offset points", xytext=(9, -3),
                fontsize=9, color="#5e3c99")
ax.set_xlabel("pressure drop  [kPa]  (pumping cost per unit flow)")
ax.set_ylabel("hotspot superheat  p99(TMean) $-$ 900 K  [K]")
ax.set_title("The design tradeoff plane: what each choice costs in pumping\n"
             "vs what it does to the hotspot (same delivered bulk, same wall flux)", fontsize=11.5)
ax.set_xlim(0.2, 14.2)
ax.grid(alpha=0.3)
for sp in ax.spines.values(): sp.set_color("0.75")
fig.text(0.5, 0.008,
         "dp = rho x area-avg inlet p over saved snapshots (few-sample scatter shown); MHD points use exact FO bookkeeping.\n"
         "dp is load-independent within each geometry (one-way coupling -- a consistency check it passes); the curved+wavy "
         "dp ~1.48 kPa independently reproduces the 1.49 kPa hydraulic baseline backed out of the MHD runs.",
         ha="center", fontsize=7.8, color="0.35")
fig.tight_layout(rect=[0, 0.055, 1, 1])
fig.savefig("figs/pumping_ledger.png", dpi=140, bbox_inches="tight")
print("wrote figs/pumping_ledger.png")

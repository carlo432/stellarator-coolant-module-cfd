#!/usr/bin/env python3
"""Hotspot FLICKER statistics — free data, nobody in the Leffler/Ferrero lineage reported this.

TPrime2Mean = resolved variance of T. sigma_T = sqrt(TPrime2Mean) on the heated wall says how
violently the hotspot flickers around its mean — the materials question (transient superheat,
fatigue) that mean-field maps hide. The MHD comparison is a single-window mechanism
diagnostic: if Lorentz drag suppresses turbulence, wall sigma_T should drop with B, but
independent MHD replicate windows are required before promoting a variance ratio.

Caveat carried on the figure: resolved-LES fluctuation only (no SGS part) at a wall-modeled
tier -> sigma_T is a LOWER bound; certify on Grace.

Outputs: figs/hotspot_flicker.png + printed table for PIPELINE_RESULTS.md.
"""
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARC = np.radians(95.0)
R_MID = 0.18

def unroll(points):
    x, y, z = points.T
    th = np.arctan2(x, 0.18 - z)
    s = th * R_MID
    t = th / ARC
    eta = y - 0.01 * np.sin(2 * np.pi * t)
    return s, eta

def wall(case, time):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam")
    rd.set_active_time_value(time)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_centers()
    s, eta = unroll(np.asarray(hw.points))
    Tm = np.asarray(hw["TMean"])
    sig = np.sqrt(np.maximum(np.asarray(hw["TPrime2Mean"]), 0.0))
    return s, eta, Tm, sig

CASES = [
    ("uniform",      "geom_pipeline",        2.2),
    ("plasma",       "geom_pipeline_plasma", 3.6),
    ("plasma B=5T",  "geom_mhd_b5",          5.0),
    ("plasma B=9.2T","geom_mhd_b92",         5.0),
]

rows = []
strips = {}
for tag, case, time in CASES:
    s, e, Tm, sig = wall(case, time)
    hot = Tm >= np.percentile(Tm, 99)              # the hotspot faces
    rows.append(dict(tag=tag,
                     Tm_p99=float(np.percentile(Tm, 99)),
                     sig_wall_mean=float(sig.mean()),
                     sig_at_hotspot=float(sig[hot].mean()),
                     sig_max=float(sig.max()),
                     excursion3=float(np.percentile(Tm, 99) + 3 * sig[hot].mean())))
    strips[tag] = (s, e, sig)
    r = rows[-1]
    print(f"{tag:14s} TMean_p99 {r['Tm_p99']:7.1f} K | sigma_T wall-mean {r['sig_wall_mean']:6.1f} "
          f"at-hotspot {r['sig_at_hotspot']:6.1f} max {r['sig_max']:6.1f} K | "
          f"3-sigma excursion ~{r['excursion3']:7.1f} K")

fig = plt.figure(figsize=(13.5, 8.2), facecolor="white")
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.15], hspace=0.42, wspace=0.25)

hi = max(np.percentile(strips[t][2], 99.5) for t in ("uniform", "plasma"))
for i, tag in enumerate(("uniform", "plasma")):
    ax = fig.add_subplot(gs[i, :])
    s, e, sig = strips[tag]
    M, xe, ye = np.histogram2d(s, e, bins=[96, 24], weights=sig)
    N, _, _ = np.histogram2d(s, e, bins=[xe, ye])
    M = (M / np.maximum(N, 1)).T
    pm = ax.pcolormesh(xe * 1e3, ye * 1e3, M, cmap="Purples", vmin=0, vmax=hi, shading="flat")
    cb = fig.colorbar(pm, ax=ax, pad=0.01)
    cb.set_label(r"$\sigma_T$ [K]")
    r = next(x for x in rows if x["tag"] == tag)
    ax.set_title(f"{tag.upper()} load -- wall temperature flicker sigma_T "
                 f"(wall-mean {r['sig_wall_mean']:.0f} K, at-hotspot {r['sig_at_hotspot']:.0f} K)",
                 fontsize=10.5, loc="left")
    ax.set_ylabel("span [mm]")
    for sp in ax.spines.values(): sp.set_color("0.75")
ax.set_xlabel("streamwise arc length s [mm]")

a = fig.add_subplot(gs[2, 0])
B = [0, 5, 9.2]
sigs = [next(r for r in rows if r["tag"] == t)["sig_at_hotspot"]
        for t in ("plasma", "plasma B=5T", "plasma B=9.2T")]
a.plot(B, sigs, "o-", color="#5e3c99", lw=2, ms=8)
for b, sv in zip(B, sigs):
    a.annotate(f"{sv:.0f} K", (b, sv), textcoords="offset points", xytext=(6, 6), fontsize=9, color="0.3")
a.set_xlabel("B [T]"); a.set_ylabel(r"hotspot $\sigma_T$ [K]")
a.set_title("Single-window MHD flicker diagnostic:\ndrag is associated with lower sigma", fontsize=10.5)
a.grid(alpha=0.3)
for sp in a.spines.values(): sp.set_color("0.75")

a = fig.add_subplot(gs[2, 1])
tags = [r["tag"] for r in rows]
exc = [r["excursion3"] for r in rows]
p99 = [r["Tm_p99"] for r in rows]
yp = np.arange(len(tags))
a.barh(yp - 0.18, np.array(p99) - 900, height=0.36, left=900, color="#c7c7d2", label="mean hotspot (TMean p99)")
a.barh(yp + 0.18, np.array(exc) - 900, height=0.36, left=900, color="#5e3c99", label="~3-sigma excursion")
a.set_yticks(yp, tags); a.invert_yaxis()
a.set_xlabel("wall T [K]")
a.set_title("what the wall actually sees:\nmean hotspot vs 3-sigma excursions", fontsize=10.5)
a.legend(fontsize=8.5, loc="lower right")
a.grid(alpha=0.3, axis="x")
for sp in a.spines.values(): sp.set_color("0.75")

fig.suptitle("Hotspot flicker: the first wall doesn't sit at its mean temperature", fontsize=13, fontweight="bold")
fig.text(0.5, 0.006,
         "sigma_T = sqrt(TPrime2Mean), resolved-LES only (no SGS) at a wall-modeled tier -> LOWER bound on flicker | "
         "windows: 0.7 s each | MHD and phase variants lack replicate variance windows | "
         "excursion estimate assumes ~Gaussian tails (indicative, not certified)",
         ha="center", fontsize=8.5, color="0.35")
fig.savefig("figs/hotspot_flicker.png", dpi=140, bbox_inches="tight")
print("wrote figs/hotspot_flicker.png")

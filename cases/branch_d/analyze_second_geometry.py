#!/usr/bin/env python3
"""Historical second-geometry contrast: load amplification across two duct models.

Four matched-local-flux runs, same pipeline and same delivered bulk/Re:
  curved+wavy uniform (geom_pipeline)       curved+wavy plasma (geom_pipeline_plasma)
  straight uniform (geom_straight_pipeline)  straight plasma (geom_straight_pipeline_plasma)

Reads TMean on heated_first_wall (p99 = hotspot). Reports:
  - within-geometry amplification A = superheat_plasma / superheat_uniform (the headline claim)
  - across-geometry ratio at matched load. This contrast changes curvature and span
    waviness together, so it cannot attribute the ratio to curvature. The separate
    nowave control in analyze_curvature_decomposition.py supplies that attribution.
Outputs figs/second_geometry_verdict.png + printed table.
"""
import numpy as np, pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_IN = 900.0
RUNS = [
    ("curved",   "uniform", "geom_pipeline",                 2.2),
    ("curved",   "plasma",  "geom_pipeline_plasma",          3.6),
    ("straight", "uniform", "geom_straight_pipeline",        1.3),
    ("straight", "plasma",  "geom_straight_pipeline_plasma", 2.7),
]

def wall(case, t):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam"); rd.set_active_time_value(t)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_data_to_point_data()
    Tm = np.asarray(hw["TMean"])
    sig = np.sqrt(np.maximum(np.asarray(hw["TPrime2Mean"]), 0.0))
    return dict(p99=float(np.percentile(Tm, 99)), mean=float(Tm.mean()),
                mx=float(Tm.max()), sig_hot=float(sig[Tm >= np.percentile(Tm, 99)].mean()))

D = {}
for geom, load, case, t in RUNS:
    D[(geom, load)] = wall(case, t)
    r = D[(geom, load)]
    print(f"{geom:9s} {load:8s} p99 {r['p99']:7.1f}  mean {r['mean']:7.1f}  max {r['mx']:7.1f}  "
          f"superheat99 {r['p99']-T_IN:6.1f}  sig_hot {r['sig_hot']:5.1f}")

def sup(g, l): return D[(g, l)]["p99"] - T_IN
A_curved = sup("curved", "plasma") / sup("curved", "uniform")
A_straight = sup("straight", "plasma") / sup("straight", "uniform")
print(f"\namplification (plasma/uniform): curved {A_curved:.3f}  straight {A_straight:.3f}")
print(f"across-geometry at matched load (curved+wavy / straight superheat): "
      f"uniform {sup('curved','uniform')/sup('straight','uniform'):.3f}  "
      f"plasma {sup('curved','plasma')/sup('straight','plasma'):.3f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4), facecolor="white",
                               gridspec_kw=dict(width_ratios=[1.15, 1]))
# left: absolute hotspot p99, grouped by geometry
geoms = ["curved", "straight"]; loads = ["uniform", "plasma"]
x = np.arange(2); w = 0.36
cols = {"uniform": "#8a8fa3", "plasma": "#7a1f1f"}
for i, load in enumerate(loads):
    vals = [D[(g, load)]["p99"] for g in geoms]
    b = ax1.bar(x + (i - 0.5) * w, np.array(vals) - T_IN, w, bottom=T_IN,
                color=cols[load], label=f"{load} load")
    for xi, v in zip(x + (i - 0.5) * w, vals):
        ax1.text(xi, v + 12, f"{v:.0f}", ha="center", fontsize=9, color="0.25")
ax1.axhline(1703, color="#7a1f1f", lw=1.2, ls="-.")
ax1.text(0.5, 1725, "FLiBe boiling 1703 K", fontsize=8.5, color="#7a1f1f", ha="center", va="bottom")
ax1.set_xticks(x, ["CURVED + WAVY", "STRAIGHT"]); ax1.set_ylabel("hotspot film T (p99) [K]")
ax1.set_ylim(900, max(2600, max(D[k]["p99"] for k in D) + 120))
ax1.set_title("Hotspot p99: combined curved+wavy duct vs straight,\nsame pipeline, same delivered bulk, same wall flux", fontsize=10.5)
ax1.legend(fontsize=9, loc="upper left"); ax1.grid(alpha=0.3, axis="y")
for sp in ax1.spines.values(): sp.set_color("0.75")
# right: amplification bars
ax2.bar([0, 1], [A_curved, A_straight], 0.5, color=["#4a6d8c", "#c07830"])
for xi, v in zip([0, 1], [A_curved, A_straight]):
    ax2.text(xi, v + 0.03, f"{v:.2f}x", ha="center", fontsize=12, fontweight="bold", color="0.2")
ax2.axhline(1.0, color="0.5", lw=1, ls=":")
ax2.set_xticks([0, 1], ["curved+wavy", "straight"])
ax2.set_ylabel("hotspot amplification  plasma / uniform  (p99 superheat)")
ax2.set_ylim(0, max(A_curved, A_straight) * 1.25)
ax2.set_title("Is the ~2x amplification geometry-robust?", fontsize=11)
ax2.grid(alpha=0.3, axis="y")
for sp in ax2.spines.values(): sp.set_color("0.75")
fig.suptitle("Second geometry: load amplification survives a combined duct change", fontsize=13, fontweight="bold")
fig.text(
    0.5,
    0.018,
    "matched cross-section/cells/inlet and LOCAL wall flux 0.5 MW/m$^2$; curved heated wall is inner "
    "surface and 9.1% smaller. Curvature and span waviness change together here; this figure cannot "
    "assign their individual effects. See the matched nowave decomposition. Wall-modeled LES tier, indicative.",
    ha="center",
    fontsize=7.4,
    color="0.35",
    wrap=True,
)
fig.tight_layout(rect=[0, 0.065, 1, 0.96])
fig.savefig("figs/second_geometry_verdict.png", dpi=140, bbox_inches="tight")
print("wrote figs/second_geometry_verdict.png")

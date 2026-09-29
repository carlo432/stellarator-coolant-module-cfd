#!/usr/bin/env python3
"""Decompose the 'curvature penalty' into CURVATURE vs SPAN-WAVINESS.

The second-geometry pass compared straight (no curvature, no waviness) against the
geom_pipeline family (curvature + waviness) and attributed the whole penalty to Dean
secondary flow. That carried a confound. geom_nowave_pipeline (curvature, NO waviness)
splits it:

    straight  ->  nowave   : CURVATURE alone
    nowave    ->  curved   : WAVINESS alone   (area-identical, flux-identical: the waviness
                             is a pure y-translation of the section, adds no wall area)

Reports superheat99 (p99 TMean - 900) per case and the two increments, for both loads.
Also reports the amplification within each geometry (must stay ~2.2-2.6x if the load-shape
result is geometry-robust across all three).
Outputs figs/curvature_decomposition.png.
"""
import numpy as np, pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_IN = 900.0
GEOMS = [("straight", "no curvature, no waviness"),
         ("nowave",   "curvature only"),
         ("curved",   "curvature + waviness")]
CASES = {
    ("straight", "uniform"): ("geom_straight_pipeline",        1.3),
    ("straight", "plasma"):  ("geom_straight_pipeline_plasma", 2.7),
    ("nowave",   "uniform"): ("geom_nowave_pipeline",          1.3),
    ("nowave",   "plasma"):  ("geom_nowave_pipeline_plasma",   2.7),
    ("curved",   "uniform"): ("geom_pipeline",                 2.2),
    ("curved",   "plasma"):  ("geom_pipeline_plasma",          3.6),
}

def wall(case, t):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam"); rd.set_active_time_value(t)
    b = rd.read()["boundary"]["heated_first_wall"]
    Tm = np.asarray(b.cell_data_to_point_data()["TMean"])
    a = b.compute_cell_sizes(length=False, volume=False)
    return float(np.percentile(Tm, 99)) - T_IN, float(np.asarray(a.cell_data["Area"]).sum())

S, A = {}, {}
for k, (c, t) in CASES.items():
    S[k], A[k] = wall(c, t)

print(f"{'geometry':10s} {'area [m2]':>11s} {'uniform sup':>12s} {'plasma sup':>11s} {'amplification':>14s}")
for g, _ in GEOMS:
    print(f"{g:10s} {A[(g,'uniform')]:11.5e} {S[(g,'uniform')]:12.1f} {S[(g,'plasma')]:11.1f} "
          f"{S[(g,'plasma')]/S[(g,'uniform')]:13.3f}x")

print("\nDECOMPOSITION (increment over the straight control):")
for load in ("uniform", "plasma"):
    s, n, c = S[("straight", load)], S[("nowave", load)], S[("curved", load)]
    tot = (c - s) / s * 100
    cur = (n - s) / s * 100
    wav = (c - n) / s * 100
    print(f"  {load:8s} total {tot:+6.1f}%  =  curvature {cur:+6.1f}%  +  waviness {wav:+6.1f}%")

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.4), facecolor="white",
                              gridspec_kw=dict(width_ratios=[1.25, 1]))
x = np.arange(3); w = 0.36
cols = {"uniform": "#8a8fa3", "plasma": "#7a1f1f"}
for i, load in enumerate(("uniform", "plasma")):
    vals = [S[(g, load)] for g, _ in GEOMS]
    ax.bar(x + (i - 0.5) * w, vals, w, color=cols[load], label=f"{load} load")
    for xi, v in zip(x + (i - 0.5) * w, vals):
        ax.text(xi, v + 12, f"{v:.0f}", ha="center", fontsize=9, color="0.25")
ax.set_xticks(x, [f"{g}\n({d})" for g, d in GEOMS], fontsize=8.5)
ax.set_ylabel("hotspot superheat  p99$(T_{mean})-900$ [K]")
ax.set_title("Decomposing the geometry penalty", fontsize=11)
ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3, axis="y")
for sp in ax.spines.values(): sp.set_color("0.75")

# waterfall of increments (curvature always +, waviness may be negative -> draw signed)
SCAT = 3.0   # window-scatter band on a superheat ratio, % (from the MHD/plasma windows)
X_TOT = 44.0
for j, load in enumerate(("uniform", "plasma")):
    s, n, c = S[("straight", load)], S[("nowave", load)], S[("curved", load)]
    cur, wav = (n - s) / s * 100, (c - n) / s * 100
    y = j * 1.0
    ax2.barh(y, cur, 0.42, color="#4a6d8c", label="curvature (Dean)" if j == 0 else None)
    ax2.text(cur / 2, y, f"{cur:+.1f}%", ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")
    lo, hi = (cur + wav, cur) if wav < 0 else (cur, cur + wav)
    ax2.barh(y, hi - lo, 0.42, left=lo, color="#c07830", label="span waviness" if j == 0 else None)
    ax2.annotate(f"{wav:+.1f}%" + ("  (within scatter)" if abs(wav) <= SCAT else ""),
                 ((lo + hi) / 2, y - 0.30), ha="center", va="center", fontsize=9, color="#8a5000")
    ax2.text(X_TOT, y, f"total {cur+wav:+.1f}%", va="center", ha="right", fontsize=9.5,
             color="0.25", fontweight="bold")
ax2.axvspan(-SCAT, SCAT, color="0.85", alpha=0.5, zorder=0)
ax2.text(0, 1.62, f"$\\pm${SCAT:.0f}% window scatter", ha="center", fontsize=8, color="0.45")
ax2.set_yticks([0, 1], ["uniform", "plasma"])
ax2.set_ylim(-0.6, 1.75)
ax2.set_xlim(-6, X_TOT + 2)
ax2.axvline(0, color="0.5", lw=1)
ax2.set_xlabel("hotspot superheat increment over the straight control [%]")
ax2.set_title("Curvature vs waviness, separated", fontsize=11)
ax2.legend(fontsize=8.5, loc="upper right", framealpha=0.95); ax2.grid(alpha=0.3, axis="x")
for sp in ax2.spines.values(): sp.set_color("0.75")

fig.suptitle("Which part of the duct shape heats the wall?", fontsize=13, fontweight="bold")
fig.text(0.5, 0.012,
         "matched cells/inlet/local wall flux | nowave vs curved is area-identical AND flux-identical (waviness is a pure\n"
         "y-translation) | straight heated wall is 9.1% larger (mid- vs inner-radius) so it absorbs 7.4% more total power "
         "yet runs coolest | wall-modeled LES tier",
         ha="center", fontsize=7.8, color="0.35")
fig.tight_layout(rect=[0, 0.075, 1, 0.95])
fig.savefig("figs/curvature_decomposition.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/curvature_decomposition.png")

#!/usr/bin/env python3
"""Sensitivity band RE-MEASURED on the clean (no-waviness) geometry.

WHY THIS EXISTS. The promoted band ("uniform-flux underestimates the hotspot 1.6-4.1x") was
computed on the curved+wavy duct. The curvature/waviness decomposition showed the span
waviness inflates the UNIFORM reference by +16.4% while leaving plasma-load hotspots
essentially untouched (-2.6%, within scatter) -- because a uniform load lets the hotspot
settle at the waviness-made weak-flushing spot, while a plasma load pins it at the flux peak.
An amplification is plasma/uniform, so a waviness-inflated denominator DEPRESSES the whole
band. The band is therefore systematically LOW, and the fix is to re-reference it on the
y_amp = 0 geometry -- measured, not scaled.

Four nowave runs (uniform + three load-shape corners at identical local-flux means), same
warm-start protocol as the wavy family. Emits the corrected band and the wavy-vs-clean
comparison figure.
"""
import json
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_IN = 900.0
# The uniform reference is taken from the DRIFT-CONTROL run (fresh 0.7 s window on an already
# developed field), not the cold-start window: the mean differs by only -0.35%, but using the
# developed value keeps every quantity in this table on the same footing as the warm-started
# plasma corners. (The variance does NOT converge as fast -- see the drift-control addendum.)
NOWAVE = [
    ("uniform", "geom_nowave_drift_ctl",        2.0, 1.00),
    ("mild",    "geom_nowave_plasma_mild",      2.7, 1.916),
    ("default", "geom_nowave_pipeline_plasma",  2.7, 3.142),  # endpoint-audited v1 peak/mean
    ("severe",  "geom_nowave_plasma_severe",    2.7, 6.686),
]

def sup99(case, t):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam"); rd.set_active_time_value(t)
    Tm = np.asarray(rd.read()["boundary"]["heated_first_wall"]
                    .cell_data_to_point_data()["TMean"])
    return float(np.percentile(Tm, 99)) - T_IN

rows = [dict(tag=t, pk=pk, sup=sup99(c, tt)) for t, c, tt, pk in NOWAVE]
u = rows[0]["sup"]
for r in rows:
    r["amp"] = r["sup"] / u

wavy = {r["tag"]: r for r in json.load(open("figs/plasma_sensitivity_results.json"))}

print(f"clean (nowave) uniform reference: {u:.1f} K   [wavy was {wavy['uniform']['superheat99']:.1f} K, "
      f"+{100*(wavy['uniform']['superheat99']/u-1):.1f}% waviness inflation]\n")
print(f"{'corner':8s} {'pk/mean':>8s} {'sup [K]':>9s} {'amp CLEAN':>10s} {'amp wavy':>9s} {'shift':>7s}")
for r in rows:
    w = wavy[r["tag"]]
    print(f"{r['tag']:8s} {r['pk']:8.2f} {r['sup']:9.1f} {r['amp']:9.3f}x {w['amp99']:8.3f}x "
          f"{100*(r['amp']/w['amp99']-1):+6.1f}%")

lo, hi = rows[1]["amp"], rows[3]["amp"]
wlo, whi = wavy["mild"]["amp99"], wavy["severe"]["amp99"]
print(f"\nBAND, clean geometry : {lo:.1f}x - {hi:.1f}x   (centre {rows[2]['amp']:.1f}x)")
print(f"BAND, wavy (promoted): {wlo:.1f}x - {whi:.1f}x   (centre {wavy['default']['amp99']:.1f}x)")
json.dump(rows, open("figs/plasma_sensitivity_nowave.json", "w"), indent=1)

fig, ax = plt.subplots(figsize=(9.4, 6.0), facecolor="white")
pk_c = [r["pk"] for r in rows]; amp_c = [r["amp"] for r in rows]
pk_w = [wavy[r["tag"]]["peak_over_mean"] for r in rows]; amp_w = [wavy[r["tag"]]["amp99"] for r in rows]
ax.plot(pk_c, amp_c, "o-", color="#1a6b52", lw=2.4, ms=9, zorder=4,
        label=f"clean geometry (y_amp=0):  {lo:.1f}x - {hi:.1f}x")
ax.plot(pk_w, amp_w, "s--", color="#7a1f1f", lw=1.6, ms=7, zorder=3, alpha=0.85,
        label=f"as-built wavy duct (promoted):  {wlo:.1f}x - {whi:.1f}x")
for r in rows:
    ax.annotate(f"{r['tag']}  {r['amp']:.2f}x", (r["pk"], r["amp"]), textcoords="offset points",
                xytext=(10, 8), fontsize=9, color="#14503e", fontweight="bold")
for r in rows[1:]:
    w = wavy[r["tag"]]
    ax.annotate(f"{w['amp99']:.2f}x", (w["peak_over_mean"], w["amp99"]), textcoords="offset points",
                xytext=(4, -15), fontsize=8.5, color="#7a1f1f")
ax.plot([1, 7], [1, 7], ":", color="0.6", lw=1.2, label="1:1 (amplification = load peaking)")
ax.set_xlim(0.7, 7.4); ax.set_ylim(0.6, 7.4)
ax.set_xlabel("delivered load peak/mean  [-]")
ax.set_ylabel("hotspot superheat amplification vs uniform  (p99)")
ax.set_title("The sensitivity band, re-measured without the waviness artifact\n"
             "(waviness inflates the UNIFORM reference, depressing every amplification)",
             fontsize=11.5)
ax.legend(fontsize=9, loc="upper left")
ax.grid(alpha=0.3)
for sp in ax.spines.values(): sp.set_color("0.75")
fig.text(0.5, 0.015,
         "same local-flux mean 0.5 MW/m$^2$ in every run; corners are the lambda / f_rad band edges. Amplification still sits\n"
         "below 1:1 (turbulent streamwise mixing shaves the peak). The correction grows with peaking: a mild load spreads like\n"
         "a uniform one and takes the same waviness penalty, while a peaked load pins its\n"
         "hotspot at the flux peak and escapes it -- so the waviness biased the reference, not the plasma cases.",
         ha="center", fontsize=7.8, color="0.35")
fig.tight_layout(rect=[0, 0.115, 1, 1])
fig.savefig("figs/plasma_sensitivity_nowave.png", dpi=140)
print("wrote figs/plasma_sensitivity_nowave.png, figs/plasma_sensitivity_nowave.json")

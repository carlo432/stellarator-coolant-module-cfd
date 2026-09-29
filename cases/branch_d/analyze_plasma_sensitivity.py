#!/usr/bin/env python3
"""Sensitivity of the Phase-2 headline (hotspot superheat vs uniform) to the lambda/f_rad bands.

Three CFD points: severe (lam 1cm, f_rad 0.3), default (3cm, 0.5), mild (5cm, 0.7), all
same total power as the uniform case. Output: hotspot amplification vs load peak/mean —
the design curve that turns the single 2.1x claim into a banded statement.
"""
import json
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASES = [
    ("uniform", "geom_pipeline", 2.2, 1.00),
    ("mild",    "geom_pipeline_plasma_mild", 5.0, None),
    ("default", "geom_pipeline_plasma", 3.6, None),
    ("severe",  "geom_pipeline_plasma_severe", 5.0, None),
]

def wall_tmean(case, time):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam")
    rd.set_active_time_value(time)
    return np.asarray(rd.read()["boundary"]["heated_first_wall"].cell_data_to_point_data()["TMean"])

def delivered_peak_over_mean(tag):
    """Peak/mean of the load the CFD actually saw (the fitted expression), not the ideal."""
    import plasma_wall_qmap as Q
    j = json.load(open(f"figs/plasma_qmap_{tag}.json")) if tag != "default" \
        else json.load(open("figs/plasma_qmap.json"))
    lam = j.get("lam", j.get("lambda_m")); frad = j.get("f_rad", 0.5)
    m = Q.build_map(lam=lam, f_rad=frad)
    grad, c, qfit, fiterr = Q.emit_expr(m["s"], m["d"], lam, frad)
    return float(qfit.max() / Q.QWALL)

T_IN = 900.0
rows = []
for tag, case, time, pk in CASES:
    Tm = wall_tmean(case, time)
    if pk is None:
        pk = delivered_peak_over_mean(tag)
    sup99 = np.percentile(Tm, 99) - T_IN
    rows.append(dict(tag=tag, peak_over_mean=pk, mean=float(Tm.mean()),
                     p99=float(np.percentile(Tm, 99)), max=float(Tm.max()),
                     superheat99=float(sup99)))
    print(f"{tag:8s} load pk/mean {pk:4.2f} | TMean mean {Tm.mean():6.1f}  p99 {np.percentile(Tm,99):6.1f}  "
          f"max {Tm.max():6.1f} | superheat99 {sup99:6.1f} K")

u = rows[0]
for r in rows:
    r["amp99"] = r["superheat99"] / u["superheat99"]

grid = json.load(open("figs/plasma_qmap_sensitivity_grid.json"))
gvals = np.asarray(grid["peak_over_mean"]).ravel()

fig, ax = plt.subplots(figsize=(8.6, 5.6), facecolor="white")
pks = [r["peak_over_mean"] for r in rows]
amps = [r["amp99"] for r in rows]
ax.plot(pks, amps, "o-", color="#7a1f1f", lw=2, ms=9, zorder=3)
for r in rows:
    ax.annotate(f"{r['tag']}\n({r['amp99']:.2f}x)", (r["peak_over_mean"], r["amp99"]),
                textcoords="offset points", xytext=(8, -4), fontsize=9.5, color="0.25")
ax.axvspan(gvals.min(), gvals.max(), color="#b3cde0", alpha=0.35, zorder=1,
           label=f"full (lambda, f_rad) band rectangle: load pk/mean {gvals.min():.1f}-{gvals.max():.1f}")
ax.plot([1, max(pks) * 1.05], [1, max(pks) * 1.05], ":", color="0.6", lw=1.2,
        label="1:1 (superheat amplification = load peaking)")
ax.set_xlabel("delivered load peak/mean  [-]")
ax.set_ylabel("wall hotspot superheat amplification vs uniform  (p99, T - 900 K)")
ax.set_title("How sensitive is the hotspot to the load-shape bands?\n"
             "3 LES runs at the band corners + center, same total power (0.5 MW/m2 mean)",
             fontsize=11.5)
ax.legend(fontsize=9, loc="upper left")
ax.grid(alpha=0.3)
for sp in ax.spines.values(): sp.set_color("0.75")
fig.text(0.5, 0.012,
         "W7-X standoff transplant; lambda 1-5 cm, f_rad 0.3-0.7 | amplification sits below 1:1 "
         "(turbulent streamwise mixing shaves the peak) | uniform-flux design assumption is "
         "non-conservative across the ENTIRE band",
         ha="center", fontsize=8.5, color="0.35")
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig("figs/plasma_sensitivity_curve.png", dpi=140)
json.dump(rows, open("figs/plasma_sensitivity_results.json", "w"), indent=1)
print("wrote figs/plasma_sensitivity_curve.png, figs/plasma_sensitivity_results.json")

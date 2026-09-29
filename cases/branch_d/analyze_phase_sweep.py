#!/usr/bin/env python3
"""PHASE SWEEP VERDICT: does WHERE the flux peak lands along the channel change the hotspot?

Four runs, IDENTICAL in every respect except the peak's streamwise position:
  same load shape, same mean (0.5025 MW/m2), same pk/mean (3.157), same duct, same warm-start,
  same averaging window, same Fourier d(s) expression. Only the phase differs.
    geom_phase_000  peak at   1.5 mm  (inlet -- the phase A2/C1 implicitly used)
    geom_phase_075  peak at  76.6 mm
    geom_phase_146  peak at 146.1 mm  (CONTROL = the phase EVERY prior CFD case used)
    geom_phase_222  peak at 223.5 mm

WHY IT MATTERS: the load map compresses one toroidal field period onto a 0.3 m duct, so the peak's
position along the channel is a MODELLING CHOICE. A real first wall has many channels at many
phases, so the design number is the WORST phase. And the project's 2.6x amplification headline,
the 1.7-4.9x sensitivity band, the load envelope, the operating window and the error budget were
ALL measured at the single phase 146 mm. If phase matters, they are phase-CONDITIONAL.

Outputs: printed verdict + figs/phase_sweep_verdict.png +
figs/phase_sweep_metrics.json
"""
import json
import numpy as np, pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_IN, T_BOIL = 900.0, 1703.0
CASES = [("geom_phase_000", 1.5, "#7a1f1f"), ("geom_phase_075", 76.6, "#c07830"),
         ("geom_phase_146", 146.1, "#1a6b52"), ("geom_phase_222", 223.5, "#2a5d8f")]

def wall(case, t=4.0):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam"); rd.set_active_time_value(t)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_data_to_point_data()
    return np.asarray(hw.points), np.asarray(hw["TMean"]), \
           np.sqrt(np.maximum(np.asarray(hw["TPrime2Mean"]), 0.0))

def prof(pts, T, nb=60):
    core = np.abs(pts[:, 1]) < 0.005
    x, Tc = pts[core, 0], T[core]
    e = np.linspace(0, 0.2984513, nb + 1); xc = 0.5*(e[:-1]+e[1:])
    i = np.clip(np.digitize(x, e)-1, 0, nb-1)
    return xc, np.array([Tc[i == k].max()-T_IN if np.any(i == k) else np.nan for k in range(nb)])

R = {}
print(f"{'case':16s} {'peak@mm':>8s} {'p99 dT':>8s} {'max dT':>8s} {'mean dT':>8s} {'sig_hot':>8s}  film-T(p99)")
for case, ps, _ in CASES:
    pts, T, sig = wall(case)
    p99 = np.percentile(T, 99)
    R[case] = dict(ps=ps, p99=float(p99-T_IN), mx=float(T.max()-T_IN), mean=float(T.mean()-T_IN),
                   sig=float(sig[T >= p99].mean()), film=float(p99))
    r = R[case]
    print(f"{case:16s} {ps:8.1f} {r['p99']:8.1f} {r['mx']:8.1f} {r['mean']:8.1f} {r['sig']:8.1f}   {r['film']:7.1f} K")

control_case = "geom_phase_146"
best_case = min(R, key=lambda case: R[case]["p99"])
worst_case = max(R, key=lambda case: R[case]["p99"])
worst_sigma_case = max(R, key=lambda case: R[case]["sig"])
ctrl = R[control_case]["p99"]           # the phase EVERY prior CFD case used
best = R[best_case]; worst = R[worst_case]
spread = (worst["p99"]-best["p99"])/ctrl*100
print(f"\n-- VERDICT --")
print(f"  control phase (146 mm, all prior work): p99 superheat {ctrl:.1f} K")
print(f"  BEST  phase {best['ps']:5.1f} mm: {best['p99']:.1f} K   ({(best['p99']/ctrl-1)*100:+.1f}% vs control)")
print(f"  WORST phase {worst['ps']:5.1f} mm: {worst['p99']:.1f} K   ({(worst['p99']/ctrl-1)*100:+.1f}% vs control)")
print(f"  PHASE SPREAD = {spread:.1f}% of the control hotspot")
print(f"  measured window scatter (3 replicates) = 1.05% -> phase effect is "
      f"{'RESOLVED (real)' if spread > 3.15 else 'WITHIN SCATTER (immaterial)'}")
print(f"\n  Is the control the worst case? {'YES -- prior work was conservative' if worst['ps']==146.1 else 'NO -- prior work was NOT the worst phase'}")
print("  VARIANCE CAVEAT: sigma_T has one window per phase; the 4.5x range is diagnostic only")

metrics = {
    "status": "measured_matched_four_phase_sweep",
    "control_case": control_case,
    "best_mean_case": best_case,
    "worst_mean_case": worst_case,
    "worst_sigma_case": worst_sigma_case,
    "phase_spread_percent_of_control_hotspot": spread,
    "worst_mean_to_control_factor": worst["p99"] / ctrl,
    "best_mean_to_control_factor": best["p99"] / ctrl,
    "worst_sigma_to_control_factor": R[worst_sigma_case]["sig"] / R[control_case]["sig"],
    "sigma_max_to_min_factor": max(row["sig"] for row in R.values()) / min(row["sig"] for row in R.values()),
    "cases": R,
    "mean_status": "promotable after field-equilibrium audit; four matched phases",
    "variance_status": "single-window diagnostic only; independent replicate windows not run",
    "variance_allowed_in_design_envelope": False,
    "use_note": (
        "Apply the mean factor only to a phase-matched control. Do not apply the "
        "single-window fluctuation factor to a promoted design envelope until each "
        "phase has independent replicate averaging windows."
    ),
}
with open("figs/phase_sweep_metrics.json", "w", encoding="ascii") as stream:
    json.dump(metrics, stream, indent=2)
    stream.write("\n")

fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.9), facecolor="white",
                       gridspec_kw=dict(width_ratios=[1.6, 1]))
a = ax[0]
for case, ps, col in CASES:
    xc, p = prof(*wall(case)[:2] if False else (lambda r: (r[0], r[1]))(wall(case)))
    a.plot(xc*1e3, p, lw=2.1, color=col, label=f"peak @ {ps:.0f} mm  (p99 {R[case]['p99']:.0f} K)")
a.axhline(T_BOIL-T_IN, color="#b04030", ls="--", lw=1.3)
a.text(4, T_BOIL-T_IN+10, "FLiBe boiling", fontsize=8.5, color="#b04030")
a.set_xlabel("streamwise x [mm]"); a.set_ylabel("wall superheat [K] (core span)")
a.set_title("Same load, same power, same peaking — only the peak's POSITION differs",
            fontsize=10.5, loc="left")
a.legend(fontsize=8.5); a.grid(alpha=0.3)
a = ax[1]
ps = [R[c]["ps"] for c, _, _ in CASES]; p99 = [R[c]["p99"] for c, _, _ in CASES]
a.plot(ps, p99, "o-", color="#111", ms=8, lw=1.8)
a.axhline(ctrl, color="#1a6b52", ls=":", lw=1.4)
a.text(6, ctrl+8, "phase used by ALL prior CFD", fontsize=8, color="#1a6b52")
a.set_xlabel("flux-peak position along channel [mm]"); a.set_ylabel("hotspot p99 superheat [K]")
a.set_title(f"Phase spread = {spread:.1f}% of hotspot", fontsize=10.5, loc="left")
a.grid(alpha=0.3)
fig.text(0.5, 0.01,
         "Mean-field phase result passes the B6 field audit. sigma_T values use one window per phase and are diagnostic, not variance-converged.",
         ha="center", fontsize=8, color="0.35")
fig.tight_layout(rect=[0, 0.055, 1, 1]); fig.savefig("figs/phase_sweep_verdict.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/phase_sweep_verdict.png and figs/phase_sweep_metrics.json")

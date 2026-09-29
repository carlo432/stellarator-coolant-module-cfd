#!/usr/bin/env python3
"""Operating envelope: hotspot film temperature vs total wall load.

The thermal scalar is PASSIVE (one-way coupled, constant properties), so superheat scales
exactly linearly with load: T(k) = 900 + k*S1. One CFD point at k=2 verifies linearity
survives the numerics (limiters, inletOutlet). The envelope is drawn as a WEDGE:
  upper edge = this tier's wall-modeled LES result;
  lower edge = the same superheat divided by 3.7 using a transferred Ferrero-channel h.
This is a bounded design bracket, not a confidence interval or same-case correction.

Material bands (film-temperature limits, cited loosely at this tier):
  alloy-corrosion film band 973-1073 K (salt-facing Ni-alloy long-life practice)
  FLiBe boiling 1703 K (hard fluid limit)
"""
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_IN = 900.0
H_CORR = 3.7               # bounded transfer ratio to Ferrero's different ARC channel
PHASE_MEAN_WORST = 1.0237193962   # matched sweep: worst mean / 146 mm control

def wall_stats(case, time):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam")
    rd.set_active_time_value(time)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_data_to_point_data()
    Tm = np.asarray(hw["TMean"])
    sig = np.sqrt(np.maximum(np.asarray(hw["TPrime2Mean"]), 0.0))
    hot = Tm >= np.percentile(Tm, 99)
    return float(np.percentile(Tm, 99)), float(sig[hot].mean())

# linearity verification
p99_1, sig_1 = wall_stats("geom_pipeline_plasma", 3.6)
p99_2, sig_2 = wall_stats("geom_pipeline_plasma_k2", 5.0)
ratio = (p99_2 - T_IN) / (p99_1 - T_IN)
print(f"linearity check: superheat(k=2)/superheat(k=1) = {ratio:.3f} (exact = 2.000)")
p99_u, sig_u = wall_stats("geom_pipeline", 2.2)

S_pl_control = p99_1 - T_IN  # prior CFD used one phase (146 mm)
S_pl = PHASE_MEAN_WORST * S_pl_control
S_un = p99_u - T_IN
S_fl = PHASE_MEAN_WORST * S_pl_control + 3 * sig_1

k = np.linspace(0.1, 3.2, 200)
load = 0.5 * k               # mean surface load [MW/m2]
fig, ax = plt.subplots(figsize=(9.6, 6.2), facecolor="white")
for S, col, lab in [(S_un, "#8a8fa3", "uniform load (p99)"),
                    (S_pl, "#7a1f1f", "plasma-shaped p99 (worst mean phase)"),
                    (S_fl, "#b05a00", "plasma +3-sigma (worst mean; control-phase sigma)")]:
    up = T_IN + k * S
    lo = T_IN + k * S / H_CORR
    ax.fill_between(load, lo, up, alpha=0.16, color=col)
    ax.plot(load, up, "-", color=col, lw=2, label=lab + " -- current-model edge")
    ax.plot(load, lo, "--", color=col, lw=1.2, label=lab + " -- Ferrero transfer edge")
ax.errorbar([1.0], [p99_2], yerr=[0.025 * (p99_2 - T_IN)], fmt="o", color="#1a355e", ms=8,
            capsize=4, zorder=5, label=f"k=2 CFD verification (ratio {ratio:.3f} vs 2.000)")
ax.axhspan(973, 1073, color="#c7b299", alpha=0.4)
ax.text(1.55, 1020, "alloy-corrosion film band", fontsize=8.5, color="0.3", va="center", ha="right")
ax.axhline(1703, color="#7a1f1f", lw=1.4, ls="-.")
ax.text(1.55, 1723, "FLiBe boiling 1703 K", fontsize=8.5, color="#7a1f1f", ha="right")
ax.set_xlim(0.05, 1.6); ax.set_ylim(900, 3400)
ax.set_xlabel("mean surface load [MW/m$^2$]  (design point = 0.5)")
ax.set_ylabel("hotspot film temperature [K]")
ax.set_title("Operating envelope: hotspot film T vs wall load\n"
             "(matched phase sweep supplies the 1.024x worst-mean correction)",
             fontsize=11.5)
ax.axvline(0.5, color="0.7", lw=1, ls=":")
ax.legend(fontsize=7.5, loc="lower right", framealpha=0.95)
ax.grid(alpha=0.3)
for sp in ax.spines.values(): sp.set_color("0.75")
fig.text(0.5, 0.005,
         "volumetric NWL scales with the surface load (k applies to all sources) | flicker sigma scales linearly (verified: "
         f"sigma ratio {sig_2/sig_1:.2f}) | four-phase mean spread -2.6% to +2.4% | "
         "phase sigma used one window/case, so no worst-phase variance correction is promoted | "
         "Ferrero edge is a bounded cross-geometry transfer | film-T limits indicative",
         ha="center", fontsize=8, color="0.35")
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig("figs/load_envelope.png", dpi=140, bbox_inches="tight")
print("wrote figs/load_envelope.png")

# allowable-load crossings
for lab, S in [("uniform", S_un), ("plasma", S_pl), ("plasma+3sig", S_fl)]:
    for lim, lname in [(1023, "alloy-band mid"), (1703, "FLiBe boiling")]:
        k_up = (lim - T_IN) / S
        k_lo = (lim - T_IN) / (S / H_CORR)
        print(f"{lab:12s} crosses {lname:14s} at load {0.5*k_up:.2f}-{0.5*k_lo:.2f} MW/m2 (current-model .. Ferrero-transfer)")

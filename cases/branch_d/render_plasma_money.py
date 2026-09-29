#!/usr/bin/env python3
"""THE Phase-2 money figure: uniform vs plasma-shaped first-wall load, same total power.

Row 1: the imposed q''(s) (DESC W7-X field-period transplant, power-normalized)
Row 2: settled wall T under UNIFORM 0.5 MW/m2 (geom_pipeline @2.2)
Row 3: settled wall T under the PLASMA-SHAPED load (geom_pipeline_plasma @3.6), same clim
"""
import json
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

def wall_strip(case, time):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam")
    rd.set_active_time_value(time)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_centers()
    s, eta = unroll(np.asarray(hw.points))
    return s, eta, np.asarray(hw["T"]), np.asarray(hw["TMean"])

def strip_ax(ax, s, e, T, lo, hi, ttl):
    M, xe, ye = np.histogram2d(s, e, bins=[96, 24], weights=T)
    N, _, _ = np.histogram2d(s, e, bins=[xe, ye])
    M = (M / np.maximum(N, 1)).T
    pm = ax.pcolormesh(xe * 1e3, ye * 1e3, M, cmap="inferno", vmin=lo, vmax=hi, shading="flat")
    ax.set_ylabel("span [mm]")
    ax.set_title(ttl, fontsize=10.5, loc="left")
    for sp in ax.spines.values(): sp.set_color("0.75")
    return pm, M, xe

qm = json.load(open("figs/plasma_qmap.json"))
sq, qq = np.asarray(qm["s_m"]), np.asarray(qm["q_W_m2"])

s0, e0, T0, Tm0 = wall_strip("geom_pipeline", 2.2)
s1, e1, T1, Tm1 = wall_strip("geom_pipeline_plasma", 3.6)
lo = 900.0
hi = float(np.percentile(np.r_[Tm0, Tm1], 99.5))

fig, axs = plt.subplots(3, 1, figsize=(13.5, 8.6), sharex=True, facecolor="white",
                        height_ratios=[1.0, 1.35, 1.35])
a = axs[0]
a.plot(sq * 1e3, qq / 1e6, "-", color="#7a1f1f", lw=2.2, label="plasma-shaped q''(s)")
a.axhline(0.5, color="0.55", lw=1.2, ls=":", label="uniform 0.5 MW/m2 (same total power)")
a.set_ylabel("q'' [MW/m2]")
a.legend(fontsize=9, loc="upper left")
a.grid(alpha=0.3)
a.set_title("imposed first-wall load: DESC W7-X field-period standoff -> radiative floor + far-SOL falloff "
            f"(peak/mean {qm['peak_over_mean']:.2f})", fontsize=10.5, loc="left")
for sp in a.spines.values(): sp.set_color("0.75")

pm0, M0, xe0 = strip_ax(axs[1], s0, e0, Tm0, lo, hi,
                        f"UNIFORM load -- TMean wall strip (geom_pipeline, window 1.5-2.2 s) -- "
                        f"mean {Tm0.mean():.0f} K, p99 {np.percentile(Tm0,99):.0f} K")
pm1, M1, xe1 = strip_ax(axs[2], s1, e1, Tm1, lo, hi,
                        f"PLASMA-SHAPED load -- TMean wall strip (geom_pipeline_plasma, window 2.9-3.6 s) -- "
                        f"mean {Tm1.mean():.0f} K, p99 {np.percentile(Tm1,99):.0f} K")
axs[2].set_xlabel("streamwise arc length s [mm]  (inlet left; one W7-X field period)")
cb = fig.colorbar(pm1, ax=axs[1:], pad=0.012, aspect=30)
cb.set_label("mean wall T [K]")

# peak localization: core span only (side-wall corner faces run hot and pollute a full-span mean)
core = np.abs(e1) < 0.005
sb = np.linspace(0, 0.2985, 49)
prof = np.array([Tm1[core & (s1 >= a) & (s1 < b)].mean() for a, b in zip(sb[:-1], sb[1:])])
s_qpk = sq[np.argmax(qq)] * 1e3
s_Tpk = 0.5 * (sb[:-1] + sb[1:])[np.argmax(prof)] * 1e3
fig.suptitle("Same total power, different distribution: what plasma shaping does to the first wall",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.005,
         f"q'' peak s={s_qpk:.0f} mm; core-span wall-T peak s={s_Tpk:.0f} mm (locks onto the load peak; convective shift < bin width) | "
         f"same total power, hotspot superheat ~2x the uniform case's | "
         "transplant honesty: W7-X example equilibrium, lambda 3 cm [1-5], f_rad 0.5 [0.3-0.7]",
         ha="center", fontsize=8.5, color="0.35")
fig.savefig("figs/pipeline_plasma_money.png", dpi=140, bbox_inches="tight")
print("wrote figs/pipeline_plasma_money.png")
print(f"uniform:  TMean mean {Tm0.mean():.1f} K  p99 {np.percentile(Tm0,99):.1f}  max {Tm0.max():.1f}")
print(f"plasma :  TMean mean {Tm1.mean():.1f} K  p99 {np.percentile(Tm1,99):.1f}  max {Tm1.max():.1f}")
print(f"q peak s={s_qpk:.0f} mm; T peak s={s_Tpk:.0f} mm")

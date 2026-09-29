#!/usr/bin/env python3
"""Pipeline result maps (Tasks 1.1/1.2): span-averaged (s,d) field maps + wall-BC comparison.

Unrolls the curved duct into (s = streamwise arc length at mid radius, d = depth from the
heated first wall) and span-averages, which is the honest 2D reduction for this geometry
(a flat plane slice would exit the duct because of the y-wobble).

Fig 1  pipeline_tritium_map.png : span-averaged CtritMean(s,d) + depth-zone means
       -> the Ferrero-analog "does tritium accumulate in recirculation zones?" verdict.
Fig 2  pipeline_firstwall_T_comparison.png : heated-wall T(s, span) strips,
       geom_heat 0.8 (naive molecular-only fixedGradient, over-delivers where y+>11.5)
       vs geom_pipeline 2.2 (exprMixed local-D_eff BC, exact 0.5 MW/m2 + OpenMC zones).
"""
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARC = np.radians(95.0)
R_MID, R_FW, W, H = 0.18, 0.165, 0.030, 0.016

def unroll(points):
    x, y, z = points.T
    r = np.hypot(x, 0.18 - z)
    th = np.arctan2(x, 0.18 - z)
    s = th * R_MID
    d = r - R_FW
    t = th / ARC
    eta = y - 0.01 * np.sin(2 * np.pi * t)
    return s, d, eta

def load(case, time):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam")
    rd.set_active_time_value(time)
    return rd.read()

def cellmap(m, field):
    vol = m["internalMesh"] if "internalMesh" in m.keys() else m.combine()
    cc = vol.cell_centers()
    return unroll(np.asarray(cc.points)), np.asarray(cc[field])

def binned(s, d, v, ns=96, nd=24):
    H2, xe, ye = np.histogram2d(s, d, bins=[ns, nd], weights=v)
    N, _, _ = np.histogram2d(s, d, bins=[xe, ye])
    return (H2 / np.maximum(N, 1)).T, xe, ye

# ---------- Fig 1: tritium ----------
pipe = load("geom_pipeline", 2.2)
(s, d, eta), C = cellmap(pipe, "CtritMean")
M, xe, ye = binned(s, d, C)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 4.6), width_ratios=[3.2, 1],
                             facecolor="white")
pm = a1.pcolormesh(xe * 1e3, ye * 1e3, M * 1e9, cmap="Blues", shading="flat")
cb = fig.colorbar(pm, ax=a1, pad=0.015)
cb.set_label("span-averaged tritium concentration  [nmol m$^{-3}$]")
a1.set_xlabel("streamwise arc length s [mm]  (inlet left, outlet right)")
a1.set_ylabel("depth from first wall d [mm]")
a1.set_title("Tritium builds monotonically along the flow path -- no recirculation trapping",
             fontsize=11.5)
for sp in a1.spines.values(): sp.set_color("0.75")
# depth-zone means (from the run's volFieldValue FOs, window 1.5-2.2)
zC = [7.490e-8, 6.817e-8, 6.841e-8, 6.784e-8, 6.798e-8, 7.338e-8]
zc = (np.arange(6) + 0.5) * 5.0
a2.barh(zc, np.array(zC) * 1e9, height=3.6, color="#4a7fb5", edgecolor="white", linewidth=2)
a2.set_ylim(0, 30); a2.invert_yaxis()
a2.set_xlabel("zone mean C$_T$ [nmol m$^{-3}$]")
a2.set_ylabel("depth zone band [mm]")
a2.set_title("wall-adjacent zones hold ~10% more\n(slower, older fluid)", fontsize=10)
for sp in a2.spines.values(): sp.set_color("0.75")
fig.suptitle("Tritium transport in the stellarator coolant slice -- OpenMC Li-6/Li-7 sources, "
             "LES scalar (window 1.5-2.2 s; mass balance closes to 0.02%)", fontsize=12.5)
fig.text(0.5, 0.012,
         "geom_pipeline: NWL 1.0 MW/m2, Li-6 90%, slice TBR 0.037 (3 cm channel; 0.67 for 0.5 m column) | "
         "outlet mean 136 nmol/m3 = generation/flow exactly | verdict: convective washout, no accumulation hotspots at this curvature",
         ha="center", fontsize=8.5, color="0.35")
fig.tight_layout(rect=[0, 0.045, 1, 0.94])
fig.savefig("figs/pipeline_tritium_map.png", dpi=140)
print("wrote figs/pipeline_tritium_map.png")
plt.close(fig)

# ---------- Fig 2: first-wall T comparison ----------
def wall_strip(case, time, field="T"):
    m = load(case, time)
    hw = m["boundary"]["heated_first_wall"]
    cc = hw.cell_centers()
    (s, d, eta) = unroll(np.asarray(cc.points))
    return s, eta, np.asarray(cc[field])

s0, e0, T0 = wall_strip("geom_heat", 0.8)
s1, e1, T1 = wall_strip("geom_pipeline", 2.2)
lo = 900
hi = float(np.percentile(np.r_[T0, T1], 99))
fig, axs = plt.subplots(2, 1, figsize=(13.5, 6.0), sharex=True, facecolor="white")
for ax, (s_, e_, T_, ttl) in zip(axs, [
    (s0, e0, T0, f"geom_heat @0.8: TRANSIENT (BC on 0.2 s = 1.3 FT), naive BC -- mean {T0.mean():.0f} K (NOT converged)"),
    (s1, e1, T1, f"geom_pipeline @2.2: settled (9 FT soak), exact-flux BC + OpenMC zones -- mean {T1.mean():.0f} K"),
]):
    ns, ne = 96, 24
    M2, xe2, ye2 = np.histogram2d(s_, e_, bins=[ns, ne], weights=T_)
    N2, _, _ = np.histogram2d(s_, e_, bins=[xe2, ye2])
    M2 = (M2 / np.maximum(N2, 1)).T
    pm = ax.pcolormesh(xe2 * 1e3, ye2 * 1e3, M2, cmap="inferno", vmin=lo, vmax=hi, shading="flat")
    ax.set_ylabel("span [mm]")
    ax.set_title(ttl, fontsize=10.5, loc="left")
    for sp in ax.spines.values(): sp.set_color("0.75")
cb = fig.colorbar(pm, ax=axs, pad=0.012, aspect=28)
cb.set_label("instantaneous wall T [K]")
axs[1].set_xlabel("streamwise arc length s [mm]")
fig.suptitle("Heated first-wall temperature: geom_heat's published numbers were a 0.2 s warming transient; "
             "the pipeline case is thermally settled with exactly 0.5 MW/m$^2$ delivered", fontsize=12)
fig.text(0.5, 0.005, "wall dT ~ 340 K at settled state -> effective h ~ 1.5 kW/m2K vs Ferrero-channel ~5.5 kW/m2K: "
         "the wall-modelled near-wall closure remains bounded; C1/C9 thermal certification is void pending exact scalar balance. "
         "Old BC also over-delivered flux up to ~35x on the y+>11.5 face tail (exprMixed removes this).",
         ha="center", fontsize=8.5, color="0.35")
fig.savefig("figs/pipeline_firstwall_T_comparison.png", dpi=140, bbox_inches="tight")
print("wrote figs/pipeline_firstwall_T_comparison.png")
print(f"wall T: old mean {T0.mean():.1f} max {T0.max():.1f} | new mean {T1.mean():.1f} max {T1.max():.1f}")

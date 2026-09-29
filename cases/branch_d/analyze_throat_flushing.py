#!/usr/bin/env python3
"""FLUSHING-PRESCRIPTION CFD TEST -- verdict.

Baseline  = geom_straight_pipeline_plasma (no throat), t=2.7
Throat    = geom_throat_pipeline_plasma  (1.33x pinch at the flux peak), t=4.0
Both: identical heated wall (z=+0.015 plane, same area), identical plasma flux q''(s),
same delivered bulk. ONLY difference = near-wall flow accelerated ~1.33x over a 34 mm strip.

The prescription assumes the DEVELOPED-duct law dT_wall ~ U^-0.70 applies LOCALLY, so a 1.33x
speed-up should cut the local superheat by 1.33^0.70 = 1.219x (-18%). This test asks whether a
REAL contraction delivers that, or whether throat re-development / downstream separation breaks
local similarity.

Outputs: printed verdict + figs/throat_flushing_verdict.png
"""
import numpy as np, pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_IN = 900.0
T_BOIL = 1703.0
EXP = 0.70                         # the measured flushing-law exponent (OUR value, not D-B's 0.80)
SPEEDUP = 1.33
X_PEAK = 0.1461                    # flux-peak location (from plasma_qmap argmax)
STRIP = 0.034                      # throat FWHM
PRED_LOCAL_FACTOR = SPEEDUP ** EXP  # predicted local dT reduction factor (=1.219)

def wall(case, t):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam"); rd.set_active_time_value(t)
    hw = rd.read()["boundary"]["heated_first_wall"].cell_data_to_point_data()
    pts = np.asarray(hw.points)
    Tm = np.asarray(hw["TMean"])
    sig = np.sqrt(np.maximum(np.asarray(hw["TPrime2Mean"]), 0.0))
    return pts, Tm, sig

def scal(Tm, sig):
    p99 = np.percentile(Tm, 99)
    return dict(p99=float(p99), mean=float(Tm.mean()), mx=float(Tm.max()),
                sig_hot=float(sig[Tm >= p99].mean()))

def profile(pts, Tm, nb=60):
    """core-span (|y|<5mm) streamwise superheat profile, binned in x."""
    core = np.abs(pts[:, 1]) < 0.005
    x = pts[core, 0]; T = Tm[core]
    edges = np.linspace(x.min(), x.max(), nb + 1)
    xc = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(x, edges) - 1, 0, nb - 1)
    prof = np.array([T[idx == i].max() - T_IN if np.any(idx == i) else np.nan for i in range(nb)])
    return xc, prof

pb, Tb, sb = wall("geom_straight_pipeline_plasma", 2.7)
pt, Tt, st = wall("geom_throat_pipeline_plasma", 4.0)
B, Th = scal(Tb, sb), scal(Tt, st)
xb, prb = profile(pb, Tb); xt, prt = profile(pt, Tt)

print(f"{'':22s} {'p99 dT':>8s} {'max dT':>8s} {'mean dT':>8s} {'sig_hot':>8s}  film-T(p99)")
for name, r in (("baseline (no throat)", B), ("throat 1.33x", Th)):
    print(f"{name:22s} {r['p99']-T_IN:8.1f} {r['mx']-T_IN:8.1f} {r['mean']-T_IN:8.1f} "
          f"{r['sig_hot']:8.1f}   {r['p99']:7.1f} K")

# --- local reduction at the throat strip (|x-xpeak| < strip/2) vs the same strip in baseline ---
def strip_peak(xc, prof):
    m = np.abs(xc - X_PEAK) < STRIP / 2
    return np.nanmax(prof[m])
loc_b = strip_peak(xb, prb); loc_t = strip_peak(xt, prt)
local_factor = loc_b / loc_t
# --- downstream separation check: hottest point in the expansion, x in (peak+strip/2, peak+3*strip) ---
def down_peak(xc, prof):
    m = (xc > X_PEAK + STRIP / 2) & (xc < X_PEAK + 3 * STRIP)
    return np.nanmax(prof[m]) if np.any(m) else np.nan
dn_b = down_peak(xb, prb); dn_t = down_peak(xt, prt)

print(f"\n-- LOCAL superheat at the flux-peak strip --")
print(f"  baseline {loc_b:7.1f} K   throat {loc_t:7.1f} K   -> local reduction {local_factor:.3f}x")
print(f"  PREDICTED by dT~U^-0.70 (developed-duct law applied locally): {PRED_LOCAL_FACTOR:.3f}x")
verdict = ("MATCHES prediction (local similarity HOLDS)" if abs(local_factor-PRED_LOCAL_FACTOR)<0.06
           else "UNDER-performs (local similarity BREAKS -- re-development)" if local_factor<PRED_LOCAL_FACTOR
           else "OVER-performs (contraction beats the developed law)")
print(f"  VERDICT: {verdict}")
print(f"\n-- downstream expansion (separation check) --")
print(f"  baseline {dn_b:7.1f} K   throat {dn_t:7.1f} K   ({'no new hot spot' if dn_t<=loc_t+5 else 'HOT SPOT appears downstream'})")
print(f"\n-- BOILING (current-model tier edge) --")
print(f"  baseline p99 film-T {B['p99']:7.1f} K vs boiling {T_BOIL}  -> {'ABOVE (boils)' if B['p99']>T_BOIL else 'below'}")
print(f"  throat   p99 film-T {Th['p99']:7.1f} K vs boiling {T_BOIL}  -> {'ABOVE (still boils)' if Th['p99']>T_BOIL else 'BELOW (boiling avoided)'}")

# near-heated-wall streamwise velocity, to PROVE the throat accelerated the flow it should have
def nearwall_U(case, t):
    open(f"{case}/case.foam", "w").close()
    rd = pv.OpenFOAMReader(f"{case}/case.foam"); rd.set_active_time_value(t)
    c = rd.read()["internalMesh"].cell_centers()
    P = np.asarray(c.points); U = np.asarray(c["UMean"])
    near = (P[:, 2] > P[:, 2].max() - 0.0007) & (np.abs(P[:, 1]) < 0.003)
    x = P[near, 0]; ux = U[near, 0]
    edges = np.linspace(0, 0.2984513, 40); xc = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(x, edges) - 1, 0, 38)
    return xc, np.array([ux[idx == i].mean() if np.any(idx == i) else np.nan for i in range(39)])
xub, uub = nearwall_U("geom_straight_pipeline_plasma", 2.7)
xut, uut = nearwall_U("geom_throat_pipeline_plasma", 4.0)
u_ratio = uut / uub

fig, ax = plt.subplots(2, 1, figsize=(10.6, 7.4), facecolor="white", sharex=True,
                       gridspec_kw=dict(height_ratios=[2.1, 1.0], hspace=0.13))
a = ax[0]
a.plot(xb*1e3, prb, lw=2.2, color="#7a1f1f", label=f"baseline, no throat (p99 dT {B['p99']-T_IN:.0f} K)")
a.plot(xt*1e3, prt, lw=2.4, color="#1a6b52", label=f"throat 1.33x (p99 dT {Th['p99']-T_IN:.0f} K)")
a.axvspan((X_PEAK-STRIP/2)*1e3, (X_PEAK+STRIP/2)*1e3, color="#1a6b52", alpha=0.10)
a.text(X_PEAK*1e3, a.get_ylim()[1]*0.97, "throat strip", ha="center", fontsize=8.5, color="#1a6b52")
a.axhline(T_BOIL-T_IN, color="#b04030", ls="--", lw=1.4)
a.text(6, T_BOIL-T_IN+8, f"FLiBe boiling ({T_BOIL:.0f} K film)", fontsize=8.5, color="#b04030")
a.set_ylabel("wall superheat  T$-$T$_{in}$  [K]")
a.set_title(f"Flushing prescription CFD test: a real 1.33x throat barely cuts the hotspot\n"
            f"local reduction {local_factor:.2f}x MEASURED vs {PRED_LOCAL_FACTOR:.2f}x predicted (dT~U$^{{-0.70}}$) "
            f"-> local similarity BREAKS; boiling NOT avoided", fontsize=10.5, loc="left")
a.legend(fontsize=9); a.grid(alpha=0.3)
a = ax[1]
a.plot(xut*1e3, u_ratio, lw=2.0, color="#2a5d8f")
a.axhline(SPEEDUP, color="0.5", ls="--", lw=1.2); a.text(6, SPEEDUP+0.01, f"target {SPEEDUP}x", fontsize=8, color="0.4")
a.axhline(1.0, color="0.7", ls=":", lw=1)
a.axvspan((X_PEAK-STRIP/2)*1e3, (X_PEAK+STRIP/2)*1e3, color="#1a6b52", alpha=0.10)
a.set_ylabel("near-wall U$_x$\n throat / baseline"); a.set_xlabel("streamwise x [mm]")
a.set_ylim(0.8, 1.45); a.grid(alpha=0.3)
a.text(X_PEAK*1e3, 1.34, f"flow DID speed up {u_ratio[np.nanargmin(np.abs(xut-X_PEAK))]:.2f}x here\n"
       f"...yet the wall barely cooled", ha="center", fontsize=8.2, color="#2a5d8f")
fig.savefig("figs/throat_flushing_verdict.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/throat_flushing_verdict.png")

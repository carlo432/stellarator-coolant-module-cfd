#!/usr/bin/env python3
"""THE PRESCRIPTION: how much local coolant acceleration does the plasma-shaped hotspot need?

Closing the project's arc. plasma -> neutrons -> coolant -> hotspot (2.6x amplification,
straddles FLiBe boiling) -> conduction CANNOT fix it (spreading length 7 mm vs a 53 mm peak)
-> therefore the remedy is FLOW-SIDE. This turns that sentence into a number.

INPUTS, all measured by this project, none assumed:
  * the D-branch flushing-speed law, wall-resolved CHT velocity sweep:  dT_wall ~ U^-0.70
    (fitted over Re 5k-20k; Dittus-Boelter would give -0.80, so this is the measured, not the
    textbook, exponent)
  * the plasma load map q''(s) (DESC W7-X standoff transplant)
  * the wall's conduction filter A_m (analytic; it still helps a few %, so include it)
  * this tier's effective h ~ 1.5 kW/m2K, and Ferrero's different-channel h ~ 5.5 kW/m2K -- the
    envelope wedge, carried through so the prescription is a RANGE, not a point.

WHAT IS ASKED: the local velocity multiplier m(s) = U_local/U_bulk needed to hold the
coolant-side film temperature at or below a target (alloy band top 1073 K; FLiBe boiling
1703 K), given the load actually delivered at s.

  film superheat  dT(s) = dT_ref(s) * m(s)^-0.70      (flushing law, local)
  require dT(s) <= T_target - T_bulk   =>   m(s) >= (dT_ref(s) / dT_allow)^(1/0.70)

Then convert the required multiplier into a duct-area contraction (mass conservation,
incompressible): A_throat/A_duct = 1/m.

HONESTY. The flushing law was fitted on wall-resolved CHT of a straight channel with a
UNIFORM load, and is applied here pointwise to a peaked load -- i.e. a local-similarity
assumption. It ignores the streamwise history the LES shows (turbulent mixing shaves the
peak, which is why amplification sits below 1:1). It therefore OVERSTATES the required
acceleration. The contraction also costs pumping (~m^2 on dynamic head) and that is reported.
This is a sizing prescription to hand an advisor, not a certified design.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_BULK, T_BOIL, ALLOY_HI = 900.0, 1703.0, 1073.0
N_LAW = 0.70                       # measured exponent, dT ~ U^-n
H_LES, H_CORR = 1.5e3, 5.5e3
KAPPA, THICK = 23.6, 3e-3          # Hastelloy-N, 3 mm
U_BULK = 1.99

q = json.load(open("figs/plasma_qmap.json"))
s = np.asarray(q["s_m"]); qs = np.asarray(q["q_W_m2"])
L = s[-1] - s[0]; N = len(s)
m_ = np.arange(N); k = m_ * np.pi / L
C = np.array([np.trapezoid(qs * np.cos(mm * np.pi * (s - s[0]) / L), s) for mm in m_])
c = 2.0 / L * C; c[0] = np.trapezoid(qs, s) / L

def film_superheat(h):
    """with the wall's conduction filter included (it is the honest reference)."""
    den = KAPPA * k * np.sinh(k * THICK) + h * np.cosh(k * THICK); den[0] = h
    return sum(c[i] / den[i] * np.cos(k[i] * (s - s[0])) for i in range(1, N)) + c[0] / h

print(f"flushing law exponent n = {N_LAW} (measured, wall-resolved CHT sweep); U_bulk = {U_BULK} m/s\n")
print(f"{'target':22s} {'h':>6s} {'peak film [K]':>14s} {'max m needed':>13s} {'throat area':>12s} {'pumping ~m^2':>13s}")
out = {}
for htag, h in (("LES", H_LES), ("corr", H_CORR)):
    dT = film_superheat(h)
    out[htag] = dT
    for ttag, T in (("alloy band top 1073 K", ALLOY_HI), ("FLiBe boiling 1703 K", T_BOIL)):
        allow = T - T_BULK
        mreq = np.maximum(dT / allow, 1.0) ** (1.0 / N_LAW)
        mmax = mreq.max()
        print(f"{ttag:22s} {htag:>6s} {T_BULK+dT.max():13.0f} {mmax:12.2f}x {1/mmax:11.2f} {mmax**2:12.1f}x")

# fraction of the wall that needs help, and the length over which
print()
for htag, h in (("LES", H_LES), ("corr", H_CORR)):
    dT = out[htag]
    for ttag, T in (("alloy", ALLOY_HI), ("boiling", T_BOIL)):
        need = dT > (T - T_BULK)
        frac = 100 * need.mean()
        span = (s[need].max() - s[need].min()) * 1e3 if need.any() else 0.0
        print(f"h={htag:4s} vs {ttag:8s}: {frac:5.1f}% of the wall needs acceleration, over a {span:5.0f} mm span")

fig, (a1, a2) = plt.subplots(2, 1, figsize=(10.2, 7.6), facecolor="white", sharex=True,
                             gridspec_kw=dict(height_ratios=[1, 1], hspace=0.16))
for hkey, htag, h, col in (("LES", "current-model edge", H_LES, "#7a1f1f"),
                           ("corr", "Ferrero transfer edge", H_CORR, "#1a6b52")):
    dT = out[hkey]
    a1.plot(s * 1e3, T_BULK + dT, "-", lw=2, color=col, label=f"film T, {htag}")
a1.axhline(ALLOY_HI, color="#8a6d3b", ls="--", lw=1.3); a1.axhline(T_BOIL, color="#7a1f1f", ls="-.", lw=1.3)
a1.text(2, ALLOY_HI + 25, "alloy band top 1073 K", fontsize=8.5, color="#8a6d3b")
a1.text(2, T_BOIL + 25, "FLiBe boiling 1703 K", fontsize=8.5, color="#7a1f1f")
a1.set_ylabel("coolant-side film T [K]"); a1.legend(fontsize=9, loc="upper right"); a1.grid(alpha=0.3)
a1.set_title("The prescription: local flushing needed to hold the plasma-shaped hotspot", fontsize=12, loc="left")
for sp in a1.spines.values(): sp.set_color("0.75")

for htag, col in (("LES", "#7a1f1f"), ("corr", "#1a6b52")):
    dT = out[htag]
    for T, ls, tt in ((ALLOY_HI, "-", "alloy"), (T_BOIL, ":", "boiling")):
        mreq = np.maximum(dT / (T - T_BULK), 1.0) ** (1.0 / N_LAW)
        a2.plot(s * 1e3, mreq, ls, lw=1.8, color=col,
                label=f"{'current model' if htag=='LES' else 'Ferrero transfer'} -> {tt}")
a2.axhline(1.0, color="0.5", lw=1)
a2.set_xlabel("streamwise arc length s [mm]")
a2.set_ylabel("required local speed-up  $m = U_{local}/U_{bulk}$")
a2.legend(fontsize=8.5, ncol=2, loc="upper right"); a2.grid(alpha=0.3)
for sp in a2.spines.values(): sp.set_color("0.75")

fig.suptitle("Conduction cannot fix it, so the flow must: sizing the local flushing remedy",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.005,
         "flushing law dT ~ U^-0.70 measured on this project's wall-resolved CHT sweep, applied pointwise (local-similarity: "
         "ignores streamwise mixing, so it OVERSTATES the requirement) | conduction filter included | contraction costs ~m^2 "
         "in dynamic head",
         ha="center", fontsize=7.8, color="0.35")
fig.subplots_adjust(top=0.90, bottom=0.10)
fig.savefig("figs/flushing_remedy.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/flushing_remedy.png")

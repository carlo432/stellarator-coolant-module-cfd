#!/usr/bin/env python3
"""The OPERATING WINDOW: what combination of inlet temperature and local flushing keeps the
plasma-shaped first-wall hotspot inside its material limits?

The remedy analysis showed two very different problems hiding under one "the alloy band is
exceeded" bullet:
  * staying below FLiBe BOILING is cheap (a ~1.3x local speed-up over ~34 mm, even at the
    current-model edge);
  * staying inside the ALLOY FILM BAND is not a flushing problem at all -- at T_in = 900 K it
    demands a 12x speed-up (142x dynamic head) at the current-model edge. It is an INLET
    TEMPERATURE problem: with a 900 K inlet and a 1073 K limit there are only 173 K of film
    budget to spend.

So plot the actual design space. Two levers, both real:
  m    = local coolant speed-up at the flux peak (flushing law, dT ~ U^-0.70)
  T_in = coolant inlet temperature, bounded below by FLiBe freezing (~732 K)

Constraint at the hotspot:   T_in + dT_peak(h) * m^-0.70  <=  T_limit
  =>  m >= ( dT_peak(h) / (T_limit - T_in) )^(1/0.70)

dT_peak(h) is the CONDUCTION-FILTERED film superheat from analyze_conduction_credit's exact
wall solution, evaluated at both edges of the envelope wedge (h = 1.5 and 5.5 kW/m2K).

Caveats carried on the figure: passive-scalar superheat is h-linear and property-independent
here, so lowering T_in does not (in this model) change dT_peak; in reality a colder inlet means
a MORE viscous FLiBe and a slightly worse h, so the T_in lever is optimistic at the cold end.
Pumping cost of a contraction scales ~m^2 in dynamic head. Local-similarity overstates m.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_BOIL, ALLOY_HI, T_FREEZE = 1703.0, 1073.0, 732.0
N_LAW = 0.70
H_LES, H_CORR = 1.5e3, 5.5e3
KAPPA, THICK, T_REF = 23.6, 3e-3, 900.0
PHASE_WORST = 1.024  # measured four-phase sweep; control was 2.4% below worst

q = json.load(open("figs/plasma_qmap.json"))
s = np.asarray(q["s_m"]); qs = np.asarray(q["q_W_m2"])
L = s[-1] - s[0]; N = len(s)
mm_ = np.arange(N); k = mm_ * np.pi / L
C = np.array([np.trapezoid(qs * np.cos(j * np.pi * (s - s[0]) / L), s) for j in mm_])
c = 2.0 / L * C; c[0] = np.trapezoid(qs, s) / L

def dT_peak(h):
    den = KAPPA * k * np.sinh(k * THICK) + h * np.cosh(k * THICK); den[0] = h
    th = sum(c[i] / den[i] * np.cos(k[i] * (s - s[0])) for i in range(1, N)) + c[0] / h
    return th.max()

dP_control = {"LES": dT_peak(H_LES), "corr": dT_peak(H_CORR)}
dP = {name: PHASE_WORST * value for name, value in dP_control.items()}
print(f"conduction-filtered peak film superheat: current model {dP['LES']:.0f} K, Ferrero transfer {dP['corr']:.0f} K\n")
print("phase correction: curves use measured worst phase (+2.4% vs the prior control); "
      "four-phase control-relative range is -2.6% to +2.4%\n")

Tin = np.linspace(T_FREEZE + 8, 940, 400)
def m_req(dt, Tlim):
    allow = np.maximum(Tlim - Tin, 1e-6)
    return np.maximum(dt / allow, 1.0) ** (1.0 / N_LAW)

print(f"{'edge':6s} {'limit':10s} {'m at T_in=900':>14s} {'T_in for m=1':>14s} {'feasible above freeze?':>24s}")
for htag in ("LES", "corr"):
    for ltag, Tl in (("alloy 1073", ALLOY_HI), ("boiling 1703", T_BOIL)):
        m900 = max(dP[htag] / (Tl - 900.0), 1.0) ** (1 / N_LAW)
        Tin1 = Tl - dP[htag]                       # T_in that needs no speed-up
        ok = "YES" if Tin1 > T_FREEZE else f"NO (needs T_in={Tin1:.0f} K < freeze)"
        print(f"{htag:6s} {ltag:10s} {m900:13.2f}x {Tin1:13.0f} K {ok:>24s}")

fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8), facecolor="white", sharey=True)
for ax, (ltag, Tl, col) in zip(axes, (("alloy film band top (1073 K)", ALLOY_HI, "#8a6d3b"),
                                      ("FLiBe boiling (1703 K)", T_BOIL, "#7a1f1f"))):
    for htag, ls, lab, al in (("LES", "-", "current-model edge (this tier)", 0.16),
                              ("corr", "--", "Ferrero-channel transfer edge (bounded)", 0.10)):
        m = m_req(dP[htag], Tl)
        ax.plot(Tin, m, ls, lw=2.2, color=col, label=lab)
        # INFEASIBLE is BELOW the requirement curve: you have less speed-up than the limit needs
        ax.fill_between(Tin, 0.8, m, color=col, alpha=al)
    ax.axvline(T_FREEZE, color="#2a5d8f", lw=1.6)
    ax.text(T_FREEZE + 4, 19, "FLiBe freezing 732 K", rotation=90, fontsize=8.5, color="#2a5d8f", va="top")
    ax.axvline(900, color="0.5", ls=":", lw=1.2)
    ax.text(901, 24, "design inlet\n900 K", fontsize=8, color="0.4")
    ax.axhline(1.0, color="0.4", lw=1)
    ax.axhline(2.0, color="0.75", ls=":", lw=1)
    ax.text(742, 2.08, "m = 2 (4x dynamic head)", fontsize=7.5, color="0.45")
    ax.set_yscale("log"); ax.set_ylim(0.8, 30)
    ax.set_xlabel("coolant inlet temperature $T_{in}$ [K]")
    ax.set_title(f"hold the film below {ltag}", fontsize=11)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8.5, loc="center left")
    for sp in ax.spines.values(): sp.set_color("0.75")
    ax.text(0.5, 0.955, "unshaded = FEASIBLE (speed-up exceeds what the limit demands)",
            transform=ax.transAxes, ha="center", fontsize=7.5, color="0.45")
    ax.text(0.5, 0.06, "shaded = limit VIOLATED at the hotspot", transform=ax.transAxes,
            ha="center", fontsize=7.5, color="0.35")
axes[0].set_ylabel("required local speed-up  $m = U_{local}/U_{bulk}$")

fig.suptitle("Operating window: worst measured load/channel phase (+2.4% vs prior control)",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.005,
         "flushing law dT ~ U$^{-0.70}$ (measured, wall-resolved CHT sweep) with the wall's exact conduction filter | "
         "four-phase spread is 5.0% (-2.6% to +2.4% vs control) | Ferrero edge is a cross-geometry design anchor | "
         "a colder inlet raises FLiBe viscosity, so the T$_{in}$ lever is optimistic at the cold end",
         ha="center", fontsize=7.8, color="0.35")
fig.subplots_adjust(top=0.86, bottom=0.16)
fig.savefig("figs/operating_window.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/operating_window.png")

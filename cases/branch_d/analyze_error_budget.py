#!/usr/bin/env python3
"""THE ERROR BUDGET for the headline number: the plasma-shaped first-wall hotspot temperature.

Motivation (earned the hard way). Every claim in this project that failed audit failed at a
caveat that was ASSERTED rather than MEASURED: "same total power" (wall areas differed 9.1%),
"curvature causes the penalty" (half of it was span waviness), "within window scatter" (the
scatter was 3x smaller than assumed, so the MHD verdict was wrong). This table exists so no
remaining number rests on an unmeasured "~".

Each contribution is tagged by how it is KNOWN:
  MEASURED   - quantified in this project from runs/replicates
  DERIVED    - exact solution of a stated model (conduction filter)
  ESTIMATED  - literature correlation applied to our conditions (Sieder-Tate)
  BOUNDED    - scenario span between the current model and a transferred channel anchor

Sign convention: a POSITIVE contribution means the reported LES hotspot is TOO HOT by that
much under the stated scenario (i.e. applying that term lowers the wall-temperature estimate).

The 3.7 heat-transfer edge is transferred from a different Ferrero channel. Therefore the
combined extrema form a model-to-anchor DESIGN WEDGE, not a confidence interval or a guarantee
that the physical truth lies between them.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = 2097.0     # LES plasma hotspot p99, curved+wavy duct, K
T_IN = 900.0
SUP = BASE - T_IN

# name, low%, high%, kind, note   (percent of superheat; + lowers the scenario estimate)
ITEMS = [
    # This is a scenario span between two different modeled channels, not a confidence bound.
    # Keep zero as the current-model edge and 73% as the transferred Ferrero edge.
    ("wall-model / Ferrero h gap\n(cross-geometry transfer)", 0.0, 73.0, "BOUNDED",
     "project eff. h 1.5 vs Ferrero-channel 5.5 kW/m2K; cross-geometry transfer edge = -73%.\n"
     "This is a bounded assumption, not a same-case measured correction."),
    ("constant-property nu\n(FLiBe mu is Arrhenius)", 10.0, 20.0, "ESTIMATED",
     "Sieder-Tate (mu_b/mu_w)^0.14 at the settled wall -> Nu +17%, superheat -15%.\n"
     "OVERLAPS the transferred h gap above -- shown for attribution, excluded from the total."),
    ("solid-wall conduction\nspreading", 3.0, 7.0, "DERIVED",
     "exact linear wall solution; 3 mm Ni alloy, spreading length 6.9 mm vs 53 mm peak."),
    ("span waviness\n(modelling artifact)", -2.6, 2.6, "MEASURED",
     "wavy duct's plasma hotspot is 2.6% COOLER than clean; 2.5 sd -> not resolved. Sign unknown."),
    ("MHD Lorentz drag\n@ 9.2 T (if claimed)", -6.7, 0.0, "MEASURED",
     "raises the mean hotspot 6.7% (6.4 sd). NEGATIVE: the LES is too COOL if B is claimed."),
    ("averaging-window scatter\n(1 sd)", -1.05, 1.05, "MEASURED",
     "3 independent windows on the plasma baseline: sup99 = 1201.8 +- 12.6 K."),
    ("load/channel phase\n(4 positions)", -2.6, 2.4, "MEASURED",
     "same load, power and peaking at 4 phases; control-relative range -2.6% to +2.4%."),
]

print(f"HEADLINE: LES plasma hotspot p99 = {BASE:.0f} K  (superheat {SUP:.0f} K over a {T_IN:.0f} K inlet)\n")
print(f"{'contribution':40s} {'range [% of superheat]':>24s}  {'kind':10s}")
for n, lo, hi, kind, _ in ITEMS:
    lbl = n.replace("\n", " ")
    print(f"{lbl:40s} {lo:+10.1f} .. {hi:+7.1f}   {kind:10s}")

# Combine the transferred h edge, conduction, MHD, waviness, scatter, and phase.
# The property term is excluded because it overlaps the transferred h gap.
lo_tot = 0.0 + 3.0 - 6.7 - 2.6 - 1.05 - 2.6       # hottest scenario edge
hi_tot = 73.0 + 7.0 + 0.0 + 2.6 + 1.05 + 2.4      # coolest scenario edge
T_hot = T_IN + SUP * (1 - lo_tot / 100)
T_cool = T_IN + SUP * (1 - hi_tot / 100)
print(f"\ncombined correction to the superheat: {lo_tot:.1f}% .. {hi_tot:.1f}%")
print(f"=> model-to-anchor design wedge = [{T_cool:.0f}, {T_hot:.0f}] K   (LES reports {BASE:.0f} K)")
boil_in = T_cool < 1703 < T_hot
print(f"   FLiBe boiling 1703 K sits {'INSIDE' if boil_in else 'OUTSIDE'} that wedge -> "
      f"{'the boiling question is STILL OPEN across the stated scenarios' if boil_in else 'boiling is settled'}")
print(f"   alloy band top 1073 K: {'exceeded even at the COOLEST bound' if T_cool > 1073 else 'reachable at the coolest bound'}")
print("\nproperty term (10-20%) is EXCLUDED from the total: it overlaps the transferred h gap.")
print("The remaining 3.7x model-to-anchor difference is not decomposed by this cross-geometry")
print("comparison and must not be assigned to near-wall closure without a matched test.")

fig, ax = plt.subplots(figsize=(11.4, 6.2), facecolor="white")
KC = {"MEASURED": "#1a6b52", "DERIVED": "#2a5d8f", "ESTIMATED": "#c07830", "BOUNDED": "#7a1f1f"}
y = np.arange(len(ITEMS))[::-1]
for yi, (n, lo, hi, kind, _) in zip(y, ITEMS):
    ax.barh(yi, hi - lo, left=lo, height=0.55, color=KC[kind], alpha=0.85)
    ax.plot([lo, hi], [yi, yi], color="0.2", lw=0.8)
    ax.text(hi + 1.2, yi, f"{lo:+.1f} .. {hi:+.1f}%", va="center", fontsize=8.5, color="0.3")
ax.axvline(0, color="0.4", lw=1.2)
ax.set_yticks(y, [n for n, *_ in ITEMS], fontsize=8.5)
ax.set_xlabel("scenario contribution to hotspot superheat [%; positive lowers the estimate]")
ax.set_title("Error budget for the plasma-shaped first-wall hotspot\n"
             f"LES reports {BASE:.0f} K; the stated scenarios span a design wedge of "
             f"[{T_cool:.0f}, {T_hot:.0f}] K", fontsize=12)
ax.grid(alpha=0.3, axis="x")
for sp in ax.spines.values(): sp.set_color("0.75")
hs = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.85) for c in KC.values()]
ax.legend(hs, KC.keys(), fontsize=8.5, loc="lower right", title="how it is known", title_fontsize=8.5)
fig.text(0.5, 0.005,
         "the property term is PART OF the h wedge (Ferrero's model uses real FLiBe), so it is shown for attribution "
         "and excluded from the total | MHD is negative: the LES is too COOL if B >= 5 T is claimed | "
         "the Ferrero edge is a cross-geometry transfer, so this is a design wedge, not a confidence interval",
         ha="center", fontsize=7.8, color="0.35")
fig.subplots_adjust(left=0.22, bottom=0.14, top=0.88)
fig.savefig("figs/error_budget.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/error_budget.png")

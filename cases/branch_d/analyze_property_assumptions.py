#!/usr/bin/env python3
"""Re-test two constitutive assumptions against the hotspot temperatures we now actually have.

Both were justified in scoping_buoyancy_mhd.md using a MODEST wall dT. The plasma-shaped load
raised the hotspot far above that, so the justifications are stale even where the conclusions
survive. No new runs: this is a properties/dimensionless-group audit of the banked results.

1. BUOYANCY. Scoping quoted Ri <= 5e-4 -> neglect. Recompute Gr/Re^2 with the REAL wall-to-bulk
   dT (up to ~1223 K at the plasma hotspot).
2. CONSTANT PROPERTIES. The runs use a single nu = 3.09278e-6 m2/s. FLiBe viscosity is strongly
   Arrhenius, mu(T) = 1.16e-4 * exp(3755/T) Pa.s, so mu varies by ~10x between bulk and hotspot
   film. Quantify with the Sieder-Tate liquid correction Nu/Nu_cp = (mu_b/mu_w)^0.14.

IMPORTANT (avoid double-counting): the load-envelope wedge already carries an h-correction of
3.7x down to Ferrero's channel correlation, and that correlation is fitted to real, variable-
property FLiBe. So the Sieder-Tate factor below is NOT an additional independent correction to
stack on the wedge. It overlaps the model-to-anchor difference, but the cross-geometry
comparison cannot quantify how much of that difference comes from properties versus geometry,
closure, or operating state.
"""
import math

RHO, NU_CASE = 1940.0, 3.09278e-6
G, BETA, DH, U = 9.81, 2.0e-4, 0.0207, 1.99
T_BULK, T_BOIL = 900.6, 1703.0

def mu(T):  return 1.16e-4 * math.exp(3755.0 / T)          # Pa.s
def nu(T):  return mu(T) / RHO

T_ref = 3755.0 / math.log(NU_CASE * RHO / 1.16e-4)
Re = U * DH / NU_CASE
print(f"case nu = {NU_CASE:.6g} m2/s  <=>  mu = {NU_CASE*RHO*1e3:.3f} mPa.s  <=>  T_ref = {T_ref:.0f} K")
print(f"(so the constant-property runs are evaluated near the mean fluid temperature, not the wall)")
print(f"Re_Dh = {Re:.0f}\n")

print("1. BUOYANCY  (Gr/Re^2; forced convection dominant while << 0.1)")
print(f"   {'case':24s} {'wall dT [K]':>11s} {'Gr':>11s} {'Gr/Re^2':>9s}  {'margin to 0.1':>13s}")
for tag, dT in (("scoping assumption", 50.0), ("uniform hotspot", 536.3),
                ("plasma hotspot", 1223.0)):
    Gr = G * BETA * dT * DH**3 / NU_CASE**2
    ri = Gr / Re**2
    print(f"   {tag:24s} {dT:11.0f} {Gr:11.3e} {ri:9.4f}  {0.1/ri:12.0f}x")
print("   VERDICT: forced convection dominates the tested range (0.013 << 0.1), but the margin")
print("   is ~8x, not the ~200x the")
print("   scoping number implied. scoping_buoyancy_mhd.md's 'Ri <= 5e-4' should be restated.\n")

print("2. CONSTANT PROPERTIES  (Sieder-Tate  Nu/Nu_cp = (mu_b/mu_w)^0.14)")
print(f"   {'wall T [K]':>11s} {'mu_w/mu_b':>10s} {'Nu factor':>10s} {'superheat factor':>17s}   note")
for tag, Tw in (("settled mean wall", 1240.0), ("uniform hotspot", 1436.3),
                ("alloy band mid", 1023.0), ("plasma hotspot", 2097.0)):
    r = mu(Tw) / mu(T_BULK)
    f = (1.0 / r) ** 0.14
    note = "ABOVE FLiBe boiling -- extrapolation invalid" if Tw > T_BOIL else \
           ("outside Sieder-Tate range (mu_b/mu_w > 9.75)" if 1/r > 9.75 else "in range")
    print(f"   {Tw:11.0f} {r:10.3f} {f:9.3f}x {1/f:16.3f}x   {note}")
print("\n   Reading: at the settled mean wall, variable properties would raise Nu ~20% and so")
print("   LOWER the predicted superheat ~17%. Our constant-property wall T is therefore biased")
print("   HIGH under this correlation estimate. This does not make every project wall temperature")
print("   a physical bound; the later matched A3 cases provide the transferable bracket.")
print("   Do NOT stack this on the envelope wedge: Ferrero's different-channel model (the transfer edge)")
print("   already uses real variable-property FLiBe. The factor overlaps the 3.7x difference,")
print("   whose decomposition requires a matched geometry/closure comparison.")

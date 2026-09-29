#!/usr/bin/env python3
"""C6: FLiBe property-realism sensitivity note. Quantifies the constant-property assumption (T-dependent
viscosity/Pr over the wall superheat and the full operating range) and the buoyancy significance (Gr/Re^2)
for the channel, to bound the error in the constant-property RANS-CHT / LES baseline.

Correlations (Romatoski & Hu 2017 review; Cantor viscosity), FLiBe (LiF-BeF2 2:1):
  rho [kg/m3] = 2415.6 - 0.49072*T
  mu  [Pa s]  = 1.16e-4 * exp(3755/T)
  cp  [J/kgK] ~ 2386 (weak T-dep, taken constant)
  k   [W/mK]  ~ 1.1  (weak T-dep, taken constant)
"""
import math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

CP, K = 2386.0, 1.1
def rho(T): return 2415.6 - 0.49072*T
def mu(T):  return 1.16e-4*math.exp(3755.0/T)
def nu(T):  return mu(T)/rho(T)
def Pr(T):  return mu(T)*CP/K

# project reference: momentum LES used nu=3.09278e-6 -> solve for the implied reference T
def find_T_for_nu(target, lo=700, hi=1000):
    for _ in range(100):
        mid=0.5*(lo+hi)
        if nu(mid)>target: lo=mid
        else: hi=mid
    return 0.5*(lo+hi)
Tref = find_T_for_nu(3.09278e-6)

print(f"Project nu=3.09278e-6 m^2/s  ->  implied reference T = {Tref:.1f} K")
print(f"  at Tref: rho={rho(Tref):.1f}  mu={mu(Tref):.5f}  nu={nu(Tref):.4e}  Pr={Pr(Tref):.2f}")
print()
print(f"{'T[K]':>6} {'rho':>8} {'mu[Pa.s]':>10} {'nu[m2/s]':>11} {'Pr':>7}")
for T in (730, 800, 850, 900, 950, 960):
    print(f"{T:6d} {rho(T):8.1f} {mu(T):10.5f} {nu(T):11.4e} {Pr(T):7.2f}")

# --- constant-property error over the wall superheat (~28 K from the CHT baseline) ---
dT_sup = 28.0
Tw = Tref - dT_sup    # near-wall is COOLER than bulk in a wall-cooled... here wall is hotter; use +/- band
T_lo, T_hi = Tref, Tref + dT_sup
print()
print(f"Over the ~{dT_sup:.0f} K wall superheat ({Tref:.0f}->{Tref+dT_sup:.0f} K):")
print(f"  nu: {nu(Tref):.3e} -> {nu(Tref+dT_sup):.3e}  ({100*(nu(Tref+dT_sup)/nu(Tref)-1):+.1f}%)")
print(f"  Pr: {Pr(Tref):.2f} -> {Pr(Tref+dT_sup):.2f}  ({100*(Pr(Tref+dT_sup)/Pr(Tref)-1):+.1f}%)")
print(f"Over the full operating range 730->950 K:")
print(f"  nu: {nu(730):.3e} -> {nu(950):.3e}  (x{nu(730)/nu(950):.2f})")
print(f"  Pr: {Pr(730):.1f} -> {Pr(950):.1f}  (x{Pr(730)/Pr(950):.2f})")

# --- buoyancy significance: Gr/Re^2 (Richardson number) ---
g = 9.81
beta = 0.49072/rho(Tref)         # thermal expansion -1/rho drho/dT
Dh = 0.0208696                   # C2 curved hydraulic diameter
Re = 10000.0
Gr = g*beta*dT_sup*Dh**3/nu(Tref)**2
Ri = Gr/Re**2
print()
print(f"Buoyancy: beta={beta:.3e} 1/K  Gr={Gr:.3e}  Re={Re:.0f}  Ri=Gr/Re^2={Ri:.3e}")
print(f"  -> {'forced convection DOMINATES; buoyancy negligible' if Ri<0.1 else 'MIXED convection; buoyancy matters'} (threshold Ri~0.1)")

# --- plot ---
Ts=[730+5*i for i in range(47)]
fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.5))
a1.plot(Ts,[nu(T)*1e6 for T in Ts],'-',color='tab:blue',label=r'$\nu$ [$\times10^{-6}$ m$^2$/s]')
a1.axvline(Tref,ls=':',color='gray'); a1.axvspan(Tref,Tref+dT_sup,alpha=0.15,color='tab:red',label='wall superheat band')
a1.set_xlabel('T [K]'); a1.set_ylabel(r'$\nu\;[\times10^{-6}$ m$^2$/s]'); a1.set_title('FLiBe kinematic viscosity'); a1.grid(alpha=0.3); a1.legend(fontsize=8)
a2.plot(Ts,[Pr(T) for T in Ts],'-',color='tab:green')
a2.axvline(Tref,ls=':',color='gray'); a2.axvspan(Tref,Tref+dT_sup,alpha=0.15,color='tab:red')
a2.set_xlabel('T [K]'); a2.set_ylabel('Pr'); a2.set_title(f'FLiBe Prandtl (Pr_ref={Pr(Tref):.1f})'); a2.grid(alpha=0.3)
out=Path(__file__).resolve().parent/"figs"/"c6_flibe_property_realism.png"
fig.tight_layout(); fig.savefig(out,dpi=140); print("wrote",out)

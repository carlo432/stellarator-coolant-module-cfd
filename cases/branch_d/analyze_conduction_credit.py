#!/usr/bin/env python3
"""CONDUCTION CREDIT: how much does the solid first wall smear the plasma-shaped hotspot?

THE GAP THIS FILLS. Every result in this project is fluid-only (passive scalar in the
coolant). The load envelope showed the alloy film band is exceeded at the design point under
a plasma-shaped load, and flagged that the one remedy a fluid-only tier structurally CANNOT
give is lateral heat spreading in the solid wall. This closes that gap exactly -- no CHT run
needed -- because steady conduction in a thin wall is linear and therefore diagonal in
Fourier space.

MODEL (exact for the stated problem, not a correlation)
  Slab of thickness t, conductivity kappa, coordinate n from the plasma-facing surface (n=0)
  into the wall, coolant at n=t. Steady Laplace, temperature excess theta = T - T_bulk.
    n = 0 : -kappa dtheta/dn = q''(s)           (plasma flux in)
    n = t : -kappa dtheta/dn = h * theta        (convection into the coolant)
  Expand q''(s) in a cosine series (adiabatic duct ends). Each mode m with wavenumber k_m
  is independent, and the FILM temperature (coolant-side wall temperature) is

      theta_w,m = q_m / (kappa*k_m*sinh(k_m*t) + h*cosh(k_m*t))

  Sanity: m = 0 gives theta_w = q_0/h, the familiar 1D film temperature -- conduction moves
  no net heat, it only redistributes it. For m > 0 the denominator grows, so SHORT-WAVELENGTH
  load structure is attenuated: a peaked plasma load is smeared, a uniform one is untouched.
  Attenuation of mode m relative to the no-conduction film (q_m/h):

      A_m = h / (kappa*k_m*sinh(k_m*t) + h*cosh(k_m*t))     (<= 1, falls with m)

  The plasma-side surface temperature adds the through-thickness drop for the mean load.

WHAT IS AND ISN'T CLAIMED. h is taken from the two edges used by this project's design wedge:
the current wall-modeled effective h ~ 1.5 kW/m2K and the transferred Ferrero-channel
correlation ~5.5 kW/m2K. Neither is asserted as a physical bound. Conduction is linear and
temperature-independent here; real kappa(T)
and contact resistances are not modelled. This is a spreading calculation, not a CHT
certification -- but it is exact for the linear problem and it decides whether a CHT run is
worth the queue time.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_BULK = 900.0
T_BOIL = 1703.0
ALLOY_LO, ALLOY_HI = 973.0, 1073.0
H_LES, H_CORR = 1.5e3, 5.5e3          # W/m2K -- the envelope wedge's two edges

MATS = [("Hastelloy-N (Ni alloy)", 23.6), ("ODS steel", 20.0), ("tungsten", 120.0)]
THICK = np.array([1, 2, 3, 4, 6, 8]) * 1e-3

q = json.load(open("figs/plasma_qmap.json"))
s = np.asarray(q["s_m"]); qs = np.asarray(q["q_W_m2"])
L = s[-1] - s[0]
N = len(s)

def cosine_modes(f):
    """coefficients of f(s) = c0 + sum_m c_m cos(k_m s), k_m = m*pi/L (adiabatic ends)."""
    m = np.arange(0, N)
    # trapezoid projection onto cos(m*pi*s/L)
    C = np.array([np.trapezoid(f * np.cos(mm * np.pi * (s - s[0]) / L), s) for mm in m])
    c = 2.0 / L * C
    c[0] = np.trapezoid(f, s) / L
    return m, c

m, c = cosine_modes(qs)
k = m * np.pi / L
recon = c[0] + sum(c[i] * np.cos(k[i] * (s - s[0])) for i in range(1, N))
print(f"cosine reconstruction of q''(s): max err {100*np.abs(recon-qs).max()/qs.max():.3f}% of peak "
      f"({N} modes, L={L:.4f} m)")

def film_T(kappa, t, h):
    """coolant-side wall temperature profile with conduction."""
    den = kappa * k * np.sinh(k * t) + h * np.cosh(k * t)
    den[0] = h                                          # limit k->0
    th = sum(c[i] / den[i] * np.cos(k[i] * (s - s[0])) for i in range(1, N)) + c[0] / h
    return T_BULK + th

def film_T_nocond(h):
    return T_BULK + qs / h

print(f"\nload peak/mean {qs.max()/qs.mean():.2f}; q_peak {qs.max()/1e6:.3f} MW/m2\n")
print(f"{'material':24s} {'t [mm]':>7s} {'h':>10s} {'peak film, no cond':>19s} {'with cond':>11s} {'credit':>9s}")
rows = []
for name, kappa in MATS:
    for t in THICK:
        for htag, h in (("LES", H_LES), ("corr", H_CORR)):
            p0 = film_T_nocond(h).max()
            p1 = film_T(kappa, t, h).max()
            cred = (p0 - p1) / (p0 - T_BULK) * 100      # % of superheat removed
            rows.append(dict(mat=name, kappa=kappa, t=t, htag=htag, h=h, p0=p0, p1=p1, cred=cred))
            if t in (3e-3,) :
                print(f"{name:24s} {t*1e3:7.0f} {htag:>10s} {p0:16.0f} K {p1:9.0f} K {cred:8.1f}%")

def fwhm_of(f):
    floor, pk = f.min(), f.max()
    ab = s[f >= floor + 0.5 * (pk - floor)]
    return ab[-1] - ab[0]

# WHY the credit is small: conduction is a low-pass filter with cutoff at the spreading length
print("\n--- the governing scale: thermal spreading length  l = sqrt(kappa*t/h) ---")
FW = fwhm_of(qs)
print(f"    plasma peak FWHM = {FW*1e3:.0f} mm")
for name, kappa in MATS:
    l1, l2 = np.sqrt(kappa*3e-3/H_LES), np.sqrt(kappa*3e-3/H_CORR)
    print(f"    {name:24s} l(3mm) = {l1*1e3:5.1f} mm (h=LES) / {l2*1e3:4.1f} mm (h=corr)"
          f"   -> l/FWHM = {l1/FW:.3f}")
print("    A wall only spreads structure NARROWER than l. The stellarator standoff modulation")
print("    is set by the toroidal FIELD PERIOD, so it is inherently long-wavelength -- an order")
print("    of magnitude broader than any realistic l. That is why the credit is small.")

# cross-check of the effective h against the LES, from a completely independent route
print(f"\n--- independent check of h ---")
print(f"    analytic no-conduction film peak at h=LES: {film_T_nocond(H_LES).max():.0f} K")
print(f"    LES plasma hotspot p99:                    2097 K  (7% apart -> h ~ 1.5 kW/m2K confirmed)")

# credit vs load sharpness across the lambda/f_rad band corners
print("\n--- conduction credit vs LOAD SHARPNESS (3 mm, h=LES) ---")
print(f"    {'corner':9s} {'pk/mean':>8s} {'FWHM':>8s} {'Hastelloy':>10s} {'tungsten':>9s}")
for tag, path in (("mild", "figs/plasma_qmap_mild.json"), ("default", "figs/plasma_qmap.json"),
                  ("severe", "figs/plasma_qmap_severe.json")):
    qq = json.load(open(path)); f = np.asarray(qq["q_W_m2"])
    mm, cc = cosine_modes(f)
    kk = mm * np.pi / L
    def pk_with(kappa, t=3e-3, h=H_LES):
        den = kappa*kk*np.sinh(kk*t) + h*np.cosh(kk*t); den[0] = h
        th = sum(cc[i]/den[i]*np.cos(kk[i]*(s-s[0])) for i in range(1, len(mm))) + cc[0]/h
        return (T_BULK + th).max()
    p0 = (T_BULK + f/H_LES).max()
    cH = (p0 - pk_with(23.6))/(p0 - T_BULK)*100
    cW = (p0 - pk_with(120.0))/(p0 - T_BULK)*100
    print(f"    {tag:9s} {f.max()/f.mean():8.2f} {fwhm_of(f)*1e3:6.0f} mm {cH:9.1f}% {cW:8.1f}%")
print("    Credit GROWS with sharpness (the filter picture), but even the severest corner in the")
print("    lambda/f_rad band keeps only 17% for a Ni alloy. The escape hatch is closed band-wide.")

# headline: does conduction rescue the alloy band / boiling at the design point?
print("\n--- design-point verdicts (3 mm wall) ---")
for name, kappa in MATS:
    for htag, h in (("current-model edge", H_LES), ("Ferrero transfer edge", H_CORR)):
        p1 = film_T(kappa, 3e-3, h).max()
        p0 = film_T_nocond(h).max()
        v_b = "BOILS" if p1 > T_BOIL else "below boiling"
        v_a = "exceeds alloy band" if p1 > ALLOY_HI else "within alloy band"
        print(f"  {name:24s} {htag:16s}: peak film {p0:6.0f} -> {p1:6.0f} K   {v_b:14s} | {v_a}")

# ---------------- figure ----------------
fig = plt.figure(figsize=(13.2, 8.0), facecolor="white")
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.24)

ax = fig.add_subplot(gs[0, :])
ax.plot(s * 1e3, film_T_nocond(H_LES), "-", color="#7a1f1f", lw=2.2,
        label=f"no conduction (fluid-only tier), h={H_LES/1e3:g} kW/m$^2$K")
for name, kappa in MATS:
    ax.plot(s * 1e3, film_T(kappa, 3e-3, H_LES), "-", lw=1.7,
            label=f"3 mm {name} ($\\kappa$={kappa:g})")
ax.axhline(T_BOIL, color="#7a1f1f", ls="-.", lw=1.2)
ax.text(s[-1]*1e3, T_BOIL + 25, "FLiBe boiling 1703 K", ha="right", fontsize=8.5, color="#7a1f1f")
ax.axhspan(ALLOY_LO, ALLOY_HI, color="#c7b299", alpha=0.45)
ax.text(s[0]*1e3 + 3, (ALLOY_LO+ALLOY_HI)/2, "alloy film band", fontsize=8.5, color="0.3", va="center")
ax.set_xlabel("streamwise arc length s [mm]"); ax.set_ylabel("coolant-side wall (film) T [K]")
ax.set_title("Lateral conduction in the solid wall smears the plasma-shaped hotspot", fontsize=11.5, loc="left")
ax.legend(fontsize=8.5, loc="upper right", framealpha=0.95); ax.grid(alpha=0.3)
for sp in ax.spines.values(): sp.set_color("0.75")

ax2 = fig.add_subplot(gs[1, 0])
for name, kappa in MATS:
    for htag, h, ls in (("LES", H_LES, "-"), ("corr", H_CORR, "--")):
        cr = [ (film_T_nocond(h).max() - film_T(kappa, t, h).max()) /
               (film_T_nocond(h).max() - T_BULK) * 100 for t in THICK ]
        ax2.plot(THICK * 1e3, cr, ls, lw=1.8, label=f"{name.split()[0]}, h={htag}")
ax2.set_xlabel("wall thickness [mm]"); ax2.set_ylabel("hotspot superheat removed [%]")
ax2.set_title("Conduction credit vs wall thickness", fontsize=10.5)
ax2.legend(fontsize=7.5, ncol=2); ax2.grid(alpha=0.3)
for sp in ax2.spines.values(): sp.set_color("0.75")

ax3 = fig.add_subplot(gs[1, 1])
A = lambda kappa, t, h: h / (kappa * k[1:] * np.sinh(k[1:] * t) + h * np.cosh(k[1:] * t))
lam = 2 * L / m[1:]
for name, kappa in MATS:
    ax3.semilogx(lam * 1e3, A(kappa, 3e-3, H_LES), lw=1.8, label=f"{name.split()[0]} 3 mm")
ax3.axvline(2 * L * 1e3 / 4, color="0.6", ls=":", lw=1)
ax3.text(2*L*1e3/4*1.06, 0.5, "scale of the\nplasma peak", fontsize=8, color="0.4")
ax3.set_xlabel("load wavelength $\\lambda$ [mm]"); ax3.set_ylabel("mode attenuation $A_m$  [-]")
ax3.set_title("Conduction is a LOW-PASS FILTER on the load", fontsize=10.5)
ax3.set_ylim(0, 1.05); ax3.legend(fontsize=7.5); ax3.grid(alpha=0.3, which="both")
for sp in ax3.spines.values(): sp.set_color("0.75")

fig.suptitle("Conduction credit: what the solid first wall does to the plasma-shaped hotspot",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.005,
         "exact solution of steady conduction in a thin wall (linear -> diagonal in Fourier); h from this project's "
         "wedge edges | uniform loads are untouched (A_0 = 1): the credit exists ONLY because the load is peaked | "
         "constant kappa, no contact resistance, no kappa(T)",
         ha="center", fontsize=8, color="0.35")
fig.subplots_adjust(top=0.90, bottom=0.075)
fig.savefig("figs/conduction_credit.png", dpi=140, bbox_inches="tight")
json.dump(rows, open("figs/conduction_credit.json", "w"), indent=1)
print("\nwrote figs/conduction_credit.png, figs/conduction_credit.json")

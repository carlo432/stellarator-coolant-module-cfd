#!/usr/bin/env python3
"""PLASMA LOAD MODEL, SECOND GENERATION — a physically defined first-wall load.

WHAT WAS WRONG WITH v1 (plasma_wall_qmap.py)
  v1: q(s) = QWALL * (f_rad + (1-f_rad) * exp(-(d(s)-d_min)/lambda) / mean(...))
      - the exponential term was labelled "Eich-style convective". Eich's lambda_q is the
        DIVERTOR parallel heat-flux width (millimetres). A main-chamber first wall does not
        receive parallel flux: measured here from the DESC W7-X equilibrium, the field-line
        incidence angle on a cylindrical wall is 3.4 deg on average (|B.Rhat|/|B| <= 0.118,
        and exactly 0 at the outboard-most point, where the flux-surface normal IS R-hat).
        Parallel deposition is suppressed by sin(alpha) ~ 0.06. So the plan task "acquire Gao
        et al. lambda_q" was chasing the WRONG QUANTITY.
      - the radiative part was a FLAT FLOOR (f_rad, uniform in s). But radiation is a
        line-of-sight quantity: a wall element that sits 3 cm from the plasma sees a much
        larger solid angle than one 27 cm away. The radiative load must modulate too.
      - the volumetric NEUTRON source was uniform in s, while the surface load modulated.
        Neutrons travel straight and are unaffected by B, so the same standoff variation that
        modulates the surface load MUST modulate the neutron wall load. v1 was internally
        inconsistent.

WHAT v2 DOES
  A single line-of-sight emission integral over the plasma volume, evaluated at each wall
  point x_w on the cylinder R = R_wall, Z = 0, toroidal angle phi_w:

      q(phi_w) = Integral  eps(rho) * max(cos theta, 0) / (4 pi r^2)  dV

  with r = |x - x_w|, cos theta = (x - x_w).n_hat / r, n_hat the inward wall normal, and
  dV = sqrt(g) drho dtheta dzeta from DESC. Three emissivities, three physical channels:

    RADIATION   eps ~ exp(-(1-rho)/w)          mantle line radiation (edge-peaked)
    NEUTRONS    eps ~ (1-rho^2)^3              fusion source (core-peaked, n^2<sigma v>)
    (bremsstrahlung eps ~ (1-rho^2)^1.5 carried as a radiation sensitivity)

  plus the CROSS-FIELD channel, which is what the exponential legitimately models:

    q_perp(s) ~ exp(-(d(s) - d_min)/lambda_far)   filamentary/blob transport to the wall
                lambda_far = 3 cm [1-5], a MAIN-CHAMBER far-SOL width -- correctly cited as
                such, not as a divertor lambda_q

  and the small parallel remnant q_par ~ sin(alpha(zeta)) * exp(-d/lambda_far), computed
  from the real |B.Rhat|/|B| rather than assumed away.

OUTPUTS: figs/plasma_qmap_v2.{png,json} — the surface load q''(s) with its three channels
separated, AND the neutron wall-load modulation NWL(s) that v1 wrongly held uniform.

FOLLOW-UP AUDIT STATUS (2026-07-09): diagnostic only. Field-line incidence, cross-field
transport labeling, and neutron wall-load modulation survive. Radiative/total peaking
magnitudes and any v2 hotspot extrapolation are withdrawn until the line-of-sight radiation
integral gets a near-field/far-field or analytic near-field treatment.

HONESTY: still a TRANSPLANT of a published configuration onto an idealized duct; emissivity
profiles are shape assumptions (bands carried), not a transport solution; the wall is a
cylinder, not real CAD. What v2 removes is the mislabelled physics and the internal
inconsistency, not the transplant.
"""
from __future__ import annotations
import json, warnings
import numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QWALL = 0.5e6
CLEARANCE = 0.03
LAM_FAR = 0.03                 # main-chamber far-SOL width [m], band 1-5 cm
F_RAD = 0.5                    # radiated fraction of the exhaust power, band 0.3-0.7
S_DUCT = np.radians(95.0) * 0.18

NR, NT, NZ_FULL = 40, 24, 96   # plasma volume grid (full torus for the LOS integral)
RHO_CLUSTER = 3.0              # radial clustering toward rho=1. REQUIRED: the mantle emitter sits
                               # closest to the wall and the 1/r^2 kernel is near-field sensitive.
                               # A uniform grid gives peak/mean 4.47 at NR=10, 3.43 at NR=40 -- still
                               # falling. Edge-clustered: 3.145 (NR=20) -> 3.133 (NR=40), converged.
NW = 192                       # wall points per field period (periodic grid; means are
                               # exact simple averages on it)

from desc.examples import get
from desc.grid import Grid

eq = get("W7-X"); NFP = eq.NFP

# ---- plasma volume, cartesian, with Jacobian ----
u = np.linspace(0, 1, NR)
rho = 0.05 + 0.95 * (1 - (1 - u) ** RHO_CLUSTER)     # clustered toward the edge
drho = np.gradient(rho)
th = np.linspace(0, 2*np.pi, NT, endpoint=False)
ze = np.linspace(0, 2*np.pi, NZ_FULL, endpoint=False)
RH, TH, ZE = np.meshgrid(rho, th, ze, indexing="ij")
DR, _, _ = np.meshgrid(drho, th, ze, indexing="ij")
nodes = np.c_[RH.ravel(), TH.ravel(), ZE.ravel()]
o = eq.compute(["R", "Z", "sqrt(g)"], grid=Grid(nodes, sort=False))
Rp = np.asarray(o["R"]); Zp = np.asarray(o["Z"]); Jg = np.abs(np.asarray(o["sqrt(g)"]))
PHI = ZE.ravel()
Xp = np.c_[Rp*np.cos(PHI), Rp*np.sin(PHI), Zp]
dV = Jg * DR.ravel() * (2*np.pi/NT) * (2*np.pi/NZ_FULL)
rr = RH.ravel()

EPS = {
    "radiation (mantle)":  np.exp(-(1.0 - rr) / 0.10),
    "radiation (brems)":   np.maximum(1.0 - rr**2, 0.0) ** 1.5,
    "neutrons (fusion)":   np.maximum(1.0 - rr**2, 0.0) ** 3.0,
}

# ---- wall geometry: cylinder at R_wall, Z=0 ----
zw = np.linspace(0, 2*np.pi/NFP, NW, endpoint=False)   # PERIODIC. The closest-approach point
# sits exactly on zeta=0; an endpoint-INCLUSIVE grid duplicates it, inflating every channel mean
# and deflating peak/mean with a slow 1/N error (5.10 -> 5.56 for the transport channel as
# NW goes 49 -> 1537). With endpoint=False the periodic average is converged at NW=49 (5.5766).
gth = np.linspace(0, 2*np.pi, 257)
G_TH, G_ZE = np.meshgrid(gth, zw, indexing="ij")
gn = np.c_[np.ones(G_TH.size), G_TH.ravel(), G_ZE.ravel()]
go = eq.compute(["R", "Z", "B", "|B|"], grid=Grid(gn, sort=False))
Rl = np.asarray(go["R"]).reshape(257, NW)
Bv = np.asarray(go["B"]).reshape(257, NW, 3)
Bm = np.asarray(go["|B|"]).reshape(257, NW)
iout = np.argmax(Rl, axis=0); jj = np.arange(NW)
R_out = Rl[iout, jj]
sin_alpha = np.abs(Bv[iout, jj, 0]) / Bm[iout, jj]      # |B.Rhat|/|B| at the outboard LCFS
R_WALL = R_out.max() + CLEARANCE
d = R_WALL - R_out                                       # standoff [m]

def los(eps):
    """line-of-sight emission integral at each wall point (unnormalised)."""
    out = np.empty(NW)
    for i, phi in enumerate(zw):
        xw = np.array([R_WALL*np.cos(phi), R_WALL*np.sin(phi), 0.0])
        nhat = np.array([-np.cos(phi), -np.sin(phi), 0.0])    # inward
        v = Xp - xw
        r2 = np.einsum("ij,ij->i", v, v)
        r = np.sqrt(r2)
        cosw = np.maximum((v @ nhat) / r, 0.0)
        out[i] = np.sum(eps * dV * cosw / (4*np.pi*r2))
    return out

chan = {k: los(v) for k, v in EPS.items()}
print(f"NFP={NFP}  R_wall={R_WALL:.4f} m  standoff d: {d.min()*100:.1f}-{d.max()*100:.1f} cm")
print(f"field-line incidence sin(alpha): {sin_alpha.min():.4f}-{sin_alpha.max():.4f} "
      f"(mean {np.degrees(np.arcsin(sin_alpha)).mean():.2f} deg) -> parallel deposition suppressed\n")

# ---- assemble the SURFACE load ----
g_perp = np.exp(-(d - d.min())/LAM_FAR)                  # cross-field (blob) transport
q_par  = sin_alpha * g_perp                              # small parallel remnant
q_rad  = chan["radiation (mantle)"]
def nrm(a): return a / a.mean()
# exhaust power splits f_rad radiated / (1-f_rad) transported; the transported part reaches the
# wall by cross-field transport, with a small parallel remnant riding on the same profile
q_transport = nrm(g_perp) * (1 - sin_alpha.mean()) + nrm(q_par) * sin_alpha.mean()
q_surf = QWALL * (F_RAD * nrm(q_rad) + (1 - F_RAD) * nrm(q_transport))

nwl = nrm(chan["neutrons (fusion)"])
print(f"{'channel':26s} {'peak/mean':>10s} {'min/mean':>9s}")
for k, v in (("radiation (mantle, LOS)", nrm(q_rad)), ("radiation (brems, LOS)", nrm(chan['radiation (brems)'])),
             ("cross-field exp(-d/lam)", nrm(g_perp)), ("parallel remnant", nrm(q_par)),
             ("NEUTRON wall load", nwl)):
    print(f"{k:26s} {v.max():10.3f} {v.min():9.3f}")
print(f"\nTOTAL surface load q'': peak/mean {q_surf.max()/q_surf.mean():.3f}  "
      f"(DIAGNOSTIC ONLY; radiative/total peaking withdrawn)")
print(f"mean = {q_surf.mean()/1e6:.4f} MW/m2 (target 0.5; v1 endpoint audit found 0.5026)")
print("WARNING: do not promote v2 total/radiative peaking or hotspot extrapolation before "
      "near-field/far-field radiation treatment.")
print(f"\n*** v1 INCONSISTENCY QUANTIFIED: the neutron wall load is NOT uniform -- it modulates")
print(f"    {nwl.min():.3f}..{nwl.max():.3f} x mean (a {100*(nwl.max()-nwl.min()):.0f}% swing) because neutrons")
print(f"    travel straight from a core-peaked source and the standoff swings {d.min()*100:.0f}-{d.max()*100:.0f} cm. ***")

s = np.linspace(0, S_DUCT, NW, endpoint=False)
json.dump(dict(source="DESC W7-X example equilibrium (v2 line-of-sight model)", NFP=int(NFP),
               R_wall_m=float(R_WALL), lambda_far_m=LAM_FAR, lambda_far_band_m=[0.01, 0.05],
               f_rad=F_RAD, f_rad_band=[0.3, 0.7],
               s_m=s.tolist(), d_m=d.tolist(), q_W_m2=q_surf.tolist(),
               nwl_rel=nwl.tolist(), sin_alpha=sin_alpha.tolist(),
               q_rad_rel=nrm(q_rad).tolist(), q_transport_rel=nrm(q_transport).tolist(),
               peak_over_mean=float(q_surf.max()/q_surf.mean()),
               notes=("radiative + neutron loads from a 3D line-of-sight emission integral over the DESC "
                      "plasma volume; cross-field far-SOL exponential correctly labelled (NOT Eich divertor "
                      "lambda_q -- field-line incidence on the wall is 3.4 deg, parallel deposition suppressed); "
                      "neutron wall load modulates, fixing v1's internal inconsistency")),
          open("figs/plasma_qmap_v2.json", "w"), indent=1)

fig, axes = plt.subplots(3, 1, figsize=(10.4, 9.6), facecolor="white", sharex=True,
                         gridspec_kw=dict(hspace=0.18))
a = axes[0]
a.plot(s*1e3, d*100, color="#2a5d8f", lw=2)
a.set_ylabel("standoff d [cm]"); a.grid(alpha=0.3)
a.set_title("Plasma load model v2: line-of-sight radiation + cross-field transport", fontsize=12, loc="left")
a2 = a.twinx(); a2.plot(s*1e3, np.degrees(np.arcsin(sin_alpha)), color="#8a6d3b", lw=1.4, ls="--")
a2.set_ylabel("field-line incidence [deg]", color="#8a6d3b", fontsize=9)
a2.tick_params(labelsize=8, colors="#8a6d3b")

a = axes[1]
a.plot(s*1e3, nrm(q_rad), lw=2, color="#7a1f1f", label="radiation (LOS, mantle)")
a.plot(s*1e3, nrm(chan["radiation (brems)"]), lw=1.2, ls=":", color="#7a1f1f", label="radiation (LOS, brems)")
a.plot(s*1e3, nrm(g_perp), lw=2, color="#1a6b52", label="cross-field  $e^{-d/\\lambda_{far}}$")
a.axhline(1.0, color="0.6", lw=1, ls=":")
a.plot(s*1e3, np.full_like(s, 1.0), lw=1.4, color="0.55", ls="--", label="v1's FLAT radiative floor")
a.set_ylabel("channel / its mean"); a.legend(fontsize=8, ncol=2); a.grid(alpha=0.3)

a = axes[2]
a.plot(s*1e3, q_surf/1e6, lw=2.4, color="#111", label=f"v2 total q'' (pk/mean {q_surf.max()/q_surf.mean():.2f})")
qv1 = json.load(open("figs/plasma_qmap.json"))
a.plot(np.asarray(qv1["s_m"])*1e3, np.asarray(qv1["q_W_m2"])/1e6, lw=1.6, ls="--", color="#7a1f1f",
       label=f"v1 (pk/mean {qv1['peak_over_mean']:.2f})")
a.plot(s*1e3, nwl*0.5, lw=1.8, color="#c07830", label="NEUTRON wall load (rel., x0.5) -- v1 held this UNIFORM")
a.set_xlabel("streamwise arc length s [mm]"); a.set_ylabel("q'' [MW/m$^2$]")
a.legend(fontsize=8.5); a.grid(alpha=0.3)
for ax in list(axes):
    for sp in ax.spines.values(): sp.set_color("0.75")
fig.text(0.5, 0.004,
         "LOS integral: q(x_w) = int eps(rho) cos(theta)/(4 pi r^2) dV over the DESC plasma volume | radiation and neutrons "
         "both modulate with standoff; only the cross-field channel needs lambda | field lines graze the wall (3.4 deg), so "
         "Eich divertor lambda_q is the WRONG quantity for a main-chamber load",
         ha="center", fontsize=7.6, color="0.35")
fig.subplots_adjust(bottom=0.075, top=0.955)
fig.savefig("figs/plasma_qmap_v2.png", dpi=140, bbox_inches="tight")
print("\nwrote figs/plasma_qmap_v2.{png,json}")

#!/usr/bin/env python3
"""Task 1.1/1.2 of STELLARATOR_PIPELINE_PLAN.md: OpenMC depth-resolved heating + tritium source.

Model: 1D FLiBe column, 14.1 MeV monodirectional planar source, reflective lateral faces.
The first `--channel-depth` metres represent OUR duct (radial depth, first wall at x=0);
the remaining thickness is backing FLiBe standing in for the blanket continuum behind the
channel (gives back-scatter; a pure 3 cm slab with vacuum behind would under-deposit).
Honesty ledger: 1D slab vs curved slice, no W/Inconel first-wall stack in front -- stated
in the plan; production fidelity is a Grace-era full-geometry OpenMC job.

Normalization: neutron wall load NWL [W/m2] incident at x=0 => source rate per area
S'' = NWL / E_src. Per-zone results:
  q'''_i  = H_i [eV/src] * EV_TO_J * S'' / dx    [W/m3]
  S_T,i   = T_i [t/src]  * S'' / (dx * N_A)      [mol/m3/s]
The 0.5 MW/m2 plasma SURFACE flux is a separate wall BC in the CFD case, NOT included here.

Outputs (for make_geom_pipeline_case.py):
  figs/openmc_depth_source.csv       per-bin profile over the whole column
  figs/openmc_depth_source.json      manifest: per-zone q''', S_T, TBRs, provenance
  figs/openmc_depth_source.png       profile plot with channel zones marked

Usage: openmc_depth_source.py [--nwl 1.0e6] [--particles 20000] [--batches 10] [--no-run]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from datetime import date
from pathlib import Path

import numpy as np

XS_PATH = os.environ.get("OPENMC_CROSS_SECTIONS", "/root/openmc_data/combined_neutron_photon_all.xml")
OPENMC_EXE = os.environ.get("OPENMC_EXE", "/root/openmc-src/build/bin/openmc")
EV_TO_J = 1.602176634e-19
N_A = 6.02214076e23
E_SRC_EV = 14.1e6

ap = argparse.ArgumentParser()
ap.add_argument("--case", default="openmc_pipeline_column")
ap.add_argument("--nwl", type=float, default=1.0e6, help="neutron wall load [W/m2] (plan default 1.0 MW/m2)")
ap.add_argument("--channel-depth", type=float, default=0.030, help="duct radial depth [m]")
ap.add_argument("--zones", type=int, default=6, help="depth zones inside the channel")
ap.add_argument("--slab-thickness", type=float, default=0.50, help="total FLiBe column incl. backing [m]")
ap.add_argument("--li6", type=float, default=0.90, help="Li-6 enrichment (ARC-like, matches C3)")
ap.add_argument("--particles", type=int, default=20000)
ap.add_argument("--batches", type=int, default=10)
ap.add_argument("--rhocp", type=float, default=4.65e6, help="rho*cp of FLiBe [J/m3K]")
ap.add_argument("--no-run", action="store_true", help="postprocess an existing statepoint only")
args = ap.parse_args()

case = Path(args.case)
case.mkdir(exist_ok=True)
dx = args.channel_depth / args.zones          # zone width; the whole column is binned at this pitch
nbins = int(round(args.slab_thickness / dx))
assert abs(nbins * dx - args.slab_thickness) < 1e-9, "slab thickness must be a multiple of zone width"

os.environ["OPENMC_CROSS_SECTIONS"] = XS_PATH
import openmc  # noqa: E402

if not args.no_run:
    flibe = openmc.Material(name="flibe_2LiF_BeF2_li6_enriched")
    flibe.set_density("g/cm3", 1.94)
    flibe.add_nuclide("Li6", 2.0 * args.li6, "ao")
    flibe.add_nuclide("Li7", 2.0 * (1.0 - args.li6), "ao")
    flibe.add_nuclide("Be9", 1.0, "ao")
    flibe.add_nuclide("F19", 4.0, "ao")
    openmc.Materials([flibe]).export_to_xml(case / "materials.xml")

    hw = 0.008 * 100  # cm; lateral extent is irrelevant (reflective), keep C3-like proportions
    L = args.slab_thickness * 100  # OpenMC works in cm
    x0 = openmc.XPlane(x0=0.0, boundary_type="vacuum")
    x1 = openmc.XPlane(x0=L, boundary_type="vacuum")
    ymin = openmc.YPlane(y0=-hw, boundary_type="reflective")
    ymax = openmc.YPlane(y0=hw, boundary_type="reflective")
    zmin = openmc.ZPlane(z0=-hw, boundary_type="reflective")
    zmax = openmc.ZPlane(z0=hw, boundary_type="reflective")
    slab = openmc.Cell(fill=flibe, region=+x0 & -x1 & +ymin & -ymax & +zmin & -zmax)
    openmc.Geometry([slab]).export_to_xml(case / "geometry.xml")

    st = openmc.Settings()
    st.run_mode = "fixed source"
    st.particles = args.particles
    st.batches = args.batches
    st.inactive = 0
    st.seed = 170705
    st.source = openmc.IndependentSource(
        space=openmc.stats.Box((1e-6, -hw, -hw), (2e-6, hw, hw)),
        angle=openmc.stats.Monodirectional(reference_uvw=(1.0, 0.0, 0.0)),
        energy=openmc.stats.Discrete([E_SRC_EV], [1.0]),
    )
    st.export_to_xml(case / "settings.xml")

    mesh = openmc.RegularMesh()
    mesh.dimension = (nbins, 1, 1)
    mesh.lower_left = (0.0, -hw, -hw)
    mesh.upper_right = (L, hw, hw)
    tally = openmc.Tally(name="depth_source")
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ["heating-local", "H3-production", "flux"]
    openmc.Tallies([tally]).export_to_xml(case / "tallies.xml")

    print(f"=== OpenMC run: {nbins} bins x {dx*1e3:.1f} mm, {args.particles}x{args.batches} particles ===")
    subprocess.run([OPENMC_EXE], cwd=case, check=True,
                   stdout=open(case / "openmc.log", "w"), stderr=subprocess.STDOUT)

sp_file = sorted(case.glob("statepoint.*.h5"))[-1]
sp = openmc.StatePoint(sp_file)
t = sp.get_tally(name="depth_source")
heat = t.get_values(scores=["heating-local"]).reshape(nbins)      # eV / source
heat_sd = t.get_values(scores=["heating-local"], value="std_dev").reshape(nbins)
trit = t.get_values(scores=["(n,Xt)"]).reshape(nbins)             # tritons / source ("H3-production" canonicalized)
trit_sd = t.get_values(scores=["(n,Xt)"], value="std_dev").reshape(nbins)

src_per_area = args.nwl / (E_SRC_EV * EV_TO_J)                    # n / m2 / s
qvol = heat * EV_TO_J * src_per_area / dx                          # W/m3
qvol_sd = heat_sd * EV_TO_J * src_per_area / dx
strit = trit * src_per_area / (dx * N_A)                           # mol/m3/s
strit_sd = trit_sd * src_per_area / (dx * N_A)

nz = args.zones
tbr_column = float(trit.sum())
tbr_channel = float(trit[:nz].sum())
edep_column = float(heat.sum()) / E_SRC_EV                         # fraction of source energy deposited
edep_channel = float(heat[:nz].sum()) / E_SRC_EV

# --- duct zone geometry (matches gen_curved_stellarator_slice_c2.py / geom_heat parameters) ---
GEOM = dict(axis_point=[0.0, 0.0, 0.18], axis_dir=[0.0, 1.0, 0.0],
            r_first_wall=0.165, r_back_wall=0.195, arc_rad=np.radians(95.0), height=0.016)
zones = []
for i in range(nz):
    r_in = GEOM["r_first_wall"] + i * dx
    r_out = r_in + dx
    vol = 0.5 * GEOM["arc_rad"] * (r_out**2 - r_in**2) * GEOM["height"]   # analytic annular band
    zones.append(dict(
        zone=i, depth0_m=i * dx, depth1_m=(i + 1) * dx, r_inner_m=r_in, r_outer_m=r_out,
        qvol_W_m3=float(qvol[i]), qvol_relerr=float(qvol_sd[i] / qvol[i]) if qvol[i] > 0 else None,
        S_T_K_per_s=float(qvol[i] / args.rhocp),
        strit_mol_m3_s=float(strit[i]), strit_relerr=float(strit_sd[i] / strit[i]) if strit[i] > 0 else None,
        volume_analytic_m3=float(vol), power_W=float(qvol[i] * vol),
        trit_gen_mol_s=float(strit[i] * vol),
    ))

monotone = all(zones[i]["qvol_W_m3"] >= zones[i + 1]["qvol_W_m3"] for i in range(nz - 1))
total_power = sum(z["power_W"] for z in zones)
total_trit = sum(z["trit_gen_mol_s"] for z in zones)

figs = Path("figs"); figs.mkdir(exist_ok=True)
with open(figs / "openmc_depth_source.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["x0_m", "x1_m", "x_mid_m", "qvol_W_m3", "qvol_sd_W_m3", "strit_mol_m3_s", "strit_sd_mol_m3_s", "in_channel"])
    for i in range(nbins):
        w.writerow([i * dx, (i + 1) * dx, (i + 0.5) * dx, qvol[i], qvol_sd[i], strit[i], strit_sd[i], int(i < nz)])

manifest = dict(
    date=str(date.today()), script="openmc_depth_source.py", statepoint=str(sp_file),
    cross_sections=XS_PATH, openmc_version=openmc.__version__,
    particles=args.particles, batches=args.batches, li6_fraction=args.li6,
    neutron_wall_load_W_m2=args.nwl, source_energy_eV=E_SRC_EV, src_per_area_n_m2_s=src_per_area,
    slab_thickness_m=args.slab_thickness, channel_depth_m=args.channel_depth, zone_width_m=dx,
    rhocp_J_m3K=args.rhocp, geometry=GEOM, zones=zones,
    tbr_column=tbr_column, tbr_channel_slice=tbr_channel,
    edep_fraction_column=edep_column, edep_fraction_channel=edep_channel,
    qvol_monotone_decreasing_in_channel=monotone,
    total_channel_power_W=total_power, total_channel_tritium_mol_s=total_trit,
    notes="1D FLiBe column, no solid FW stack, neutron-only transport with local heating (KERMA); "
          "backing FLiBe behind channel provides back-scatter. Surface 0.5 MW/m2 flux is separate.",
)
(figs / "openmc_depth_source.json").write_text(json.dumps(manifest, indent=2))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
xm = (np.arange(nbins) + 0.5) * dx
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5), facecolor="white")
for ax, y, ysd, lab, ttl in [
    (axs[0], qvol, qvol_sd, "q''' [W/m3]", "Volumetric heating"),
    (axs[1], strit, strit_sd, "S_T [mol/m3/s]", "Tritium generation"),
]:
    ax.errorbar(xm, y, yerr=ysd, fmt=".-", ms=3, lw=0.8, ecolor="0.7")
    ax.axvspan(0, args.channel_depth, color="tab:orange", alpha=0.15, label="our channel (6 zones)")
    ax.set_yscale("log"); ax.set_xlabel("depth from first wall [m]"); ax.set_ylabel(lab)
    ax.set_title(ttl); ax.legend(); ax.grid(alpha=0.3)
fig.suptitle(f"OpenMC FLiBe depth source -- NWL {args.nwl/1e6:.2f} MW/m2, Li-6 {args.li6:.0%}, "
             f"{args.particles*args.batches:,} particles", fontsize=12)
fig.tight_layout()
fig.savefig(figs / "openmc_depth_source.png", dpi=130)

print(f"deposited energy fraction: column {edep_column:.4f}, channel(3cm) {edep_channel:.4f}")
print(f"TBR: column {tbr_column:.4f}, channel slice {tbr_channel:.4f}")
print(f"channel q''' by zone [W/m3]: " + " ".join(f"{z['qvol_W_m3']:.4g}" for z in zones))
print(f"  monotone decreasing: {monotone}")
print(f"channel S_T by zone [K/s]:  " + " ".join(f"{z['S_T_K_per_s']:.4g}" for z in zones))
print(f"total channel power {total_power:.2f} W, tritium {total_trit:.4g} mol/s")
print("wrote figs/openmc_depth_source.{csv,json,png}")

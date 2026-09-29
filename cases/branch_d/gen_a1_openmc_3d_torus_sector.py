#!/usr/bin/env python3
"""Build, run, and postprocess the A1 native-OpenMC 3D toroidal blanket model.

The CSG model is a full outboard annular blanket ring. A non-axisymmetric DESC
W7-X fusion-source mesh supplies the five-field-period modulation. Results are
folded to one field period and compared with the existing v2 line-of-sight
neutron-wall-load proxy and the earlier 1D OpenMC column.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/revengence-mpl-a1")

import numpy as np


WORK = Path(__file__).resolve().parent
REPO = WORK.parents[1]
DEFAULT_CASE = WORK / "a1_openmc_3d_torus_sector"
DEFAULT_RESULTS = REPO / "results" / "branch_d" / "a1_openmc_3d_torus_sector"
V2_PATH = WORK / "figs" / "plasma_qmap_v2.json"
PROXY_PATH = WORK / "figs" / "openmc_depth_source.json"
XS_PATH = os.environ.get("OPENMC_CROSS_SECTIONS", "/root/openmc_data/combined_neutron_photon_all.xml")
OPENMC_EXE = os.environ.get("OPENMC_EXE", "/root/openmc-src/build/bin/openmc")
EV_TO_J = 1.602176634e-19
E_SOURCE_EV = 14.1e6
NFP = 5
REFERENCE_NWL_W_M2 = 1.0e6


def import_openmc():
    os.environ["OPENMC_CROSS_SECTIONS"] = XS_PATH
    import openmc  # noqa: PLC0415

    return openmc


def make_steel_proxy(openmc, name: str):
    material = openmc.Material(name=name)
    material.set_density("g/cm3", 7.8)
    material.add_nuclide("Fe56", 0.88, "ao")
    material.add_nuclide("Cr52", 0.09, "ao")
    material.add_nuclide("Ni58", 0.03, "ao")
    return material


def make_flibe(openmc, li6_fraction: float):
    material = openmc.Material(name="flibe_2LiF_BeF2_li6_enriched")
    material.set_density("g/cm3", 1.94)
    material.add_nuclide("Li6", 2.0 * li6_fraction, "ao")
    material.add_nuclide("Li7", 2.0 * (1.0 - li6_fraction), "ao")
    material.add_nuclide("Be9", 1.0, "ao")
    material.add_nuclide("F19", 4.0, "ao")
    return material


def desc_fusion_source(
    *,
    n_rho: int,
    n_theta: int,
    n_zeta: int,
    n_r_mesh: int,
    n_phi_mesh: int,
    n_z_mesh: int,
    wall_radius_m: float,
) -> dict[str, np.ndarray | float | str]:
    """Histogram a core-peaked DESC fusion source onto a cylindrical mesh."""
    from desc.examples import get
    from desc.grid import Grid

    equilibrium = get("W7-X")
    rho = np.linspace(0.025, 0.995, n_rho)
    drho = np.gradient(rho)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)
    rh, th, ze = np.meshgrid(rho, theta, zeta, indexing="ij")
    dr, _, _ = np.meshgrid(drho, theta, zeta, indexing="ij")
    nodes = np.c_[rh.ravel(), th.ravel(), ze.ravel()]
    output = equilibrium.compute(["R", "Z", "sqrt(g)"], grid=Grid(nodes, sort=False))
    radius_m = np.asarray(output["R"], dtype=float)
    z_m = np.asarray(output["Z"], dtype=float)
    jacobian = np.abs(np.asarray(output["sqrt(g)"], dtype=float))
    phi = ze.ravel()
    rho_flat = rh.ravel()
    d_volume = jacobian * dr.ravel() * (2.0 * np.pi / n_theta) * (2.0 * np.pi / n_zeta)
    emissivity = np.maximum(1.0 - rho_flat**2, 0.0) ** 3.0
    weights = emissivity * d_volume

    radial_margin_m = 0.03
    z_margin_m = 0.03
    r_min_m = float(radius_m.min() - radial_margin_m)
    r_max_m = float(min(radius_m.max() + radial_margin_m, wall_radius_m - 0.005))
    z_bound_m = float(np.max(np.abs(z_m)) + z_margin_m)
    r_edges_cm = np.linspace(r_min_m * 100.0, r_max_m * 100.0, n_r_mesh + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi_mesh + 1)
    z_edges_cm = np.linspace(-z_bound_m * 100.0, z_bound_m * 100.0, n_z_mesh + 1)
    samples = np.c_[radius_m * 100.0, np.mod(phi, 2.0 * np.pi), z_m * 100.0]
    histogram, _ = np.histogramdd(
        samples,
        bins=(r_edges_cm, phi_edges, z_edges_cm),
        weights=weights,
    )
    captured = float(histogram.sum() / weights.sum())
    if captured < 0.999:
        raise ValueError(f"DESC source mesh captured only {captured:.6f} of source weight")
    if np.count_nonzero(histogram) == 0:
        raise ValueError("DESC source histogram is empty")

    return {
        "source": "DESC W7-X example equilibrium; emissivity proportional to (1-rho^2)^3",
        "r_edges_cm": r_edges_cm,
        "phi_edges_rad": phi_edges,
        "z_edges_cm": z_edges_cm,
        "strengths": histogram,
        "capture_fraction": captured,
        "nonzero_bins": int(np.count_nonzero(histogram)),
        "sample_points": int(len(weights)),
        "plasma_R_min_m": float(radius_m.min()),
        "plasma_R_max_m": float(radius_m.max()),
        "plasma_abs_Z_max_m": float(np.max(np.abs(z_m))),
    }


def cylindrical_cell_volumes(mesh) -> np.ndarray:
    r0 = np.asarray(mesh.r_grid[:-1], dtype=float)
    r1 = np.asarray(mesh.r_grid[1:], dtype=float)
    dphi = np.diff(np.asarray(mesh.phi_grid, dtype=float))
    dz = np.diff(np.asarray(mesh.z_grid, dtype=float))
    radial_area = 0.5 * (r1**2 - r0**2)
    return radial_area[:, None, None] * dphi[None, :, None] * dz[None, None, :]


def build_model(
    case: Path,
    *,
    source_data: dict,
    wall_radius_m: float,
    first_wall_m: float,
    blanket_m: float,
    back_wall_m: float,
    li6_fraction: float,
    n_phi_tally: int,
    n_blanket_depth: int,
    particles: int,
    batches: int,
    seed: int,
) -> dict:
    openmc = import_openmc()
    steel_front = make_steel_proxy(openmc, "first_wall_steel_proxy")
    flibe = make_flibe(openmc, li6_fraction)
    steel_back = make_steel_proxy(openmc, "back_wall_steel_proxy")
    openmc.Materials([steel_front, flibe, steel_back]).export_to_xml(case / "materials.xml")

    r_inner_cm = float(source_data["r_edges_cm"][0] - 2.0)
    r_wall_cm = wall_radius_m * 100.0
    r_front_outer_cm = (wall_radius_m + first_wall_m) * 100.0
    r_blanket_outer_cm = (wall_radius_m + first_wall_m + blanket_m) * 100.0
    r_back_outer_cm = (wall_radius_m + first_wall_m + blanket_m + back_wall_m) * 100.0
    z_half_cm = float(max(abs(source_data["z_edges_cm"][0]), abs(source_data["z_edges_cm"][-1])) + 2.0)

    radial_inner = openmc.ZCylinder(r=r_inner_cm, boundary_type="vacuum", name="inner_radial_vacuum")
    wall_inner = openmc.ZCylinder(r=r_wall_cm, name="plasma_facing_wall")
    wall_outer = openmc.ZCylinder(r=r_front_outer_cm, name="first_wall_outer")
    blanket_outer = openmc.ZCylinder(r=r_blanket_outer_cm, name="blanket_outer")
    back_outer = openmc.ZCylinder(r=r_back_outer_cm, boundary_type="vacuum", name="outer_radial_vacuum")
    z_min = openmc.ZPlane(z0=-z_half_cm, boundary_type="vacuum", name="lower_vacuum")
    z_max = openmc.ZPlane(z0=z_half_cm, boundary_type="vacuum", name="upper_vacuum")
    z_region = +z_min & -z_max

    cavity = openmc.Cell(name="plasma_and_standoff_void", region=+radial_inner & -wall_inner & z_region)
    first_wall = openmc.Cell(name="first_wall", fill=steel_front, region=+wall_inner & -wall_outer & z_region)
    blanket = openmc.Cell(name="flibe_blanket", fill=flibe, region=+wall_outer & -blanket_outer & z_region)
    back_wall = openmc.Cell(name="back_wall", fill=steel_back, region=+blanket_outer & -back_outer & z_region)
    cells = [cavity, first_wall, blanket, back_wall]
    openmc.Geometry(cells).export_to_xml(case / "geometry.xml")

    source_mesh = openmc.CylindricalMesh(
        r_grid=source_data["r_edges_cm"],
        phi_grid=source_data["phi_edges_rad"],
        z_grid=source_data["z_edges_cm"],
        name="desc_fusion_source_mesh",
    )
    source = openmc.IndependentSource(
        space=openmc.stats.MeshSpatial(
            source_mesh,
            strengths=np.asarray(source_data["strengths"]).ravel(order="F"),
            volume_normalized=False,
        ),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([E_SOURCE_EV], [1.0]),
        constraints={"domains": [cavity], "rejection_strategy": "resample"},
    )
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = 0
    settings.seed = seed
    settings.source = source
    settings.photon_transport = False
    settings.export_to_xml(case / "settings.xml")

    phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi_tally + 1)
    wall_mesh = openmc.CylindricalMesh(
        r_grid=[r_wall_cm, r_front_outer_cm],
        phi_grid=phi_edges,
        z_grid=[-z_half_cm, z_half_cm],
        name="plasma_facing_wall_current_mesh",
    )
    wall_surface_filter = openmc.MeshSurfaceFilter(wall_mesh)
    wall_current = openmc.Tally(name="wall_current_by_toroidal_bin")
    wall_current.filters = [wall_surface_filter]
    wall_current.scores = ["current"]

    wall_uncollided = openmc.Tally(name="uncollided_wall_current_by_toroidal_bin")
    wall_uncollided.filters = [
        wall_surface_filter,
        openmc.CollisionFilter([0]),
    ]
    wall_uncollided.scores = ["current"]

    blanket_depth_edges = np.linspace(r_front_outer_cm, r_blanket_outer_cm, n_blanket_depth + 1)
    heat_r_edges = np.r_[r_wall_cm, blanket_depth_edges, r_back_outer_cm]
    heat_mesh = openmc.CylindricalMesh(
        r_grid=heat_r_edges,
        phi_grid=phi_edges,
        z_grid=[-z_half_cm, z_half_cm],
        name="blanket_heating_mesh",
    )
    mesh_tally = openmc.Tally(name="heating_and_tritium_by_toroidal_depth_bin")
    mesh_tally.filters = [openmc.MeshFilter(heat_mesh)]
    mesh_tally.scores = ["heating-local", "H3-production", "flux"]

    cell_tally = openmc.Tally(name="heating_and_tritium_by_region")
    cell_tally.filters = [openmc.CellFilter([first_wall, blanket, back_wall])]
    cell_tally.scores = ["heating-local", "H3-production", "flux"]
    openmc.Tallies([wall_current, wall_uncollided, mesh_tally, cell_tally]).export_to_xml(
        case / "tallies.xml"
    )

    return {
        "geometry_kind": "full 360-degree outboard annular blanket ring; five field periods folded in postprocessing",
        "length_unit": "OpenMC CSG dimensions are cm; manifest values below are reported in m",
        "wall_radius_m": wall_radius_m,
        "first_wall_thickness_m": first_wall_m,
        "blanket_thickness_m": blanket_m,
        "back_wall_thickness_m": back_wall_m,
        "inner_vacuum_radius_m": r_inner_cm / 100.0,
        "z_half_height_m": z_half_cm / 100.0,
        "n_phi_tally": n_phi_tally,
        "n_blanket_depth": n_blanket_depth,
        "cell_ids": {cell.name: cell.id for cell in cells},
        "source_mesh_id": source_mesh.id,
        "wall_mesh_id": wall_mesh.id,
        "heat_mesh_id": heat_mesh.id,
        "heat_r_edges_cm": heat_r_edges.tolist(),
    }


def latest_statepoint(case: Path) -> Path:
    files = sorted(case.glob("statepoint.*.h5"))
    if not files:
        raise FileNotFoundError(f"no OpenMC statepoint in {case}")
    return files[-1]


def score_values(tally, score: str) -> tuple[np.ndarray, np.ndarray]:
    mean = np.asarray(tally.get_values(scores=[score]), dtype=float).ravel()
    sd = np.asarray(tally.get_values(scores=[score], value="std_dev"), dtype=float).ravel()
    return mean, sd


def fold_period(values: np.ndarray, uncertainties: np.ndarray, nfp: int) -> tuple[np.ndarray, np.ndarray]:
    if len(values) % nfp:
        raise ValueError("toroidal bin count must be divisible by NFP")
    n_phase = len(values) // nfp
    stacked = np.vstack([values[period * n_phase : (period + 1) * n_phase] for period in range(nfp)])
    stacked_sd = np.vstack(
        [uncertainties[period * n_phase : (period + 1) * n_phase] for period in range(nfp)]
    )
    return stacked.mean(axis=0), np.sqrt(np.sum(stacked_sd**2, axis=0)) / nfp


def read_wall_current(tally, n_phi: int) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    mesh_filter = tally.filters[0]
    mean = np.asarray(tally.mean, dtype=float).reshape(-1)
    sd = np.asarray(tally.std_dev, dtype=float).reshape(-1)
    direction_totals: dict[str, float] = {}
    direction_arrays: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for direction in ("x-min in", "x-min out"):
        values = np.zeros(n_phi)
        errors = np.zeros(n_phi)
        for idx, bin_value in enumerate(mesh_filter.bins):
            radial_idx, phi_idx, z_idx, surface = bin_value
            if radial_idx == 1 and z_idx == 1 and surface == direction:
                values[phi_idx - 1] += mean[idx]
                errors[phi_idx - 1] = np.hypot(errors[phi_idx - 1], sd[idx])
        direction_arrays[direction] = (values, errors)
        direction_totals[direction] = float(values.sum())

    incoming, incoming_sd = direction_arrays["x-min in"]
    if incoming.sum() <= 0.0:
        raise ValueError(f"no incoming radial wall current; direction totals={direction_totals}")
    return incoming, incoming_sd, direction_totals


def relative_uncertainty(mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    return np.divide(sd, np.abs(mean), out=np.zeros_like(sd), where=np.abs(mean) > 0.0)


def postprocess(case: Path, results: Path, manifest: dict, source_data: dict) -> dict:
    openmc = import_openmc()
    statepoint = latest_statepoint(case)
    sp = openmc.StatePoint(statepoint)
    n_phi = int(manifest["n_phi_tally"])
    n_phase = n_phi // NFP

    wall_tally = sp.get_tally(name="wall_current_by_toroidal_bin")
    current, current_sd, direction_totals = read_wall_current(wall_tally, n_phi)
    uncollided_tally = sp.get_tally(name="uncollided_wall_current_by_toroidal_bin")
    uncollided, uncollided_sd, uncollided_direction_totals = read_wall_current(
        uncollided_tally, n_phi
    )
    z_height_cm = 2.0 * manifest["z_half_height_m"] * 100.0
    dphi = 2.0 * np.pi / n_phi
    wall_bin_area_cm2 = manifest["wall_radius_m"] * 100.0 * dphi * z_height_cm
    current_density = current / wall_bin_area_cm2
    current_density_sd = current_sd / wall_bin_area_cm2
    folded_current, folded_current_sd = fold_period(current_density, current_density_sd, NFP)
    current_rel = folded_current / folded_current.mean()
    current_rel_sd = folded_current_sd / folded_current.mean()
    uncollided_density = uncollided / wall_bin_area_cm2
    uncollided_density_sd = uncollided_sd / wall_bin_area_cm2
    folded_uncollided, folded_uncollided_sd = fold_period(
        uncollided_density, uncollided_density_sd, NFP
    )
    uncollided_rel = folded_uncollided / folded_uncollided.mean()
    uncollided_rel_sd = folded_uncollided_sd / folded_uncollided.mean()

    mesh_tally = sp.get_tally(name="heating_and_tritium_by_toroidal_depth_bin")
    mesh_filter = mesh_tally.filters[0]
    mesh = mesh_filter.mesh
    dimensions = tuple(int(value) for value in mesh.dimension)
    heat, heat_sd = score_values(mesh_tally, "heating-local")
    trit, trit_sd = score_values(mesh_tally, "(n,Xt)")
    heat = heat.reshape(dimensions, order="F")
    heat_sd = heat_sd.reshape(dimensions, order="F")
    trit = trit.reshape(dimensions, order="F")
    trit_sd = trit_sd.reshape(dimensions, order="F")
    volumes_cm3 = cylindrical_cell_volumes(mesh)

    incident_fraction = float(current.sum())
    full_wall_area_m2 = 2.0 * np.pi * manifest["wall_radius_m"] * (2.0 * manifest["z_half_height_m"])
    source_rate_s = REFERENCE_NWL_W_M2 * full_wall_area_m2 / (
        incident_fraction * E_SOURCE_EV * EV_TO_J
    )
    heat_w = heat * EV_TO_J * source_rate_s
    qvol_w_m3 = heat_w / (volumes_cm3 * 1.0e-6)
    heat_phi = heat.sum(axis=(0, 2))
    heat_phi_sd = np.sqrt(np.sum(heat_sd**2, axis=(0, 2)))
    folded_heat, folded_heat_sd = fold_period(heat_phi, heat_phi_sd, NFP)
    heat_rel = folded_heat / folded_heat.mean()
    heat_rel_sd = folded_heat_sd / folded_heat.mean()
    trit_phi = trit.sum(axis=(0, 2))
    trit_phi_sd = np.sqrt(np.sum(trit_sd**2, axis=(0, 2)))
    folded_trit, folded_trit_sd = fold_period(trit_phi, trit_phi_sd, NFP)

    cell_tally = sp.get_tally(name="heating_and_tritium_by_region")
    cell_filter = cell_tally.filters[0]
    cell_ids = [int(value) for value in cell_filter.bins]
    cell_heat, cell_heat_sd = score_values(cell_tally, "heating-local")
    cell_trit, cell_trit_sd = score_values(cell_tally, "(n,Xt)")
    region_by_id = {int(value): name for name, value in manifest["cell_ids"].items()}
    region_rows = []
    for idx, cell_id in enumerate(cell_ids):
        region_rows.append(
            {
                "region": region_by_id.get(cell_id, str(cell_id)),
                "cell_id": cell_id,
                "heating_eV_per_source": float(cell_heat[idx]),
                "heating_relative_uncertainty": float(relative_uncertainty(cell_heat, cell_heat_sd)[idx]),
                "tritons_per_source": float(cell_trit[idx]),
                "tritium_relative_uncertainty": float(relative_uncertainty(cell_trit, cell_trit_sd)[idx]),
            }
        )
    blanket_row = next(row for row in region_rows if row["region"] == "flibe_blanket")
    tbr_outboard = float(blanket_row["tritons_per_source"])
    tbr_per_incident = tbr_outboard / incident_fraction

    v2 = json.loads(V2_PATH.read_text())
    v2_nwl = np.asarray(v2["nwl_rel"], dtype=float)
    v2_d = np.asarray(v2["d_m"], dtype=float)
    phase_fraction = (np.arange(n_phase) + 0.5) / n_phase
    s_phase = phase_fraction * float(np.asarray(v2["s_m"], dtype=float)[-1])
    if len(v2_nwl) % n_phase == 0:
        samples_per_bin = len(v2_nwl) // n_phase
        los_rel = v2_nwl.reshape(n_phase, samples_per_bin).mean(axis=1)
        standoff = v2_d.reshape(n_phase, samples_per_bin).mean(axis=1)
    else:
        dense_fraction = (np.arange(len(v2_nwl)) + 0.5) / len(v2_nwl)
        los_rel = np.interp(phase_fraction, dense_fraction, v2_nwl, period=1.0)
        standoff = np.interp(phase_fraction, dense_fraction, v2_d, period=1.0)
    los_rel /= los_rel.mean()
    correlation = float(np.corrcoef(uncollided_rel, los_rel)[0, 1])
    rmse = float(np.sqrt(np.mean((uncollided_rel - los_rel) ** 2)))

    proxy = json.loads(PROXY_PATH.read_text()) if PROXY_PATH.exists() else {}
    proxy_tbr = proxy.get("tbr_column")
    proxy_channel_tbr = proxy.get("tbr_channel_slice")
    depth_heat = heat.sum(axis=(1, 2))
    depth_heat_sd = np.sqrt(np.sum(heat_sd**2, axis=(1, 2)))
    depth_trit = trit.sum(axis=(1, 2))
    depth_trit_sd = np.sqrt(np.sum(trit_sd**2, axis=(1, 2)))
    r_edges_m = np.asarray(mesh.r_grid, dtype=float) / 100.0
    depth_heat_rel_unc = relative_uncertainty(depth_heat, depth_heat_sd)
    depth_trit_rel_unc = relative_uncertainty(depth_trit, depth_trit_sd)

    leakage = None
    for row in sp.global_tallies:
        if row["name"].decode() == "leakage":
            leakage = {"mean": float(row["mean"]), "std_dev": float(row["std_dev"])}
            break

    resolution_check = None
    prior_path = results / "a1_openmc_3d_summary_nphi80.json"
    if prior_path.exists():
        prior = json.loads(prior_path.read_text())
        resolution_check = {
            "coarse_toroidal_bins_full_torus": int(prior["geometry"]["n_phi_tally"]),
            "fine_toroidal_bins_full_torus": n_phi,
            "uncollided_min_change_percent": float(
                100.0
                * (
                    uncollided_rel.min()
                    / prior["openmc_uncollided_wall_current_min_over_mean"]
                    - 1.0
                )
            ),
            "uncollided_max_change_percent": float(
                100.0
                * (
                    uncollided_rel.max()
                    / prior["openmc_uncollided_wall_current_max_over_mean"]
                    - 1.0
                )
            ),
            "outboard_tbr_change_percent": float(
                100.0 * (tbr_outboard / prior["tbr_outboard_per_fusion_source"] - 1.0)
            ),
        }

    results.mkdir(parents=True, exist_ok=True)
    phase_rows = []
    for idx in range(n_phase):
        phase_rows.append(
            {
                "phase_bin": idx,
                "phase_fraction": float(phase_fraction[idx]),
                "phase_angle_rad": float(phase_fraction[idx] * 2.0 * np.pi / NFP),
                "mapped_s_m": float(s_phase[idx]),
                "standoff_m": float(standoff[idx]),
                "openmc_wall_current_rel": float(current_rel[idx]),
                "openmc_wall_current_rel_sd": float(current_rel_sd[idx]),
                "openmc_uncollided_wall_current_rel": float(uncollided_rel[idx]),
                "openmc_uncollided_wall_current_rel_sd": float(uncollided_rel_sd[idx]),
                "v2_los_nwl_rel": float(los_rel[idx]),
                "openmc_heating_rel": float(heat_rel[idx]),
                "openmc_heating_rel_sd": float(heat_rel_sd[idx]),
                "tritons_per_source_folded_bin": float(folded_trit[idx]),
                "tritons_per_source_folded_bin_sd": float(folded_trit_sd[idx]),
            }
        )
    phase_csv = results / "a1_toroidal_wall_load.csv"
    with phase_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(phase_rows[0]))
        writer.writeheader()
        writer.writerows(phase_rows)

    depth_rows = []
    for idx in range(len(depth_heat)):
        depth_rows.append(
            {
                "radial_bin": idx,
                "r0_m": float(r_edges_m[idx]),
                "r1_m": float(r_edges_m[idx + 1]),
                "depth0_m": float(r_edges_m[idx] - manifest["wall_radius_m"]),
                "depth1_m": float(r_edges_m[idx + 1] - manifest["wall_radius_m"]),
                "heating_eV_per_source": float(depth_heat[idx]),
                "heating_relative_uncertainty": float(relative_uncertainty(depth_heat, depth_heat_sd)[idx]),
                "normalized_mean_qvol_W_m3": float(
                    np.sum(heat_w[idx]) / (np.sum(volumes_cm3[idx]) * 1.0e-6)
                ),
                "tritons_per_source": float(depth_trit[idx]),
                "tritium_relative_uncertainty": float(relative_uncertainty(depth_trit, depth_trit_sd)[idx]),
            }
        )
    depth_csv = results / "a1_radial_heating_tbr.csv"
    with depth_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(depth_rows[0]))
        writer.writeheader()
        writer.writerows(depth_rows)

    region_csv = results / "a1_region_heating_tbr.csv"
    with region_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(region_rows[0]))
        writer.writeheader()
        writer.writerows(region_rows)

    summary = {
        "status": "A1 native-OpenMC CSG 3D outboard annular blanket result; not reactor CAD neutronics",
        "case": str(case.relative_to(REPO)),
        "statepoint": str(statepoint.relative_to(REPO)),
        "source": source_data["source"],
        "geometry": manifest,
        "particles": manifest["particles"],
        "batches": manifest["batches"],
        "source_histogram_capture_fraction": source_data["capture_fraction"],
        "source_histogram_nonzero_bins": source_data["nonzero_bins"],
        "incident_outboard_wall_current_per_source": incident_fraction,
        "wall_current_direction_totals": direction_totals,
        "uncollided_outboard_wall_current_per_source": float(uncollided.sum()),
        "uncollided_wall_current_direction_totals": uncollided_direction_totals,
        "uncollided_fraction_of_total_incoming_current": float(uncollided.sum() / current.sum()),
        "reference_neutron_wall_load_W_m2": REFERENCE_NWL_W_M2,
        "normalizing_source_rate_per_s": source_rate_s,
        "openmc_wall_current_min_over_mean": float(current_rel.min()),
        "openmc_wall_current_max_over_mean": float(current_rel.max()),
        "openmc_wall_current_peak_to_trough": float(current_rel.max() / current_rel.min()),
        "openmc_uncollided_wall_current_min_over_mean": float(uncollided_rel.min()),
        "openmc_uncollided_wall_current_max_over_mean": float(uncollided_rel.max()),
        "openmc_uncollided_wall_current_peak_to_trough": float(
            uncollided_rel.max() / uncollided_rel.min()
        ),
        "v2_los_min_over_mean_on_a1_bins": float(los_rel.min()),
        "v2_los_max_over_mean_on_a1_bins": float(los_rel.max()),
        "v2_los_exact_grid_min_over_mean": float(v2_nwl.min() / v2_nwl.mean()),
        "v2_los_exact_grid_max_over_mean": float(v2_nwl.max() / v2_nwl.mean()),
        "openmc_vs_v2_pearson_r": correlation,
        "openmc_vs_v2_normalized_rmse": rmse,
        "openmc_to_exact_v2_peak_departure_ratio": float(
            (uncollided_rel.max() - 1.0) / (v2_nwl.max() / v2_nwl.mean() - 1.0)
        ),
        "openmc_to_exact_v2_peak_to_trough_ratio": float(
            (uncollided_rel.max() / uncollided_rel.min())
            / ((v2_nwl.max() / v2_nwl.mean()) / (v2_nwl.min() / v2_nwl.mean()))
        ),
        "wall_load_verdict": (
            "OpenMC confirms the toroidal phase and existence of standoff-driven neutron-wall-load "
            "modulation but does not confirm the full 0.69-1.69 LOS amplitude."
        ),
        "tbr_outboard_per_fusion_source": tbr_outboard,
        "tbr_outboard_per_neutron_entering_wall": tbr_per_incident,
        "tbr_1d_proxy_column_per_incident_source": proxy_tbr,
        "tbr_1d_proxy_first_3cm_per_incident_source": proxy_channel_tbr,
        "conditional_tbr_change_vs_1d_proxy_percent": (
            float(100.0 * (tbr_per_incident / proxy_tbr - 1.0)) if proxy_tbr else None
        ),
        "total_heating_eV_per_source": float(cell_heat.sum()),
        "energy_deposition_fraction": float(cell_heat.sum() / E_SOURCE_EV),
        "leakage_fraction": leakage,
        "openmc_runtime_seconds": {key: float(value) for key, value in sp.runtime.items()},
        "max_wall_current_relative_uncertainty": float(relative_uncertainty(current, current_sd).max()),
        "max_uncollided_wall_current_relative_uncertainty": float(
            relative_uncertainty(uncollided, uncollided_sd).max()
        ),
        "max_nonzero_mesh_heating_relative_uncertainty": float(
            relative_uncertainty(heat[heat > 0.0], heat_sd[heat > 0.0]).max()
        ),
        "p95_nonzero_mesh_heating_relative_uncertainty": float(
            np.quantile(relative_uncertainty(heat[heat > 0.0], heat_sd[heat > 0.0]), 0.95)
        ),
        "max_radial_bin_heating_relative_uncertainty": float(depth_heat_rel_unc.max()),
        "max_flibe_radial_bin_tritium_relative_uncertainty": float(
            depth_trit_rel_unc[1:-1].max()
        ),
        "toroidal_resolution_check": resolution_check,
        "phase_csv": str(phase_csv.relative_to(REPO)),
        "depth_csv": str(depth_csv.relative_to(REPO)),
        "region_csv": str(region_csv.relative_to(REPO)),
        "limitations": [
            "outboard annular CSG blanket, not a non-axisymmetric CAD first wall",
            "DESC fusion source is histogrammed onto a finite cylindrical source mesh",
            "source emissivity uses the assumed core-peaked (1-rho^2)^3 profile",
            "top, bottom, and inboard boundaries are vacuum; TBR is an outboard-sector metric",
            "90% Li-6 and proxy steel are retained for comparison with the earlier C3/1D models",
            "no photon transport, depletion, activation, magnet geometry, or thermal feedback",
        ],
    }
    summary_path = results / "a1_openmc_3d_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    np.savez_compressed(
        results / "a1_mesh_fields.npz",
        qvol_W_m3=qvol_w_m3,
        heating_eV_per_source=heat,
        tritons_per_source=trit,
        r_edges_cm=np.asarray(mesh.r_grid),
        phi_edges_rad=np.asarray(mesh.phi_grid),
        z_edges_cm=np.asarray(mesh.z_grid),
    )
    figure = results / "a1_openmc_3d_toroidal_neutronics.png"
    plot_result(figure, summary, phase_rows, depth_rows, source_data, manifest)
    print(json.dumps(summary, indent=2))
    print(f"wrote {summary_path.relative_to(REPO)}")
    print(f"wrote {figure.relative_to(REPO)}")
    return summary


def plot_result(
    output: Path,
    summary: dict,
    phase_rows: list[dict],
    depth_rows: list[dict],
    source_data: dict,
    manifest: dict,
) -> None:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(14.5, 9.0), constrained_layout=True)
    ax3d = fig.add_subplot(2, 2, 1, projection="3d")
    period = 2.0 * np.pi / NFP
    phi = np.linspace(0.0, period, 80)
    z = np.linspace(-manifest["z_half_height_m"], manifest["z_half_height_m"], 12)
    pp, zz = np.meshgrid(phi, z)
    radii = [
        (manifest["wall_radius_m"], "#8c969f", 0.42),
        (manifest["wall_radius_m"] + manifest["first_wall_thickness_m"], "#19856f", 0.25),
        (
            manifest["wall_radius_m"]
            + manifest["first_wall_thickness_m"]
            + manifest["blanket_thickness_m"],
            "#19856f",
            0.18,
        ),
    ]
    for radius, color, alpha in radii:
        ax3d.plot_surface(radius * np.cos(pp), radius * np.sin(pp), zz, color=color, alpha=alpha, linewidth=0)
    ax3d.set_title("Native CSG outboard blanket sector")
    ax3d.set_xlabel("x [m]")
    ax3d.set_ylabel("y [m]")
    ax3d.set_zlabel("z [m]")
    ax3d.view_init(elev=24, azim=-58)
    ax3d.set_box_aspect((1.0, 0.8, 0.45))

    ax = fig.add_subplot(2, 2, 2)
    x = np.asarray([row["phase_fraction"] for row in phase_rows])
    current = np.asarray([row["openmc_wall_current_rel"] for row in phase_rows])
    current_sd = np.asarray([row["openmc_wall_current_rel_sd"] for row in phase_rows])
    uncollided = np.asarray([row["openmc_uncollided_wall_current_rel"] for row in phase_rows])
    uncollided_sd = np.asarray(
        [row["openmc_uncollided_wall_current_rel_sd"] for row in phase_rows]
    )
    los = np.asarray([row["v2_los_nwl_rel"] for row in phase_rows])
    ax.errorbar(x, current, yerr=current_sd, marker="o", color="#1f6aa5", lw=1.7, capsize=2, label="OpenMC total")
    ax.errorbar(
        x,
        uncollided,
        yerr=uncollided_sd,
        marker="^",
        color="#19856f",
        lw=2,
        capsize=2,
        label="OpenMC zero-collision",
    )
    ax.plot(x, los, marker="s", color="#c07830", lw=1.7, ls="--", label="v2 LOS proxy")
    ax.axhline(1.0, color="0.55", lw=1, ls=":")
    ax.set_xlabel("fraction of one field period")
    ax.set_ylabel("neutron wall current / mean")
    ax.set_title("Toroidal neutron wall-load modulation")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    ax.text(
        0.02,
        0.04,
        f"total {summary['openmc_wall_current_min_over_mean']:.2f}-{summary['openmc_wall_current_max_over_mean']:.2f}x\n"
        f"zero-collision {summary['openmc_uncollided_wall_current_min_over_mean']:.2f}-{summary['openmc_uncollided_wall_current_max_over_mean']:.2f}x\n"
        f"LOS {summary['v2_los_min_over_mean_on_a1_bins']:.2f}-{summary['v2_los_max_over_mean_on_a1_bins']:.2f}x\n"
        f"r = {summary['openmc_vs_v2_pearson_r']:.2f}",
        transform=ax.transAxes,
        va="bottom",
    )

    ax = fig.add_subplot(2, 2, 3)
    source = np.asarray(source_data["strengths"], dtype=float).sum(axis=2)
    source /= source.max()
    phi_edges = np.asarray(source_data["phi_edges_rad"])
    r_edges = np.asarray(source_data["r_edges_cm"]) / 100.0
    mesh = ax.pcolormesh(phi_edges, r_edges, source, shading="auto", cmap="magma")
    ax.set_xlabel("toroidal angle [rad]")
    ax.set_ylabel("major radius R [m]")
    ax.set_title("DESC fusion-source mesh, z-integrated")
    fig.colorbar(mesh, ax=ax, label="relative source strength")

    ax = fig.add_subplot(2, 2, 4)
    depth_mid = np.asarray([(row["depth0_m"] + row["depth1_m"]) / 2.0 for row in depth_rows])
    qvol = np.asarray([row["normalized_mean_qvol_W_m3"] for row in depth_rows]) / 1.0e6
    trit = np.asarray([row["tritons_per_source"] for row in depth_rows])
    ax.step(depth_mid, qvol, where="mid", color="#b23a2b", lw=2, label="heating")
    ax.set_xlabel("depth from plasma-facing wall [m]")
    ax.set_ylabel("q''' at 1 MW/m2 mean NWL [MW/m3]", color="#b23a2b")
    ax.tick_params(axis="y", labelcolor="#b23a2b")
    ax.grid(alpha=0.25)
    twin = ax.twinx()
    twin.step(depth_mid, trit, where="mid", color="#287a68", lw=2, label="tritium")
    twin.set_ylabel("tritons per source by radial bin", color="#287a68")
    twin.tick_params(axis="y", labelcolor="#287a68")
    ax.set_title("Radial heating and tritium production")
    ax.text(
        0.98,
        0.94,
        f"TBR/source = {summary['tbr_outboard_per_fusion_source']:.3f}\n"
        f"TBR/incident = {summary['tbr_outboard_per_neutron_entering_wall']:.3f}\n"
        f"1D proxy = {summary['tbr_1d_proxy_column_per_incident_source']:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
    )

    fig.suptitle("A1 native-OpenMC 3D toroidal neutronics: wall load, heating, and TBR", fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, default=DEFAULT_CASE)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--particles", type=int, default=20_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--n-phi-tally", type=int, default=80)
    parser.add_argument("--n-blanket-depth", type=int, default=10)
    parser.add_argument("--n-rho-source", type=int, default=20)
    parser.add_argument("--n-theta-source", type=int, default=16)
    parser.add_argument("--n-zeta-source", type=int, default=80)
    parser.add_argument("--n-r-source-mesh", type=int, default=24)
    parser.add_argument("--n-phi-source-mesh", type=int, default=80)
    parser.add_argument("--n-z-source-mesh", type=int, default=24)
    parser.add_argument("--first-wall-m", type=float, default=0.005)
    parser.add_argument("--blanket-m", type=float, default=0.50)
    parser.add_argument("--back-wall-m", type=float, default=0.005)
    parser.add_argument("--li6", type=float, default=0.90)
    parser.add_argument("--seed", type=int, default=170711)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-run", action="store_true")
    parser.add_argument("--build-only", action="store_true")
    args = parser.parse_args()
    if args.n_phi_tally % NFP:
        raise ValueError("n-phi-tally must be divisible by five field periods")

    v2 = json.loads(V2_PATH.read_text())
    wall_radius_m = float(v2["R_wall_m"])
    case = args.case if args.case.is_absolute() else WORK / args.case
    results = args.results if args.results.is_absolute() else REPO / args.results
    manifest_path = case / "a1_model_manifest.json"
    source_path = case / "a1_desc_source_mesh.npz"

    if args.no_run:
        manifest = json.loads(manifest_path.read_text())
        source_npz = np.load(source_path)
        source_data = {
            "source": manifest["source_description"],
            "r_edges_cm": source_npz["r_edges_cm"],
            "phi_edges_rad": source_npz["phi_edges_rad"],
            "z_edges_cm": source_npz["z_edges_cm"],
            "strengths": source_npz["strengths"],
            "capture_fraction": manifest["source_capture_fraction"],
            "nonzero_bins": manifest["source_nonzero_bins"],
        }
        postprocess(case, results, manifest, source_data)
        return

    if case.exists():
        if not args.force:
            raise FileExistsError(f"{case} exists; pass --force to rebuild this generated case")
        shutil.rmtree(case)
    case.mkdir(parents=True)
    source_data = desc_fusion_source(
        n_rho=args.n_rho_source,
        n_theta=args.n_theta_source,
        n_zeta=args.n_zeta_source,
        n_r_mesh=args.n_r_source_mesh,
        n_phi_mesh=args.n_phi_source_mesh,
        n_z_mesh=args.n_z_source_mesh,
        wall_radius_m=wall_radius_m,
    )
    model = build_model(
        case,
        source_data=source_data,
        wall_radius_m=wall_radius_m,
        first_wall_m=args.first_wall_m,
        blanket_m=args.blanket_m,
        back_wall_m=args.back_wall_m,
        li6_fraction=args.li6,
        n_phi_tally=args.n_phi_tally,
        n_blanket_depth=args.n_blanket_depth,
        particles=args.particles,
        batches=args.batches,
        seed=args.seed,
    )
    manifest = {
        **model,
        "source_description": source_data["source"],
        "source_capture_fraction": source_data["capture_fraction"],
        "source_nonzero_bins": source_data["nonzero_bins"],
        "source_sample_points": source_data["sample_points"],
        "source_mesh_dimensions": [args.n_r_source_mesh, args.n_phi_source_mesh, args.n_z_source_mesh],
        "plasma_R_min_m": source_data["plasma_R_min_m"],
        "plasma_R_max_m": source_data["plasma_R_max_m"],
        "plasma_abs_Z_max_m": source_data["plasma_abs_Z_max_m"],
        "standoff_range_m": [float(min(v2["d_m"])), float(max(v2["d_m"]))],
        "li6_fraction": args.li6,
        "particles": args.particles,
        "batches": args.batches,
        "seed": args.seed,
        "source_energy_eV": E_SOURCE_EV,
        "cross_sections": XS_PATH,
        "openmc_executable": OPENMC_EXE,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    np.savez_compressed(
        source_path,
        r_edges_cm=source_data["r_edges_cm"],
        phi_edges_rad=source_data["phi_edges_rad"],
        z_edges_cm=source_data["z_edges_cm"],
        strengths=source_data["strengths"],
    )
    print(
        f"built {case.relative_to(REPO)} with {source_data['nonzero_bins']} nonzero source bins; "
        f"captured {source_data['capture_fraction']:.8f} of DESC source weight"
    )
    if args.build_only:
        print("build-only requested; run OpenMC in the case directory, then rerun with --no-run")
        return
    subprocess.run(
        [OPENMC_EXE],
        cwd=case,
        check=True,
        stdout=(case / "openmc.log").open("w"),
        stderr=subprocess.STDOUT,
        env={**os.environ, "OPENMC_CROSS_SECTIONS": XS_PATH},
    )
    postprocess(case, results, manifest, source_data)


if __name__ == "__main__":
    main()

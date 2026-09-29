#!/usr/bin/env python3
"""Generate and run a tiny C3 OpenMC material-stack heating tally.

This is the next rung above the all-FLiBe slab smoke case: a simple CSG stack
with steel-proxy wall regions and one FLiBe region. It is still not reactor
neutronics or a blanket model. The goal is only to produce a local
material/depth source-placement artifact that can replace a literature-lumped
fluid/solid split in the C3 CHT comparison.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
from pathlib import Path

import numpy as np


XS_PATH = os.environ.get("OPENMC_CROSS_SECTIONS", "/root/openmc_data/combined_neutron_photon_all.xml")
OPENMC_EXE_CANDIDATES = [
    "/root/openmc-src/build/bin/openmc",
    "/root/miniforge3/envs/openmc-env/bin/openmc",
    "openmc",
]
EV_TO_J = 1.602176634e-19
C3_REFERENCE_POWER_W = 427.0320391
C3_REFERENCE_HEAT_FLUX_W_M2 = 100000.0


REGION_NAMES = ("front_wall", "flibe", "back_wall")
REGION_MATERIAL_CLASSES = {
    "front_wall": "solid",
    "flibe": "fluid",
    "back_wall": "solid",
}


def pick_openmc_exec() -> str:
    for candidate in OPENMC_EXE_CANDIDATES:
        if candidate == "openmc" or Path(candidate).exists():
            return candidate
    return "openmc"


def import_openmc():
    os.environ["OPENMC_CROSS_SECTIONS"] = XS_PATH
    import openmc  # noqa: PLC0415

    return openmc


def make_steel_proxy(openmc, name: str):
    """Small ENDF-safe steel proxy using nuclides confirmed in local data."""
    mat = openmc.Material(name=name)
    mat.set_density("g/cm3", 7.8)
    mat.add_nuclide("Fe56", 0.88, "ao")
    mat.add_nuclide("Cr52", 0.09, "ao")
    mat.add_nuclide("Ni58", 0.03, "ao")
    return mat


def make_flibe(openmc, li6_fraction: float):
    flibe = openmc.Material(name="flibe_2LiF_BeF2_li6_enriched")
    flibe.set_density("g/cm3", 1.94)
    # 2LiF-BeF2 atom ratios: Li:Be:F = 2:1:4.
    flibe.add_nuclide("Li6", 2.0 * li6_fraction, "ao")
    flibe.add_nuclide("Li7", 2.0 * (1.0 - li6_fraction), "ao")
    flibe.add_nuclide("Be9", 1.0, "ao")
    flibe.add_nuclide("F19", 4.0, "ao")
    return flibe


def build_model(
    case: Path,
    *,
    front_wall_thickness: float,
    flibe_thickness: float,
    back_wall_thickness: float,
    half_width: float,
    half_height: float,
    bins: int,
    particles: int,
    batches: int,
    li6_fraction: float,
) -> dict:
    openmc = import_openmc()

    front_steel = make_steel_proxy(openmc, "front_wall_steel_proxy")
    flibe = make_flibe(openmc, li6_fraction)
    back_steel = make_steel_proxy(openmc, "back_wall_steel_proxy")
    openmc.Materials([front_steel, flibe, back_steel]).export_to_xml(case / "materials.xml")

    x0_value = 0.0
    x1_value = front_wall_thickness
    x2_value = x1_value + flibe_thickness
    x3_value = x2_value + back_wall_thickness

    x0 = openmc.XPlane(x0=x0_value, boundary_type="vacuum")
    x1 = openmc.XPlane(x0=x1_value)
    x2 = openmc.XPlane(x0=x2_value)
    x3 = openmc.XPlane(x0=x3_value, boundary_type="vacuum")
    ymin = openmc.YPlane(y0=-half_width, boundary_type="reflective")
    ymax = openmc.YPlane(y0=half_width, boundary_type="reflective")
    zmin = openmc.ZPlane(z0=-half_height, boundary_type="reflective")
    zmax = openmc.ZPlane(z0=half_height, boundary_type="reflective")
    yz_region = +ymin & -ymax & +zmin & -zmax

    front_cell = openmc.Cell(
        name="front_wall",
        fill=front_steel,
        region=+x0 & -x1 & yz_region,
    )
    flibe_cell = openmc.Cell(
        name="flibe",
        fill=flibe,
        region=+x1 & -x2 & yz_region,
    )
    back_cell = openmc.Cell(
        name="back_wall",
        fill=back_steel,
        region=+x2 & -x3 & yz_region,
    )
    cells = [front_cell, flibe_cell, back_cell]
    openmc.Geometry(cells).export_to_xml(case / "geometry.xml")

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = 0
    settings.seed = 170706
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            (1.0e-8, -half_width, -half_height),
            (2.0e-8, half_width, half_height),
            only_fissionable=False,
        ),
        angle=openmc.stats.Monodirectional(reference_uvw=(1.0, 0.0, 0.0)),
        energy=openmc.stats.Discrete([14.1e6], [1.0]),
    )
    settings.export_to_xml(case / "settings.xml")

    total_thickness = x3_value
    mesh = openmc.RegularMesh()
    mesh.dimension = (bins, 1, 1)
    mesh.lower_left = (0.0, -half_width, -half_height)
    mesh.upper_right = (total_thickness, half_width, half_height)
    mesh_filter = openmc.MeshFilter(mesh)
    cell_filter = openmc.CellFilter(cells)

    heat_cell = openmc.Tally(name="heating_local_by_region")
    heat_cell.filters = [cell_filter]
    heat_cell.scores = ["heating-local"]

    flux_cell = openmc.Tally(name="flux_by_region")
    flux_cell.filters = [cell_filter]
    flux_cell.scores = ["flux"]

    heat_depth = openmc.Tally(name="heating_local_by_depth")
    heat_depth.filters = [mesh_filter]
    heat_depth.scores = ["heating-local"]

    heat_cell_depth = openmc.Tally(name="heating_local_by_region_depth")
    heat_cell_depth.filters = [cell_filter, mesh_filter]
    heat_cell_depth.scores = ["heating-local"]

    openmc.Tallies([heat_cell, flux_cell, heat_depth, heat_cell_depth]).export_to_xml(
        case / "tallies.xml"
    )

    region_bounds = {
        "front_wall": [x0_value, x1_value],
        "flibe": [x1_value, x2_value],
        "back_wall": [x2_value, x3_value],
    }
    cell_ids = {cell.name: cell.id for cell in cells}
    material_ids = {
        "front_wall": front_steel.id,
        "flibe": flibe.id,
        "back_wall": back_steel.id,
    }
    return {
        "total_thickness_m": total_thickness,
        "region_bounds_m": region_bounds,
        "cell_ids": cell_ids,
        "material_ids": material_ids,
    }


def latest_statepoint(case: Path) -> Path:
    statepoints = sorted(case.glob("statepoint.*.h5"))
    if not statepoints:
        raise FileNotFoundError(f"no statepoint files in {case}")
    return statepoints[-1]


def relative_uncertainty(mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    return np.divide(sd, np.abs(mean), out=np.zeros_like(sd), where=np.abs(mean) > 0.0)


def normalized_positive(values: np.ndarray) -> tuple[np.ndarray, float]:
    positive = np.maximum(values, 0.0)
    total = float(np.sum(positive))
    if total <= 0.0:
        raise ValueError("OpenMC tally produced no positive heating values")
    return positive / total, total


def write_region_csv(
    path: Path,
    *,
    region_bounds: dict,
    cell_mean: np.ndarray,
    cell_sd: np.ndarray,
    cell_shape: np.ndarray,
    target_power: float,
) -> list[dict]:
    rows = []
    rel_unc = relative_uncertainty(cell_mean, cell_sd)
    for i, region in enumerate(REGION_NAMES):
        x0, x1 = region_bounds[region]
        power = float(cell_shape[i] * target_power)
        row = {
            "region": region,
            "material_class": REGION_MATERIAL_CLASSES[region],
            "x0_m": x0,
            "x1_m": x1,
            "thickness_m": x1 - x0,
            "heating_mean_eV_per_source": float(cell_mean[i]),
            "relative_uncertainty": float(rel_unc[i]),
            "normalized_fraction": float(cell_shape[i]),
            "normalized_power_W": power,
        }
        rows.append(row)

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_depth_csv(
    path: Path,
    *,
    total_thickness: float,
    depth_mean: np.ndarray,
    depth_sd: np.ndarray,
    depth_shape: np.ndarray,
    target_heat_flux: float,
) -> list[dict]:
    bins = len(depth_mean)
    dx = total_thickness / bins
    rel_unc = relative_uncertainty(depth_mean, depth_sd)
    rows = []
    for i in range(bins):
        rows.append(
            {
                "bin": i,
                "x0_m": i * dx,
                "x1_m": (i + 1) * dx,
                "x_mid_m": (i + 0.5) * dx,
                "heating_mean_eV_per_source": float(depth_mean[i]),
                "relative_uncertainty": float(rel_unc[i]),
                "normalized_shape": float(depth_shape[i]),
                "normalized_qvol_W_m3": float(depth_shape[i] * target_heat_flux / dx),
            }
        )
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_region_depth_csv(
    path: Path,
    *,
    total_thickness: float,
    region_depth_mean: np.ndarray,
    region_depth_sd: np.ndarray,
    total_heating_eV: float,
    target_heat_flux: float,
) -> list[dict]:
    bins = region_depth_mean.shape[1]
    dx = total_thickness / bins
    rel_unc = relative_uncertainty(region_depth_mean, region_depth_sd)
    rows = []
    for r_idx, region in enumerate(REGION_NAMES):
        for b_idx in range(bins):
            mean = float(region_depth_mean[r_idx, b_idx])
            shape = max(mean, 0.0) / total_heating_eV if total_heating_eV > 0.0 else 0.0
            rows.append(
                {
                    "region": region,
                    "material_class": REGION_MATERIAL_CLASSES[region],
                    "bin": b_idx,
                    "x0_m": b_idx * dx,
                    "x1_m": (b_idx + 1) * dx,
                    "x_mid_m": (b_idx + 0.5) * dx,
                    "heating_mean_eV_per_source": mean,
                    "relative_uncertainty": float(rel_unc[r_idx, b_idx]),
                    "normalized_shape": shape,
                    "normalized_qvol_W_m3": shape * target_heat_flux / dx,
                }
            )
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def plot_outputs(out_prefix: Path, region_rows: list[dict], depth_rows: list[dict]) -> str | None:
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        x_mid = np.array([row["x_mid_m"] for row in depth_rows])
        qvol = np.array([row["normalized_qvol_W_m3"] for row in depth_rows]) / 1.0e6
        names = [row["region"] for row in region_rows]
        powers = [row["normalized_power_W"] for row in region_rows]
        colors = ["#5c6670", "#236a8c", "#8a6f48"]

        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
        axes[0].bar(names, powers, color=colors)
        axes[0].set_ylabel("Normalized power (W)")
        axes[0].set_title("C3 material-stack heating split")
        axes[0].tick_params(axis="x", rotation=18)

        axes[1].step(x_mid, qvol, where="mid", color="#236a8c", linewidth=2.0)
        axes[1].set_xlabel("Depth into material stack (m)")
        axes[1].set_ylabel("Normalized qvol (MW/m3)")
        axes[1].set_title("Depth heating shape")
        axes[1].grid(True, alpha=0.25)

        fig.tight_layout()
        png = out_prefix.with_suffix(".png")
        fig.savefig(png, dpi=180)
        plt.close(fig)
        return str(png)
    except Exception:
        return None


def postprocess(
    case: Path,
    *,
    total_thickness: float,
    region_bounds: dict,
    bins: int,
    target_power: float,
    target_heat_flux: float,
    out_prefix: Path,
) -> dict:
    openmc = import_openmc()
    sp = openmc.StatePoint(latest_statepoint(case))

    cell_tally = sp.get_tally(name="heating_local_by_region")
    cell_mean = np.asarray(cell_tally.mean).reshape((len(REGION_NAMES),))
    cell_sd = np.asarray(cell_tally.std_dev).reshape((len(REGION_NAMES),))
    cell_shape, cell_total = normalized_positive(cell_mean)

    depth_tally = sp.get_tally(name="heating_local_by_depth")
    depth_mean = np.asarray(depth_tally.mean).reshape((bins,))
    depth_sd = np.asarray(depth_tally.std_dev).reshape((bins,))
    depth_shape, depth_total = normalized_positive(depth_mean)

    region_depth_tally = sp.get_tally(name="heating_local_by_region_depth")
    region_depth_mean = np.asarray(region_depth_tally.mean).reshape((len(REGION_NAMES), bins))
    region_depth_sd = np.asarray(region_depth_tally.std_dev).reshape((len(REGION_NAMES), bins))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    region_csv = out_prefix.with_name(out_prefix.name + "_region_split.csv")
    depth_csv = out_prefix.with_name(out_prefix.name + "_depth_profile.csv")
    region_depth_csv = out_prefix.with_name(out_prefix.name + "_region_depth.csv")

    region_rows = write_region_csv(
        region_csv,
        region_bounds=region_bounds,
        cell_mean=cell_mean,
        cell_sd=cell_sd,
        cell_shape=cell_shape,
        target_power=target_power,
    )
    depth_rows = write_depth_csv(
        depth_csv,
        total_thickness=total_thickness,
        depth_mean=depth_mean,
        depth_sd=depth_sd,
        depth_shape=depth_shape,
        target_heat_flux=target_heat_flux,
    )
    write_region_depth_csv(
        region_depth_csv,
        total_thickness=total_thickness,
        region_depth_mean=region_depth_mean,
        region_depth_sd=region_depth_sd,
        total_heating_eV=cell_total,
        target_heat_flux=target_heat_flux,
    )

    fluid_power = sum(row["normalized_power_W"] for row in region_rows if row["material_class"] == "fluid")
    solid_power = sum(row["normalized_power_W"] for row in region_rows if row["material_class"] == "solid")
    fluid_fraction = fluid_power / target_power
    solid_fraction = solid_power / target_power

    summary = {
        "case": str(case),
        "statepoint": str(latest_statepoint(case)),
        "status": "C3 material-stack OpenMC source-placement tally; not reactor neutronics",
        "total_thickness_m": total_thickness,
        "bins": bins,
        "target_power_W": target_power,
        "target_heat_flux_W_m2": target_heat_flux,
        "reference_area_m2": target_power / target_heat_flux,
        "cell_total_heating_eV_per_source": cell_total,
        "depth_total_heating_eV_per_source": depth_total,
        "cell_total_heating_J_per_source": cell_total * EV_TO_J,
        "fluid_fraction": fluid_fraction,
        "solid_fraction": solid_fraction,
        "fluid_power_W": fluid_power,
        "solid_power_W": solid_power,
        "front_wall_power_W": region_rows[0]["normalized_power_W"],
        "flibe_power_W": region_rows[1]["normalized_power_W"],
        "back_wall_power_W": region_rows[2]["normalized_power_W"],
        "max_region_relative_uncertainty": max(row["relative_uncertainty"] for row in region_rows),
        "max_depth_relative_uncertainty": max(row["relative_uncertainty"] for row in depth_rows),
        "region_csv": str(region_csv),
        "depth_csv": str(depth_csv),
        "region_depth_csv": str(region_depth_csv),
    }
    png = plot_outputs(out_prefix, region_rows, depth_rows)
    if png:
        summary["png"] = png

    json_path = out_prefix.with_suffix(".json")
    json_path.write_text(json.dumps(summary, indent=2))

    md_path = out_prefix.with_suffix(".md")
    md_path.write_text(
        "\n".join(
            [
                "# C3 OpenMC Material-Stack Heating Split",
                "",
                "Status: local source-placement tally, not reactor neutronics.",
                "",
                "## Model",
                "",
                "- geometry: steel-proxy front wall, FLiBe layer, steel-proxy back wall",
                "- source: fixed `14.1 MeV` monodirectional neutron source",
                "- tally score: `heating-local`",
                f"- total material-stack thickness: `{total_thickness:g} m`",
                f"- depth bins: `{bins}`",
                f"- target comparison power: `{target_power:.6g} W`",
                f"- target heat-flux scale: `{target_heat_flux:.6g} W/m2`",
                "",
                "## Tallied Split",
                "",
                "| Region | Class | Fraction | Power, W | Rel. uncertainty |",
                "|---|---:|---:|---:|---:|",
            ]
            + [
                "| {region} | {material_class} | `{normalized_fraction:.5f}` | "
                "`{normalized_power_W:.5f}` | `{relative_uncertainty:.4f}` |".format(**row)
                for row in region_rows
            ]
            + [
                "",
                "Aggregated CHT powers:",
                "",
                f"- fluid/FLiBe: `{fluid_power:.5f} W` (`{fluid_fraction:.5f}`)",
                f"- solid/wall: `{solid_power:.5f} W` (`{solid_fraction:.5f}`)",
                "",
                "## Interpretation",
                "",
                "This replaces the first C3 literature-lumped split with a tiny local",
                "OpenMC material-stack split. It is useful for source-placement plumbing",
                "and sensitivity direction only. It does not include Be multiplier zones,",
                "true stellarator blanket geometry, gamma transport validation, tritium,",
                "MHD, depletion, or coupled thermal feedback.",
                "",
                f"- region split CSV: `{region_csv}`",
                f"- depth profile CSV: `{depth_csv}`",
                f"- region-depth CSV: `{region_depth_csv}`",
                f"- JSON: `{json_path}`",
            ]
        )
        + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="c3_openmc_material_stack")
    parser.add_argument("--front-wall-thickness", type=float, default=0.005)
    parser.add_argument("--flibe-thickness", type=float, default=0.020)
    parser.add_argument("--back-wall-thickness", type=float, default=0.005)
    parser.add_argument("--half-width", type=float, default=0.008)
    parser.add_argument("--half-height", type=float, default=0.008)
    parser.add_argument("--bins", type=int, default=30)
    parser.add_argument("--particles", type=int, default=4000)
    parser.add_argument("--batches", type=int, default=8)
    parser.add_argument("--li6-fraction", type=float, default=0.90)
    parser.add_argument("--target-power", type=float, default=C3_REFERENCE_POWER_W)
    parser.add_argument("--target-heat-flux", type=float, default=C3_REFERENCE_HEAT_FLUX_W_M2)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    work = Path(__file__).resolve().parent
    case = work / args.case
    if case.exists() and args.overwrite:
        shutil.rmtree(case)
    case.mkdir(parents=True, exist_ok=True)

    geometry_meta = build_model(
        case,
        front_wall_thickness=args.front_wall_thickness,
        flibe_thickness=args.flibe_thickness,
        back_wall_thickness=args.back_wall_thickness,
        half_width=args.half_width,
        half_height=args.half_height,
        bins=args.bins,
        particles=args.particles,
        batches=args.batches,
        li6_fraction=args.li6_fraction,
    )

    manifest = {
        "case": str(case),
        "cross_sections": XS_PATH,
        "openmc_exec": pick_openmc_exec(),
        "front_wall_thickness_m": args.front_wall_thickness,
        "flibe_thickness_m": args.flibe_thickness,
        "back_wall_thickness_m": args.back_wall_thickness,
        "half_width_m": args.half_width,
        "half_height_m": args.half_height,
        "bins": args.bins,
        "particles": args.particles,
        "batches": args.batches,
        "li6_fraction": args.li6_fraction,
        "target_power_W": args.target_power,
        "target_heat_flux_W_m2": args.target_heat_flux,
        "run_requested": args.run,
        **geometry_meta,
    }
    (case / "c3_openmc_material_stack_manifest.json").write_text(json.dumps(manifest, indent=2))

    if args.run:
        openmc = import_openmc()
        openmc.run(cwd=str(case), openmc_exec=pick_openmc_exec(), output=True)
        summary = postprocess(
            case,
            total_thickness=geometry_meta["total_thickness_m"],
            region_bounds=geometry_meta["region_bounds_m"],
            bins=args.bins,
            target_power=args.target_power,
            target_heat_flux=args.target_heat_flux,
            out_prefix=work / "figs/c3_openmc_material_stack_heating_split",
        )
        manifest["summary"] = summary
        (case / "c3_openmc_material_stack_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

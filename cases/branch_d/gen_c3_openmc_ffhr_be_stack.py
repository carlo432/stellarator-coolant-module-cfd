#!/usr/bin/env python3
"""Generate and run a tiny C3 FFHR-like Be/FLiBe material-stack tally.

This is a local source-placement artifact, not a reactor neutronics model.  It
extends the earlier steel/FLiBe/steel stack with the Yamanishi-style breeding
zone scaffold: front FLiBe, Be multiplier proxy, and rear FLiBe. The tallied
region powers are aggregated back to the two CHT source classes used by the
current OpenFOAM model: fluid/FLiBe and solid/wall+Be.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
from pathlib import Path

import numpy as np

from gen_c3_openmc_material_stack import (
    C3_REFERENCE_HEAT_FLUX_W_M2,
    C3_REFERENCE_POWER_W,
    EV_TO_J,
    XS_PATH,
    import_openmc,
    latest_statepoint,
    make_flibe,
    make_steel_proxy,
    pick_openmc_exec,
    relative_uncertainty,
)


REGION_NAMES = (
    "front_wall",
    "front_flibe",
    "be_multiplier_proxy",
    "rear_flibe",
    "back_wall",
)
REGION_MATERIAL_CLASSES = {
    "front_wall": "solid",
    "front_flibe": "fluid",
    "be_multiplier_proxy": "solid",
    "rear_flibe": "fluid",
    "back_wall": "solid",
}


def make_be_proxy(openmc, name: str):
    mat = openmc.Material(name=name)
    mat.set_density("g/cm3", 1.85)
    mat.add_nuclide("Be9", 1.0, "ao")
    return mat


def positive_shape(values: np.ndarray) -> tuple[np.ndarray, float]:
    positive = np.maximum(values, 0.0)
    total = float(np.sum(positive))
    if total <= 0.0:
        raise ValueError("OpenMC tally produced no positive heating values")
    return positive / total, total


def build_model(
    case: Path,
    *,
    front_wall_thickness: float,
    front_flibe_thickness: float,
    be_thickness: float,
    rear_flibe_thickness: float,
    back_wall_thickness: float,
    half_width: float,
    half_height: float,
    bins: int,
    particles: int,
    batches: int,
    li6_fraction: float,
) -> dict:
    openmc = import_openmc()

    materials_by_region = {
        "front_wall": make_steel_proxy(openmc, "front_wall_steel_proxy"),
        "front_flibe": make_flibe(openmc, li6_fraction),
        "be_multiplier_proxy": make_be_proxy(openmc, "be_multiplier_proxy"),
        "rear_flibe": make_flibe(openmc, li6_fraction),
        "back_wall": make_steel_proxy(openmc, "back_wall_steel_proxy"),
    }
    materials_by_region["front_flibe"].name = "front_flibe_2LiF_BeF2_li6_enriched"
    materials_by_region["rear_flibe"].name = "rear_flibe_2LiF_BeF2_li6_enriched"
    openmc.Materials(list(materials_by_region.values())).export_to_xml(case / "materials.xml")

    thicknesses = {
        "front_wall": front_wall_thickness,
        "front_flibe": front_flibe_thickness,
        "be_multiplier_proxy": be_thickness,
        "rear_flibe": rear_flibe_thickness,
        "back_wall": back_wall_thickness,
    }
    x_values = [0.0]
    for region in REGION_NAMES:
        x_values.append(x_values[-1] + thicknesses[region])

    planes = []
    for i, x_val in enumerate(x_values):
        boundary = "vacuum" if i in (0, len(x_values) - 1) else "transmission"
        planes.append(openmc.XPlane(x0=x_val, boundary_type=boundary))
    ymin = openmc.YPlane(y0=-half_width, boundary_type="reflective")
    ymax = openmc.YPlane(y0=half_width, boundary_type="reflective")
    zmin = openmc.ZPlane(z0=-half_height, boundary_type="reflective")
    zmax = openmc.ZPlane(z0=half_height, boundary_type="reflective")
    yz_region = +ymin & -ymax & +zmin & -zmax

    cells = []
    region_bounds = {}
    for idx, region in enumerate(REGION_NAMES):
        region_bounds[region] = [x_values[idx], x_values[idx + 1]]
        cells.append(
            openmc.Cell(
                name=region,
                fill=materials_by_region[region],
                region=+planes[idx] & -planes[idx + 1] & yz_region,
            )
        )
    openmc.Geometry(cells).export_to_xml(case / "geometry.xml")

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = 0
    settings.seed = 170707
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

    total_thickness = x_values[-1]
    mesh = openmc.RegularMesh()
    mesh.dimension = (bins, 1, 1)
    mesh.lower_left = (0.0, -half_width, -half_height)
    mesh.upper_right = (total_thickness, half_width, half_height)
    mesh_filter = openmc.MeshFilter(mesh)
    cell_filter = openmc.CellFilter(cells)

    heat_cell = openmc.Tally(name="heating_local_by_region")
    heat_cell.filters = [cell_filter]
    heat_cell.scores = ["heating-local"]

    heat_depth = openmc.Tally(name="heating_local_by_depth")
    heat_depth.filters = [mesh_filter]
    heat_depth.scores = ["heating-local"]

    heat_cell_depth = openmc.Tally(name="heating_local_by_region_depth")
    heat_cell_depth.filters = [cell_filter, mesh_filter]
    heat_cell_depth.scores = ["heating-local"]

    openmc.Tallies([heat_cell, heat_depth, heat_cell_depth]).export_to_xml(case / "tallies.xml")

    return {
        "total_thickness_m": total_thickness,
        "region_bounds_m": region_bounds,
        "cell_ids": {cell.name: cell.id for cell in cells},
        "material_ids": {name: mat.id for name, mat in materials_by_region.items()},
        "thicknesses_m": thicknesses,
    }


def write_region_outputs(
    out_prefix: Path,
    *,
    region_bounds: dict,
    cell_mean: np.ndarray,
    cell_sd: np.ndarray,
    cell_shape: np.ndarray,
    target_power: float,
) -> tuple[list[dict], list[dict], Path, Path]:
    rel_unc = relative_uncertainty(cell_mean, cell_sd)
    rows = []
    for idx, region in enumerate(REGION_NAMES):
        x0, x1 = region_bounds[region]
        rows.append(
            {
                "region": region,
                "material_class": REGION_MATERIAL_CLASSES[region],
                "x0_m": x0,
                "x1_m": x1,
                "thickness_m": x1 - x0,
                "heating_mean_eV_per_source": float(cell_mean[idx]),
                "relative_uncertainty": float(rel_unc[idx]),
                "normalized_fraction": float(cell_shape[idx]),
                "normalized_power_W": float(cell_shape[idx] * target_power),
            }
        )

    class_rows = []
    for klass in ("fluid", "solid"):
        selected = [row for row in rows if row["material_class"] == klass]
        class_rows.append(
            {
                "material_class": klass,
                "normalized_fraction": sum(row["normalized_fraction"] for row in selected),
                "normalized_power_W": sum(row["normalized_power_W"] for row in selected),
                "regions": "+".join(row["region"] for row in selected),
            }
        )

    region_csv = out_prefix.with_name(out_prefix.name + "_region_split.csv")
    with region_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    class_csv = out_prefix.with_name(out_prefix.name + "_class_split.csv")
    with class_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(class_rows[0].keys()))
        writer.writeheader()
        writer.writerows(class_rows)

    return rows, class_rows, region_csv, class_csv


def write_depth_outputs(
    out_prefix: Path,
    *,
    total_thickness: float,
    depth_mean: np.ndarray,
    depth_sd: np.ndarray,
    depth_shape: np.ndarray,
    target_heat_flux: float,
) -> tuple[list[dict], Path]:
    bins = len(depth_mean)
    dx = total_thickness / bins
    rel_unc = relative_uncertainty(depth_mean, depth_sd)
    rows = []
    for idx in range(bins):
        rows.append(
            {
                "bin": idx,
                "x0_m": idx * dx,
                "x1_m": (idx + 1) * dx,
                "x_mid_m": (idx + 0.5) * dx,
                "heating_mean_eV_per_source": float(depth_mean[idx]),
                "relative_uncertainty": float(rel_unc[idx]),
                "normalized_shape": float(depth_shape[idx]),
                "normalized_qvol_W_m3": float(depth_shape[idx] * target_heat_flux / dx),
            }
        )

    depth_csv = out_prefix.with_name(out_prefix.name + "_depth_profile.csv")
    with depth_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows, depth_csv


def write_region_depth_output(
    out_prefix: Path,
    *,
    total_thickness: float,
    region_depth_mean: np.ndarray,
    region_depth_sd: np.ndarray,
    total_heating_eV: float,
    target_heat_flux: float,
) -> Path:
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

    region_depth_csv = out_prefix.with_name(out_prefix.name + "_region_depth.csv")
    with region_depth_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return region_depth_csv


def plot_outputs(out_prefix: Path, region_rows: list[dict], depth_rows: list[dict]) -> str | None:
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        x_mid = np.array([row["x_mid_m"] for row in depth_rows])
        qvol = np.array([row["normalized_qvol_W_m3"] for row in depth_rows]) / 1.0e6
        names = [row["region"] for row in region_rows]
        powers = [row["normalized_power_W"] for row in region_rows]
        colors = ["#5c6670", "#236a8c", "#b99b47", "#2d7d9a", "#7a6650"]

        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
        axes[0].bar(names, powers, color=colors)
        axes[0].set_ylabel("Normalized power (W)")
        axes[0].set_title("C3 FFHR-like stack heating split")
        axes[0].tick_params(axis="x", rotation=24)

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
    cell_shape, cell_total = positive_shape(cell_mean)

    depth_tally = sp.get_tally(name="heating_local_by_depth")
    depth_mean = np.asarray(depth_tally.mean).reshape((bins,))
    depth_sd = np.asarray(depth_tally.std_dev).reshape((bins,))
    depth_shape, depth_total = positive_shape(depth_mean)

    region_depth_tally = sp.get_tally(name="heating_local_by_region_depth")
    region_depth_mean = np.asarray(region_depth_tally.mean).reshape((len(REGION_NAMES), bins))
    region_depth_sd = np.asarray(region_depth_tally.std_dev).reshape((len(REGION_NAMES), bins))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    region_rows, class_rows, region_csv, class_csv = write_region_outputs(
        out_prefix,
        region_bounds=region_bounds,
        cell_mean=cell_mean,
        cell_sd=cell_sd,
        cell_shape=cell_shape,
        target_power=target_power,
    )
    depth_rows, depth_csv = write_depth_outputs(
        out_prefix,
        total_thickness=total_thickness,
        depth_mean=depth_mean,
        depth_sd=depth_sd,
        depth_shape=depth_shape,
        target_heat_flux=target_heat_flux,
    )
    region_depth_csv = write_region_depth_output(
        out_prefix,
        total_thickness=total_thickness,
        region_depth_mean=region_depth_mean,
        region_depth_sd=region_depth_sd,
        total_heating_eV=cell_total,
        target_heat_flux=target_heat_flux,
    )

    fluid = next(row for row in class_rows if row["material_class"] == "fluid")
    solid = next(row for row in class_rows if row["material_class"] == "solid")
    summary = {
        "case": str(case),
        "statepoint": str(latest_statepoint(case)),
        "status": "C3 FFHR-like Be/FLiBe stack OpenMC source-placement tally; not reactor neutronics",
        "total_thickness_m": total_thickness,
        "bins": bins,
        "target_power_W": target_power,
        "target_heat_flux_W_m2": target_heat_flux,
        "reference_area_m2": target_power / target_heat_flux,
        "cell_total_heating_eV_per_source": cell_total,
        "depth_total_heating_eV_per_source": depth_total,
        "cell_total_heating_J_per_source": cell_total * EV_TO_J,
        "fluid_fraction": fluid["normalized_fraction"],
        "solid_fraction": solid["normalized_fraction"],
        "fluid_power_W": fluid["normalized_power_W"],
        "solid_power_W": solid["normalized_power_W"],
        "max_region_relative_uncertainty": max(row["relative_uncertainty"] for row in region_rows),
        "max_depth_relative_uncertainty": max(row["relative_uncertainty"] for row in depth_rows),
        "region_csv": str(region_csv),
        "class_csv": str(class_csv),
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
                "# C3 OpenMC FFHR-Like Be/FLiBe Stack Heating Split",
                "",
                "Status: local source-placement tally, not reactor neutronics.",
                "",
                "## Model",
                "",
                "- geometry: steel-proxy wall, front FLiBe, Be multiplier proxy, rear FLiBe, steel-proxy wall",
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
                f"- fluid/FLiBe: `{fluid['normalized_power_W']:.5f} W` (`{fluid['normalized_fraction']:.5f}`)",
                f"- solid/wall+Be: `{solid['normalized_power_W']:.5f} W` (`{solid['normalized_fraction']:.5f}`)",
                "",
                "## Interpretation",
                "",
                "This is the first C3 stack that includes the Be layer from the FFHR",
                "literature anchor. It is still only a CSG source-placement toy: the Be",
                "layer is treated as a solid proxy, FLiBe-in-pebble gaps are not resolved,",
                "and the result is aggregated back into the two-region CHT model.",
                "",
                f"- region split CSV: `{region_csv}`",
                f"- class split CSV: `{class_csv}`",
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
    parser.add_argument("--case", default="c3_openmc_ffhr_be_stack")
    parser.add_argument("--front-wall-thickness", type=float, default=0.005)
    parser.add_argument("--front-flibe-thickness", type=float, default=0.020)
    parser.add_argument("--be-thickness", type=float, default=0.120)
    parser.add_argument("--rear-flibe-thickness", type=float, default=0.140)
    parser.add_argument("--back-wall-thickness", type=float, default=0.005)
    parser.add_argument("--half-width", type=float, default=0.008)
    parser.add_argument("--half-height", type=float, default=0.008)
    parser.add_argument("--bins", type=int, default=58)
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
        front_flibe_thickness=args.front_flibe_thickness,
        be_thickness=args.be_thickness,
        rear_flibe_thickness=args.rear_flibe_thickness,
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
        "front_flibe_thickness_m": args.front_flibe_thickness,
        "be_thickness_m": args.be_thickness,
        "rear_flibe_thickness_m": args.rear_flibe_thickness,
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
    manifest_path = case / "c3_openmc_ffhr_be_stack_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

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
            out_prefix=work / "figs/c3_openmc_ffhr_be_stack_heating_split",
        )
        manifest["summary"] = summary
        manifest_path.write_text(json.dumps(manifest, indent=2))

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

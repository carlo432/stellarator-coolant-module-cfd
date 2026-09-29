#!/usr/bin/env python3
"""Generate and optionally run a tiny C3 OpenMC FLiBe slab heating tally.

This is a C3 smoke/taste artifact, not a reactor neutronics model.  The tally
shape is normalized to a chosen surface heat-flux scale before being handed to
the thermal side.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
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


def pick_openmc_exec() -> str:
    for candidate in OPENMC_EXE_CANDIDATES:
        if candidate == "openmc" or Path(candidate).exists():
            return candidate
    return "openmc"


def import_openmc():
    os.environ["OPENMC_CROSS_SECTIONS"] = XS_PATH
    import openmc  # noqa: PLC0415

    return openmc


def build_model(
    case: Path,
    *,
    slab_thickness: float,
    half_width: float,
    half_height: float,
    bins: int,
    particles: int,
    batches: int,
    li6_fraction: float,
) -> None:
    openmc = import_openmc()

    flibe = openmc.Material(name="flibe_2LiF_BeF2_li6_enriched")
    flibe.set_density("g/cm3", 1.94)
    # 2LiF-BeF2 atom ratios: Li:Be:F = 2:1:4.
    flibe.add_nuclide("Li6", 2.0 * li6_fraction, "ao")
    flibe.add_nuclide("Li7", 2.0 * (1.0 - li6_fraction), "ao")
    flibe.add_nuclide("Be9", 1.0, "ao")
    flibe.add_nuclide("F19", 4.0, "ao")
    openmc.Materials([flibe]).export_to_xml(case / "materials.xml")

    x0 = openmc.XPlane(x0=0.0, boundary_type="vacuum")
    x1 = openmc.XPlane(x0=slab_thickness, boundary_type="vacuum")
    ymin = openmc.YPlane(y0=-half_width, boundary_type="reflective")
    ymax = openmc.YPlane(y0=half_width, boundary_type="reflective")
    zmin = openmc.ZPlane(z0=-half_height, boundary_type="reflective")
    zmax = openmc.ZPlane(z0=half_height, boundary_type="reflective")
    slab = openmc.Cell(fill=flibe, region=+x0 & -x1 & +ymin & -ymax & +zmin & -zmax)
    openmc.Geometry([slab]).export_to_xml(case / "geometry.xml")

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = 0
    settings.seed = 170705
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

    mesh = openmc.RegularMesh()
    mesh.dimension = (bins, 1, 1)
    mesh.lower_left = (0.0, -half_width, -half_height)
    mesh.upper_right = (slab_thickness, half_width, half_height)
    mesh_filter = openmc.MeshFilter(mesh)

    heat = openmc.Tally(name="heating_local_by_depth")
    heat.filters = [mesh_filter]
    heat.scores = ["heating-local"]

    flux = openmc.Tally(name="flux_by_depth")
    flux.filters = [mesh_filter]
    flux.scores = ["flux"]

    openmc.Tallies([heat, flux]).export_to_xml(case / "tallies.xml")


def latest_statepoint(case: Path) -> Path:
    statepoints = sorted(case.glob("statepoint.*.h5"))
    if not statepoints:
        raise FileNotFoundError(f"no statepoint files in {case}")
    return statepoints[-1]


def postprocess(
    case: Path,
    *,
    slab_thickness: float,
    bins: int,
    target_heat_flux: float,
    out_prefix: Path,
) -> dict:
    openmc = import_openmc()
    sp = openmc.StatePoint(latest_statepoint(case))
    heat_tally = sp.get_tally(name="heating_local_by_depth")
    heat_mean = np.asarray(heat_tally.mean).reshape((bins,))
    heat_sd = np.asarray(heat_tally.std_dev).reshape((bins,))

    source_kind = "heating-local"
    shape_raw = np.maximum(heat_mean, 0.0)
    if not np.any(shape_raw > 0.0):
        flux_tally = sp.get_tally(name="flux_by_depth")
        shape_raw = np.maximum(np.asarray(flux_tally.mean).reshape((bins,)), 0.0)
        heat_sd = np.zeros_like(shape_raw)
        source_kind = "flux-proxy"

    total = float(np.sum(shape_raw))
    if total <= 0.0:
        # Last-resort deterministic attenuation shape so downstream scripts
        # still have a sane file. This should not trigger in a successful run.
        mids = np.linspace(0.5 / bins, 1.0 - 0.5 / bins, bins)
        shape_raw = np.exp(-3.0 * mids)
        total = float(np.sum(shape_raw))
        source_kind = "deterministic-fallback"

    shape = shape_raw / total
    dx = slab_thickness / bins
    qvol = shape * target_heat_flux / dx
    mids = (np.arange(bins) + 0.5) * dx
    rel_unc = np.divide(heat_sd, np.abs(heat_mean), out=np.zeros_like(heat_sd), where=np.abs(heat_mean) > 0)
    eV_per_source = float(np.sum(np.maximum(heat_mean, 0.0)))
    joule_per_source = eV_per_source * EV_TO_J
    peak_to_mean = float(np.max(qvol) / np.mean(qvol))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_prefix.with_suffix(".csv")
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "bin",
                "x0_m",
                "x1_m",
                "x_mid_m",
                "heating_mean_eV_per_source",
                "relative_uncertainty",
                "normalized_shape",
                "normalized_qvol_W_m3",
            ]
        )
        for i in range(bins):
            writer.writerow(
                [
                    i,
                    i * dx,
                    (i + 1) * dx,
                    mids[i],
                    float(heat_mean[i]),
                    float(rel_unc[i]),
                    float(shape[i]),
                    float(qvol[i]),
                ]
            )

    summary = {
        "case": str(case),
        "statepoint": str(latest_statepoint(case)),
        "source_kind": source_kind,
        "slab_thickness_m": slab_thickness,
        "bins": bins,
        "target_heat_flux_W_m2": target_heat_flux,
        "bin_width_m": dx,
        "openmc_total_heating_eV_per_source": eV_per_source,
        "openmc_total_heating_J_per_source": joule_per_source,
        "normalized_qvol_mean_W_m3": float(np.mean(qvol)),
        "normalized_qvol_peak_W_m3": float(np.max(qvol)),
        "normalized_qvol_peak_to_mean": peak_to_mean,
        "front_half_power_fraction": float(np.sum(shape[: max(1, bins // 2)])),
        "max_relative_uncertainty": float(np.max(rel_unc)),
        "csv": str(csv_path),
    }
    json_path = out_prefix.with_suffix(".json")
    json_path.write_text(json.dumps(summary, indent=2))

    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        ax.step(mids, qvol / 1.0e6, where="mid", color="#255f85", linewidth=2.0)
        ax.set_xlabel("Depth into FLiBe slab (m)")
        ax.set_ylabel("Normalized volumetric source (MW/m^3)")
        ax.set_title("C3 OpenMC FLiBe slab heating shape")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(out_prefix.with_suffix(".png"), dpi=180)
        plt.close(fig)
        summary["png"] = str(out_prefix.with_suffix(".png"))
        json_path.write_text(json.dumps(summary, indent=2))
    except Exception as exc:  # pragma: no cover - plotting is best-effort
        summary["plot_error"] = str(exc)
        json_path.write_text(json.dumps(summary, indent=2))

    md_path = out_prefix.with_suffix(".md")
    md_path.write_text(
        "\n".join(
            [
                "# C3 OpenMC FLiBe Slab Heating Profile",
                "",
                "Status: smoke/taste, not reactor neutronics.",
                "",
                f"- OpenMC source kind used for shape: `{source_kind}`",
                f"- slab thickness: `{slab_thickness:g} m`",
                f"- mesh bins: `{bins}`",
                f"- target thermal normalization: `{target_heat_flux:g} W/m^2`",
                f"- total OpenMC heating: `{eV_per_source:.6e} eV/source`",
                f"- normalized mean q''': `{np.mean(qvol):.6e} W/m^3`",
                f"- normalized peak q''': `{np.max(qvol):.6e} W/m^3`",
                f"- peak/mean: `{peak_to_mean:.3f}`",
                f"- front-half power fraction: `{summary['front_half_power_fraction']:.3f}`",
                f"- max relative uncertainty: `{summary['max_relative_uncertainty']:.3f}`",
                "",
                "Interpretation: this creates a one-way neutronics-to-thermal",
                "heating-shape file. Absolute reactor normalization is not claimed;",
                "the profile is normalized to the same `100 kW/m^2` scale used by",
                "the local C2 thermal comparison.",
                "",
                f"- CSV: `{csv_path}`",
                f"- JSON: `{json_path}`",
            ]
        )
        + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="c3_openmc_flibe_slab")
    parser.add_argument("--slab-thickness", type=float, default=0.20)
    parser.add_argument("--half-width", type=float, default=0.008)
    parser.add_argument("--half-height", type=float, default=0.008)
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument("--particles", type=int, default=5000)
    parser.add_argument("--batches", type=int, default=10)
    parser.add_argument("--li6-fraction", type=float, default=0.90)
    parser.add_argument("--target-heat-flux", type=float, default=100000.0)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    work = Path(__file__).resolve().parent
    case = work / args.case
    if case.exists() and args.overwrite:
        shutil.rmtree(case)
    case.mkdir(parents=True, exist_ok=True)

    build_model(
        case,
        slab_thickness=args.slab_thickness,
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
        "slab_thickness_m": args.slab_thickness,
        "half_width_m": args.half_width,
        "half_height_m": args.half_height,
        "bins": args.bins,
        "particles": args.particles,
        "batches": args.batches,
        "li6_fraction": args.li6_fraction,
        "target_heat_flux_W_m2": args.target_heat_flux,
        "run_requested": args.run,
    }
    (case / "c3_openmc_manifest.json").write_text(json.dumps(manifest, indent=2))

    if args.run:
        openmc = import_openmc()
        openmc.run(cwd=str(case), openmc_exec=pick_openmc_exec(), output=True)
        postprocess(
            case,
            slab_thickness=args.slab_thickness,
            bins=args.bins,
            target_heat_flux=args.target_heat_flux,
            out_prefix=work / "figs/c3_openmc_flibe_slab_heating_profile",
        )

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

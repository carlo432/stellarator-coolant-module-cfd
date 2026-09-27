#!/usr/bin/env python3
"""Create Branch D v27 RANS-vs-LES mean-field and resolved-TKE figures."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RANS_TIME = ROOT / "cases/baseline_module/openfoam_cases/v22_isothermal_rans_flibe_re10000_bent_recirculation/1000"
LES_TIME = ROOT / "cases/branch_d/v27_thermal_les_flibe_re10000_bent_recirculation_wale_100kw/0.500073"
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"


def _internal_block(text: str) -> tuple[str, int, str]:
    match = re.search(
        r"internalField\s+nonuniform\s+List<(?P<kind>\w+)>\s+"
        r"(?P<count>\d+)\s*\((?P<body>.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError("Could not find nonuniform internalField block")
    return match.group("kind"), int(match.group("count")), match.group("body")


def read_scalar_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "scalar":
        raise ValueError(f"{path} is {kind}, not scalar")
    data = np.fromstring(body, sep=" ")
    if data.size != count:
        raise ValueError(f"{path}: expected {count}, got {data.size}")
    return data


def read_vector_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "vector":
        raise ValueError(f"{path} is {kind}, not vector")
    rows = re.findall(r"\(([^()]+)\)", body)
    data = np.array([[float(x) for x in row.split()] for row in rows], dtype=float)
    if data.shape != (count, 3):
        raise ValueError(f"{path}: expected {(count, 3)}, got {data.shape}")
    return data


def read_symmtensor_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "symmTensor":
        raise ValueError(f"{path} is {kind}, not symmTensor")
    rows = re.findall(r"\(([^()]+)\)", body)
    data = np.array([[float(x) for x in row.split()] for row in rows], dtype=float)
    if data.shape != (count, 6):
        raise ValueError(f"{path}: expected {(count, 6)}, got {data.shape}")
    return data


def nearest_midplane_indices(centers: np.ndarray, z_value: float = 0.02) -> tuple[np.ndarray, float]:
    z = centers[:, 2]
    unique_z = np.unique(np.round(z, 9))
    actual = float(unique_z[np.argmin(np.abs(unique_z - z_value))])
    indices = np.where(np.abs(z - actual) <= 1.0e-9)[0]
    if indices.size == 0:
        raise ValueError(f"No cells found for z={actual}")
    return indices, actual


def zone_mask(centers: np.ndarray, x_min: float, x_max: float, y_max: float) -> np.ndarray:
    return (
        (centers[:, 0] >= x_min)
        & (centers[:, 0] <= x_max)
        & (centers[:, 1] <= y_max)
    )


def stats(values: np.ndarray) -> tuple[float, float, float]:
    return float(np.mean(values)), float(np.max(values)), float(np.min(values))


def write_metrics(
    path: Path,
    centers: np.ndarray,
    rans_u: np.ndarray,
    les_u: np.ndarray,
    tke: np.ndarray,
    tprime_rms: np.ndarray,
) -> None:
    speed_rans = np.linalg.norm(rans_u, axis=1)
    speed_les = np.linalg.norm(les_u, axis=1)
    speed_diff = speed_les - speed_rans
    vector_diff = np.linalg.norm(les_u - rans_u, axis=1)
    ux_diff = les_u[:, 0] - rans_u[:, 0]

    zones = {
        "inlet_duct": (centers[:, 0] < 0.10),
        "chamber": (centers[:, 0] >= 0.10) & (centers[:, 0] <= 0.30),
        "near_heated_wall": zone_mask(centers, 0.10, 0.30, -0.02),
        "hotspot_zone": zone_mask(centers, 0.22, 0.30, -0.02),
        "downstream_hotspot_band": zone_mask(centers, 0.25, 0.27, -0.02),
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "zone",
                "cell_count",
                "rans_mean_Ux_m_s",
                "les_mean_Ux_m_s",
                "mean_Ux_difference_m_s",
                "rans_mean_speed_m_s",
                "les_mean_speed_m_s",
                "mean_speed_difference_m_s",
                "mean_vector_difference_m_s",
                "max_vector_difference_m_s",
                "mean_resolved_tke_m2_s2",
                "max_resolved_tke_m2_s2",
                "mean_Tprime_rms_K",
                "max_Tprime_rms_K",
            ]
        )
        for name, mask in zones.items():
            if not np.any(mask):
                continue
            writer.writerow(
                [
                    name,
                    int(np.count_nonzero(mask)),
                    float(np.mean(rans_u[mask, 0])),
                    float(np.mean(les_u[mask, 0])),
                    float(np.mean(ux_diff[mask])),
                    float(np.mean(speed_rans[mask])),
                    float(np.mean(speed_les[mask])),
                    float(np.mean(speed_diff[mask])),
                    float(np.mean(vector_diff[mask])),
                    float(np.max(vector_diff[mask])),
                    float(np.mean(tke[mask])),
                    float(np.max(tke[mask])),
                    float(np.mean(tprime_rms[mask])),
                    float(np.max(tprime_rms[mask])),
                ]
            )


def tricontour_panel(ax, triangulation, values, title, cmap, label, levels=36, vmin=None, vmax=None):
    contour = ax.tricontourf(triangulation, values, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.tricontour(triangulation, values, levels=min(12, levels), colors="k", linewidths=0.18, alpha=0.22)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(True, color="white", linewidth=0.3, alpha=0.25)
    cbar = plt.colorbar(contour, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label(label, fontsize=9)


def annotate_hotspot(ax):
    ax.plot([0.26], [-0.035], marker="x", color="white", markersize=7, markeredgewidth=1.6)
    ax.text(0.252, -0.030, "hotspot", color="white", fontsize=8, va="bottom", ha="right")


def make_mean_field_figure(centers, indices, rans_u, les_u):
    xy = centers[indices, :2]
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1])

    rans_speed = np.linalg.norm(rans_u[indices], axis=1)
    les_speed = np.linalg.norm(les_u[indices], axis=1)
    vec_diff = np.linalg.norm((les_u - rans_u)[indices], axis=1)

    vmax = max(float(np.percentile(rans_speed, 99)), float(np.percentile(les_speed, 99)))

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.2), dpi=180)
    tricontour_panel(axes[0], tri, rans_speed, "RANS v22 |U|", "viridis", "|U| [m/s]", vmin=0, vmax=vmax)
    tricontour_panel(axes[1], tri, les_speed, "LES v27 time-mean |UMean|", "viridis", "|UMean| [m/s]", vmin=0, vmax=vmax)
    tricontour_panel(axes[2], tri, vec_diff, "Mean-field vector difference", "magma", "|UMean - U_RANS| [m/s]")
    for ax in axes:
        annotate_hotspot(ax)
    fig.suptitle("Branch D RANS-vs-LES Mean Velocity Comparison, z=0.02 m", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "LES is coarse, wall-modeled, short-sample, and illustrative; comparison uses v27 fieldAverage at t=0.500073.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    out = OUT_DIR / "v27_rans_vs_les_mean_velocity_midplane.png"
    fig.savefig(out, bbox_inches="tight")
    advisor = ADVISOR_DIR / "22_v27_rans_vs_les_mean_velocity_midplane.png"
    fig.savefig(advisor, bbox_inches="tight")
    plt.close(fig)
    return out, advisor


def make_tke_figure(centers, indices, tke, tprime_rms):
    xy = centers[indices, :2]
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1])

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), dpi=180)
    tricontour_panel(
        axes[0],
        tri,
        tke[indices],
        "Resolved TKE from UPrime2Mean",
        "plasma",
        "0.5 tr(UPrime2Mean) [m2/s2]",
    )
    tricontour_panel(
        axes[1],
        tri,
        tprime_rms[indices],
        "Thermal fluctuation RMS",
        "inferno",
        "sqrt(TPrime2Mean) [K]",
    )
    for ax in axes:
        annotate_hotspot(ax)
    fig.suptitle("Branch D v27 LES Resolved Turbulence and Thermal Fluctuation, z=0.02 m", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "Coarse wall-modeled thermal LES smoke test: short, developing, not grid-converged, and not wall-surface sampling.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.92))
    out = OUT_DIR / "v27_resolved_tke_temperature_fluctuation_midplane.png"
    fig.savefig(out, bbox_inches="tight")
    advisor = ADVISOR_DIR / "23_v27_resolved_tke_temperature_fluctuation_midplane.png"
    fig.savefig(advisor, bbox_inches="tight")
    plt.close(fig)
    return out, advisor


def write_slice_csv(path, centers, indices, rans_u, les_u, tke, tprime_rms):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "cell",
                "x",
                "y",
                "z",
                "rans_Ux",
                "rans_Uy",
                "rans_Uz",
                "rans_Umag",
                "les_UMean_x",
                "les_UMean_y",
                "les_UMean_z",
                "les_UMean_mag",
                "mean_vector_difference",
                "resolved_tke",
                "Tprime_rms",
            ]
        )
        for idx in indices:
            ru = rans_u[idx]
            lu = les_u[idx]
            writer.writerow(
                [
                    int(idx),
                    *centers[idx],
                    *ru,
                    float(np.linalg.norm(ru)),
                    *lu,
                    float(np.linalg.norm(lu)),
                    float(np.linalg.norm(lu - ru)),
                    float(tke[idx]),
                    float(tprime_rms[idx]),
                ]
            )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ADVISOR_DIR.mkdir(parents=True, exist_ok=True)

    rans_c = read_vector_field(RANS_TIME / "C")
    rans_u = read_vector_field(RANS_TIME / "U")
    les_c = read_vector_field(LES_TIME / "C")
    les_u = read_vector_field(LES_TIME / "UMean")
    uprime = read_symmtensor_field(LES_TIME / "UPrime2Mean")
    tprime2 = read_scalar_field(LES_TIME / "TPrime2Mean")

    if rans_c.shape != les_c.shape:
        raise SystemExit(f"RANS/LES cell count mismatch: {rans_c.shape} vs {les_c.shape}")
    max_center_delta = float(np.max(np.linalg.norm(rans_c - les_c, axis=1)))
    if max_center_delta > 1.0e-10:
        raise SystemExit(f"RANS/LES cell centers are not aligned; max delta={max_center_delta}")

    # OpenFOAM symmTensor order is xx, xy, xz, yy, yz, zz.
    tke = 0.5 * (uprime[:, 0] + uprime[:, 3] + uprime[:, 5])
    tprime_rms = np.sqrt(np.maximum(tprime2, 0.0))

    indices, actual_z = nearest_midplane_indices(rans_c)
    write_metrics(OUT_DIR / "v27_rans_vs_les_mean_tke_metrics.csv", rans_c, rans_u, les_u, tke, tprime_rms)
    write_slice_csv(OUT_DIR / "v27_rans_vs_les_midplane_values.csv", rans_c, indices, rans_u, les_u, tke, tprime_rms)

    mean_out, mean_advisor = make_mean_field_figure(rans_c, indices, rans_u, les_u)
    tke_out, tke_advisor = make_tke_figure(rans_c, indices, tke, tprime_rms)

    print(f"actual_z={actual_z}")
    print(f"max_center_delta={max_center_delta}")
    for path in [
        OUT_DIR / "v27_rans_vs_les_mean_tke_metrics.csv",
        OUT_DIR / "v27_rans_vs_les_midplane_values.csv",
        mean_out,
        mean_advisor,
        tke_out,
        tke_advisor,
    ]:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()

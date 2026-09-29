#!/usr/bin/env python3
"""C5 sampled mean/TKE comparison from retained LES snapshots.

The C5 BFS run did not retain OpenFOAM fieldAverage UMean/UPrime2Mean fields.
This script computes a matched, caveated snapshot-average from retained
post-0.5 s velocity fields for BFS and outlet-offset v26.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "cases" / "branch_d"
FIGS = WORK / "figs"


@dataclass(frozen=True)
class CaseSpec:
    name: str
    label: str
    case_dir: Path
    center_time: str
    z_target: float
    near_wall_y_max: float
    hotspot_x_min: float
    hotspot_x_max: float


CASES = [
    CaseSpec(
        name="bfs",
        label="BFS",
        case_dir=WORK / "c5_bfs_les_smoke_wale",
        center_time="0.799989",
        z_target=0.02,
        near_wall_y_max=-0.005,
        hotspot_x_min=0.10,
        hotspot_x_max=0.24,
    ),
    CaseSpec(
        name="outlet_offset",
        label="outlet-offset",
        case_dir=WORK / "v26_les_smoke_flibe_re10000_bent_recirculation_wale",
        center_time="0.800034",
        z_target=0.02,
        near_wall_y_max=-0.02,
        hotspot_x_min=0.15,
        hotspot_x_max=0.28,
    ),
]


def internal_block(text: str) -> tuple[str, int, str]:
    match = re.search(
        r"internalField\s+nonuniform\s+List<(?P<kind>\w+)>\s+"
        r"(?P<count>\d+)\s*\((?P<body>.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError("Could not find nonuniform internalField block")
    return match.group("kind"), int(match.group("count")), match.group("body")


def read_vector_field(path: Path) -> np.ndarray:
    kind, count, body = internal_block(path.read_text(encoding="ascii"))
    if kind != "vector":
        raise ValueError(f"{path} is {kind}, not vector")
    rows = re.findall(r"\(([^()]+)\)", body)
    data = np.array([[float(x) for x in row.split()] for row in rows], dtype=float)
    if data.shape != (count, 3):
        raise ValueError(f"{path}: expected {(count, 3)}, got {data.shape}")
    return data


def available_times(case_dir: Path, min_time: float = 0.5) -> list[str]:
    out: list[tuple[float, str]] = []
    for path in case_dir.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value >= min_time and (path / "U").exists():
            out.append((value, path.name))
    return [name for _, name in sorted(out)]


def nearest_midplane_indices(centers: np.ndarray, z_target: float) -> tuple[np.ndarray, float]:
    z_values = np.unique(np.round(centers[:, 2], 9))
    actual = float(z_values[np.argmin(np.abs(z_values - z_target))])
    indices = np.where(np.abs(centers[:, 2] - actual) <= 1.0e-9)[0]
    if indices.size == 0:
        raise RuntimeError(f"no cells at z={actual}")
    return indices, actual


def compute_sampled_stats(spec: CaseSpec) -> dict[str, object]:
    centers = read_vector_field(spec.case_dir / spec.center_time / "C")
    times = available_times(spec.case_dir)
    if len(times) < 2:
        raise RuntimeError(f"{spec.name}: need at least 2 retained U snapshots, got {times}")

    fields = [read_vector_field(spec.case_dir / time / "U") for time in times]
    stack = np.stack(fields, axis=0)
    mean_u = np.mean(stack, axis=0)
    fluct = stack - mean_u[None, :, :]
    var = np.mean(fluct * fluct, axis=0)
    tke = 0.5 * np.sum(var, axis=1)
    speed = np.linalg.norm(mean_u, axis=1)

    indices, actual_z = nearest_midplane_indices(centers, spec.z_target)
    return {
        "spec": spec,
        "centers": centers,
        "times": times,
        "mean_u": mean_u,
        "speed": speed,
        "tke": tke,
        "indices": indices,
        "actual_z": actual_z,
    }


def zone_masks(spec: CaseSpec, centers: np.ndarray) -> dict[str, np.ndarray]:
    x = centers[:, 0]
    y = centers[:, 1]
    return {
        "all_fluid": np.ones(centers.shape[0], dtype=bool),
        "inlet_upstream": x < 0.10,
        "chamber": (x >= 0.10) & (x <= 0.30),
        "near_heated_wall": (x >= 0.10) & (x <= 0.30) & (y <= spec.near_wall_y_max),
        "hotspot_band": (x >= spec.hotspot_x_min)
        & (x <= spec.hotspot_x_max)
        & (y <= spec.near_wall_y_max),
        "downstream": (x >= 0.24) & (x <= 0.32),
    }


def write_metrics(results: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "geometry",
                "snapshot_count",
                "snapshot_times",
                "zone",
                "cell_count",
                "mean_speed_m_s",
                "max_speed_m_s",
                "mean_Ux_m_s",
                "min_Ux_m_s",
                "max_Ux_m_s",
                "mean_sampled_tke_m2_s2",
                "max_sampled_tke_m2_s2",
            ]
        )
        for result in results:
            spec: CaseSpec = result["spec"]  # type: ignore[assignment]
            centers = result["centers"]  # type: ignore[assignment]
            mean_u = result["mean_u"]  # type: ignore[assignment]
            speed = result["speed"]  # type: ignore[assignment]
            tke = result["tke"]  # type: ignore[assignment]
            times = result["times"]  # type: ignore[assignment]
            for zone, mask in zone_masks(spec, centers).items():
                if not np.any(mask):
                    continue
                writer.writerow(
                    [
                        spec.label,
                        len(times),
                        " ".join(times),
                        zone,
                        int(np.count_nonzero(mask)),
                        float(np.mean(speed[mask])),
                        float(np.max(speed[mask])),
                        float(np.mean(mean_u[mask, 0])),
                        float(np.min(mean_u[mask, 0])),
                        float(np.max(mean_u[mask, 0])),
                        float(np.mean(tke[mask])),
                        float(np.max(tke[mask])),
                    ]
                )


def write_nearwall_profile(results: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "geometry",
                "x_center_m",
                "x_min_m",
                "x_max_m",
                "cell_count",
                "mean_Ux_m_s",
                "mean_speed_m_s",
                "reverse_flow_fraction",
                "mean_sampled_tke_m2_s2",
            ]
        )
        for result in results:
            spec: CaseSpec = result["spec"]  # type: ignore[assignment]
            centers = result["centers"]  # type: ignore[assignment]
            mean_u = result["mean_u"]  # type: ignore[assignment]
            speed = result["speed"]  # type: ignore[assignment]
            tke = result["tke"]  # type: ignore[assignment]
            x = centers[:, 0]
            y = centers[:, 1]
            near_wall = (x >= 0.10) & (x <= 0.32) & (y <= spec.near_wall_y_max)
            bins = np.linspace(0.10, 0.32, 23)
            for x_min, x_max in zip(bins[:-1], bins[1:]):
                mask = near_wall & (x >= x_min) & (x < x_max)
                if not np.any(mask):
                    continue
                writer.writerow(
                    [
                        spec.label,
                        0.5 * (x_min + x_max),
                        x_min,
                        x_max,
                        int(np.count_nonzero(mask)),
                        float(np.mean(mean_u[mask, 0])),
                        float(np.mean(speed[mask])),
                        float(np.mean(mean_u[mask, 0] < 0.0)),
                        float(np.mean(tke[mask])),
                    ]
                )


def make_nearwall_profile_figure(profile_path: Path, out_path: Path) -> None:
    rows = list(csv.DictReader(profile_path.open(encoding="ascii")))
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.4), dpi=180, sharex=True)
    colors = {"BFS": "tab:red", "outlet-offset": "tab:blue"}
    for geometry in ("BFS", "outlet-offset"):
        subset = [row for row in rows if row["geometry"] == geometry]
        x = np.array([float(row["x_center_m"]) for row in subset])
        ux = np.array([float(row["mean_Ux_m_s"]) for row in subset])
        reverse_fraction = np.array([float(row["reverse_flow_fraction"]) for row in subset])
        tke = np.array([float(row["mean_sampled_tke_m2_s2"]) for row in subset])
        axes[0].plot(x, ux, "o-", ms=3.2, lw=1.4, color=colors[geometry], label=geometry)
        axes[1].plot(x, reverse_fraction, "o-", ms=3.2, lw=1.4, color=colors[geometry], label=geometry)
        axes[2].plot(x, tke, "o-", ms=3.2, lw=1.4, color=colors[geometry], label=geometry)
    axes[0].axhline(0, color="0.25", lw=0.8)
    axes[0].set_ylabel("near-wall mean Ux [m/s]")
    axes[1].set_ylabel("reverse-flow fraction")
    axes[2].set_ylabel("sampled TKE proxy [m2/s2]")
    axes[2].set_xlabel("x [m]")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(
        "C5 sampled near-wall recirculation profile, retained t>=0.5 s snapshots",
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.012,
        "Geometry-specific near-wall bands; sampled statistics are not long fieldAverage outputs.",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.95))
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def tricontour(ax, centers: np.ndarray, indices: np.ndarray, values: np.ndarray, title: str, cmap: str, label: str) -> None:
    xy = centers[indices, :2]
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1])
    vals = values[indices]
    contour = ax.tricontourf(tri, vals, levels=34, cmap=cmap)
    ax.tricontour(tri, vals, levels=10, colors="k", linewidths=0.15, alpha=0.2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(True, color="white", linewidth=0.25, alpha=0.25)
    cbar = plt.colorbar(contour, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(label, fontsize=8.5)


def make_figure(results: list[dict[str, object]], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2), dpi=180)
    for row, result in enumerate(results):
        spec: CaseSpec = result["spec"]  # type: ignore[assignment]
        centers = result["centers"]  # type: ignore[assignment]
        indices = result["indices"]  # type: ignore[assignment]
        speed = result["speed"]  # type: ignore[assignment]
        tke = result["tke"]  # type: ignore[assignment]
        times = result["times"]  # type: ignore[assignment]
        actual_z = result["actual_z"]  # type: ignore[assignment]
        tricontour(
            axes[row, 0],
            centers,
            indices,
            speed,
            f"{spec.label}: sampled mean |U|, z={actual_z:.4f} m",
            "viridis",
            "sampled mean |U| [m/s]",
        )
        tricontour(
            axes[row, 1],
            centers,
            indices,
            tke,
            f"{spec.label}: sampled resolved TKE, {len(times)} snapshots",
            "plasma",
            "0.5 sum var(U) [m2/s2]",
        )
    fig.suptitle(
        "C5 Geometry LES: Sampled Mean Velocity and Resolved Fluctuation Proxy",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.012,
        "Snapshot statistics from retained t>=0.5 s fields; not OpenFOAM fieldAverage, not grid/statistically converged.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.95))
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_summary(results: list[dict[str, object]], metrics_path: Path, fig_path: Path, path: Path) -> None:
    def metric_lookup(label: str, zone: str) -> dict[str, str]:
        with metrics_path.open(encoding="ascii") as handle:
            for row in csv.DictReader(handle):
                if row["geometry"] == label and row["zone"] == zone:
                    return row
        raise KeyError((label, zone))

    bfs_hot = metric_lookup("BFS", "hotspot_band")
    off_hot = metric_lookup("outlet-offset", "hotspot_band")
    bfs_near = metric_lookup("BFS", "near_heated_wall")
    off_near = metric_lookup("outlet-offset", "near_heated_wall")

    text = f"""# C5 Sampled Mean/TKE Comparison

Date: 2026-06-26

## Status

The BFS C5 run did not retain OpenFOAM `fieldAverage` outputs
(`UMean/UPrime2Mean`).  To avoid wasting the completed run, this postprocess
uses retained instantaneous `U` snapshots after `t >= 0.5 s` and computes a
matched sampled mean velocity and fluctuation proxy for BFS and outlet-offset.

This is a useful C5 comparison artifact, but it is weaker than a long
fieldAverage LES comparison.

## Inputs

| Geometry | Snapshots | Source |
|---|---:|---|
| BFS | {len(results[0]['times'])} | `cases/branch_d/c5_bfs_les_smoke_wale/` |
| outlet-offset | {len(results[1]['times'])} | `cases/branch_d/v26_les_smoke_flibe_re10000_bent_recirculation_wale/` |

## Key Metrics

| Geometry/zone | Mean sampled speed [m/s] | Mean sampled TKE [m2/s2] | Max sampled TKE [m2/s2] |
|---|---:|---:|---:|
| BFS near heated wall | {float(bfs_near['mean_speed_m_s']):.5g} | {float(bfs_near['mean_sampled_tke_m2_s2']):.5g} | {float(bfs_near['max_sampled_tke_m2_s2']):.5g} |
| outlet-offset near heated wall | {float(off_near['mean_speed_m_s']):.5g} | {float(off_near['mean_sampled_tke_m2_s2']):.5g} | {float(off_near['max_sampled_tke_m2_s2']):.5g} |
| BFS hotspot band | {float(bfs_hot['mean_speed_m_s']):.5g} | {float(bfs_hot['mean_sampled_tke_m2_s2']):.5g} | {float(bfs_hot['max_sampled_tke_m2_s2']):.5g} |
| outlet-offset hotspot band | {float(off_hot['mean_speed_m_s']):.5g} | {float(off_hot['mean_sampled_tke_m2_s2']):.5g} | {float(off_hot['max_sampled_tke_m2_s2']):.5g} |

## Artifacts

- `{metrics_path.relative_to(ROOT)}`
- `{fig_path.relative_to(ROOT)}`
- `cases/branch_d/figs/c5_sampled_nearwall_profile.csv`
- `cases/branch_d/figs/c5_sampled_nearwall_profile.png`

## Claim Boundary

Allowed:

- C5 now has BFS and outlet-offset probe-level and sampled-field comparison
  artifacts.
- The sampled field comparison supports the qualitative statement that both
  simplified recirculation geometries carry resolved velocity fluctuations in
  their geometry-specific low-flow/recovery regions.

Not allowed:

- Do not call this a statistically converged LES mean/TKE comparison.
- Do not treat the sampled TKE proxy as equivalent to a long `UPrime2Mean`
  fieldAverage.
- Do not claim thermal or wall-temperature behavior from the C5 BFS run.

## Remaining Work

Full C5 closure would require rerunning or extending BFS and outlet-offset with
explicit `fieldAverage` over a longer averaging window, then promoting the same
near-wall recirculation profile and resolved-TKE metric from those averaged
fields.
"""
    path.write_text(text, encoding="ascii")


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    results = [compute_sampled_stats(spec) for spec in CASES]
    metrics_path = FIGS / "c5_sampled_mean_tke_metrics.csv"
    fig_path = FIGS / "c5_bfs_offset_sampled_mean_tke.png"
    profile_path = FIGS / "c5_sampled_nearwall_profile.csv"
    profile_fig_path = FIGS / "c5_sampled_nearwall_profile.png"
    summary_path = FIGS / "c5_sampled_mean_tke_summary.md"
    write_metrics(results, metrics_path)
    write_nearwall_profile(results, profile_path)
    make_figure(results, fig_path)
    make_nearwall_profile_figure(profile_path, profile_fig_path)
    write_summary(results, metrics_path, fig_path, summary_path)
    print(f"wrote {metrics_path}")
    print(f"wrote {fig_path}")
    print(f"wrote {profile_path}")
    print(f"wrote {profile_fig_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()

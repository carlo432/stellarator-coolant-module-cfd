#!/usr/bin/env python3
"""Extract a v30 backward-facing-step near-wall reattachment proxy."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RANS_TIME = ROOT / "cases/baseline_module/openfoam_cases/v30_isothermal_rans_flibe_re10000_bfs/1000"
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"

STEP_X_M = 0.10
STEP_HEIGHT_M = 0.02
HEATED_X_MIN = 0.10
HEATED_X_MAX = 0.30
NEAR_WALL_Y_MAX = -0.015
BIN_COUNT = 40


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


def read_vector_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "vector":
        raise ValueError(f"{path} is {kind}, not vector")
    rows = re.findall(r"\(([^()]+)\)", body)
    data = np.array([[float(x) for x in row.split()] for row in rows], dtype=float)
    if data.shape != (count, 3):
        raise ValueError(f"{path}: expected {(count, 3)}, got {data.shape}")
    return data


def binned_rows(centers: np.ndarray, velocity: np.ndarray) -> list[dict[str, float | int]]:
    mask = (
        (centers[:, 0] >= HEATED_X_MIN)
        & (centers[:, 0] <= HEATED_X_MAX)
        & (centers[:, 1] <= NEAR_WALL_Y_MAX)
    )
    if not np.any(mask):
        raise ValueError("No BFS near-wall cells found")
    x = centers[mask, 0]
    u = velocity[mask]
    speed = np.linalg.norm(u, axis=1)
    edges = np.linspace(HEATED_X_MIN, HEATED_X_MAX, BIN_COUNT + 1)
    bin_id = np.digitize(x, edges, right=False) - 1
    bin_id = np.clip(bin_id, 0, BIN_COUNT - 1)

    rows: list[dict[str, float | int]] = []
    for i in range(BIN_COUNT):
        local = bin_id == i
        if not np.any(local):
            continue
        x_center = float((edges[i] + edges[i + 1]) * 0.5)
        rows.append(
            {
                "bin": i,
                "x_center_m": x_center,
                "x_minus_step_over_h": (x_center - STEP_X_M) / STEP_HEIGHT_M,
                "cell_count": int(np.count_nonzero(local)),
                "mean_Ux_m_s": float(np.mean(u[local, 0])),
                "mean_Uy_m_s": float(np.mean(u[local, 1])),
                "mean_speed_m_s": float(np.mean(speed[local])),
                "reverse_flow_fraction": float(np.mean(u[local, 0] < 0.0)),
            }
        )
    return rows


def first_negative_to_positive_crossing(rows: list[dict[str, float | int]]) -> float | None:
    previous_x = None
    previous_u = None
    saw_negative = False
    for row in rows:
        x = float(row["x_center_m"])
        u = float(row["mean_Ux_m_s"])
        if previous_u is not None and saw_negative and previous_u < 0.0 <= u:
            slope = (u - previous_u) / max(x - float(previous_x), 1.0e-12)
            return float(previous_x) - previous_u / slope
        if u < 0.0:
            saw_negative = True
        previous_x = x
        previous_u = u
    return None


def reverse_fraction_first_below(rows: list[dict[str, float | int]], threshold: float = 0.25) -> float | None:
    saw_reverse = False
    for row in rows:
        value = float(row["reverse_flow_fraction"])
        if value >= threshold:
            saw_reverse = True
        if saw_reverse and value < threshold:
            return float(row["x_center_m"])
    return None


def write_profile(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, float | int]]) -> None:
    crossing = first_negative_to_positive_crossing(rows)
    reverse_below = reverse_fraction_first_below(rows)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value", "note"])
        writer.writerow(["step_x_m", STEP_X_M, "BFS step location"])
        writer.writerow(["step_height_m", STEP_HEIGHT_M, "BFS step height h"])
        writer.writerow(["near_wall_band_y_max_m", NEAR_WALL_Y_MAX, "cells y <= this value"])
        writer.writerow(["near_wall_cell_count", sum(int(row["cell_count"]) for row in rows), "sampled cells"])
        writer.writerow(
            [
                "mean_Ux_reattachment_proxy_x_m",
                "" if crossing is None else crossing,
                "first negative-to-positive crossing of binned near-wall mean Ux; not Cf",
            ]
        )
        writer.writerow(
            [
                "mean_Ux_reattachment_proxy_x_over_h",
                "" if crossing is None else (crossing - STEP_X_M) / STEP_HEIGHT_M,
                "(x_r - x_step) / h",
            ]
        )
        writer.writerow(
            [
                "reverse_fraction_first_below_0p25_x_m",
                "" if reverse_below is None else reverse_below,
                "first bin after reverse region where Ux<0 fraction drops below 0.25",
            ]
        )
        writer.writerow(
            [
                "reverse_fraction_first_below_0p25_x_over_h",
                "" if reverse_below is None else (reverse_below - STEP_X_M) / STEP_HEIGHT_M,
                "(x - x_step) / h",
            ]
        )
        writer.writerow(["mean_reverse_flow_fraction", np.mean([float(row["reverse_flow_fraction"]) for row in rows]), "across bins"])
        writer.writerow(["min_binned_mean_Ux_m_s", min(float(row["mean_Ux_m_s"]) for row in rows), "near-wall band"])


def make_plot(path: Path, advisor_path: Path, rows: list[dict[str, float | int]]) -> None:
    xh = np.array([float(row["x_minus_step_over_h"]) for row in rows])
    mean_ux = np.array([float(row["mean_Ux_m_s"]) for row in rows])
    rev = np.array([float(row["reverse_flow_fraction"]) for row in rows])
    crossing = first_negative_to_positive_crossing(rows)
    crossing_xh = None if crossing is None else (crossing - STEP_X_M) / STEP_HEIGHT_M

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), dpi=180)

    ax = axes[0]
    ax.plot(xh, mean_ux, color="#f08c00", linewidth=2.1)
    ax.axhline(0.0, color="#343a40", linestyle="--", linewidth=1.0)
    if crossing_xh is not None:
        ax.axvline(crossing_xh, color="#e03131", linestyle=":", linewidth=1.5)
        ax.text(crossing_xh + 0.12, max(mean_ux) * 0.82, f"x_r/h ~ {crossing_xh:.2f}", fontsize=8.5)
    ax.fill_between(xh, np.minimum(mean_ux, 0.0), 0.0, color="#f08c00", alpha=0.2)
    ax.set_title("Near-wall mean Ux proxy", fontsize=11, fontweight="bold")
    ax.set_xlabel("(x - x_step) / h")
    ax.set_ylabel("Mean Ux [m/s]")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ax.plot(xh, rev, color="#d9480f", linewidth=2.1)
    ax.axhline(0.25, color="#495057", linestyle="--", linewidth=1.0, label="0.25 threshold")
    ax.set_title("Reverse-flow fraction", fontsize=11, fontweight="bold")
    ax.set_xlabel("(x - x_step) / h")
    ax.set_ylabel("Fraction of near-wall cells with Ux < 0")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, frameon=True)

    fig.suptitle("v30 BFS Reattachment Proxy From RANS Near-Wall Velocity", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        f"h={STEP_HEIGHT_M:.3f} m, wall band y<={NEAR_WALL_Y_MAX:.3f} m; proxy uses velocity sign, not wall shear/Cf.",
        ha="center",
        fontsize=8.8,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    advisor_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(advisor_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    centers = read_vector_field(RANS_TIME / "C")
    velocity = read_vector_field(RANS_TIME / "U")
    rows = binned_rows(centers, velocity)
    profile_path = OUT_DIR / "v30_bfs_reattachment_proxy_profile.csv"
    summary_path = OUT_DIR / "v30_bfs_reattachment_proxy_summary.csv"
    fig_path = OUT_DIR / "v30_bfs_reattachment_proxy.png"
    advisor_path = ADVISOR_DIR / "27_v30_bfs_reattachment_proxy.png"
    write_profile(profile_path, rows)
    write_summary(summary_path, rows)
    make_plot(fig_path, advisor_path, rows)
    print(profile_path)
    print(summary_path)
    print(fig_path)
    print(advisor_path)


if __name__ == "__main__":
    main()

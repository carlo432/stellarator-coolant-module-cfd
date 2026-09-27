#!/usr/bin/env python3
"""Compare v22 RANS and v27 LES mean near-wall recirculation metrics."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RANS_TIME = ROOT / "cases/baseline_module/openfoam_cases/v22_isothermal_rans_flibe_re10000_bent_recirculation/1000"
LES_TIME = ROOT / "cases/branch_d/v27_thermal_les_flibe_re10000_bent_recirculation_wale_100kw/0.500073"
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"

HEATED_X_MIN = 0.10
HEATED_X_MAX = 0.30
NEAR_WALL_Y_MAX = -0.035
BIN_COUNT = 36


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


def binned_nearwall_rows(centers: np.ndarray, rans_u: np.ndarray, les_u: np.ndarray) -> list[dict[str, float | int]]:
    mask = (
        (centers[:, 0] >= HEATED_X_MIN)
        & (centers[:, 0] <= HEATED_X_MAX)
        & (centers[:, 1] <= NEAR_WALL_Y_MAX)
    )
    if not np.any(mask):
        raise ValueError("No near-wall cells found")

    x = centers[mask, 0]
    rans = rans_u[mask]
    les = les_u[mask]
    edges = np.linspace(HEATED_X_MIN, HEATED_X_MAX, BIN_COUNT + 1)
    bin_id = np.digitize(x, edges, right=False) - 1
    bin_id = np.clip(bin_id, 0, BIN_COUNT - 1)

    rows: list[dict[str, float | int]] = []
    for i in range(BIN_COUNT):
        local = bin_id == i
        if not np.any(local):
            continue
        rans_speed = np.linalg.norm(rans[local], axis=1)
        les_speed = np.linalg.norm(les[local], axis=1)
        rows.append(
            {
                "bin": i,
                "x_center_m": float((edges[i] + edges[i + 1]) * 0.5),
                "cell_count": int(np.count_nonzero(local)),
                "rans_mean_Ux_m_s": float(np.mean(rans[local, 0])),
                "les_mean_Ux_m_s": float(np.mean(les[local, 0])),
                "rans_mean_speed_m_s": float(np.mean(rans_speed)),
                "les_mean_speed_m_s": float(np.mean(les_speed)),
                "rans_reverse_flow_fraction": float(np.mean(rans[local, 0] < 0.0)),
                "les_reverse_flow_fraction": float(np.mean(les[local, 0] < 0.0)),
            }
        )
    return rows


def zero_crossing_proxy(rows: list[dict[str, float | int]], key: str) -> float | None:
    """Return first negative-to-positive binned mean-Ux crossing, if present."""
    previous_x = None
    previous_u = None
    saw_negative = False
    for row in rows:
        x = float(row["x_center_m"])
        u = float(row[key])
        if previous_u is not None and saw_negative and previous_u < 0.0 <= u:
            slope = (u - previous_u) / max(x - float(previous_x), 1.0e-12)
            return float(previous_x) - previous_u / slope
        if u < 0.0:
            saw_negative = True
        previous_x = x
        previous_u = u
    return None


def reverse_fraction_end_proxy(rows: list[dict[str, float | int]], key: str, threshold: float = 0.25) -> float | None:
    """Return first x after a reverse-flow region where reverse fraction drops below threshold."""
    saw_reverse = False
    for row in rows:
        value = float(row[key])
        if value >= threshold:
            saw_reverse = True
        if saw_reverse and value < threshold:
            return float(row["x_center_m"])
    return None


def reverse_fraction_last_high_proxy(
    rows: list[dict[str, float | int]], key: str, threshold: float = 0.25
) -> float | None:
    """Return the last x where reverse fraction is at or above threshold."""
    xs = [float(row["x_center_m"]) for row in rows if float(row[key]) >= threshold]
    return max(xs) if xs else None


def write_profile(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, float | int]]) -> None:
    rans_zero = zero_crossing_proxy(rows, "rans_mean_Ux_m_s")
    les_zero = zero_crossing_proxy(rows, "les_mean_Ux_m_s")
    rans_reverse_end = reverse_fraction_end_proxy(rows, "rans_reverse_flow_fraction")
    les_reverse_end = reverse_fraction_end_proxy(rows, "les_reverse_flow_fraction")
    rans_last_high = reverse_fraction_last_high_proxy(rows, "rans_reverse_flow_fraction")
    les_last_high = reverse_fraction_last_high_proxy(rows, "les_reverse_flow_fraction")

    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "rans_v22", "les_v27_UMean", "note"])
        writer.writerow(
            [
                "near_wall_cell_count",
                sum(int(row["cell_count"]) for row in rows),
                sum(int(row["cell_count"]) for row in rows),
                f"x={HEATED_X_MIN}-{HEATED_X_MAX} m, y<={NEAR_WALL_Y_MAX} m",
            ]
        )
        writer.writerow(
            [
                "mean_Ux_first_negative_to_positive_x_m",
                "" if rans_zero is None else rans_zero,
                "" if les_zero is None else les_zero,
                "near-wall binned mean-Ux first-crossing proxy; not wall shear or sustained reattachment",
            ]
        )
        writer.writerow(
            [
                "reverse_fraction_first_below_0p25_after_reverse_x_m",
                "" if rans_reverse_end is None else rans_reverse_end,
                "" if les_reverse_end is None else les_reverse_end,
                "first local bin after reverse-flow region where Ux<0 fraction drops below 0.25",
            ]
        )
        writer.writerow(
            [
                "reverse_fraction_last_ge_0p25_x_m",
                "" if rans_last_high is None else rans_last_high,
                "" if les_last_high is None else les_last_high,
                "last bin where Ux<0 fraction remains at or above 0.25",
            ]
        )
        writer.writerow(
            [
                "heated_wall_mean_reverse_fraction",
                float(np.mean([float(row["rans_reverse_flow_fraction"]) for row in rows])),
                float(np.mean([float(row["les_reverse_flow_fraction"]) for row in rows])),
                "mean across x bins in the near-wall heated-wall band",
            ]
        )
        writer.writerow(
            [
                "heated_wall_min_mean_Ux_m_s",
                min(float(row["rans_mean_Ux_m_s"]) for row in rows),
                min(float(row["les_mean_Ux_m_s"]) for row in rows),
                "minimum binned mean Ux in near-wall band",
            ]
        )


def make_plot(path: Path, advisor_path: Path, rows: list[dict[str, float | int]]) -> None:
    x = np.array([float(row["x_center_m"]) for row in rows])
    rans_ux = np.array([float(row["rans_mean_Ux_m_s"]) for row in rows])
    les_ux = np.array([float(row["les_mean_Ux_m_s"]) for row in rows])
    rans_rev = np.array([float(row["rans_reverse_flow_fraction"]) for row in rows])
    les_rev = np.array([float(row["les_reverse_flow_fraction"]) for row in rows])

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), dpi=180)

    ax = axes[0]
    ax.plot(x, rans_ux, color="#1c7ed6", linewidth=2.0, label="RANS v22 mean Ux")
    ax.plot(x, les_ux, color="#e03131", linewidth=2.0, label="LES v27 UMean_x")
    ax.axhline(0.0, color="#343a40", linewidth=1.0, linestyle="--")
    ax.fill_between(x, np.minimum(rans_ux, 0.0), 0.0, color="#1c7ed6", alpha=0.16)
    ax.fill_between(x, np.minimum(les_ux, 0.0), 0.0, color="#e03131", alpha=0.16)
    ax.set_title("Near-wall mean streamwise velocity", fontsize=11, fontweight="bold")
    ax.set_xlabel("x along heated wall [m]")
    ax.set_ylabel("Mean Ux [m/s]")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, frameon=True)

    ax = axes[1]
    ax.plot(x, rans_rev, color="#1c7ed6", linewidth=2.0, label="RANS v22 Ux<0 fraction")
    ax.plot(x, les_rev, color="#e03131", linewidth=2.0, label="LES v27 UMean_x<0 fraction")
    ax.axhline(0.25, color="#495057", linewidth=1.0, linestyle="--", label="0.25 threshold")
    ax.set_title("Reverse-flow fraction in near-wall band", fontsize=11, fontweight="bold")
    ax.set_xlabel("x along heated wall [m]")
    ax.set_ylabel("Fraction of near-wall cells with Ux < 0")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, frameon=True)

    fig.suptitle("v27 RANS-vs-LES Near-Wall Recirculation Profile", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        f"Outlet-offset heated wall band: x={HEATED_X_MIN:.2f}-{HEATED_X_MAX:.2f} m, y<={NEAR_WALL_Y_MAX:.3f} m; reattachment is a velocity proxy, not Cf.",
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
    rans_u = read_vector_field(RANS_TIME / "U")
    les_u = read_vector_field(LES_TIME / "UMean")
    rows = binned_nearwall_rows(centers, rans_u, les_u)

    profile_path = OUT_DIR / "v27_rans_vs_les_nearwall_recirculation_profile.csv"
    summary_path = OUT_DIR / "v27_rans_vs_les_nearwall_recirculation_summary.csv"
    fig_path = OUT_DIR / "v27_rans_vs_les_nearwall_recirculation_profile.png"
    advisor_path = ADVISOR_DIR / "26_v27_rans_vs_les_nearwall_recirculation_profile.png"

    write_profile(profile_path, rows)
    write_summary(summary_path, rows)
    make_plot(fig_path, advisor_path, rows)

    print(profile_path)
    print(summary_path)
    print(fig_path)
    print(advisor_path)


if __name__ == "__main__":
    main()

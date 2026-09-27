#!/usr/bin/env python3
"""Plot Branch D peak wall dT versus near-wall flushing speed."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"


GEOMETRIES = [
    {
        "geometry": "90-deg bend",
        "case": "v28/v29",
        "nearwall_speed_m_s": 0.86,
        "peak_wall_dT_K": 59.9,
        "note": "weak/local separation, reattaches; no hotspot amplification",
    },
    {
        "geometry": "straight",
        "case": "v13/v14",
        "nearwall_speed_m_s": 0.47,
        "peak_wall_dT_K": 63.52,
        "note": "baseline streamwise development",
    },
    {
        "geometry": "outlet-offset",
        "case": "v22/v23",
        "nearwall_speed_m_s": 0.39,
        "peak_wall_dT_K": 77.00,
        "note": "broad slow near-wall recirculation; hotspot-prone",
    },
    {
        "geometry": "backward-facing step",
        "case": "v30/v31",
        "nearwall_speed_m_s": 0.296,
        "peak_wall_dT_K": 77.40,
        "note": "forward-predicted slow near-wall flow; hot peak confirmed",
    },
]


def write_csv(path: Path, pearson_r: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "geometry",
                "case",
                "nearwall_mean_speed_m_s",
                "peak_wall_deltaT_K",
                "pearson_r_all_four",
                "note",
            ]
        )
        for row in GEOMETRIES:
            writer.writerow(
                [
                    row["geometry"],
                    row["case"],
                    row["nearwall_speed_m_s"],
                    row["peak_wall_dT_K"],
                    pearson_r,
                    row["note"],
                ]
            )


def make_plot(path: Path, advisor_path: Path, pearson_r: float) -> None:
    speeds = np.array([row["nearwall_speed_m_s"] for row in GEOMETRIES], dtype=float)
    peaks = np.array([row["peak_wall_dT_K"] for row in GEOMETRIES], dtype=float)

    colors = {
        "90-deg bend": "#2f9e44",
        "straight": "#1c7ed6",
        "outlet-offset": "#e03131",
        "backward-facing step": "#f08c00",
    }

    fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=180)
    for row in GEOMETRIES:
        ax.scatter(
            row["nearwall_speed_m_s"],
            row["peak_wall_dT_K"],
            s=120,
            color=colors[row["geometry"]],
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
        )
        display_geometry = "BFS" if row["geometry"] == "backward-facing step" else row["geometry"]
        label = f"{display_geometry}\n{row['case']}"
        x_offset = 0.018
        ha = "left"
        ax.text(
            row["nearwall_speed_m_s"] + x_offset,
            row["peak_wall_dT_K"] + 0.8,
            label,
            fontsize=9,
            ha=ha,
            va="bottom",
        )

    coeff = np.polyfit(speeds, peaks, 1)
    x_line = np.linspace(min(speeds) - 0.03, max(speeds) + 0.03, 100)
    y_line = np.polyval(coeff, x_line)
    ax.plot(x_line, y_line, color="#495057", linestyle="--", linewidth=1.2, label="linear fit, n=4")

    ax.set_title("Branch D: Peak Wall dT Tracks Near-Wall Flushing Speed", fontsize=13, fontweight="bold")
    ax.set_xlabel("Near-wall mean speed over heated-wall region [m/s]")
    ax.set_ylabel("Peak heated-wall temperature rise [K]")
    ax.grid(True, alpha=0.28)
    ax.set_xlim(0.25, 0.92)
    ax.set_ylim(56, 81.0)
    ax.legend(loc="upper right", frameon=True)
    ax.text(
        0.46,
        57.2,
        f"Pearson r = {pearson_r:.2f} (n=4)\n"
        "BFS confirmed the slow-flow / hot-peak prediction.\n"
        "Illustrative trend, not a validated correlation.",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#ced4da", "alpha": 0.95},
    )

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    advisor_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(advisor_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    speeds = np.array([row["nearwall_speed_m_s"] for row in GEOMETRIES], dtype=float)
    peaks = np.array([row["peak_wall_dT_K"] for row in GEOMETRIES], dtype=float)
    pearson_r = float(np.corrcoef(speeds, peaks)[0, 1])

    csv_path = OUT_DIR / "branch_d_nearwall_flushing_law.csv"
    fig_path = OUT_DIR / "branch_d_nearwall_flushing_law.png"
    advisor_path = ADVISOR_DIR / "24_branch_d_nearwall_flushing_law.png"

    write_csv(csv_path, pearson_r)
    make_plot(fig_path, advisor_path, pearson_r)

    print(f"pearson_r={pearson_r:.6f}")
    print(csv_path)
    print(fig_path)
    print(advisor_path)


if __name__ == "__main__":
    main()

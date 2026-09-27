#!/usr/bin/env python3
"""Plot Branch D bent-vs-straight heat-flux sensitivity."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(".")
OUTDIR = ROOT / "results/branch_d"


ROWS = [
    {
        "case": "v14",
        "geometry": "straight",
        "heat_flux_kw_m2": 100.0,
        "wall_avg_dt_K": 43.80,
        "wall_peak_dt_K": 63.52,
    },
    {
        "case": "v16",
        "geometry": "straight",
        "heat_flux_kw_m2": 500.0,
        "wall_avg_dt_K": 218.89,
        "wall_peak_dt_K": 316.61,
    },
    {
        "case": "v23",
        "geometry": "outlet_offset",
        "heat_flux_kw_m2": 100.0,
        "wall_avg_dt_K": 44.16,
        "wall_peak_dt_K": 77.00,
    },
    {
        "case": "v24",
        "geometry": "outlet_offset",
        "heat_flux_kw_m2": 500.0,
        "wall_avg_dt_K": 220.83,
        "wall_peak_dt_K": 385.16,
    },
]


def enriched_rows() -> list[dict[str, float | str]]:
    straight_by_flux = {
        row["heat_flux_kw_m2"]: row for row in ROWS if row["geometry"] == "straight"
    }
    out: list[dict[str, float | str]] = []
    for row in ROWS:
        straight = straight_by_flux[row["heat_flux_kw_m2"]]
        peak_gap = row["wall_peak_dt_K"] - straight["wall_peak_dt_K"]
        avg_gap = row["wall_avg_dt_K"] - straight["wall_avg_dt_K"]
        peak_ratio = row["wall_peak_dt_K"] / straight["wall_peak_dt_K"]
        avg_ratio = row["wall_avg_dt_K"] / straight["wall_avg_dt_K"]
        out.append(
            {
                **row,
                "avg_gap_vs_straight_K": avg_gap,
                "peak_gap_vs_straight_K": peak_gap,
                "avg_ratio_vs_straight": avg_ratio,
                "peak_ratio_vs_straight": peak_ratio,
            }
        )
    return out


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case",
        "geometry",
        "heat_flux_kw_m2",
        "wall_avg_dt_K",
        "wall_peak_dt_K",
        "avg_gap_vs_straight_K",
        "peak_gap_vs_straight_K",
        "avg_ratio_vs_straight",
        "peak_ratio_vs_straight",
    ]
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def values(geometry: str, key: str) -> list[float]:
    rows = sorted((row for row in ROWS if row["geometry"] == geometry), key=lambda item: item["heat_flux_kw_m2"])
    return [float(row[key]) for row in rows]


def plot(path: Path) -> None:
    heat_flux = sorted({float(row["heat_flux_kw_m2"]) for row in ROWS})
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), dpi=160, sharex=True)

    styles = {
        "straight": {"label": "straight", "color": "#2459a6", "marker": "o"},
        "outlet_offset": {"label": "outlet offset", "color": "#b3261e", "marker": "s"},
    }

    for geometry, style in styles.items():
        axes[0].plot(
            heat_flux,
            values(geometry, "wall_avg_dt_K"),
            linewidth=2.0,
            marker=style["marker"],
            color=style["color"],
            label=style["label"],
        )
        axes[1].plot(
            heat_flux,
            values(geometry, "wall_peak_dt_K"),
            linewidth=2.0,
            marker=style["marker"],
            color=style["color"],
            label=style["label"],
        )

    axes[0].set_title("Wall-average response")
    axes[1].set_title("Wall-peak response")
    for ax in axes:
        ax.set_xlabel("Heat flux [kW/m^2]")
        ax.set_ylabel("Temperature rise T - Tin [K]")
        ax.grid(True, alpha=0.25)
        ax.legend()

    axes[1].annotate(
        "+68.6 K peak gap",
        xy=(500.0, 385.16),
        xytext=(280.0, 360.0),
        arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.0},
        fontsize=9,
    )
    fig.suptitle("Straight vs outlet-offset Branch D heat-flux sensitivity")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    rows = enriched_rows()
    csv_path = OUTDIR / "branch_d_heat_flux_sensitivity.csv"
    png_path = OUTDIR / "branch_d_heat_flux_sensitivity.png"
    write_csv(csv_path, rows)
    plot(png_path)
    print(csv_path)
    print(png_path)


if __name__ == "__main__":
    main()

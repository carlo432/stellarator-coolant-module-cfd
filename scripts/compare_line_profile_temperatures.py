#!/usr/bin/env python3
"""Overlay temperature line profiles from generated CSV files."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


def read_csv(path: Path, axis: str) -> tuple[list[float], list[float]]:
    x_values: list[float] = []
    t_values: list[float] = []
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_values.append(float(row[axis]))
            t_values.append(float(row["T"]))
    return x_values, t_values


def plot_pair(out_path: Path, title: str, axis: str, axis_label: str, low_csv: Path, high_csv: Path) -> None:
    x_low, t_low = read_csv(low_csv, axis)
    x_high, t_high = read_csv(high_csv, axis)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    ax.plot(x_low, t_low, color="#3066be", linewidth=2.0, label="10 kW/m^2")
    ax.plot(x_high, t_high, color="#b3261e", linewidth=2.0, label="100 kW/m^2")
    ax.set_title(title)
    ax.set_xlabel(axis_label)
    ax.set_ylabel("Temperature-like scalar T [K]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    root = Path("results/line_profiles")
    v6 = "v6_scalar_temperature_water_re10000_heatflux"
    v7 = "v7_scalar_temperature_water_re10000_heatflux100kw"

    plot_pair(
        root / "comparison_vertical_chamber_midplane_temperature.png",
        "Vertical chamber midplane",
        "y",
        "y at x=0.20 m, z=0.02 m [m]",
        root / f"{v6}_vertical_chamber_midplane.csv",
        root / f"{v7}_vertical_chamber_midplane.csv",
    )
    plot_pair(
        root / "comparison_streamwise_near_heated_wall_temperature.png",
        "Streamwise near heated wall",
        "x",
        "x at y=-0.035 m, z=0.02 m [m]",
        root / f"{v6}_streamwise_near_heated_wall.csv",
        root / f"{v7}_streamwise_near_heated_wall.csv",
    )


if __name__ == "__main__":
    main()

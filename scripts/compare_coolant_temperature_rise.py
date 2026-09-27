#!/usr/bin/env python3
"""Compare 100 kW/m^2 temperature rise profiles for water and FLiBe-like cases."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


def read_delta_t(path: Path, axis: str, inlet_temperature: float) -> tuple[list[float], list[float]]:
    axis_values: list[float] = []
    delta_t_values: list[float] = []
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            axis_values.append(float(row[axis]))
            delta_t_values.append(float(row["T"]) - inlet_temperature)
    return axis_values, delta_t_values


def plot_comparison(
    out_path: Path,
    title: str,
    axis: str,
    axis_label: str,
    water_csv: Path,
    flibe_csv: Path,
) -> None:
    x_water, dt_water = read_delta_t(water_csv, axis, inlet_temperature=300.0)
    x_flibe, dt_flibe = read_delta_t(flibe_csv, axis, inlet_temperature=800.0)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    ax.plot(x_water, dt_water, color="#3066be", linewidth=2.0, label="water, Re 10k")
    ax.plot(x_flibe, dt_flibe, color="#b3261e", linewidth=2.0, label="FLiBe-like, Re 10k")
    ax.set_title(title)
    ax.set_xlabel(axis_label)
    ax.set_ylabel("Temperature rise T - Tin [K]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    root = Path("results/line_profiles")
    water = "v7_scalar_temperature_water_re10000_heatflux100kw"
    flibe = "v9_scalar_temperature_flibe_re10000_heatflux100kw"

    plot_comparison(
        root / "comparison_100kw_water_flibe_vertical_deltaT.png",
        "100 kW/m^2 vertical temperature rise",
        "y",
        "y at x=0.20 m, z=0.02 m [m]",
        root / f"{water}_vertical_chamber_midplane.csv",
        root / f"{flibe}_vertical_chamber_midplane.csv",
    )
    plot_comparison(
        root / "comparison_100kw_water_flibe_streamwise_deltaT.png",
        "100 kW/m^2 streamwise temperature rise",
        "x",
        "x at y=-0.035 m, z=0.02 m [m]",
        root / f"{water}_streamwise_near_heated_wall.csv",
        root / f"{flibe}_streamwise_near_heated_wall.csv",
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compare FLiBe-like 100 kW/m^2 temperature-rise profiles across refinements."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


def read_delta_t(path: Path, axis: str) -> tuple[list[float], list[float]]:
    axis_values: list[float] = []
    delta_t_values: list[float] = []
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            axis_values.append(float(row[axis]))
            delta_t_values.append(float(row["T"]) - 800.0)
    return axis_values, delta_t_values


def plot_comparison(out_path: Path, title: str, axis: str, axis_label: str, cases: list[tuple[str, Path, str]]) -> None:
    styles = ["-", "-", "-"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    for i, (label, csv_path, color) in enumerate(cases):
        if not csv_path.exists():
            continue
        x_values, dt_values = read_delta_t(csv_path, axis)
        ax.plot(x_values, dt_values, color=color, linewidth=2.0, linestyle=styles[i], label=label)
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
    coarse = "v9_scalar_temperature_flibe_re10000_heatflux100kw"
    refined = "v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005"
    thermal_refined = "v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025"

    vertical_cases = [
        ("coarse v9", root / f"{coarse}_vertical_chamber_midplane.csv", "#7b3f98"),
        ("refined v12", root / f"{refined}_vertical_chamber_midplane.csv", "#0f766e"),
        ("thermal refined v14", root / f"{thermal_refined}_vertical_chamber_midplane.csv", "#b45309"),
    ]
    streamwise_cases = [
        ("coarse v9", root / f"{coarse}_streamwise_near_heated_wall.csv", "#7b3f98"),
        ("refined v12", root / f"{refined}_streamwise_near_heated_wall.csv", "#0f766e"),
        ("thermal refined v14", root / f"{thermal_refined}_streamwise_near_heated_wall.csv", "#b45309"),
    ]

    plot_comparison(
        root / "comparison_flibe_refinement_vertical_deltaT.png",
        "FLiBe-like 100 kW/m^2 vertical refinement check",
        "y",
        "y at x=0.20 m, z=0.02 m [m]",
        vertical_cases,
    )
    plot_comparison(
        root / "comparison_flibe_refinement_streamwise_deltaT.png",
        "FLiBe-like 100 kW/m^2 streamwise refinement check",
        "x",
        "x at y=-0.035 m, z=0.02 m [m]",
        streamwise_cases,
    )


if __name__ == "__main__":
    main()

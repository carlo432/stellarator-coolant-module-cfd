#!/usr/bin/env python3
"""Overlay heated-wall temperature-rise profiles across mesh refinements."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


def read_profile(path: Path) -> tuple[list[float], list[float], list[float]]:
    x_values: list[float] = []
    mean_values: list[float] = []
    max_values: list[float] = []
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_values.append(float(row["x"]))
            mean_values.append(float(row["dT_mean"]))
            max_values.append(float(row["dT_max"]))
    return x_values, mean_values, max_values


def main() -> None:
    root = Path("results/wall_profiles")
    cases = [
        ("coarse v9", root / "v9_flibe_coarse_100kw_heated_wall_surface_profile.csv", "#7b3f98"),
        ("refined v12", root / "v12_flibe_refined_100kw_heated_wall_surface_profile.csv", "#0f766e"),
        ("thermal refined v14", root / "v14_flibe_thermal_refined_100kw_heated_wall_surface_profile.csv", "#b45309"),
    ]

    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    for label, csv_path, color in cases:
        if not csv_path.exists():
            continue
        x_values, mean_values, max_values = read_profile(csv_path)
        ax.plot(x_values, mean_values, color=color, linewidth=2.0, label=f"{label} mean")
        ax.plot(x_values, max_values, color=color, linewidth=1.4, linestyle="--", label=f"{label} max")
    ax.set_title("FLiBe-like heated-wall refinement check")
    ax.set_xlabel("x along heated wall [m]")
    ax.set_ylabel("Wall temperature rise T - Tin [K]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    out_path = root / "comparison_flibe_heated_wall_refinement_deltaT.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()

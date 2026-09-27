#!/usr/bin/env python3
"""Create a compact mesh-sensitivity summary table and plot."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


CASES = [
    {
        "case": "v9 coarse",
        "cells": 3792,
        "heated_wall_faces": 160,
        "pressure_drop_pa": 1860.0,
        "outlet_avg_t": 800.123,
        "heated_wall_avg_t": 1084.43,
        "max_t": 1119.34,
        "wall_surface_avg_dt": 284.43,
        "wall_surface_max_dt": 319.34,
    },
    {
        "case": "v12 refined",
        "cells": 28976,
        "heated_wall_faces": 1280,
        "pressure_drop_pa": 1858.0,
        "outlet_avg_t": 800.195,
        "heated_wall_avg_t": 869.171,
        "max_t": 893.437,
        "wall_surface_avg_dt": 69.17,
        "wall_surface_max_dt": 93.44,
    },
    {
        "case": "v14 thermal refined",
        "cells": 73344,
        "heated_wall_faces": 3216,
        "pressure_drop_pa": 1774.0,
        "outlet_avg_t": 800.200,
        "heated_wall_avg_t": 843.798,
        "max_t": 863.519,
        "wall_surface_avg_dt": 43.80,
        "wall_surface_max_dt": 63.52,
    },
]


def write_csv(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(CASES[0].keys()))
        writer.writeheader()
        writer.writerows(CASES)


def plot_summary(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = [case["case"] for case in CASES]
    wall_avg = [case["wall_surface_avg_dt"] for case in CASES]
    wall_max = [case["wall_surface_max_dt"] for case in CASES]
    pressure = [case["pressure_drop_pa"] for case in CASES]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 4.2), dpi=160)

    x = range(len(labels))
    ax1.plot(x, wall_avg, color="#0f766e", marker="o", linewidth=2.0, label="wall avg")
    ax1.plot(x, wall_max, color="#b3261e", marker="o", linewidth=2.0, linestyle="--", label="wall max")
    ax1.set_xticks(list(x), labels, rotation=15, ha="right")
    ax1.set_ylabel("Heated-wall T - Tin [K]")
    ax1.set_title("Wall temperature trend")
    ax1.grid(True, alpha=0.25)
    ax1.legend()

    ax2.plot(x, pressure, color="#2459a6", marker="o", linewidth=2.0)
    ax2.set_xticks(list(x), labels, rotation=15, ha="right")
    ax2.set_ylabel("Pressure drop [Pa]")
    ax2.set_title("RANS pressure drop trend")
    ax2.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    outdir = Path("results/mesh_sensitivity")
    write_csv(outdir / "flibe_100kw_mesh_sensitivity_summary.csv")
    plot_summary(outdir / "flibe_100kw_mesh_sensitivity_summary.png")
    print(outdir / "flibe_100kw_mesh_sensitivity_summary.csv")
    print(outdir / "flibe_100kw_mesh_sensitivity_summary.png")


if __name__ == "__main__":
    main()

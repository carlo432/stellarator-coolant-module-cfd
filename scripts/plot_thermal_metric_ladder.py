#!/usr/bin/env python3
"""Build the thermal-metric ladder for Branch D/CHT claim control.

This is a reporting artifact, not a new CFD result. It combines previously computed
CHT grid data with the existing passive-scalar GCI correction so the report can
say which thermal metric is trustworthy and which is only diagnostic.
"""

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

OUT_DIR.mkdir(parents=True, exist_ok=True)
ADVISOR_DIR.mkdir(parents=True, exist_ok=True)


CHT_GRID_ROWS = [
    {
        "case": "cht_coarse",
        "model": "k-epsilon wall-function CHT",
        "cells": 4560,
        "first_cell_y_plus": 16.7,
        "interface_dT_K": 23.94,
        "wall_function_band": "below 30; buffer layer",
        "status": "wall-function invalid region",
    },
    {
        "case": "cht_medium",
        "model": "k-epsilon wall-function CHT",
        "cells": 20480,
        "first_cell_y_plus": 9.8,
        "interface_dT_K": 43.55,
        "wall_function_band": "below 30; buffer layer",
        "status": "wall-function invalid region",
    },
    {
        "case": "cht_fine",
        "model": "k-epsilon wall-function CHT",
        "cells": 84864,
        "first_cell_y_plus": 6.2,
        "interface_dT_K": 28.61,
        "wall_function_band": "below 30; buffer layer",
        "status": "wall-function invalid region",
    },
]


LADDER_ROWS = [
    {
        "rank": "4",
        "metric": "Passive scalar peak wall dT",
        "representative_dT_K": 77.00,
        "values_K": "94.43 / 77.00 / 105.38",
        "grid_status": "non-monotonic; formal GCI invalid",
        "physics_status": "fixed-gradient passive scalar; no turbulent thermal diffusivity",
        "trust_verdict": "same-mesh hotspot diagnostic only",
        "recommended_use": "locate hotspot-prone regions; do not quote as converged wall temperature",
        "correlation_check": "about 3-4x above 25-26 K duct-correlation film dT",
    },
    {
        "rank": "3",
        "metric": "Passive scalar wall-average dT",
        "representative_dT_K": 44.16,
        "values_K": "73.74 / 44.16 / 32.58",
        "grid_status": "monotonic but still mesh-sensitive",
        "physics_status": "energy-balance-bound average from passive scalar",
        "trust_verdict": "more stable than peak but weak geometry discriminator",
        "recommended_use": "global heat-removal context and apparent Nu/St diagnostics",
        "correlation_check": "not a direct wall-film validation metric",
    },
    {
        "rank": "2",
        "metric": "CHT interface dT, k-epsilon wall functions",
        "representative_dT_K": 43.55,
        "values_K": "23.94 / 43.55 / 28.61",
        "grid_status": "non-monotonic because y+ = 16.7 / 9.8 / 6.2 sits in buffer layer",
        "physics_status": "physical CHT with energy balance and solid conduction, but wall-function invalid",
        "trust_verdict": "physical but wall-function scattered; use wall-resolved result instead",
        "recommended_use": "show CHT path is physical; do not quote a converged k-epsilon digit",
        "correlation_check": "fine case 28.61 K agrees; sequence still wall-function scattered",
    },
    {
        "rank": "1",
        "metric": "Wall-resolved CHT interface dT, kOmegaSST",
        "representative_dT_K": "27-28",
        "values_K": "28.16 / 27.21",
        "grid_status": "two wall-resolved meshes 3.4% apart; y+ avg 0.57 / 0.59",
        "physics_status": "solid conduction exact; kOmegaSST integrates to the wall",
        "trust_verdict": "grid-confirmed and correlation-validated wall-temperature result",
        "recommended_use": "quote 27-28 K interface superheat at 100 kW/m2",
        "correlation_check": "within 6-10% of Dittus-Boelter/Gnielinski film dT",
    },
]


def write_csvs() -> None:
    with (OUT_DIR / "cht_wall_function_grid_sensitivity.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(CHT_GRID_ROWS[0]))
        writer.writeheader()
        writer.writerows(CHT_GRID_ROWS)

    with (OUT_DIR / "thermal_metric_ladder.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(LADDER_ROWS[0]))
        writer.writeheader()
        writer.writerows(LADDER_ROWS)


def plot_ladder() -> None:
    labels = [
        "Passive\npeak",
        "Passive\nwall avg",
        "k-eps CHT\ninterface",
        "wall-resolved\nCHT target",
    ]
    y = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    fig.patch.set_facecolor("white")

    # Passive peak and passive wall-average mesh sequences.
    sequences = [
        [94.43, 77.00, 105.38],
        [73.74, 44.16, 32.58],
        [23.94, 43.55, 28.61],
    ]
    colors = ["#c92a2a", "#1c7ed6", "#2f9e44"]
    markers = ["o", "s", "^"]
    for idx, vals in enumerate(sequences):
        xs = np.array(vals, dtype=float)
        jit = np.linspace(-0.13, 0.13, len(xs))
        ax.scatter(xs, np.full_like(xs, y[idx]) + jit, s=95, color=colors[idx], marker=markers[idx], zorder=3)
        ax.plot(xs, np.full_like(xs, y[idx]), color=colors[idx], lw=2, alpha=0.45)

    # Wall-resolved kOmegaSST CHT result and refined confirmation.
    ax.scatter([28.16, 27.21], [y[3] - 0.04, y[3] + 0.04], s=130, facecolor="white", edgecolor="#495057", marker="D", lw=2, zorder=3)
    ax.text(32.0, y[3], "27-28 K\n3.4% apart", va="center", ha="left", fontsize=10, color="#343a40")

    ax.axvspan(26.0, 31.0, color="#2f9e44", alpha=0.10, label="wall-resolved CHT result region")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("wall/interface superheat dT [K]")
    ax.set_title("Thermal Metric Ladder: Which Wall-Heating Number To Trust")
    ax.grid(axis="x", alpha=0.25)
    ax.set_xlim(15, 112)

    notes = [
        "same-mesh hotspot diagnostic only",
        "energy-bound average; weak geometry discriminator",
        "physical CHT but k-eps wall functions invalid at y+ 6-17",
        "grid-confirmed; matches correlations within 6-10%",
    ]
    for idx, note in enumerate(notes):
        ax.text(
            111,
            y[idx],
            note,
            va="center",
            ha="right",
            fontsize=9,
            color="#495057",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.6},
        )

    fig.text(
        0.125,
        0.035,
        "Sequences are coarse / medium / fine where available. k-eps CHT is not grid-converged; wall-resolved CHT is 27-28 K and correlation-validated.",
        fontsize=9,
        color="#495057",
    )

    fig.tight_layout(rect=[0, 0.07, 1, 1])
    out = OUT_DIR / "thermal_metric_ladder.png"
    advisor = ADVISOR_DIR / "30_thermal_metric_ladder.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)


def main() -> None:
    write_csvs()
    plot_ladder()
    print(OUT_DIR / "thermal_metric_ladder.csv")
    print(OUT_DIR / "cht_wall_function_grid_sensitivity.csv")
    print(OUT_DIR / "thermal_metric_ladder.png")
    print(ADVISOR_DIR / "30_thermal_metric_ladder.png")


if __name__ == "__main__":
    main()

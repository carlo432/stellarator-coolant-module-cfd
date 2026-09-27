#!/usr/bin/env python3
"""Create report-ready 100 vs 500 kW/m^2 heat-flux sensitivity plots."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(".")
TIN = 800.0


CASES = [
    {
        "label": "100 kW/m^2 v14",
        "heat_flux": 100.0,
        "case_dir": ROOT
        / "cases/baseline_module/openfoam_cases/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025",
        "line_prefix": "v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025",
        "wall_csv": ROOT / "results/wall_profiles/v14_flibe_thermal_refined_100kw_heated_wall_surface_profile.csv",
        "color": "#0f766e",
    },
    {
        "label": "500 kW/m^2 v16",
        "heat_flux": 500.0,
        "case_dir": ROOT
        / "cases/baseline_module/openfoam_cases/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025",
        "line_prefix": "v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025",
        "wall_csv": ROOT / "results/wall_profiles/v16_flibe_thermal_refined_500kw_heated_wall_surface_profile.csv",
        "color": "#b3261e",
    },
]


def read_outlet_average(case_dir: Path) -> float:
    text = (case_dir / "2/T").read_text(encoding="ascii")
    match = re.search(
        r"outlet\s*\{.*?value\s+nonuniform\s+List<scalar>\s+\d+\s*\((?P<body>.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError(f"Could not find nonuniform outlet values in {case_dir / '2/T'}")
    values = [float(x) for x in match.group("body").split()]
    return sum(values) / len(values)


def read_wall_stats(path: Path) -> tuple[float, float]:
    weighted_sum = 0.0
    count_sum = 0.0
    max_dt = None
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            count = float(row["count"])
            weighted_sum += count * float(row["T_mean"])
            count_sum += count
            row_max_dt = float(row["dT_max"])
            max_dt = row_max_dt if max_dt is None else max(max_dt, row_max_dt)
    if count_sum == 0.0 or max_dt is None:
        raise ValueError(f"No wall-profile rows found in {path}")
    wall_avg_t = weighted_sum / count_sum
    return wall_avg_t, TIN + max_dt


def read_line_max_dt(path: Path) -> float:
    max_dt = None
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dt = float(row["T"]) - TIN
            max_dt = dt if max_dt is None else max(max_dt, dt)
    if max_dt is None:
        raise ValueError(f"No line-profile rows found in {path}")
    return max_dt


def write_summary_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case",
        "heat_flux_kw_m2",
        "outlet_avg_t",
        "outlet_avg_dt",
        "wall_surface_avg_t",
        "wall_surface_avg_dt",
        "wall_surface_max_t",
        "wall_surface_max_dt",
        "vertical_line_max_dt",
        "streamwise_line_max_dt",
    ]
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_summary(path: Path, rows: list[dict[str, float | str]]) -> None:
    heat_flux = [float(row["heat_flux_kw_m2"]) for row in rows]
    wall_avg = [float(row["wall_surface_avg_dt"]) for row in rows]
    wall_max = [float(row["wall_surface_max_dt"]) for row in rows]
    outlet = [float(row["outlet_avg_dt"]) for row in rows]

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    ax.plot(heat_flux, wall_avg, marker="o", linewidth=2.0, color="#0f766e", label="wall avg")
    ax.plot(heat_flux, wall_max, marker="s", linewidth=2.0, color="#b3261e", label="wall max")
    ax.plot(heat_flux, outlet, marker="^", linewidth=2.0, color="#2459a6", label="outlet avg")
    ax.set_title("FLiBe-like thermal-refined heat-flux sensitivity")
    ax.set_xlabel("Imposed heat flux [kW/m^2]")
    ax.set_ylabel("Temperature rise T - Tin [K]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def read_wall_profile(path: Path) -> tuple[list[float], list[float], list[float]]:
    x_values: list[float] = []
    mean_dt: list[float] = []
    max_dt: list[float] = []
    with path.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_values.append(float(row["x"]))
            mean_dt.append(float(row["dT_mean"]))
            max_dt.append(float(row["dT_max"]))
    return x_values, mean_dt, max_dt


def plot_wall_overlay(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    for case in CASES:
        x_values, mean_dt, max_dt = read_wall_profile(case["wall_csv"])
        ax.plot(x_values, mean_dt, color=case["color"], linewidth=2.0, label=f"{case['label']} mean")
        ax.plot(x_values, max_dt, color=case["color"], linewidth=1.4, linestyle="--", label=f"{case['label']} max")
    ax.set_title("Thermal-refined heated-wall heat-flux sensitivity")
    ax.set_xlabel("x along heated wall [m]")
    ax.set_ylabel("Wall temperature rise T - Tin [K]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    outdir = ROOT / "results/heat_flux_sensitivity"
    rows: list[dict[str, float | str]] = []
    for case in CASES:
        line_root = ROOT / "results/line_profiles"
        outlet_avg_t = read_outlet_average(case["case_dir"])
        wall_avg_t, wall_max_t = read_wall_stats(case["wall_csv"])
        vertical_max_dt = read_line_max_dt(line_root / f"{case['line_prefix']}_vertical_chamber_midplane.csv")
        streamwise_max_dt = read_line_max_dt(line_root / f"{case['line_prefix']}_streamwise_near_heated_wall.csv")
        rows.append(
            {
                "case": case["label"],
                "heat_flux_kw_m2": case["heat_flux"],
                "outlet_avg_t": outlet_avg_t,
                "outlet_avg_dt": outlet_avg_t - TIN,
                "wall_surface_avg_t": wall_avg_t,
                "wall_surface_avg_dt": wall_avg_t - TIN,
                "wall_surface_max_t": wall_max_t,
                "wall_surface_max_dt": wall_max_t - TIN,
                "vertical_line_max_dt": vertical_max_dt,
                "streamwise_line_max_dt": streamwise_max_dt,
            }
        )

    csv_path = outdir / "flibe_thermal_refined_heat_flux_summary.csv"
    summary_png = outdir / "flibe_thermal_refined_heat_flux_summary.png"
    wall_png = outdir / "comparison_flibe_thermal_refined_100kw_500kw_wall_deltaT.png"
    write_summary_csv(csv_path, rows)
    plot_summary(summary_png, rows)
    plot_wall_overlay(wall_png)
    print(csv_path)
    print(summary_png)
    print(wall_png)


if __name__ == "__main__":
    main()

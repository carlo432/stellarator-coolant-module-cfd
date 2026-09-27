#!/usr/bin/env python3
"""Generate the AH CHT correlation-validation table and figure."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
ADVISOR_DIR.mkdir(parents=True, exist_ok=True)


Q_WALL = 100_000.0
RE = 10_000.0
PR = 14.4
DH = 0.02667
K_FLUID = 1.0
FRICTION_FACTOR = 0.0315

CHT_DT = {
    "wall-resolved CHT, 50k": 28.16,
    "wall-resolved CHT, 98k": 27.21,
}
PASSIVE_PEAKS = {
    "passive peak, medium": 77.00,
    "passive peak, fine": 105.38,
}


def nu_from_dt(dt_k: float) -> float:
    return Q_WALL * DH / (dt_k * K_FLUID)


def dt_from_nu(nu: float) -> float:
    return Q_WALL * DH / (nu * K_FLUID)


def dittus_boelter() -> float:
    return 0.023 * RE**0.8 * PR**0.4


def gnielinski() -> float:
    numerator = (FRICTION_FACTOR / 8.0) * (RE - 1000.0) * PR
    denominator = 1.0 + 12.7 * math.sqrt(FRICTION_FACTOR / 8.0) * (PR ** (2.0 / 3.0) - 1.0)
    return numerator / denominator


def build_rows() -> list[dict[str, str]]:
    nu_db = dittus_boelter()
    nu_gn = gnielinski()
    dt_db = dt_from_nu(nu_db)
    dt_gn = dt_from_nu(nu_gn)

    rows: list[dict[str, str]] = []

    for label, dt in CHT_DT.items():
        nu = nu_from_dt(dt)
        rows.append(
            {
                "metric": label,
                "deltaT_K": f"{dt:.2f}",
                "Nu": f"{nu:.1f}",
                "comparison": "CHT result",
                "difference_vs_gnielinski_pct": f"{100.0 * (dt - dt_gn) / dt_gn:.1f}",
                "interpretation": "wall-resolved CHT agrees with duct correlations",
            }
        )

    rows.extend(
        [
            {
                "metric": "Dittus-Boelter",
                "deltaT_K": f"{dt_db:.2f}",
                "Nu": f"{nu_db:.1f}",
                "comparison": "textbook turbulent duct correlation",
                "difference_vs_gnielinski_pct": f"{100.0 * (dt_db - dt_gn) / dt_gn:.1f}",
                "interpretation": "independent film-temperature scale",
            },
            {
                "metric": "Gnielinski",
                "deltaT_K": f"{dt_gn:.2f}",
                "Nu": f"{nu_gn:.1f}",
                "comparison": "textbook turbulent duct correlation",
                "difference_vs_gnielinski_pct": "0.0",
                "interpretation": "independent film-temperature scale",
            },
        ]
    )

    for label, dt in PASSIVE_PEAKS.items():
        rows.append(
            {
                "metric": label,
                "deltaT_K": f"{dt:.2f}",
                "Nu": f"{nu_from_dt(dt):.1f}",
                "comparison": "retracted passive-scalar peak metric",
                "difference_vs_gnielinski_pct": f"{100.0 * (dt - dt_gn) / dt_gn:.1f}",
                "interpretation": "about 3-4x above the correlation film-temperature scale",
            }
        )

    return rows


def write_csv(rows: list[dict[str, str]]) -> Path:
    out = OUT_DIR / "cht_correlation_validation.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return out


def plot(rows: list[dict[str, str]]) -> tuple[Path, Path]:
    plot_rows = [row for row in rows if "passive peak" not in row["metric"]] + [
        row for row in rows if "passive peak" in row["metric"]
    ]
    labels = [row["metric"] for row in plot_rows]
    values = [float(row["deltaT_K"]) for row in plot_rows]
    colors = ["#2f9e44", "#2f9e44", "#1c7ed6", "#1c7ed6", "#c92a2a", "#c92a2a"]

    fig, ax = plt.subplots(figsize=(11.5, 5.9))
    fig.patch.set_facecolor("white")
    bars = ax.barh(labels, values, color=colors, alpha=0.86)
    ax.invert_yaxis()
    ax.set_xlabel("film/interface superheat dT [K]")
    ax.set_title("CHT Wall Temperature Validation Against Turbulent-Duct Correlations")
    ax.grid(axis="x", alpha=0.25)
    ax.set_xlim(0, 115)

    for bar, value in zip(bars, values):
        ax.text(value + 1.0, bar.get_y() + bar.get_height() / 2.0, f"{value:.1f} K", va="center", fontsize=9)

    ax.text(
        0.01,
        -0.12,
        "Re=10,000, Pr=14.4, Dh=0.02667 m, k=1 W/m/K, q''=100 kW/m2. Passive peaks are shown only as the retracted metric.",
        transform=ax.transAxes,
        fontsize=9,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.07, 1, 1])

    out = OUT_DIR / "cht_correlation_validation.png"
    advisor = ADVISOR_DIR / "31_cht_correlation_validation.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)
    return out, advisor


def main() -> None:
    rows = build_rows()
    csv_out = write_csv(rows)
    fig_out, advisor_out = plot(rows)
    print(csv_out)
    print(fig_out)
    print(advisor_out)


if __name__ == "__main__":
    main()

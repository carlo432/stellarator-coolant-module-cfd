#!/usr/bin/env python3
"""Generate the AI FLiBe thermal-margin table and figure."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
ADVISOR_DIR.mkdir(parents=True, exist_ok=True)


T_INLET = 800.0
MELT_K = 732.0
FLIBE_LIQUID_UPPER_K = 1703.0
RAFM_UPPER_K = 873.0
PINT_PEAK_LOWER_K = 873.0
PINT_PEAK_UPPER_K = 973.0

BASE_BULK_RISE_100_K = 0.37
BASE_BULK_RISE_500_K = 1.85


def case_row(
    label: str,
    heat_flux_kw_m2: float,
    bulk_rise_reference: float,
    interface_low: float,
    interface_high: float,
    solid_base_low: float,
    solid_base_high: float,
    status: str,
    interpretation: str,
) -> dict[str, str]:
    superheat_low = interface_low - T_INLET
    superheat_high = interface_high - T_INLET
    base_drop_low = solid_base_low - interface_low
    base_drop_high = solid_base_high - interface_high

    return {
        "case": label,
        "heat_flux_kW_m2": f"{heat_flux_kw_m2:.0f}",
        "status": status,
        "inlet_temperature_K": f"{T_INLET:.2f}",
        "bulk_rise_reference_K": f"{bulk_rise_reference:.2f}",
        "interface_superheat_K_low": f"{superheat_low:.2f}",
        "interface_superheat_K_high": f"{superheat_high:.2f}",
        "interface_temperature_K_low": f"{interface_low:.2f}",
        "interface_temperature_K_high": f"{interface_high:.2f}",
        "solid_base_drop_K_low": f"{base_drop_low:.2f}",
        "solid_base_drop_K_high": f"{base_drop_high:.2f}",
        "solid_base_temperature_K_low": f"{solid_base_low:.2f}",
        "solid_base_temperature_K_high": f"{solid_base_high:.2f}",
        "interface_margin_above_melt_K_low": f"{interface_low - MELT_K:.2f}",
        "solid_base_margin_above_melt_K_low": f"{solid_base_low - MELT_K:.2f}",
        "solid_base_margin_to_raFM_873K_K": f"{RAFM_UPPER_K - solid_base_high:.2f}",
        "solid_base_margin_to_pint_973K_K": f"{PINT_PEAK_UPPER_K - solid_base_high:.2f}",
        "interpretation": interpretation,
    }


def build_rows() -> list[dict[str, str]]:
    return [
        case_row(
            "wall-resolved CHT bridge",
            100.0,
            BASE_BULK_RISE_100_K,
            827.212171,
            828.160944,
            840.545483,
            841.494288,
            "run",
            "thermally plausible at 800 K inlet; structural margin to 873 K is real but not large",
        ),
        case_row(
            "wall-resolved CHT high-heat-flux run",
            500.0,
            BASE_BULK_RISE_500_K,
            940.804663,
            940.804663,
            1007.47153,
            1007.47153,
            "run",
            "computed warning case; solid base exceeds 973 K and RAFM upper-window guidance",
        ),
    ]


def write_csv(rows: list[dict[str, str]]) -> Path:
    out = OUT_DIR / "flibe_thermal_margin_ai.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return out


def plot(rows: list[dict[str, str]]) -> tuple[Path, Path]:
    labels = ["100 kW/m2\nCHT run", "500 kW/m2\nCHT run"]
    interface_mid = [
        (
            float(row["interface_temperature_K_low"])
            + float(row["interface_temperature_K_high"])
        )
        / 2.0
        for row in rows
    ]
    interface_err = [
        (
            float(row["interface_temperature_K_high"])
            - float(row["interface_temperature_K_low"])
        )
        / 2.0
        for row in rows
    ]
    base_mid = [
        (
            float(row["solid_base_temperature_K_low"])
            + float(row["solid_base_temperature_K_high"])
        )
        / 2.0
        for row in rows
    ]
    base_err = [
        (
            float(row["solid_base_temperature_K_high"])
            - float(row["solid_base_temperature_K_low"])
        )
        / 2.0
        for row in rows
    ]

    x = range(len(rows))
    fig, ax = plt.subplots(figsize=(10.8, 6.0))
    fig.patch.set_facecolor("white")

    ax.axhspan(MELT_K, FLIBE_LIQUID_UPPER_K, color="#d8f3dc", alpha=0.35, label="reported FLiBe liquid window")
    ax.axhspan(PINT_PEAK_LOWER_K, PINT_PEAK_UPPER_K, color="#ffe8cc", alpha=0.55, label="Pint blanket peak operating band")
    ax.axhline(T_INLET, color="#495057", linestyle=":", linewidth=1.6, label="inlet 800 K")
    ax.axhline(RAFM_UPPER_K, color="#c92a2a", linestyle="--", linewidth=1.6, label="RAFM upper-window guide, 873 K")
    ax.axhline(MELT_K, color="#2b8a3e", linestyle="-.", linewidth=1.2, label="FLiBe melt/liquidus guide, 732 K")

    ax.errorbar(
        [i - 0.08 for i in x],
        interface_mid,
        yerr=interface_err,
        fmt="o",
        markersize=8,
        capsize=5,
        color="#1971c2",
        label="coolant-facing interface",
    )
    ax.errorbar(
        [i + 0.08 for i in x],
        base_mid,
        yerr=base_err,
        fmt="s",
        markersize=7,
        capsize=5,
        color="#e67700",
        label="heated solid base",
    )

    for i, row in enumerate(rows):
        ax.text(
            i - 0.08,
            interface_mid[i] + 12,
            f"{interface_mid[i]:.0f} K",
            ha="center",
            fontsize=9,
            color="#1971c2",
        )
        ax.text(
            i + 0.08,
            base_mid[i] + 12,
            f"{base_mid[i]:.0f} K",
            ha="center",
            fontsize=9,
            color="#e67700",
        )

    ax.set_xticks(list(x), labels)
    ax.set_ylabel("temperature [K]")
    ax.set_ylim(700, 1045)
    ax.set_title("AI: FLiBe Thermal Margin from Wall-Resolved CHT Bridge")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper left", fontsize=8.5, ncol=2)
    ax.text(
        0.0,
        -0.16,
        "Both rows are wall-resolved CHT runs. The 500 kW/m2 case confirms linear heat-flux scaling "
        "for this fixed-geometry bridge; solid-base temperatures include q*t/kappa through the 4 mm solid.",
        transform=ax.transAxes,
        fontsize=9,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.07, 1, 1])

    out = OUT_DIR / "flibe_thermal_margin_ai.png"
    advisor = ADVISOR_DIR / "32_flibe_thermal_margin.png"
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

#!/usr/bin/env python3
"""Generate the AP wall-temperature versus pumping-power tradeoff."""

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


RHO = 1940.0
U0 = 1.1598
RE0 = 10_000.0
INLET_AREA = 0.0008
Q0 = U0 * INLET_AREA
DP0_STRAIGHT_PA = 1774.444086
PUMP0_W = DP0_STRAIGHT_PA * Q0
DT_INTERFACE_100KW_K = 28.16
SOLID_DROP_100KW_K = 13.333333
TIN_K = 800.0
PRESSURE_EXPONENT = 1.8
HEAT_TRANSFER_EXPONENT = -0.7047
AN_SWEEP_DT_BY_RE = {
    5000: 45.25,
    7500: 34.39,
    10000: 28.16,
    15000: 21.09,
    20000: 17.03,
}

GEOMETRY_MARKERS = [
    {
        "geometry": "90-deg bend",
        "case_pair": "v28/v29",
        "pressure_drop_Pa": 1770.038928,
        "pressure_source": "patch-average inlet pressure from v28 postProcessing",
        "nearwall_speed_m_s": 0.86,
        "passive_peak_dT_K": 59.9,
        "passive_mean_dT_K": 43.51,
    },
    {
        "geometry": "straight",
        "case_pair": "v13/v14",
        "pressure_drop_Pa": 1774.444086,
        "pressure_source": "patch-average inlet minus outlet pressure from v13 postProcessing",
        "nearwall_speed_m_s": 0.47,
        "passive_peak_dT_K": 63.52,
        "passive_mean_dT_K": 43.80,
    },
    {
        "geometry": "outlet-offset",
        "case_pair": "v22/v23",
        "pressure_drop_Pa": 1983.0,
        "pressure_source": "documented Branch D v22 result",
        "nearwall_speed_m_s": 0.39,
        "passive_peak_dT_K": 77.00,
        "passive_mean_dT_K": 44.16,
    },
    {
        "geometry": "backward-facing step",
        "case_pair": "v30/v31",
        "pressure_drop_Pa": None,
        "pressure_source": "unresolved: patch-average pressure sign failed check",
        "nearwall_speed_m_s": 0.296,
        "passive_peak_dT_K": 77.4,
        "passive_mean_dT_K": 46.37,
    },
]


def build_tradeoff_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for re in [5000, 7500, 10000, 15000, 20000]:
        ratio = re / RE0
        u = U0 * ratio
        q_vol = Q0 * ratio
        dp = DP0_STRAIGHT_PA * ratio**PRESSURE_EXPONENT
        pump = dp * q_vol
        interface_dt = AN_SWEEP_DT_BY_RE[re]
        base_dt = interface_dt + SOLID_DROP_100KW_K
        rows.append(
            {
                "Re": f"{re:.0f}",
                "U_m_s": f"{u:.4f}",
                "Q_m3_s": f"{q_vol:.8f}",
                "pressure_drop_Pa": f"{dp:.2f}",
                "pump_power_W": f"{pump:.4f}",
                "pump_power_relative_to_Re10000": f"{pump / PUMP0_W:.3f}",
                "interface_superheat_K": f"{interface_dt:.2f}",
                "solid_base_superheat_K": f"{base_dt:.2f}",
                "interface_temperature_K": f"{TIN_K + interface_dt:.2f}",
                "solid_base_temperature_K": f"{TIN_K + base_dt:.2f}",
                "law_note": "interface dT uses actual AN wall-resolved CHT sweep row; fit exponent -0.7047; dP as U^1.8; pump as U^2.8",
            }
        )
    return rows


def build_geometry_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in GEOMETRY_MARKERS:
        dp = item["pressure_drop_Pa"]
        pump = None if dp is None else dp * Q0
        rows.append(
            {
                "geometry": item["geometry"],
                "case_pair": item["case_pair"],
                "pressure_drop_Pa": "" if dp is None else f"{dp:.2f}",
                "pump_power_W_at_U0": "" if pump is None else f"{pump:.4f}",
                "pump_power_relative_to_straight": "" if pump is None else f"{pump / PUMP0_W:.3f}",
                "nearwall_speed_m_s": f"{item['nearwall_speed_m_s']:.3f}",
                "passive_peak_dT_K_same_mesh": f"{item['passive_peak_dT_K']:.2f}",
                "passive_mean_dT_K": f"{item['passive_mean_dT_K']:.2f}",
                "pressure_source": item["pressure_source"],
                "thermal_status": "passive-scalar diagnostic; not trusted physical wall-temperature scale",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(tradeoff_rows: list[dict[str, str]], geometry_rows: list[dict[str, str]]) -> tuple[Path, Path]:
    pump = [float(row["pump_power_W"]) for row in tradeoff_rows]
    interface_dt = [float(row["interface_superheat_K"]) for row in tradeoff_rows]
    base_dt = [float(row["solid_base_superheat_K"]) for row in tradeoff_rows]
    reynolds = [row["Re"] for row in tradeoff_rows]

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6), gridspec_kw={"width_ratios": [1.18, 1.0]})
    fig.patch.set_facecolor("white")

    ax = axes[0]
    ax.plot(pump, interface_dt, marker="o", color="#1971c2", label="interface superheat")
    ax.plot(pump, base_dt, marker="s", color="#e67700", label="solid-base superheat")
    for x, y, re in zip(pump, interface_dt, reynolds):
        ax.text(x, y + 1.2, f"Re {re}", fontsize=7.5, ha="center")
    ax.set_xscale("log")
    ax.set_xlabel("pumping power [W], straight-channel scaling")
    ax.set_ylabel("temperature rise above inlet [K]")
    ax.set_title("Fitted CHT Cooling Tradeoff at 100 kW/m2")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    ax2 = axes[1]
    valid = [row for row in geometry_rows if row["pump_power_W_at_U0"]]
    missing = [row for row in geometry_rows if not row["pump_power_W_at_U0"]]
    geom_pump = [float(row["pump_power_W_at_U0"]) for row in valid]
    passive_peak = [float(row["passive_peak_dT_K_same_mesh"]) for row in valid]
    labels = [row["geometry"] for row in valid]
    ax2.scatter(geom_pump, passive_peak, s=72, color="#5c940d")
    for x, y, label in zip(geom_pump, passive_peak, labels):
        if label == "90-deg bend":
            offset, ha, va = (8, -13), "left", "center"
        elif label == "straight":
            offset, ha, va = (8, 10), "left", "center"
        else:
            offset, ha, va = (8, 8), "left", "center"
        ax2.annotate(label, xy=(x, y), xytext=offset, textcoords="offset points", fontsize=8, ha=ha, va=va)
    if missing:
        ax2.text(
            0.04,
            0.90,
            "BFS pressure marker omitted:\nsign check unresolved.",
            transform=ax2.transAxes,
            fontsize=7.8,
            color="#c92a2a",
            va="top",
        )
    if geom_pump:
        ax2.set_xlim(min(geom_pump) - 0.03, max(geom_pump) + 0.055)
    ax2.set_ylim(58.2, 79.0)
    ax2.set_xlabel("pumping power at Re=10,000 [W]")
    ax2.set_ylabel("passive peak dT [K]")
    ax2.set_title("Existing Geometry Markers")
    ax2.grid(True, alpha=0.25)

    fig.text(
        0.015,
        0.02,
        "AP now uses the AN wall-resolved CHT fit for the design curve. Geometry markers retain passive peak dT only as same-mesh diagnostics.",
        fontsize=9,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])

    out = OUT_DIR / "pumping_tradeoff_ap.png"
    advisor = ADVISOR_DIR / "34_pumping_tradeoff.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)
    return out, advisor


def main() -> None:
    tradeoff_rows = build_tradeoff_rows()
    geometry_rows = build_geometry_rows()
    write_csv(OUT_DIR / "pumping_tradeoff_ap.csv", tradeoff_rows)
    write_csv(OUT_DIR / "geometry_pumping_markers_ap.csv", geometry_rows)
    fig_out, advisor_out = plot(tradeoff_rows, geometry_rows)
    print(OUT_DIR / "pumping_tradeoff_ap.csv")
    print(OUT_DIR / "geometry_pumping_markers_ap.csv")
    print(fig_out)
    print(advisor_out)


if __name__ == "__main__":
    main()

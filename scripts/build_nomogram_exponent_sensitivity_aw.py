#!/usr/bin/env python3
"""Build AW sensitivity of AQ sizing to the fitted -0.7047 exponent."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/branch_d"
ADVISOR = ROOT / "results/report_figures"

OUT.mkdir(parents=True, exist_ok=True)
ADVISOR.mkdir(parents=True, exist_ok=True)

RHO = 1940.0
MU = 0.006
DH = 0.02667
PR = 14.4
K_FLUID = 1.0
TIN = 800.0
SOLID_THICKNESS = 0.004
K_SOLID = 30.0
RE0 = 10_000.0
U0 = 1.1598
DT0_CHT = 28.160944
Q0 = 100_000.0
PUMP0_W = 1.6464
PUMP_EXPONENT = 2.8
VALID_RE_MIN = 5_000.0
VALID_RE_MAX = 20_000.0

HEAT_FLUXES = [100_000, 150_000, 250_000, 350_000, 500_000]
TARGET_BASE_T = [873.0, 923.0, 973.0]


def u_from_re(reynolds: float) -> float:
    return reynolds * MU / (RHO * DH)


def pump_from_re(reynolds: float) -> float:
    return PUMP0_W * (reynolds / RE0) ** PUMP_EXPONENT


def db_re_for_allowable_dt(q: float, allowable_dt: float) -> float:
    nu_required = q * DH / (K_FLUID * allowable_dt)
    return (nu_required / (0.023 * PR**0.4)) ** (1.0 / 0.8)


def anchored_re_for_allowable_dt(q: float, allowable_dt: float, exponent: float) -> float:
    dt_at_re0 = DT0_CHT * (q / Q0)
    ratio = allowable_dt / dt_at_re0
    return RE0 * ratio ** (1.0 / exponent)


def status_for_re(raw_re: float) -> str:
    if raw_re < VALID_RE_MIN:
        return "below envelope; use Re 5000 minimum"
    if raw_re > VALID_RE_MAX:
        return "above project envelope"
    return "within project envelope"


def row_for_model(target: float, q: float, model: str, raw_re: float, allowable_dt: float, solid_drop: float) -> dict[str, str]:
    recommended_re = max(raw_re, VALID_RE_MIN)
    return {
        "target_solid_base_T_K": f"{target:.0f}",
        "heat_flux_kW_m2": f"{q / 1000:.0f}",
        "solid_drop_K": f"{solid_drop:.2f}",
        "allowable_interface_dT_K": f"{allowable_dt:.2f}",
        "model": model,
        "required_Re_raw": f"{raw_re:.0f}",
        "recommended_Re": f"{recommended_re:.0f}",
        "recommended_U_m_s": f"{u_from_re(recommended_re):.3f}",
        "recommended_pump_power_W": f"{pump_from_re(recommended_re):.3f}",
        "pump_relative_to_Re10000": f"{pump_from_re(recommended_re) / PUMP0_W:.2f}",
        "validity_status": status_for_re(raw_re),
    }


def impossible_row(target: float, q: float, model: str, allowable_dt: float, solid_drop: float) -> dict[str, str]:
    return {
        "target_solid_base_T_K": f"{target:.0f}",
        "heat_flux_kW_m2": f"{q / 1000:.0f}",
        "solid_drop_K": f"{solid_drop:.2f}",
        "allowable_interface_dT_K": f"{allowable_dt:.2f}",
        "model": model,
        "required_Re_raw": "",
        "recommended_Re": "",
        "recommended_U_m_s": "",
        "recommended_pump_power_W": "",
        "pump_relative_to_Re10000": "",
        "validity_status": "impossible: solid conduction drop alone exceeds target margin",
    }


def build_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    models = [
        ("AQ_Dittus_Boelter_full", None),
        ("same_anchor_exponent_-0.8000", -0.8),
        ("same_anchor_exponent_-0.7047", -0.7047),
    ]
    for target in TARGET_BASE_T:
        for q in HEAT_FLUXES:
            solid_drop = q * SOLID_THICKNESS / K_SOLID
            allowable_dt = target - TIN - solid_drop
            for model, exponent in models:
                if allowable_dt <= 0:
                    rows.append(impossible_row(target, q, model, allowable_dt, solid_drop))
                    continue
                if exponent is None:
                    raw_re = db_re_for_allowable_dt(q, allowable_dt)
                else:
                    raw_re = anchored_re_for_allowable_dt(q, allowable_dt, exponent)
                rows.append(row_for_model(target, q, model, raw_re, allowable_dt, solid_drop))
    return rows


def write_csv(rows: list[dict[str, str]]) -> Path:
    out = OUT / "nomogram_exponent_sensitivity_aw.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return out


def plot(rows: list[dict[str, str]]) -> tuple[Path, Path]:
    selected = [
        row
        for row in rows
        if row["target_solid_base_T_K"] == "873" and row["recommended_Re"]
    ]
    models = [
        ("AQ_Dittus_Boelter_full", "#495057", "AQ Dittus-Boelter"),
        ("same_anchor_exponent_-0.8000", "#1971c2", "same anchor, -0.800"),
        ("same_anchor_exponent_-0.7047", "#c92a2a", "same anchor, -0.705"),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.5), sharex=True)
    fig.patch.set_facecolor("white")

    for model, color, label in models:
        model_rows = [row for row in selected if row["model"] == model]
        qvals = [float(row["heat_flux_kW_m2"]) for row in model_rows]
        revals = [float(row["recommended_Re"]) for row in model_rows]
        pump = [float(row["recommended_pump_power_W"]) for row in model_rows]
        ax1.plot(qvals, revals, marker="o", label=label, color=color)
        ax2.plot(qvals, pump, marker="o", label=label, color=color)

    ax1.axhspan(VALID_RE_MIN, VALID_RE_MAX, color="#d0ebff", alpha=0.35, label="project Re envelope")
    ax1.set_yscale("log")
    ax1.set_xlabel("heat flux [kW/m2]")
    ax1.set_ylabel("recommended Re")
    ax1.set_title("873 K Solid-Base Target")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=7.5)

    ax2.set_yscale("log")
    ax2.set_xlabel("heat flux [kW/m2]")
    ax2.set_ylabel("hydraulic pumping power [W]")
    ax2.set_title("Pumping Sensitivity")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=7.5)

    fig.suptitle("AW: AQ Nomogram Sensitivity To Fitted CHT Exponent", fontsize=13, weight="bold")
    fig.text(
        0.5,
        0.028,
        "Same-anchor curves isolate exponent sensitivity around the Re=10000 wall-resolved CHT point; AQ curve shows the original Dittus-Boelter calculator.",
        ha="center",
        fontsize=8.5,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])

    out = OUT / "nomogram_exponent_sensitivity_aw.png"
    advisor = ADVISOR / "38_nomogram_exponent_sensitivity.png"
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

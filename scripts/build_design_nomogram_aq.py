#!/usr/bin/env python3
"""Build the AQ design nomogram/calculator table and figure."""

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
PUMP0_W = 1.6464
PUMP_EXPONENT = 2.8
VALID_RE_MIN = 5_000.0
VALID_RE_MAX = 20_000.0

HEAT_FLUXES = [100_000, 150_000, 250_000, 350_000, 500_000]
TARGET_BASE_T = [873.0, 923.0, 973.0]


def dittus_boelter_re_for_nu(nu: float) -> float:
    return (nu / (0.023 * PR**0.4)) ** (1.0 / 0.8)


def u_from_re(reynolds: float) -> float:
    return reynolds * MU / (RHO * DH)


def pump_from_re(reynolds: float) -> float:
    return PUMP0_W * (reynolds / RE0) ** PUMP_EXPONENT


def build_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for target in TARGET_BASE_T:
        for q in HEAT_FLUXES:
            solid_drop = q * SOLID_THICKNESS / K_SOLID
            allowable_film = target - TIN - solid_drop
            if allowable_film <= 0:
                rows.append(
                    {
                        "target_solid_base_T_K": f"{target:.0f}",
                        "heat_flux_kW_m2": f"{q / 1000:.0f}",
                        "solid_drop_K": f"{solid_drop:.2f}",
                        "allowable_interface_dT_K": f"{allowable_film:.2f}",
                        "required_Nu": "",
                        "required_Re_raw": "",
                        "recommended_Re": "",
                        "recommended_U_m_s": "",
                        "recommended_pump_power_W": "",
                        "pump_relative_to_Re10000": "",
                        "validity_status": "impossible: solid conduction drop alone exceeds target margin",
                    }
                )
                continue

            nu_required = q * DH / (K_FLUID * allowable_film)
            re_raw = dittus_boelter_re_for_nu(nu_required)
            recommended_re = max(re_raw, VALID_RE_MIN)
            status = "within project envelope"
            if re_raw < VALID_RE_MIN:
                status = "below envelope; use Re 5000 minimum for project scaling"
            elif re_raw > VALID_RE_MAX:
                status = "above project envelope; requires extrapolation or redesign"

            pump = pump_from_re(recommended_re)
            rows.append(
                {
                    "target_solid_base_T_K": f"{target:.0f}",
                    "heat_flux_kW_m2": f"{q / 1000:.0f}",
                    "solid_drop_K": f"{solid_drop:.2f}",
                    "allowable_interface_dT_K": f"{allowable_film:.2f}",
                    "required_Nu": f"{nu_required:.1f}",
                    "required_Re_raw": f"{re_raw:.0f}",
                    "recommended_Re": f"{recommended_re:.0f}",
                    "recommended_U_m_s": f"{u_from_re(recommended_re):.3f}",
                    "recommended_pump_power_W": f"{pump:.3f}",
                    "pump_relative_to_Re10000": f"{pump / PUMP0_W:.2f}",
                    "validity_status": status,
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(rows: list[dict[str, str]]) -> tuple[Path, Path]:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.5), sharex=True)
    fig.patch.set_facecolor("white")

    colors = {873.0: "#c92a2a", 923.0: "#e67700", 973.0: "#1971c2"}
    for target in TARGET_BASE_T:
        selected = [row for row in rows if float(row["target_solid_base_T_K"]) == target and row["recommended_Re"]]
        qvals = [float(row["heat_flux_kW_m2"]) for row in selected]
        revals = [float(row["recommended_Re"]) for row in selected]
        pump = [float(row["recommended_pump_power_W"]) for row in selected]
        ax1.plot(qvals, revals, marker="o", label=f"Tbase <= {target:.0f} K", color=colors[target])
        ax2.plot(qvals, pump, marker="o", label=f"Tbase <= {target:.0f} K", color=colors[target])

    ax1.axhspan(VALID_RE_MIN, VALID_RE_MAX, color="#d0ebff", alpha=0.35, label="project envelope")
    ax1.set_yscale("log")
    ax1.set_xlabel("heat flux [kW/m2]")
    ax1.set_ylabel("recommended Re")
    ax1.set_title("Required Flushing Speed")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=7.6)

    ax2.set_yscale("log")
    ax2.set_xlabel("heat flux [kW/m2]")
    ax2.set_ylabel("hydraulic pumping power [W]")
    ax2.set_title("Pumping Penalty")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=7.6)

    fig.suptitle("AQ Design Nomogram: Heat Flux + Target Solid Temperature -> Speed And Pumping", fontsize=13, weight="bold")
    fig.text(
        0.5,
        0.03,
        "Uses Dittus-Boelter Nu=0.023 Re^0.8 Pr^0.4, Tin=800 K, FLiBe-like Pr=14.4, solid drop q*t/k; hydraulic pump power scales from AP.",
        ha="center",
        fontsize=8.5,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])

    out = OUT / "design_nomogram_aq.png"
    advisor = ADVISOR / "36_design_nomogram.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)
    return out, advisor


def main() -> None:
    rows = build_rows()
    write_csv(OUT / "design_nomogram_aq.csv", rows)
    fig_out, advisor_out = plot(rows)
    print(OUT / "design_nomogram_aq.csv")
    print(fig_out)
    print(advisor_out)


if __name__ == "__main__":
    main()

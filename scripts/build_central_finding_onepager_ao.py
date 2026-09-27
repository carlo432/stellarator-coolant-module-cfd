#!/usr/bin/env python3
"""Build the AO central-finding summary figure."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
BRANCH = ROOT / "results/branch_d"
ADVISOR = ROOT / "results/report_figures"

ADVISOR.mkdir(parents=True, exist_ok=True)
BRANCH.mkdir(parents=True, exist_ok=True)

PUMP0_W = 1.6464
RE0 = 10_000.0
PUMP_EXPONENT = 2.8


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def plot() -> tuple[Path, Path]:
    flushing = read_rows(BRANCH / "branch_d_nearwall_flushing_law.csv")
    tradeoff = read_rows(BRANCH / "pumping_tradeoff_ap.csv")
    an_sweep_path = BRANCH / "cht_velocity_sweep_law_an.csv"
    an_summary_path = BRANCH / "cht_velocity_sweep_fit_summary_an.csv"

    fig = plt.figure(figsize=(12.4, 7.2))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.42], hspace=0.42, wspace=0.34)

    ax1 = fig.add_subplot(gs[0, 0])
    speeds = [float(row["nearwall_mean_speed_m_s"]) for row in flushing]
    peaks = [float(row["peak_wall_deltaT_K"]) for row in flushing]
    labels = [row["geometry"] for row in flushing]
    ax1.scatter(speeds, peaks, s=78, color="#2f9e44")
    for x, y, label in zip(speeds, peaks, labels):
        if label == "90-deg bend":
            offset = (-62, -11)
        elif label == "backward-facing step":
            offset = (8, -14)
        elif label == "outlet-offset":
            offset = (8, -8)
        else:
            offset = (8, 8)
        ax1.annotate(label, (x, y), xytext=offset, textcoords="offset points", fontsize=7.5)
    ax1.set_xlabel("near-wall mean speed [m/s]")
    ax1.set_ylabel("passive peak dT [K]")
    ax1.set_title("Velocity Mechanism")
    ax1.set_xlim(0.24, 0.95)
    ax1.set_ylim(58.7, 79.0)
    ax1.grid(True, alpha=0.25)
    ax1.text(0.03, 0.06, "passive peak = hotspot diagnostic", transform=ax1.transAxes, fontsize=7.6, color="#495057")

    ax2 = fig.add_subplot(gs[0, 1])
    names = ["passive\npeak", "CHT\nbridge", "Gnielinski"]
    values = [77.0, 28.16, 25.71]
    colors = ["#868e96", "#1971c2", "#12b886"]
    bars = ax2.bar(names, values, color=colors)
    for bar, val in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 2.0, f"{val:.1f} K", ha="center", fontsize=8)
    ax2.set_ylabel("interface / peak dT [K]")
    ax2.set_title("Which Thermal Number To Trust")
    ax2.set_ylim(0, 88)
    ax2.grid(True, axis="y", alpha=0.25)
    ax2.text(0.36, 0.82, "CHT is grid-confirmed\nand correlation-validated", transform=ax2.transAxes, fontsize=7.6, color="#1971c2")

    ax3 = fig.add_subplot(gs[0, 2])
    if an_sweep_path.exists() and an_summary_path.exists():
        an_rows = read_rows(an_sweep_path)
        an_summary = read_rows(an_summary_path)[0]
        pump = [PUMP0_W * (float(row["Re"]) / RE0) ** PUMP_EXPONENT for row in an_rows]
        interface = [float(row["dT_K"]) for row in an_rows]
        re_values = [row["Re"] for row in an_rows]
        exponent_note = f"CHT fit: dT ~ U^{float(an_summary['dT_exponent_n']):.2f}; Ppump ~ U^2.8"
    else:
        pump = [float(row["pump_power_W"]) for row in tradeoff]
        interface = [float(row["interface_superheat_K"]) for row in tradeoff]
        re_values = [row["Re"] for row in tradeoff]
        exponent_note = "design estimate: dT ~ U^-0.8; Ppump ~ U^2.8"
    ax3.plot(pump, interface, marker="o", color="#e67700")
    for x, y, re in zip(pump, interface, re_values):
        if re in {"5000", "10000", "20000"}:
            if re == "20000":
                offset = (-50, 8)
            else:
                offset = (6, 8)
            ax3.annotate(f"Re {re}", (x, y), xytext=offset, textcoords="offset points", fontsize=7.4)
    ax3.set_xscale("log")
    ax3.set_xlabel("hydraulic pumping power [W]")
    ax3.set_ylabel("interface dT [K]")
    ax3.set_title("Design Knob")
    ax3.set_xlim(min(pump) * 0.8, max(pump) * 1.45)
    ax3.set_ylim(14.2, 50.8)
    ax3.grid(True, which="both", alpha=0.25)
    ax3.text(0.03, 0.06, exponent_note, transform=ax3.transAxes, fontsize=7.6, color="#495057")

    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis("off")
    bullets = [
        "1. Geometry matters through near-wall flushing: slow, broad near-wall flow makes hotspot-prone regions.",
        "2. Passive scalar peak dT is retained only as a same-mesh hotspot diagnostic.",
        "3. Wall-resolved CHT is the trusted wall-temperature scale: 28.16/27.21 K at 100 kW/m2.",
        "4. Pumping tradeoff: Re 20k drops interface dT to 17.0 K but costs about 7x the Re 10k pump power.",
        "5. Next physics upgrades: matched-heated-length CHT, volumetric heating, and HPC-1 wall-resolved LES.",
    ]
    ax4.text(0.01, 0.95, "Five takeaways", fontsize=10, weight="bold", va="top")
    ax4.text(0.01, 0.78, "\n".join(bullets), fontsize=9, va="top", linespacing=1.55)

    fig.suptitle("Central Finding: Wall Temperature Is A Flushing-Speed Tradeoff", fontsize=15, weight="bold")
    fig.text(
        0.5,
        0.02,
        "Bridge/P4 CHT results are physical simplified-module results; Branch D passive peaks remain diagnostic context, not final wall-temperature values.",
        ha="center",
        fontsize=8.5,
        color="#495057",
    )
    fig.subplots_adjust(left=0.06, right=0.98, top=0.84, bottom=0.10, hspace=0.56, wspace=0.34)

    out = BRANCH / "central_finding_onepager_ao.png"
    advisor = ADVISOR / "35_central_finding_onepager.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)
    return out, advisor


def main() -> None:
    out, advisor = plot()
    print(out)
    print(advisor)


if __name__ == "__main__":
    main()

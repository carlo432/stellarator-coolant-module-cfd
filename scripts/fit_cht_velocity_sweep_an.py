#!/usr/bin/env python3
"""Fit the AN wall-resolved CHT velocity sweep."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "cases/branch_d/cht_wr_sweep_results.txt"
OUT = ROOT / "results/branch_d"
ADVISOR = ROOT / "results/report_figures"

OUT.mkdir(parents=True, exist_ok=True)
ADVISOR.mkdir(parents=True, exist_ok=True)

PR = 14.4
DH = 0.02667
K_FLUID = 1.0
QFLUX = 100_000.0


def read_sweep() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for line in SWEEP.read_text().splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        u, re, interface_t, dt, nu, yplus = text.split()
        rows.append(
            {
                "U_m_s": float(u),
                "Re": float(re),
                "interfaceT_K": float(interface_t),
                "dT_K": float(dt),
                "Nu_CHT": float(nu),
                "heated_yplus_avg": float(yplus),
            }
        )
    return rows


def log_fit(x: list[float], y: list[float]) -> tuple[float, float, float]:
    lx = [math.log(v) for v in x]
    ly = [math.log(v) for v in y]
    mx = sum(lx) / len(lx)
    my = sum(ly) / len(ly)
    sxx = sum((v - mx) ** 2 for v in lx)
    sxy = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    slope = sxy / sxx
    intercept = my - slope * mx
    predictions = [intercept + slope * v for v in lx]
    ss_res = sum((a - b) ** 2 for a, b in zip(ly, predictions))
    ss_tot = sum((v - my) ** 2 for v in ly)
    r2 = 1.0 - ss_res / ss_tot
    return math.exp(intercept), slope, r2


def db_nu(reynolds: float) -> float:
    return 0.023 * reynolds**0.8 * PR**0.4


def db_dt(reynolds: float) -> float:
    return QFLUX * DH / (db_nu(reynolds) * K_FLUID)


def build_outputs(rows: list[dict[str, float]]) -> tuple[list[dict[str, str]], dict[str, str]]:
    dT_c, dT_n, dT_r2 = log_fit([row["U_m_s"] for row in rows], [row["dT_K"] for row in rows])
    nu_c, nu_n, nu_r2 = log_fit([row["Re"] for row in rows], [row["Nu_CHT"] for row in rows])

    out_rows: list[dict[str, str]] = []
    for row in rows:
        fit_dt = dT_c * row["U_m_s"] ** dT_n
        fit_nu = nu_c * row["Re"] ** nu_n
        db = db_nu(row["Re"])
        db_temp = db_dt(row["Re"])
        out_rows.append(
            {
                "U_m_s": f"{row['U_m_s']:.4f}",
                "Re": f"{row['Re']:.0f}",
                "interfaceT_K": f"{row['interfaceT_K']:.5f}",
                "dT_K": f"{row['dT_K']:.2f}",
                "fit_dT_K": f"{fit_dt:.2f}",
                "Nu_CHT": f"{row['Nu_CHT']:.1f}",
                "fit_Nu": f"{fit_nu:.1f}",
                "Nu_Dittus_Boelter": f"{db:.1f}",
                "dT_Dittus_Boelter_K": f"{db_temp:.2f}",
                "dT_vs_DB_pct": f"{(row['dT_K'] / db_temp - 1.0) * 100.0:.1f}",
                "heated_yplus_avg": f"{row['heated_yplus_avg']:.3f}",
            }
        )

    summary = {
        "dT_prefactor_for_U_m_s": f"{dT_c:.4f}",
        "dT_exponent_n": f"{dT_n:.4f}",
        "dT_logfit_R2": f"{dT_r2:.5f}",
        "Nu_prefactor_for_Re": f"{nu_c:.6f}",
        "Nu_exponent_m": f"{nu_n:.4f}",
        "Nu_logfit_R2": f"{nu_r2:.5f}",
        "interpretation": "wall-resolved CHT sweep gives dT~U^n and Nu~Re^m; fitted exponent is shallower than Dittus-Boelter 0.8 but close over Re 5k-20k",
    }
    return out_rows, summary


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, summary: dict[str, str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)


def plot(rows: list[dict[str, float]], summary: dict[str, str]) -> tuple[Path, Path]:
    u = [row["U_m_s"] for row in rows]
    re = [row["Re"] for row in rows]
    dt = [row["dT_K"] for row in rows]
    nu = [row["Nu_CHT"] for row in rows]

    dT_c = float(summary["dT_prefactor_for_U_m_s"])
    dT_n = float(summary["dT_exponent_n"])
    dT_r2 = float(summary["dT_logfit_R2"])
    nu_c = float(summary["Nu_prefactor_for_Re"])
    nu_n = float(summary["Nu_exponent_m"])
    nu_r2 = float(summary["Nu_logfit_R2"])

    u_line = [min(u) * (max(u) / min(u)) ** (i / 80) for i in range(81)]
    re_line = [min(re) * (max(re) / min(re)) ** (i / 80) for i in range(81)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.4))
    fig.patch.set_facecolor("white")

    ax1.loglog(u, dt, "o", color="#1971c2", label="wall-resolved CHT sweep")
    ax1.loglog(u_line, [dT_c * val**dT_n for val in u_line], "-", color="#1971c2", alpha=0.75, label=f"fit: dT ~ U^{dT_n:.2f}")
    ax1.loglog(u_line, [db_dt(val / u[2] * re[2]) for val in u_line], "--", color="#868e96", label="Dittus-Boelter -0.8")
    for row in rows:
        ax1.annotate(f"Re {row['Re']:.0f}", (row["U_m_s"], row["dT_K"]), xytext=(6, 6), textcoords="offset points", fontsize=7.5)
    ax1.set_xlabel("U [m/s]")
    ax1.set_ylabel("interface dT [K]")
    ax1.set_title("Wall-Resolved CHT Flushing Law")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8)
    ax1.text(0.05, 0.08, f"log-fit R2 = {dT_r2:.4f}", transform=ax1.transAxes, fontsize=8, color="#495057")

    ax2.loglog(re, nu, "o", color="#e67700", label="Nu_CHT")
    ax2.loglog(re_line, [nu_c * val**nu_n for val in re_line], "-", color="#e67700", alpha=0.75, label=f"fit: Nu ~ Re^{nu_n:.2f}")
    ax2.loglog(re_line, [db_nu(val) for val in re_line], "--", color="#868e96", label="Dittus-Boelter")
    ax2.set_xlabel("Re")
    ax2.set_ylabel("Nu")
    ax2.set_title("Heat-Transfer Coefficient Scaling")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=8)
    ax2.text(0.05, 0.08, f"log-fit R2 = {nu_r2:.4f}", transform=ax2.transAxes, fontsize=8, color="#495057")

    fig.suptitle("AN Quantitative Law From Wall-Resolved CHT Sweep", fontsize=14, weight="bold")
    fig.text(
        0.5,
        0.025,
        "Sweep is straight-channel bridge CHT. It validates the CHT/correlation scale, but does not yet make Branch D geometry CHT complete.",
        ha="center",
        fontsize=8.5,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])

    out = OUT / "cht_velocity_sweep_law_an.png"
    advisor = ADVISOR / "37_cht_velocity_sweep_law.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)
    return out, advisor


def main() -> None:
    rows = read_sweep()
    out_rows, summary = build_outputs(rows)
    write_csv(OUT / "cht_velocity_sweep_law_an.csv", out_rows)
    write_summary(OUT / "cht_velocity_sweep_fit_summary_an.csv", summary)
    fig_out, advisor_out = plot(rows, summary)
    print(OUT / "cht_velocity_sweep_law_an.csv")
    print(OUT / "cht_velocity_sweep_fit_summary_an.csv")
    print(fig_out)
    print(advisor_out)


if __name__ == "__main__":
    main()

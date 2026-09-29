#!/usr/bin/env python3
"""
Overlay a C1 thermal-LES mean profile against the digitized Kawamura Pr=5 DNS.
Standalone: reads the local digitization in kawamura_pr5_dns_digitized/.

Usage:
  overlay_c1_vs_kawamura.py <profile_csv> [--label "..."] [--tag diagnostic|cert] [--out name.png]
"""
import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = Path(__file__).resolve().parent
DIG = WORK / "kawamura_pr5_dns_digitized"


def read_col(path, xk, yk):
    xs, ys = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                x, y = float(r[xk]), float(r[yk])
            except (ValueError, KeyError):
                continue
            if math.isfinite(x) and math.isfinite(y):
                xs.append(x); ys.append(y)
    return xs, ys


def kader_beta(pr):
    return (3.85 * pr ** (1 / 3) - 1.3) ** 2 + 2.12 * math.log(pr)


ap = argparse.ArgumentParser()
ap.add_argument("profile_csv")
ap.add_argument("--label", default="C1 LES")
ap.add_argument("--tag", default="diagnostic")
ap.add_argument("--out", default=None)
ap.add_argument("--pr", type=float, default=5.0)
args = ap.parse_args()

les_y, les_T = read_col(args.profile_csv, "y_plus", "T_plus")
les_yp, les_prt = read_col(args.profile_csv, "y_plus", "Pr_t_resolved")
dns_y, dns_T = read_col(DIG / "kawamura_fig2_mean_temperature.csv", "y_plus", "Theta_plus")
dns_yp, dns_prt = read_col(DIG / "kawamura_fig9_turbulent_prandtl.csv", "y_plus", "Pr_t")

import numpy as np
yk = np.logspace(0, math.log10(max(dns_y)), 150)
beta = kader_beta(args.pr)
kad = [2.12 * math.log(y) + beta for y in yk]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
ax1.semilogx(dns_y, dns_T, "o", mfc="none", color="k", label="Kawamura DNS Pr=5 (digitized)")
ax1.semilogx(yk, kad, "-", color="tab:red", lw=1, label="Kader (1981)")
ax1.semilogx([y for y in yk if y < 6], [args.pr * y for y in yk if y < 6], "--",
             color="tab:blue", lw=1, label=r"$\Theta^+=Pr\,y^+$")
ax1.semilogx(les_y, les_T, "-s", ms=3, color="tab:green", label=args.label)
ax1.set_xlabel(r"$y^+$"); ax1.set_ylabel(r"$\Theta^+ / H^+$")
ax1.set_title(f"Mean temperature [{args.tag}]")
ax1.set_ylim(0, max(max(les_T, default=50), 50) * 1.05)
ax1.grid(True, which="both", alpha=0.3); ax1.legend(fontsize=8)

if les_prt:
    ax2.semilogx(dns_yp, dns_prt, "o", mfc="none", color="k", label="Kawamura DNS Pr=5")
    ax2.semilogx(les_yp, les_prt, "-s", ms=3, color="tab:green", label=args.label)
    ax2.axhline(0.85, ls=":", color="gray")
    ax2.set_xlabel(r"$y^+$"); ax2.set_ylabel(r"$Pr_t$")
    ax2.set_title(f"Turbulent Prandtl [{args.tag}]")
    ax2.set_ylim(0, 1.6); ax2.grid(True, which="both", alpha=0.3); ax2.legend(fontsize=8)

out = args.out or f"c1_overlay_{args.tag}.png"
outp = (WORK / "kawamura_pr5_dns_digitized" / out)
fig.tight_layout(); fig.savefig(outp, dpi=130)
print("wrote", outp)

# quick numeric comparison at the log layer + centerline
if les_y and dns_y:
    def interp(xs, ys, xq):
        for i in range(1, len(xs)):
            if xs[i - 1] <= xq <= xs[i]:
                t = (xq - xs[i - 1]) / (xs[i] - xs[i - 1])
                return ys[i - 1] + t * (ys[i] - ys[i - 1])
        return ys[-1] if xq > xs[-1] else ys[0]
    print(f"centerline-ish:  LES T+({les_y[-1]:.0f})={les_T[-1]:.1f}   DNS T+(180)~{dns_T[-1]:.1f}")
    for yq in (10, 30, 50):
        print(f"  y+={yq:3d}:  LES={interp(les_y, les_T, yq):6.2f}   DNS={interp(dns_y, dns_T, yq):6.2f}")

#!/usr/bin/env python3
"""
Vision-assisted digitization of Kawamura et al. (1998) thermal-channel DNS,
Pr = 5.0, Re_tau = 180 -- the C1 certification overlay targets.

Provenance / honesty:
  - Points were read visually (AI-assisted) directly off the rendered figure crops
    under references/digitization/kawamura_1998/figures/.
  - Accuracy ~ +/- 1 Theta+ unit on the mean-temperature curve; good enough for a
    V&V *overlay* (does the LES track or miss the DNS), NOT publication-grade.
    Refine with WebPlotDigitizer if a tighter tolerance is ever needed.
  - Fig 2 (mean T+) and Fig 9 (Pr_t) were clearly legible and the Pr=5 symbol was
    unambiguous. Fig 3 is derived EXACTLY as Theta+/Pr from the Fig 2 reading
    (Kawamura plots the same data, just normalized), so the two are self-consistent.
  - Fig 6 (theta_rms+/Pr) is TOP-CLIPPED in the current crop (peak off-frame) and
    Fig 8 (-v'T'+) has a symbol-ID ambiguity for the Pr=5 series at high y+.
    Those two are deliberately NOT digitized here -- they need a full-page re-render
    + a careful WebPlotDigitizer pass before use. See README in the output dir.

Cross-check provided: Kader (1981) correlation and the conduction sublayer
asymptote Theta+ = Pr*y+, so the overlay can show DNS points bracketed by the
analytical reference the same way Kawamura's own Fig 2 does.

Writes template-matched CSVs (drop-in for references/digitization/.../digitized/)
plus a combined overlay CSV and a sanity PNG.
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PR = 5.0
RE_TAU = 180
DIGITIZER = "ai_assisted_visual_digitization_2026-06-25"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "kawamura_pr5_dns_digitized")
os.makedirs(OUTDIR, exist_ok=True)

# --- Fig 2: mean temperature Theta+ vs y+ (Pr=5), read off the log-region crop ---
FIG2 = [
    (1.0, 5.0), (1.5, 7.1), (2.0, 9.2), (3.0, 12.9), (4.0, 16.6),
    (5.0, 20.4), (6.0, 23.5), (7.0, 26.4), (8.5, 29.8), (10.0, 32.5),
    (12.0, 35.0), (15.0, 37.7), (20.0, 40.0), (25.0, 41.2), (30.0, 41.9),
    (40.0, 42.7), (50.0, 43.1), (75.0, 43.9), (100.0, 44.2),
    (140.0, 44.5), (178.0, 44.6),
]

# --- Fig 9: turbulent Prandtl number Pr_t vs y+ (Pr=5), read off the crop ---
FIG9 = [
    (1.0, 1.05), (2.0, 1.00), (3.0, 0.97), (5.0, 0.93), (7.0, 0.92),
    (10.0, 0.93), (15.0, 0.96), (20.0, 0.98), (30.0, 1.00), (40.0, 0.95),
    (50.0, 0.90), (75.0, 0.78), (100.0, 0.70), (130.0, 0.62), (160.0, 0.60),
]

PAGE3 = "references/digitization/kawamura_1998/pages/kawamura_1998_page-03.png"
PAGE5 = "references/digitization/kawamura_1998/pages/kawamura_1998_page-05.png"
NOTE = "vision-assisted approx (+/- ~1 Theta+ unit); refine with WebPlotDigitizer if publication-grade needed"


def kader_beta(pr):
    # Kader (1981): beta(Pr) = (3.85*Pr^(1/3) - 1.3)^2 + 2.12*ln(Pr)
    return (3.85 * pr ** (1.0 / 3.0) - 1.3) ** 2 + 2.12 * math.log(pr)


def kader_theta(yp, pr):
    # log-region: Theta+ = 2.12 ln(y+) + beta(Pr)
    return 2.12 * math.log(yp) + kader_beta(pr)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    print("wrote", path)


# Fig 2 -- mean temperature
write_csv(
    os.path.join(OUTDIR, "kawamura_fig2_mean_temperature.csv"),
    ["series_pr", "y_plus", "Theta_plus", "source_page_png", "digitizer", "notes"],
    [[PR, yp, th, PAGE3, DIGITIZER, NOTE] for yp, th in FIG2],
)

# Fig 3 -- near-wall, Theta+/Pr (derived exactly from Fig 2; identical data, normalized)
write_csv(
    os.path.join(OUTDIR, "kawamura_fig3_near_wall_temperature.csv"),
    ["series_pr", "y_plus", "Theta_plus_over_Pr", "source_page_png", "digitizer", "notes"],
    [[PR, yp, round(th / PR, 4), PAGE3, DIGITIZER,
      "derived = Theta+/Pr from Fig 2 digitization (self-consistent)"] for yp, th in FIG2],
)

# Fig 9 -- turbulent Prandtl number
write_csv(
    os.path.join(OUTDIR, "kawamura_fig9_turbulent_prandtl.csv"),
    ["series_pr", "y_plus", "Pr_t", "source_page_png", "digitizer", "notes"],
    [[PR, yp, prt, PAGE5, DIGITIZER, NOTE] for yp, prt in FIG9],
)

# Combined overlay-ready CSV with Kader + sublayer reference alongside the DNS
beta5 = kader_beta(PR)
with open(os.path.join(OUTDIR, "kawamura_pr5_meanT_overlay.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["y_plus", "Theta_plus_dns", "Theta_plus_kader", "Theta_plus_sublayer_Pr_y"])
    for yp, th in FIG2:
        w.writerow([yp, th, round(kader_theta(yp, PR), 3), round(PR * yp, 3)])
print("Kader beta(Pr=5) =", round(beta5, 3))

# Sanity plot -- the way Kawamura's own Fig 2 is drawn
yp2 = [p[0] for p in FIG2]
th2 = [p[1] for p in FIG2]
import numpy as np
yk = np.logspace(0, math.log10(180), 200)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
ax1.semilogx(yp2, th2, "o", mfc="none", color="k", label="Kawamura DNS Pr=5 (digitized)")
ax1.semilogx(yk, [kader_theta(y, PR) for y in yk], "-", color="tab:red", label="Kader (1981)")
ax1.semilogx(yk[yk < 6], [PR * y for y in yk[yk < 6]], "--", color="tab:blue",
             label=r"$\Theta^+=Pr\,y^+$ (conduction)")
ax1.set_xlabel(r"$y^+$"); ax1.set_ylabel(r"$\Theta^+$")
ax1.set_title(r"Fig 2 mean temperature, Pr=5, $Re_\tau=180$")
ax1.set_ylim(0, 48); ax1.grid(True, which="both", alpha=0.3); ax1.legend(fontsize=8)

yp9 = [p[0] for p in FIG9]; pr9 = [p[1] for p in FIG9]
ax2.semilogx(yp9, pr9, "o", mfc="none", color="k", label="Kawamura DNS Pr=5 (digitized)")
ax2.axhline(0.85, ls=":", color="gray", label=r"$Pr_t=0.85$ (Reynolds analogy)")
ax2.set_xlabel(r"$y^+$"); ax2.set_ylabel(r"$Pr_t$")
ax2.set_title(r"Fig 9 turbulent Prandtl, Pr=5, $Re_\tau=180$")
ax2.set_ylim(0, 1.4); ax2.grid(True, which="both", alpha=0.3); ax2.legend(fontsize=8)
fig.tight_layout()
png = os.path.join(OUTDIR, "kawamura_pr5_dns_digitized_sanity.png")
fig.savefig(png, dpi=130)
print("wrote", png)

#!/usr/bin/env python3
"""Assemble the Branch D fidelity ladder panel from existing figures."""

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]

PANELS = [
    (
        "A. RANS outlet-offset velocity",
        ROOT / "results/branch_d/v22_bent_midplane_velocity.png",
    ),
    (
        "B. Coarse LES vs URANS probe",
        ROOT / "results/branch_d/les_vs_urans_recirc_probe.png",
    ),
    (
        "C. Wall dT: straight vs offset",
        ROOT / "results/branch_d/wall_dT_bent_vs_straight.png",
    ),
    (
        "D. Heat-flux scaling",
        ROOT / "results/branch_d/branch_d_heat_flux_sensitivity.png",
    ),
]


def main() -> None:
    missing = [str(path) for _, path in PANELS if not path.exists()]
    if missing:
        raise SystemExit("Missing input figure(s):\n" + "\n".join(missing))

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=180)
    fig.patch.set_facecolor("white")

    for ax, (title, path) in zip(axes.ravel(), PANELS):
        ax.imshow(mpimg.imread(path))
        ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#333333")
            spine.set_linewidth(0.8)

    fig.suptitle(
        "Branch D Fidelity Ladder: RANS Flow, LES Probe Dynamics, Wall Response, Heat-Flux Scaling",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.012,
        "LES panel is a coarse wall-modeled smoke test: short, isothermal, not wall-resolved, not grid-converged, and not validated.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#333333",
    )
    fig.tight_layout(rect=(0.015, 0.035, 0.985, 0.955))

    out = ROOT / "results/branch_d/branch_d_fidelity_ladder_panel.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")

    advisor_out = ROOT / "results/report_figures/21_branch_d_fidelity_ladder_panel.png"
    advisor_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(advisor_out, bbox_inches="tight")

    print(out.relative_to(ROOT))
    print(advisor_out.relative_to(ROOT))


if __name__ == "__main__":
    main()

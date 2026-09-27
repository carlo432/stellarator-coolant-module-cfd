#!/usr/bin/env python3
"""Render advisor-safe status artifacts after the v2 near-field resolution."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ADVISOR = ROOT / "results" / "advisor_packet"
WORK_FIGS = ROOT / "cases" / "branch_d" / "figs"


def render_v2_status_card() -> None:
    fig = plt.figure(figsize=(14.2, 8.0), facecolor="white")
    fig.text(0.055, 0.92, "Plasma load model v2: near-field audit resolved", fontsize=25, weight="bold")
    fig.text(
        0.055,
        0.865,
        "The July 9 withdrawal was correct claim control; the July 10 cell-extent integration supplies the replacement result.",
        fontsize=14,
        color="#444444",
    )

    fig.text(0.07, 0.76, "RESOLVED RESULT", fontsize=17, weight="bold", color="#176b52")
    left = [
        "Radiative peak/mean = 1.09: approximately flat.",
        "Corrected v2 total peak/mean = 3.14 nominal.",
        "Total-load band = 2.29-3.99 over radiated fraction 0.7-0.3.",
        "Nominal v2 matches endpoint-audited v1 peak/mean = 3.142.",
        "The completed v1 hotspot band remains the thermal headline.",
    ]
    for i, line in enumerate(left):
        fig.text(0.07, 0.70 - 0.075 * i, f"• {line}", fontsize=13.2)

    fig.text(0.56, 0.76, "STILL STANDS", fontsize=17, weight="bold", color="#285c8c")
    right = [
        "Mean field-line incidence ≈ 3.5°; parallel deposition is suppressed.",
        "The exponential channel is cross-field far-SOL transport.",
        "Cross-field transport peak/mean ≈ 5.3-5.6.",
        "Neutron wall load varies 0.69-1.69× mean.",
        "A banded neutron source is optional consistency work, not a hotspot gate.",
    ]
    for i, line in enumerate(right):
        fig.text(0.56, 0.70 - 0.075 * i, f"• {line}", fontsize=13.2)

    fig.text(0.07, 0.27, "SUPERSEDED HISTORICAL CLAIMS", fontsize=15, weight="bold", color="#812c26")
    fig.text(
        0.07,
        0.215,
        "Point-source radiative peak 3.21-3.26; v2 total 3.997-4.241; 27-35% extra peaking; and the ~3.1-3.3× hotspot extrapolation.",
        fontsize=12.5,
        color="#333333",
    )
    fig.text(
        0.07,
        0.12,
        "Operational note: the v2 surface-load rerun is optional because nominal corrected v2 equals v1. See advisor figure 101 for convergence.",
        fontsize=12.5,
        weight="bold",
    )
    fig.text(0.07, 0.065, "Source: near-field cell-extent sub-quadrature audit, 2026-07-10", fontsize=10.5, color="#666666")
    fig.savefig(ADVISOR / "99_plasma_qmap_v2.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    sidecar = ADVISOR / "99_plasma_qmap_v2.json"
    data = json.loads(sidecar.read_text())
    data["resolution_status"] = {
        "status": "superseded_by_converged_nearfield_result",
        "replacement_advisor_figure": "101_plasma_qmap_v2_nearfield_corrected.png",
        "radiative_peak_over_mean": 1.0884941153992669,
        "total_peak_over_mean_nominal": 3.139150836351965,
        "total_peak_over_mean_band": [2.2918203187760966, 3.986481353927833],
        "note": "Raw point-source arrays are retained below only as withdrawn provenance.",
    }
    sidecar.write_text(json.dumps(data, indent=2) + "\n")


def render_probe_status() -> None:
    probe_path = ADVISOR / "100_radiative_nearfield_resolution_probe.json"
    probe = json.loads(probe_path.read_text())
    labels = [row["label"].replace("FWHM = ", "") for row in probe["probe_curves"]]
    values = [row["radiative_peak_over_mean"] for row in probe["probe_curves"]]
    corrected = 1.0884941153992669

    fig, (ax, note) = plt.subplots(1, 2, figsize=(12.4, 5.0), facecolor="white", gridspec_kw={"width_ratios": [1.35, 1]})
    colors = ["#aa382b", "#5d8cc1", "#4f9d69", "#7b6fb0"]
    ax.bar(np.arange(len(values)), values, color=colors)
    ax.axhline(corrected, color="#176b52", lw=2.0, ls="--", label="converged cell-extent result 1.09")
    ax.set_xticks(np.arange(len(values)), labels, rotation=17, ha="right")
    ax.set_ylabel("radiative peak / mean")
    ax.set_title("Output-smoothing probes: diagnostic upper-bound family", loc="left")
    ax.grid(axis="y", alpha=0.28)
    ax.legend(fontsize=8.8)

    note.axis("off")
    note.set_title("Superseded by source-side integration", loc="left", fontsize=12.5)
    lines = [
        "3.21, 1.41, 1.22, and 1.17 were never corrected answers.",
        "They smooth the output isotropically and only demonstrate",
        "that the point-source spike is length-scale dominated.",
        "",
        "The physical correction integrates each near plasma cell",
        "over its actual extent before evaluating wall flux.",
        "It converges to peak/mean 1.09.",
        "",
        "Use advisor figure 101 for the promoted v2 result.",
    ]
    note.text(0.0, 0.88, "\n".join(lines), va="top", fontsize=11.2, linespacing=1.45)
    fig.suptitle("Historical radiative near-field probe: retained as a diagnostic, not a model", fontsize=15)
    fig.tight_layout()
    fig.savefig(ADVISOR / "100_radiative_nearfield_resolution_probe.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    probe["status"] = "diagnostic_upper_bound_family_superseded_by_converged_result"
    probe["interpretation"] = (
        "The 3.21/1.41/1.22/1.17 output-smoothed values are diagnostic upper-bound brackets, "
        "not corrected estimates. Cell-extent source integration supersedes them."
    )
    probe["resolved_result"] = {
        "radiative_peak_over_mean": corrected,
        "advisor_figure": "101_plasma_qmap_v2_nearfield_corrected.png",
    }
    probe_path.write_text(json.dumps(probe, indent=2) + "\n")


def promote_new_results() -> None:
    shutil.copy2(
        WORK_FIGS / "plasma_qmap_v2_nearfield.png",
        ADVISOR / "101_plasma_qmap_v2_nearfield_corrected.png",
    )
    shutil.copy2(
        WORK_FIGS / "plasma_qmap_v2_nearfield.json",
        ADVISOR / "101_plasma_qmap_v2_nearfield_corrected.json",
    )
    shutil.copy2(
        WORK_FIGS / "throat_flushing_verdict.png",
        ADVISOR / "102_throat_flushing_verdict.png",
    )


def main() -> None:
    ADVISOR.mkdir(parents=True, exist_ok=True)
    render_v2_status_card()
    render_probe_status()
    promote_new_results()
    print("refreshed advisor figures 99-100 and promoted figures 101-102")


if __name__ == "__main__":
    main()

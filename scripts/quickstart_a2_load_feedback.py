#!/usr/bin/env python3
"""Regenerate a CFD-free A2 load-feedback figure from packaged CSV data."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-revengence-quickstart")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def find_package_dir() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "05_DATA" / "a2_standoff_load_feedback.csv").exists():
            return parent
        candidate = parent / "proposal"
        if (candidate / "05_DATA" / "a2_standoff_load_feedback.csv").exists():
            return candidate
    raise FileNotFoundError("Could not locate proposal/05_DATA/a2_standoff_load_feedback.csv")


def parse_args():
    package = find_package_dir()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=package / "05_DATA" / "a2_standoff_load_feedback.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=package / "03_RESULTS" / "figures" / "quickstart_a2_load_feedback.png",
    )
    return parser.parse_args()


def peak_over_mean(values: np.ndarray) -> float:
    return float(np.max(values) / np.mean(values))


def main() -> None:
    args = parse_args()
    data = np.genfromtxt(args.input, delimiter=",", names=True)
    required = {
        "phi_rad",
        "baseline_standoff_m",
        "coil_traced_standoff_m",
        "traced_3cm_qtotal_W_m2",
        "conditional_qtotal_W_m2",
        "censored_lower_qtotal_W_m2",
    }
    missing = required.difference(data.dtype.names or ())
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    phi = data["phi_rad"]
    s_mm = (phi - phi.min()) / (phi.max() - phi.min()) * 298.4513
    q_base = data["traced_3cm_qtotal_W_m2"]
    q_finite = data["conditional_qtotal_W_m2"]
    q_censored = data["censored_lower_qtotal_W_m2"]

    for name, values in (
        ("baseline", q_base),
        ("finite", q_finite),
        ("censored", q_censored),
    ):
        if not np.isclose(np.mean(values), 0.5e6, rtol=2e-3):
            raise ValueError(f"{name} mean load is not power-normalized: {np.mean(values)}")

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4), facecolor="white")
    ax = axes[0]
    ax.plot(
        s_mm,
        100 * data["baseline_standoff_m"],
        color="0.45",
        ls="--",
        lw=1.8,
        label="fixed boundary",
    )
    ax.plot(
        s_mm,
        100 * data["coil_traced_standoff_m"],
        color="#247a74",
        lw=2.2,
        label="direct field-line envelope",
    )
    ax.set_xlabel("mapped duct coordinate s [mm]")
    ax.set_ylabel("standoff [cm]")
    ax.set_title("A2 standoff reconstruction", loc="left")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1]
    series = (
        (q_base, "0.45", "--", "traced envelope, 3.0 cm"),
        (q_finite, "#a53632", "-", "finite subset 5.9 cm"),
        (q_censored, "#247a74", "-", "censored lower bound >=11.6 cm"),
    )
    for values, color, linestyle, label in series:
        ax.plot(
            s_mm,
            values / 1e6,
            color=color,
            ls=linestyle,
            lw=2.1,
            label=f"{label} (pk/mean {peak_over_mean(values):.2f})",
        )
    ax.set_xlabel("mapped duct coordinate s [mm]")
    ax.set_ylabel("surface load [MW/m2]")
    ax.set_title("Reduced-SOL width feedback", loc="left")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    fig.suptitle("Day-1 reproducibility check: packaged A2 data", x=0.06, ha="left")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"input rows: {len(data)}")
    print(
        "peak/mean: "
        f"{peak_over_mean(q_base):.3f}, "
        f"{peak_over_mean(q_finite):.3f}, "
        f"{peak_over_mean(q_censored):.3f}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot the C1 transformed Pr=5 Stage 1 shakedown Re_tau sweep."""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "c1_kawamura_stage1_shakedown.csv"
OUT = ROOT / "results" / "c1_kawamura_stage1_shakedown.png"


def as_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib.pyplot as plt

    rows = []
    with CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["run_status"] != "completed":
                continue
            rows.append(
                {
                    "case": row["case"],
                    "u_bar": as_float(row["u_bar"]),
                    "re_tau": as_float(row["measured_Re_tau"]),
                }
            )

    rows = [row for row in rows if math.isfinite(row["u_bar"]) and math.isfinite(row["re_tau"])]
    rows.sort(key=lambda row: row["u_bar"])

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot([row["u_bar"] for row in rows], [row["re_tau"] for row in rows], "o-", lw=1.8)
    ax.axhline(180.0, color="tab:red", ls="--", lw=1.2, label="Kawamura target")
    for row in rows:
        label = row["case"].replace("les_channel_c1_kawamura_transformed_pr5_re180_stage1_shakedown_", "")
        ax.annotate(
            label,
            (row["u_bar"], row["re_tau"]),
            textcoords="offset points",
            xytext=(4, 6),
            fontsize=8,
        )
    ax.set_xlabel("meanVelocityForce target u_bar [m/s]")
    ax.set_ylabel("measured Re_tau from harness")
    ax.set_title("C1 transformed Pr=5 Stage 1 shakedown calibration")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    ax.margins(x=0.12, y=0.18)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=180)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

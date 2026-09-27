#!/usr/bin/env python3
"""Plot the C1 transformed Pr=5 short Re_tau calibration sweep."""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "c1_kawamura_re180_calibration.csv"
OUT = ROOT / "results" / "c1_kawamura_re180_calibration.png"


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
    max_u_bar = max(row["u_bar"] for row in rows)
    for row in rows:
        label = row["case"].replace("les_channel_c1_kawamura_transformed_pr5_", "")
        xytext = (-58, 6) if row["u_bar"] == max_u_bar else (4, 6)
        ax.annotate(
            label,
            (row["u_bar"], row["re_tau"]),
            textcoords="offset points",
            xytext=xytext,
            fontsize=8,
        )
    ax.set_xlabel("meanVelocityForce target u_bar [m/s]")
    ax.set_ylabel("measured Re_tau from harness")
    ax.set_title("C1 transformed Pr=5 short Re_tau calibration")
    ax.margins(x=0.07, y=0.12)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=180)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

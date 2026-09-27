#!/usr/bin/env python3
"""Track written-time velocity development for the ARC 2D OpenFOAM case."""
from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def numeric_times(case: Path) -> list[tuple[float, Path]]:
    out = []
    for p in case.iterdir():
        if not p.is_dir():
            continue
        try:
            out.append((float(p.name), p))
        except ValueError:
            pass
    return sorted(out)


def latest_time(case: Path) -> Path:
    times = numeric_times(case)
    if not times:
        raise FileNotFoundError(f"no numeric time directories in {case}")
    return times[-1][1]


def read_vector_field(path: Path) -> np.ndarray:
    txt = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    if not m:
        raise ValueError(f"could not parse vector field {path}")
    n = int(m.group(1))
    vals = re.findall(r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", m.group(2))
    arr = np.array(vals, dtype=float)
    if len(arr) != n:
        raise ValueError(f"{path}: expected {n} vectors, parsed {len(arr)}")
    return arr


def field_metrics(path: Path) -> dict[str, float]:
    u = read_vector_field(path)
    mag = np.linalg.norm(u[:, :2], axis=1)
    return {
        "mean": float(np.mean(mag)),
        "max": float(np.max(mag)),
        "p95": float(np.percentile(mag, 95.0)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d")
    ap.add_argument("--rans-case", default="cases/arc_blanket_2d_rans_baseline")
    ap.add_argument("--outdir", default="results/arc_blanket_2d")
    args = ap.parse_args()

    case = Path(args.case)
    rans_case = Path(args.rans_case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for t, tdir in numeric_times(case):
        u_path = tdir / "U"
        if not u_path.exists() or t == 0.0:
            continue
        u = field_metrics(u_path)
        row = {
            "time_s": t,
            "U_mean_m_per_s": u["mean"],
            "U_p95_m_per_s": u["p95"],
            "U_max_m_per_s": u["max"],
            "UMean_mean_m_per_s": "",
            "UMean_p95_m_per_s": "",
            "UMean_max_m_per_s": "",
        }
        umean_path = tdir / "UMean"
        if umean_path.exists():
            um = field_metrics(umean_path)
            row.update(
                {
                    "UMean_mean_m_per_s": um["mean"],
                    "UMean_p95_m_per_s": um["p95"],
                    "UMean_max_m_per_s": um["max"],
                }
            )
        rows.append(row)

    if not rows:
        raise RuntimeError(f"no U fields found in {case}")

    rans_mean = None
    rans_time = None
    if rans_case.exists():
        rtdir = latest_time(rans_case)
        rans_time = rtdir.name
        rans_mean = field_metrics(rtdir / "U")["mean"]

    csv_path = outdir / "arc2d_development_metrics.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    times = np.array([float(r["time_s"]) for r in rows])
    u_mean = np.array([float(r["U_mean_m_per_s"]) for r in rows])
    u_p95 = np.array([float(r["U_p95_m_per_s"]) for r in rows])
    u_max = np.array([float(r["U_max_m_per_s"]) for r in rows])
    umean_times = np.array([float(r["time_s"]) for r in rows if r["UMean_mean_m_per_s"] != ""])
    umean_mean = np.array([float(r["UMean_mean_m_per_s"]) for r in rows if r["UMean_mean_m_per_s"] != ""])

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 6.2), sharex=True, constrained_layout=True)
    ax.plot(times, u_mean, "o-", label="instantaneous domain mean |U|", color="#0b7a75")
    ax.plot(times, u_p95, "s-", label="instantaneous p95 |U|", color="#255da8", alpha=0.85)
    if len(umean_times):
        ax.plot(umean_times, umean_mean, "D--", label="fieldAverage mean |UMean|", color="#7a3b9a")
    if rans_mean is not None:
        ax.axhline(rans_mean, color="0.25", ls=":", lw=1.5, label=f"RANS mean |U| iter {rans_time}")
    ax.set_ylabel("mean / p95 |U| [m/s]")
    ax.set_title("ARC 2D Velocity Development Tracker", fontweight="bold")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    ax2.plot(times, u_max, "^-", label="instantaneous max |U|", color="#b22821", alpha=0.85)
    ax2.set_xlabel("physical time [s]")
    ax2.set_ylabel("max |U| [m/s]")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="best", fontsize=8)
    png_path = outdir / "arc2d_development_metrics.png"
    fig.savefig(png_path, dpi=160)
    plt.close(fig)

    latest = rows[-1]
    summary = outdir / "arc2d_development_metrics_summary.md"
    summary.write_text(
        f"""# ARC 2D Development Metrics

Case: `{case}`

Latest written time: `{latest['time_s']} s`

Latest instantaneous domain mean `|U|`: `{float(latest['U_mean_m_per_s']):.6g} m/s`

Latest instantaneous p95 `|U|`: `{float(latest['U_p95_m_per_s']):.6g} m/s`

Latest instantaneous max `|U|`: `{float(latest['U_max_m_per_s']):.6g} m/s`

Latest fieldAverage mean `|UMean|`: `{latest['UMean_mean_m_per_s']}` m/s

RANS domain mean `|U|`: `{'' if rans_mean is None else f'{rans_mean:.6g}'}` m/s

Outputs:

- `{png_path}`
- `{csv_path}`

This tracks solver development only. `UMean` is only a Leffler-style averaging-window result if the case `fieldAverage` `timeStart` matches the intended nondimensional window.
"""
    )
    print(f"wrote {png_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

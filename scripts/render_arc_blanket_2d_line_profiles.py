#!/usr/bin/env python3
"""Slide-approximate line profiles for the local Leffler ARC 2D case."""
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


LINE_SETS = {
    "visible-bars": {
    # Approximate sampling bars from Leffler slide 12, mapped onto the current
    # slide-traced OpenFOAM geometry.  These are profile scaffolds, not
    # digitized thesis/CAD coordinates.
        "left_mid": ((-1.86, 0.00), (-0.82, 0.00), "s from outer wall [m]"),
        "right_mid": ((0.86, 0.00), (2.30, 0.00), "s from inner side [m]"),
    # Slide 16 plots the top/bottom line coordinates over about 14 coordinate
    # units (180-194 and -194--180), so these are short slot bars near the
    # top/bottom openings rather than long cuts through the central void.
        "top_slot": ((0.30, 1.94), (0.30, 1.80), "s downward [m]"),
        "bottom_slot": ((0.30, -1.94), (0.30, -1.80), "s upward [m]"),
    },
    "axis-matched": {
        # Diagnostic variant from Leffler slide-16 plot-axis spans.  If the
        # slide coordinates are read as centimeters, the horizontal spans are
        # about 0.25 m (left) and 0.60 m (right), not the longer visible-bar
        # cuts used above.  These are a convention test, not a confirmed CAD
        # coordinate recovery.
        "left_mid": ((-1.86, 0.00), (-1.61, 0.00), "axis-matched x [m]"),
        "right_mid": ((1.75, 0.00), (2.35, 0.00), "axis-matched x [m]"),
        "top_slot": ((0.30, 1.94), (0.30, 1.80), "axis-matched y [m]"),
        "bottom_slot": ((0.30, -1.94), (0.30, -1.80), "axis-matched y [m]"),
    },
    "slide12-detected": {
        # Same geometry spans as axis-matched, now backed by the slide-12
        # sampling-bar detector.  The detected pixel length ratios are
        # right/left ~= 2.44 and top/left ~= bottom/left ~= 0.56, consistent
        # with slide-16 coordinate spans of 60, 25, and 14 units.
        "left_mid": ((-1.86, 0.00), (-1.61, 0.00), "slide-12 detected x [m]"),
        "right_mid": ((1.75, 0.00), (2.35, 0.00), "slide-12 detected x [m]"),
        "top_slot": ((0.30, 1.94), (0.30, 1.80), "slide-12 detected y [m]"),
        "bottom_slot": ((0.30, -1.94), (0.30, -1.80), "slide-12 detected y [m]"),
    },
}


def latest_time(case: Path) -> Path:
    times = []
    for p in case.iterdir():
        if not p.is_dir():
            continue
        try:
            times.append((float(p.name), p))
        except ValueError:
            pass
    if not times:
        raise FileNotFoundError(f"no numeric time directories in {case}")
    return sorted(times)[-1][1]


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


def sample_nearest(
    xy: np.ndarray,
    field: np.ndarray,
    a: tuple[float, float],
    b: tuple[float, float],
    n: int,
    max_distance: float,
) -> dict[str, np.ndarray]:
    a_arr = np.array(a, dtype=float)
    b_arr = np.array(b, dtype=float)
    frac = np.linspace(0.0, 1.0, n)
    pts = a_arr[None, :] + frac[:, None] * (b_arr - a_arr)[None, :]
    idx = []
    nearest_dist = []
    for p in pts:
        d2 = np.sum((xy - p[None, :]) ** 2, axis=1)
        j = int(np.argmin(d2))
        idx.append(j)
        nearest_dist.append(float(np.sqrt(d2[j])))
    idx_arr = np.asarray(idx, dtype=int)
    u = field[idx_arr, :2]
    valid = np.asarray(nearest_dist) <= max_distance
    umag = np.linalg.norm(u, axis=1)
    ux = u[:, 0].copy()
    uy = u[:, 1].copy()
    umag[~valid] = np.nan
    ux[~valid] = np.nan
    uy[~valid] = np.nan
    return {
        "s": frac * float(np.linalg.norm(b_arr - a_arr)),
        "x": pts[:, 0],
        "y": pts[:, 1],
        "nearest_dist": np.asarray(nearest_dist),
        "valid": valid,
        "umag": umag,
        "ux": ux,
        "uy": uy,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d")
    ap.add_argument("--outdir", default="results/arc_blanket_2d")
    ap.add_argument("--samples", type=int, default=160)
    ap.add_argument("--max-distance", type=float, default=0.06)
    ap.add_argument("--line-set", choices=sorted(LINE_SETS), default="visible-bars")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    case = Path(args.case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tdir = latest_time(case)

    c = read_vector_field(tdir / "C")
    umean_path = tdir / "UMean"
    u = read_vector_field(umean_path if umean_path.exists() else tdir / "U")
    xy = c[:, :2]

    rows = []
    sampled = {}
    lines = LINE_SETS[args.line_set]
    for name, (a, b, xlabel) in lines.items():
        data = sample_nearest(xy, u, a, b, args.samples, args.max_distance)
        sampled[name] = (data, xlabel)
        for i in range(args.samples):
            rows.append(
                {
                    "line": name,
                    "s_m": data["s"][i],
                    "x_m": data["x"][i],
                    "y_m": data["y"][i],
                    "nearest_cell_distance_m": data["nearest_dist"][i],
                    "valid": bool(data["valid"][i]),
                    "UMean_mag_m_per_s": data["umag"][i],
                    "UMean_x_m_per_s": data["ux"][i],
                    "UMean_y_m_per_s": data["uy"][i],
                }
            )

    suffix = f"_{args.tag}" if args.tag else ""
    csv_path = outdir / f"arc2d_line_profiles_approx{suffix}.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(4, 3, figsize=(12, 10), sharex=False, constrained_layout=True)
    for row, (name, (data, xlabel)) in enumerate(sampled.items()):
        axes[row, 0].plot(data["s"], data["umag"], color="#0b7a75", lw=1.8)
        axes[row, 1].plot(data["s"], data["ux"], color="#b22821", lw=1.8)
        axes[row, 2].plot(data["s"], data["uy"], color="#255da8", lw=1.8)
        axes[row, 0].set_ylabel(name)
        for col, title in enumerate(("|UMean| [m/s]", "UMean_x [m/s]", "UMean_y [m/s]")):
            axes[row, col].set_title(title if row == 0 else "", fontsize=10, fontweight="bold")
            axes[row, col].grid(True, alpha=0.25)
            axes[row, col].set_xlabel(xlabel)
    fig.suptitle(
        f"ARC 2D Slide-Approximate Mean-Velocity Line Profiles ({args.line_set}) at t={tdir.name}s",
        fontsize=13,
        fontweight="bold",
    )
    png_path = outdir / f"arc2d_line_profiles_approx{suffix}.png"
    fig.savefig(png_path, dpi=160)
    plt.close(fig)

    max_nearest = max(float(np.max(data["nearest_dist"])) for data, _ in sampled.values())
    valid_count = sum(int(np.count_nonzero(data["valid"])) for data, _ in sampled.values())
    total_count = sum(len(data["valid"]) for data, _ in sampled.values())
    summary = outdir / f"arc2d_line_profiles_approx{suffix}_summary.md"
    summary.write_text(
        f"""# ARC 2D Approximate Line Profiles

Case: `{case}`

Latest time: `{tdir.name} s`

Inputs:

- field sampled: `{'UMean' if umean_path.exists() else 'U'}`
- line set: `{args.line_set}`
- samples per line: `{args.samples}`
- valid nearest-cell distance cutoff: `{args.max_distance:.6g} m`
- valid samples: `{valid_count} / {total_count}`
- maximum nearest-cell sampling distance: `{max_nearest:.6g} m`

Outputs:

- `{png_path}`
- `{csv_path}`

The four lines are approximate mappings of the visible black sampling bars on Leffler slide 12 onto the current slide-traced geometry. They are suitable for workflow images and profile scaffolding, not for a formal digitized comparison against Leffler's original line coordinates.
"""
    )
    print(f"wrote {png_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

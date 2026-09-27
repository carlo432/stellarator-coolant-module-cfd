#!/usr/bin/env python3
"""Quick preview plots for the local Leffler ARC 2D OpenFOAM case."""
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def latest_time(case: Path) -> Path:
    times = []
    for p in case.iterdir():
        if not p.is_dir():
            continue
        try:
            val = float(p.name)
        except ValueError:
            continue
        times.append((val, p))
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


def panel(ax, tri, values, title, label, cmap="viridis", symmetric=False) -> None:
    values = np.asarray(values)
    if symmetric:
        lim = float(np.nanpercentile(np.abs(values), 99.5))
        levels = np.linspace(-lim, lim, 31)
        extend = "both"
    else:
        lo = float(np.nanpercentile(values, 0.5))
        hi = float(np.nanpercentile(values, 99.5))
        if hi <= lo:
            hi = lo + 1e-12
        levels = np.linspace(lo, hi, 31)
        extend = "both"
    cf = ax.tricontourf(tri, np.clip(values, levels[0], levels[-1]), levels=levels, cmap=cmap, extend=extend)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    cb = plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label(label)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d")
    ap.add_argument("--outdir", default="results/arc_blanket_2d")
    args = ap.parse_args()

    case = Path(args.case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tdir = latest_time(case)
    c = read_vector_field(tdir / "C")
    u = read_vector_field(tdir / "U")
    umean = read_vector_field(tdir / "UMean") if (tdir / "UMean").exists() else u

    x = c[:, 0]
    y = c[:, 1]
    tri = mtri.Triangulation(x, y)
    umag = np.linalg.norm(u[:, :2], axis=1)
    ummag = np.linalg.norm(umean[:, :2], axis=1)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    panel(axes[0, 0], tri, umag, f"Instantaneous |U| at t={tdir.name}s", "|U| [m/s]", "turbo")
    panel(axes[0, 1], tri, u[:, 0], "Instantaneous Ux", "Ux [m/s]", "coolwarm", symmetric=True)
    panel(axes[1, 0], tri, u[:, 1], "Instantaneous Uy", "Uy [m/s]", "coolwarm", symmetric=True)
    panel(axes[1, 1], tri, ummag, "Mean |U| smoke average", "|UMean| [m/s]", "turbo")
    fig.suptitle("ARC Blanket 2D OpenFOAM Smoke Preview", fontsize=14, fontweight="bold")
    png = outdir / "arc2d_velocity_preview.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    summary = outdir / "arc2d_velocity_preview_summary.md"
    summary.write_text(
        f"""# ARC 2D Velocity Preview

Case: `{case}`

Latest time: `{tdir.name} s`

Cells plotted: `{len(x)}`

Velocity magnitude:

- instantaneous min/mean/max: `{umag.min():.6g}` / `{umag.mean():.6g}` / `{umag.max():.6g}` m/s
- mean-field min/mean/max: `{ummag.min():.6g}` / `{ummag.mean():.6g}` / `{ummag.max():.6g}` m/s

Figure:

- `{png}`

This is a staged preview field. It is not a full Leffler-equivalent long-window average unless the case controls explicitly set the matching averaging window.
"""
    )
    print(f"wrote {png}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

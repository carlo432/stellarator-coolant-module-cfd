#!/usr/bin/env python3
"""Compare ARC 2D unsteady mean velocity against the steady RANS baseline."""
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
            times.append((float(p.name), p))
        except ValueError:
            continue
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


def plot_panel(ax, tri, values, title, label, cmap, levels=None, symmetric=False) -> None:
    values = np.asarray(values)
    if levels is None:
        if symmetric:
            lim = float(np.nanpercentile(np.abs(values), 99.5))
            levels = np.linspace(-lim, lim, 31)
        else:
            lo = float(np.nanpercentile(values, 0.5))
            hi = float(np.nanpercentile(values, 99.5))
            if hi <= lo:
                hi = lo + 1e-12
            levels = np.linspace(lo, hi, 31)
    cf = ax.tricontourf(tri, np.clip(values, levels[0], levels[-1]), levels=levels, cmap=cmap, extend="both")
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    cb = plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label(label)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unsteady-case", default="cases/arc_blanket_2d")
    ap.add_argument("--rans-case", default="cases/arc_blanket_2d_rans_baseline")
    ap.add_argument("--outdir", default="results/arc_blanket_2d")
    args = ap.parse_args()

    unsteady_case = Path(args.unsteady_case)
    rans_case = Path(args.rans_case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    utdir = latest_time(unsteady_case)
    rtdir = latest_time(rans_case)
    c = read_vector_field(utdir / "C")
    u_unsteady = read_vector_field(utdir / "UMean" if (utdir / "UMean").exists() else utdir / "U")
    u_rans = read_vector_field(rtdir / "U")
    if len(u_unsteady) != len(u_rans):
        raise ValueError(f"cell-count mismatch: unsteady {len(u_unsteady)} vs RANS {len(u_rans)}")

    x = c[:, 0]
    y = c[:, 1]
    tri = mtri.Triangulation(x, y)
    mag_unsteady = np.linalg.norm(u_unsteady[:, :2], axis=1)
    mag_rans = np.linalg.norm(u_rans[:, :2], axis=1)

    mag_levels = np.linspace(
        float(np.nanpercentile(np.r_[mag_unsteady, mag_rans], 0.5)),
        float(np.nanpercentile(np.r_[mag_unsteady, mag_rans], 99.5)),
        31,
    )
    ux_lim = float(np.nanpercentile(np.abs(np.r_[u_unsteady[:, 0], u_rans[:, 0]]), 99.5))
    uy_lim = float(np.nanpercentile(np.abs(np.r_[u_unsteady[:, 1], u_rans[:, 1]]), 99.5))
    ux_levels = np.linspace(-ux_lim, ux_lim, 31)
    uy_levels = np.linspace(-uy_lim, uy_lim, 31)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), constrained_layout=True)
    plot_panel(axes[0, 0], tri, mag_unsteady, f"Unsteady mean |U| t={utdir.name}s", "|UMean| [m/s]", "turbo", mag_levels)
    plot_panel(axes[0, 1], tri, u_unsteady[:, 0], "Unsteady mean Ux", "UMean_x [m/s]", "coolwarm", ux_levels)
    plot_panel(axes[0, 2], tri, u_unsteady[:, 1], "Unsteady mean Uy", "UMean_y [m/s]", "coolwarm", uy_levels)
    plot_panel(axes[1, 0], tri, mag_rans, f"RANS |U| iter={rtdir.name}", "|U| [m/s]", "turbo", mag_levels)
    plot_panel(axes[1, 1], tri, u_rans[:, 0], "RANS Ux", "Ux [m/s]", "coolwarm", ux_levels)
    plot_panel(axes[1, 2], tri, u_rans[:, 1], "RANS Uy", "Uy [m/s]", "coolwarm", uy_levels)
    fig.suptitle("ARC 2D OpenFOAM Velocity Comparison: Unsteady Mean vs RANS", fontsize=14, fontweight="bold")
    png = outdir / "arc2d_rans_vs_unsteady_velocity.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    diff_mag = mag_unsteady - mag_rans
    diff_vec = u_unsteady[:, :2] - u_rans[:, :2]
    summary = outdir / "arc2d_rans_vs_unsteady_velocity_summary.md"
    summary.write_text(
        f"""# ARC 2D RANS vs Unsteady Velocity Comparison

Unsteady case: `{unsteady_case}` at `{utdir.name} s`

RANS case: `{rans_case}` at iteration/time `{rtdir.name}`

Cells compared: `{len(c)}`

Metrics:

- mean `|UMean|`: `{float(np.mean(mag_unsteady)):.6g} m/s`
- mean RANS `|U|`: `{float(np.mean(mag_rans)):.6g} m/s`
- mean signed `|UMean|-|U_RANS|`: `{float(np.mean(diff_mag)):.6g} m/s`
- mean absolute vector difference: `{float(np.mean(np.linalg.norm(diff_vec, axis=1))):.6g} m/s`
- max absolute vector difference: `{float(np.max(np.linalg.norm(diff_vec, axis=1))):.6g} m/s`

Figure:

- `{png}`

This is a visual baseline. The unsteady field is a staged OpenFOAM mean at `{utdir.name} s`, and the RANS case is an iteration/time `{rtdir.name}` baseline rather than a strict convergence study.
"""
    )
    print(f"wrote {png}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

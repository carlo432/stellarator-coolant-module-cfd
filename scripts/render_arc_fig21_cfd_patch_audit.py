#!/usr/bin/env python3
"""Render named inlet/outlet patch locations for the anchor-exact ARC CFD case."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PATCH_COLORS = {
    "outlet180": "#7b3294",
    "inlet_aux100": "#008837",
    "inlet_ch_a": "#d7191c",
    "inlet_main40": "#fdae61",
    "inlet_ch_b": "#2c7bb6",
}


def parse_points(geo: str) -> dict[int, tuple[float, float]]:
    return {
        int(m.group(1)): (float(m.group(2)), float(m.group(3)))
        for m in re.finditer(
            r"Point\((\d+)\)\s*=\s*\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,",
            geo,
        )
    }


def parse_lines(geo: str) -> dict[int, tuple[int, int]]:
    return {
        int(m.group(1)): (int(m.group(2)), int(m.group(3)))
        for m in re.finditer(r"Line\((\d+)\)\s*=\s*\{\s*(\d+)\s*,\s*(\d+)\s*\};", geo)
    }


def parse_loop(geo: str, loop_id: int) -> list[int]:
    m = re.search(rf"Curve Loop\({loop_id}\)\s*=\s*\{{([^}}]*)\}};", geo)
    if not m:
        raise ValueError(f"Curve Loop({loop_id}) not found")
    return [abs(int(x.strip())) for x in m.group(1).split(",") if x.strip()]


def parse_physical_surfaces(geo: str) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    for m in re.finditer(r'Physical Surface\("([^"]+)"\)\s*=\s*\{([^}]*)\};', geo):
        name = m.group(1)
        # ext[2 + (curve - 1)] maps to side surface index. Invert it.
        curves = []
        for ext_id in re.findall(r"ext\[(\d+)\]", m.group(2)):
            curves.append(int(ext_id) - 1)
        out[name] = curves
    return out


def line_xy(points: dict[int, tuple[float, float]], lines: dict[int, tuple[int, int]], lid: int) -> np.ndarray:
    a, b = lines[lid]
    return np.asarray([points[a], points[b]])


def polygon_xy(points: dict[int, tuple[float, float]], lines: dict[int, tuple[int, int]], lids: list[int]) -> np.ndarray:
    xy = [points[lines[lids[0]][0]]]
    for lid in lids:
        xy.append(points[lines[lid][1]])
    return np.asarray(xy)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_fig21_anchor_exact_cfd_candidate")
    ap.add_argument("--outdir", default="results/arc_fig21_anchor_exact_cfd_candidate")
    args = ap.parse_args()

    case = Path(args.case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    geo_path = case / "geometry" / "arc_fig21_anchor_exact_cfd.geo"
    geo = geo_path.read_text()
    points = parse_points(geo)
    lines = parse_lines(geo)
    outer = polygon_xy(points, lines, parse_loop(geo, 1))
    inner = polygon_xy(points, lines, parse_loop(geo, 2))
    patches = parse_physical_surfaces(geo)

    fig, ax = plt.subplots(figsize=(9, 10.5), constrained_layout=True)
    ax.plot(outer[:, 0], outer[:, 1], color="#17365d", lw=2.0, label="anchor-exact tank contour")
    ax.plot(inner[:, 0], inner[:, 1], color="#8b1e1e", lw=1.6, label="VV/legs/feet exclusion contour")
    ax.fill(inner[:, 0], inner[:, 1], color="#f4cccc", alpha=0.25)

    rows = [
        "# ARC Fig. 2.1 Anchor-Exact CFD Patch Audit",
        "",
        f"Case: `{case}`",
        "",
        "| Patch | Curve IDs | Length [m] | Midpoint `(x,y)` [m] |",
        "|---|---:|---:|---:|",
    ]
    for name, color in PATCH_COLORS.items():
        ids = patches.get(name, [])
        total = 0.0
        mids = []
        for lid in ids:
            xy = line_xy(points, lines, lid)
            total += float(np.linalg.norm(xy[1] - xy[0]))
            mids.append((xy[0] + xy[1]) * 0.5)
            ax.plot(xy[:, 0], xy[:, 1], color=color, lw=5.0, solid_capstyle="round")
        if mids:
            mid = np.mean(np.asarray(mids), axis=0)
            ax.scatter([mid[0]], [mid[1]], color=color, s=36, zorder=5)
            ax.annotate(name, mid, xytext=(6, 6), textcoords="offset points", fontsize=8.5, color=color)
            rows.append(f"| `{name}` | `{ids}` | `{total:.5f}` | `({mid[0]:.4f}, {mid[1]:.4f})` |")

    ax.axvline(0, color="0.55", ls="-.", lw=0.8)
    ax.axhline(0, color="0.55", ls="-.", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("local x = R - 3.3 m [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("Ferrero Fig. 2.1 Anchor-Exact CFD Patch Audit")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper right", fontsize=8)
    out = outdir / "arc_fig21_anchor_exact_cfd_patch_audit.png"
    fig.savefig(out, dpi=170)
    plt.close(fig)

    rows += [
        "",
        f"Figure: `{out}`",
        "",
        "Caveat: VV/channel inlet cuts are provisional nearest-contour cuts until the exact Ferrero/Leffler inlet locations are approved or recovered.",
        "",
    ]
    summary = outdir / "arc_fig21_anchor_exact_cfd_patch_audit_summary.md"
    summary.write_text("\n".join(rows))
    print(f"wrote {out}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

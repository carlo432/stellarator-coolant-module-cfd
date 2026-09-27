#!/usr/bin/env python3
"""Render Ferrero source geometry beside the generated dimensioned ARC outline."""

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
from PIL import Image


def parse_geo_points(geo_path: Path) -> dict[int, tuple[float, float]]:
    text = geo_path.read_text()
    points: dict[int, tuple[float, float]] = {}
    for m in re.finditer(
        r"Point\((\d+)\)\s*=\s*\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,",
        text,
    ):
        points[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return points


def parse_geo_lines(geo_path: Path) -> list[tuple[int, int, int]]:
    text = geo_path.read_text()
    lines: list[tuple[int, int, int]] = []
    for m in re.finditer(r"Line\((\d+)\)\s*=\s*\{\s*(\d+)\s*,\s*(\d+)\s*\};", text):
        lines.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
    return lines


def parse_curve_loop(geo_path: Path, loop_id: int) -> list[int]:
    text = geo_path.read_text()
    m = re.search(rf"Curve Loop\({loop_id}\)\s*=\s*\{{([^}}]*)\}};", text)
    if not m:
        raise ValueError(f"Curve Loop({loop_id}) not found in {geo_path}")
    return [abs(int(x.strip())) for x in m.group(1).split(",") if x.strip()]


def ordered_polygon(points: dict[int, tuple[float, float]], lines: list[tuple[int, int, int]], line_ids: list[int]) -> np.ndarray:
    line_by_id = {line[0]: line for line in lines}
    segment_lines = [line_by_id[lid] for lid in line_ids]
    ordered = [points[segment_lines[0][1]]]
    for _, _, b in segment_lines:
        ordered.append(points[b])
    return np.asarray(ordered)


def panel_source(ax, source: Path) -> None:
    img = Image.open(source)
    ax.imshow(img)
    ax.set_title("Ferrero thesis Fig. 2.1 source crop", fontweight="bold")
    ax.set_axis_off()


def panel_generated(ax, geo: Path) -> None:
    points = parse_geo_points(geo)
    lines = parse_geo_lines(geo)
    outer_ids = parse_curve_loop(geo, 1)
    inner_ids = parse_curve_loop(geo, 2)
    outer = ordered_polygon(points, lines, outer_ids)
    inner = ordered_polygon(points, lines, inner_ids)

    ax.plot(outer[:, 0], outer[:, 1], color="#1f4e79", lw=2.2, label="generated tank outline")
    ax.fill(outer[:, 0], outer[:, 1], color="#d9eaf7", alpha=0.25)
    ax.plot(inner[:, 0], inner[:, 1], color="#9c2f2f", lw=2.0, label="provisional VV/divertor exclusion")
    ax.fill(inner[:, 0], inner[:, 1], color="#f4cccc", alpha=0.35)

    for x, label in [(-1.779, "R=1521"), (0.0, "R=3300"), (2.503, "R=5803")]:
        ax.axvline(x, color="0.65", lw=0.8, ls="--")
        ax.text(x, 2.05, label, ha="center", va="bottom", fontsize=8)
    for y, label in [(1.938, "Z half=1938"), (0.0, "midplane"), (-1.938, "-1938")]:
        ax.axhline(y, color="0.75", lw=0.8, ls=":")
        ax.text(2.58, y, label, ha="left", va="center", fontsize=8)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-2.05, 2.8)
    ax.set_ylim(-2.15, 2.25)
    ax.set_xlabel("local x = R - 3.3 m [m]")
    ax.set_ylabel("vertical y [m]")
    ax.set_title("Generated Ferrero-dimensioned scaffold", fontweight="bold")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="lower center", fontsize=8)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d_ferrero_dimensioned_candidate")
    ap.add_argument("--outdir", default="results/arc_blanket_2d_ferrero_dimensioned_candidate")
    ap.add_argument(
        "--source",
        default="references/digitization/ferrero_arc/fig2_1_arc_geometric_measures_crop.png",
    )
    args = ap.parse_args()

    case = Path(args.case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    source = Path(args.source)
    geo = case / "geometry" / "arc_blanket_2d.geo"

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    panel_source(axes[0], source)
    panel_generated(axes[1], geo)
    fig.suptitle("Ferrero ARC Geometry Recovery Audit", fontsize=15, fontweight="bold")
    out = outdir / "ferrero_dimensioned_geometry_audit.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)

    summary = outdir / "ferrero_dimensioned_geometry_audit_summary.md"
    summary.write_text(
        f"""# Ferrero Dimensioned Geometry Audit

Case: `{case}`

Source crop: `{source}`

Generated geometry: `{geo}`

Figure:

- `{out}`

Interpretation:

- The right panel is the current generated scaffold, not a final digitized CAD recovery.
- Check the case manifest for Table 2.2 area/volume agreement; some geometry modes are silhouette-first rather than volume-matched.
- The inner VV/divertor exclusion remains the main shape-fidelity target before F1.
"""
    )
    print(f"wrote {out}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

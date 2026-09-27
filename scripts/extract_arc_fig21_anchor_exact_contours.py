#!/usr/bin/env python3
"""Extract mesh-ready contours from the anchor-exact Ferrero Fig. 2.1 JSON.

The source JSON stores the tank outline and separate VV/neck/foot polygons.  This
script raster-unions the VV pieces and extracts one exterior contour without
requiring Shapely.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
import numpy as np


def remove_consecutive_duplicates(points: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    keep = [0]
    for i in range(1, len(points)):
        if np.linalg.norm(points[i] - points[keep[-1]]) > tol:
            keep.append(i)
    out = points[keep]
    if len(out) > 1 and np.linalg.norm(out[0] - out[-1]) <= tol:
        out = out[:-1]
    return out


def signed_area(points: np.ndarray) -> float:
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def load_poly(data: dict, key: str) -> np.ndarray:
    return remove_consecutive_duplicates(np.asarray(data[key], dtype=float))


def union_mask(polys: list[np.ndarray], rr: np.ndarray, zz: np.ndarray) -> np.ndarray:
    pts = np.c_[rr.ravel(), zz.ravel()]
    mask = np.zeros(len(pts), dtype=bool)
    for poly in polys:
        mask |= MplPath(poly).contains_points(pts)
    return mask.reshape(rr.shape)


def largest_contour(mask: np.ndarray, r_grid: np.ndarray, z_grid: np.ndarray) -> np.ndarray:
    fig, ax = plt.subplots()
    cs = ax.contour(r_grid, z_grid, mask.astype(float), levels=[0.5])
    plt.close(fig)
    candidates: list[np.ndarray] = []
    if hasattr(cs, "collections"):
        for collection in cs.collections:
            for path in collection.get_paths():
                v = path.vertices
                if len(v) > 20:
                    candidates.append(remove_consecutive_duplicates(v))
    else:
        for level_segments in cs.allsegs:
            for v in level_segments:
                if len(v) > 20:
                    candidates.append(remove_consecutive_duplicates(np.asarray(v)))
    if not candidates:
        raise RuntimeError("no VV union contour extracted")
    return max(candidates, key=lambda p: abs(signed_area(p)))


def downsample_closed(points: np.ndarray, target: int) -> np.ndarray:
    if len(points) <= target:
        return points
    closed = np.vstack([points, points[0]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(seg)]
    new_s = np.linspace(0.0, s[-1], target + 1)[:-1]
    out = np.empty((target, 2))
    out[:, 0] = np.interp(new_s, s, closed[:, 0])
    out[:, 1] = np.interp(new_s, s, closed[:, 1])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", default="references/digitization/ferrero_arc/arc_fig21_outline_coords.json")
    ap.add_argument("--outdir", default="references/digitization/ferrero_arc")
    ap.add_argument("--resolution-mm", type=float, default=8.0)
    ap.add_argument("--target-inner-points", type=int, default=260)
    args = ap.parse_args()

    coords = Path(args.coords)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = json.loads(coords.read_text())

    tank = load_poly(data, "tank_inner_wall")
    vv_parts = [
        load_poly(data, "vv_outer"),
        load_poly(data, "upper_neck"),
        load_poly(data, "lower_neck"),
        load_poly(data, "upper_foot"),
        load_poly(data, "lower_foot"),
    ]

    all_pts = np.vstack([tank, *vv_parts])
    margin = 80.0
    r_min, z_min = all_pts.min(axis=0) - margin
    r_max, z_max = all_pts.max(axis=0) + margin
    r = np.arange(r_min, r_max + args.resolution_mm, args.resolution_mm)
    z = np.arange(z_min, z_max + args.resolution_mm, args.resolution_mm)
    rr, zz = np.meshgrid(r, z)
    mask = union_mask(vv_parts, rr, zz)
    inner = largest_contour(mask, rr, zz)
    inner = downsample_closed(inner, args.target_inner_points)

    if signed_area(tank) < 0:
        tank = tank[::-1]
    if signed_area(inner) > 0:
        inner = inner[::-1]

    out_json = outdir / "arc_fig21_anchor_exact_mesh_contours.json"
    payload = {
        "units": "mm",
        "source": str(coords),
        "resolution_mm": args.resolution_mm,
        "tank_inner_wall": tank.tolist(),
        "vv_exclusion_union_contour": inner.tolist(),
        "notes": [
            "VV union contour extracted by rasterizing vv_outer, necks, and feet from the anchor-exact JSON.",
            "Use for Gmsh/OpenFOAM geometry preflight; preserve original JSON for exact anchors.",
        ],
    }
    out_json.write_text(json.dumps(payload, indent=2))

    fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
    ax.plot(tank[:, 0], tank[:, 1], color="#1a355e", lw=2.0, label="tank inner wall")
    ax.plot(inner[:, 0], inner[:, 1], color="#8b1e1e", lw=1.8, label="raster-union VV exclusion contour")
    ax.fill(inner[:, 0], inner[:, 1], color="#f4cccc", alpha=0.35)
    ax.axvline(3300, color="0.55", ls="-.", lw=0.8)
    ax.axhline(0, color="0.55", ls="-.", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("R [mm]")
    ax.set_ylabel("Z [mm]")
    ax.set_title("Ferrero Fig. 2.1 Anchor-Exact Mesh Contours")
    ax.legend(fontsize=8)
    fig_path = outdir / "arc_fig21_anchor_exact_mesh_contours.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    print(f"wrote {out_json}")
    print(f"wrote {fig_path}")
    print(f"tank points={len(tank)} inner points={len(inner)}")


if __name__ == "__main__":
    main()

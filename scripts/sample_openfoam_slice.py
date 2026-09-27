#!/usr/bin/env python3
"""Create lightweight 2D contour slices from OpenFOAM ASCII vol fields.

The script uses the nearest available cell-center slab for a requested plane
and draws a triangulated contour plot. This is meant for report/slide slices,
not high-order interpolation.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
PLANE_AXES = {
    "x": ("y", "z"),
    "y": ("x", "z"),
    "z": ("x", "y"),
}


def _internal_block(text: str) -> tuple[str, int, str]:
    match = re.search(
        r"internalField\s+nonuniform\s+List<(?P<kind>\w+)>\s+"
        r"(?P<count>\d+)\s*\((?P<body>.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError("Could not find nonuniform internalField block")
    return match.group("kind"), int(match.group("count")), match.group("body")


def read_scalar_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "scalar":
        raise ValueError(f"{path} is not a scalar field")
    values = np.fromstring(body, sep=" ")
    if values.size != count:
        raise ValueError(f"{path}: expected {count} scalars, got {values.size}")
    return values


def read_vector_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "vector":
        raise ValueError(f"{path} is not a vector field")
    rows = re.findall(r"\(([^()]+)\)", body)
    data = np.array([[float(x) for x in row.split()] for row in rows], dtype=float)
    if data.shape != (count, 3):
        raise ValueError(f"{path}: expected {(count, 3)}, got {data.shape}")
    return data


def field_values(time_dir: Path, field: str, component: str) -> tuple[np.ndarray, str, str]:
    if field == "Umag":
        vectors = read_vector_field(time_dir / "U")
        return np.linalg.norm(vectors, axis=1), "|U|", "Velocity magnitude [m/s]"

    path = time_dir / field
    text = path.read_text(encoding="ascii")
    kind, _, _ = _internal_block(text)
    if kind == "scalar":
        values = read_scalar_field(path)
        label = "Temperature-like scalar T [K]" if field == "T" else field
        return values, field, label

    vectors = read_vector_field(path)
    if component == "mag":
        return np.linalg.norm(vectors, axis=1), f"|{field}|", f"{field} magnitude"

    idx = AXIS_INDEX[component]
    return vectors[:, idx], f"{field}{component}", f"{field}{component}"


def nearest_slab_indices(
    centers: np.ndarray,
    plane_axis: str,
    plane_value: float,
    half_width: float | None,
) -> tuple[np.ndarray, float, float]:
    axis = AXIS_INDEX[plane_axis]
    coords = centers[:, axis]

    if half_width is None:
        unique_coords = np.unique(np.round(coords, 9))
        nearest = float(unique_coords[np.argmin(np.abs(unique_coords - plane_value))])
        tolerance = max(1.0e-9, 1.0e-6 * max(1.0, abs(nearest)))
    else:
        nearest = plane_value
        tolerance = half_width

    indices = np.where(np.abs(coords - nearest) <= tolerance)[0]
    if indices.size == 0:
        raise ValueError(
            f"No cell centers found near {plane_axis}={plane_value} "
            f"with tolerance {tolerance}"
        )
    return indices, nearest, tolerance


def write_slice_csv(
    path: Path,
    centers: np.ndarray,
    values: np.ndarray,
    indices: np.ndarray,
    field_label: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(["cell", "x", "y", "z", field_label])
        for idx in indices:
            writer.writerow([idx, *centers[idx], values[idx]])


def plot_slice(
    path: Path,
    centers: np.ndarray,
    values: np.ndarray,
    indices: np.ndarray,
    plane_axis: str,
    actual_plane_value: float,
    field_title: str,
    colorbar_label: str,
    title: str,
    levels: int,
    cmap: str,
) -> None:
    axis_a, axis_b = PLANE_AXES[plane_axis]
    a_idx = AXIS_INDEX[axis_a]
    b_idx = AXIS_INDEX[axis_b]
    x = centers[indices, a_idx]
    y = centers[indices, b_idx]
    z = values[indices]

    triangulation = mtri.Triangulation(x, y)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=180)
    contour = ax.tricontourf(triangulation, z, levels=levels, cmap=cmap)
    ax.tricontour(triangulation, z, levels=min(12, levels), colors="k", linewidths=0.25, alpha=0.25)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"{axis_a} [m]")
    ax.set_ylabel(f"{axis_b} [m]")
    ax.set_title(f"{title}\n{field_title}, {plane_axis}={actual_plane_value:.6g} m")
    ax.grid(True, color="white", linewidth=0.35, alpha=0.35)
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label(colorbar_label)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--time", default="2")
    parser.add_argument("--field", default="T", help="Scalar/vector field name, or Umag for |U| from U")
    parser.add_argument("--component", choices=["mag", "x", "y", "z"], default="mag")
    parser.add_argument("--plane-axis", choices=["x", "y", "z"], default="z")
    parser.add_argument("--plane-value", type=float, default=0.02)
    parser.add_argument("--half-width", type=float, default=None, help="Use a finite slab instead of nearest exact center layer")
    parser.add_argument("--name", default=None)
    parser.add_argument("--title", default=None)
    parser.add_argument("--levels", type=int, default=32)
    parser.add_argument("--cmap", default=None)
    parser.add_argument("--outdir", type=Path, default=Path("results/slices"))
    args = parser.parse_args()

    time_dir = args.case / args.time
    centers = read_vector_field(time_dir / "C")
    values, field_title, colorbar_label = field_values(time_dir, args.field, args.component)
    indices, actual_plane_value, tolerance = nearest_slab_indices(
        centers,
        args.plane_axis,
        args.plane_value,
        args.half_width,
    )

    if args.name is None:
        value_tag = str(args.plane_value).replace("-", "m").replace(".", "p")
        name = f"{args.case.name}_{args.field}_{args.plane_axis}{value_tag}"
    else:
        name = args.name

    if args.cmap is None:
        cmap = "inferno" if args.field == "T" else "viridis"
    else:
        cmap = args.cmap

    title = args.title or args.case.name
    csv_path = args.outdir / f"{name}.csv"
    png_path = args.outdir / f"{name}.png"
    write_slice_csv(csv_path, centers, values, indices, field_title)
    plot_slice(
        png_path,
        centers,
        values,
        indices,
        args.plane_axis,
        actual_plane_value,
        field_title,
        colorbar_label,
        title,
        args.levels,
        cmap,
    )

    print(csv_path)
    print(png_path)
    print(f"sampled_cells={indices.size}")
    print(f"requested_{args.plane_axis}={args.plane_value}")
    print(f"actual_{args.plane_axis}={actual_plane_value}")
    print(f"tolerance={tolerance}")


if __name__ == "__main__":
    main()

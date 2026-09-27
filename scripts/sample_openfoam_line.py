#!/usr/bin/env python3
"""Nearest-cell line sampling for simple OpenFOAM ASCII fields.

This is intentionally lightweight. It reads OpenFOAM vol fields written in
ASCII format, samples the nearest cell centers along user-defined lines, and
writes CSV/PNG outputs for early report plots.
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
import numpy as np


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


def line_points(start: np.ndarray, end: np.ndarray, n: int) -> np.ndarray:
    weights = np.linspace(0.0, 1.0, n)
    return start[None, :] * (1.0 - weights[:, None]) + end[None, :] * weights[:, None]


def nearest_indices(centers: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    indices = []
    distances = []
    for p in points:
        d2 = np.sum((centers - p) ** 2, axis=1)
        idx = int(np.argmin(d2))
        indices.append(idx)
        distances.append(math.sqrt(float(d2[idx])))
    return np.array(indices, dtype=int), np.array(distances)


def write_csv(
    path: Path,
    points: np.ndarray,
    centers: np.ndarray,
    indices: np.ndarray,
    distances: np.ndarray,
    temperature: np.ndarray,
    velocity: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "sample",
                "x",
                "y",
                "z",
                "nearest_x",
                "nearest_y",
                "nearest_z",
                "nearest_distance",
                "T",
                "Ux",
                "Uy",
                "Uz",
                "Umag",
            ]
        )
        for i, idx in enumerate(indices):
            u = velocity[idx]
            writer.writerow(
                [
                    i,
                    *points[i],
                    *centers[idx],
                    distances[i],
                    temperature[idx],
                    *u,
                    float(np.linalg.norm(u)),
                ]
            )


def plot_profile(path: Path, axis_values: np.ndarray, axis_label: str, t_values: np.ndarray, u_values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax1 = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    ax1.plot(axis_values, t_values, color="#b3261e", linewidth=2.0, label="T")
    ax1.set_xlabel(axis_label)
    ax1.set_ylabel("Temperature-like scalar T [K]", color="#b3261e")
    ax1.tick_params(axis="y", labelcolor="#b3261e")
    ax1.grid(True, alpha=0.25)

    ax2 = ax1.twinx()
    ax2.plot(axis_values, u_values, color="#2459a6", linewidth=1.8, linestyle="--", label="|U|")
    ax2.set_ylabel("Velocity magnitude [m/s]", color="#2459a6")
    ax2.tick_params(axis="y", labelcolor="#2459a6")

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--time", default="2")
    parser.add_argument("--outdir", type=Path, default=Path("results/line_profiles"))
    args = parser.parse_args()

    time_dir = args.case / args.time
    centers = read_vector_field(time_dir / "C")
    temperature = read_scalar_field(time_dir / "T")
    velocity = read_vector_field(time_dir / "U")

    lines = {
        "vertical_chamber_midplane": {
            "start": np.array([0.20, -0.038, 0.020]),
            "end": np.array([0.20, 0.038, 0.020]),
            "axis": "y",
            "label": "y at x=0.20 m, z=0.02 m [m]",
        },
        "streamwise_near_heated_wall": {
            "start": np.array([0.105, -0.035, 0.020]),
            "end": np.array([0.295, -0.035, 0.020]),
            "axis": "x",
            "label": "x at y=-0.035 m, z=0.02 m [m]",
        },
    }

    for name, spec in lines.items():
        points = line_points(spec["start"], spec["end"], 80)
        indices, distances = nearest_indices(centers, points)
        sampled_t = temperature[indices]
        sampled_u = np.linalg.norm(velocity[indices], axis=1)
        axis = points[:, 1] if spec["axis"] == "y" else points[:, 0]

        csv_path = args.outdir / f"{args.case.name}_{name}.csv"
        png_path = args.outdir / f"{args.case.name}_{name}.png"
        write_csv(csv_path, points, centers, indices, distances, temperature, velocity)
        plot_profile(png_path, axis, spec["label"], sampled_t, sampled_u)
        print(csv_path)
        print(png_path)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Extract x-wise wall-temperature profiles from OpenFOAM VTP patch exports."""

from __future__ import annotations

import argparse
import base64
import csv
import os
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


VTK_DTYPES = {
    "Float32": np.dtype("<f4"),
    "Float64": np.dtype("<f8"),
    "Int32": np.dtype("<i4"),
    "UInt64": np.dtype("<u8"),
}


def decode_data_array(element: ET.Element) -> np.ndarray:
    vtk_type = element.attrib["type"]
    dtype = VTK_DTYPES[vtk_type]
    encoded = "".join((element.text or "").split())
    blob = base64.b64decode(encoded)
    if len(blob) < 8:
        raise ValueError(f"DataArray {element.attrib.get('Name')} is missing UInt64 payload header")

    payload_length = struct.unpack("<Q", blob[:8])[0]
    payload = blob[8 : 8 + payload_length]
    values = np.frombuffer(payload, dtype=dtype).copy()

    components = int(element.attrib.get("NumberOfComponents", "1"))
    if components > 1:
        values = values.reshape((-1, components))
    return values


def find_data_array(parent: ET.Element, name: str) -> ET.Element:
    for element in parent.iter("DataArray"):
        if element.attrib.get("Name") == name:
            return element
    raise ValueError(f"Could not find DataArray named {name}")


def read_vtp_wall_cells(path: Path) -> tuple[np.ndarray, np.ndarray]:
    root = ET.parse(path).getroot()
    piece = root.find("./PolyData/Piece")
    if piece is None:
        raise ValueError(f"{path} does not contain a PolyData Piece")

    points_parent = piece.find("Points")
    polys_parent = piece.find("Polys")
    cell_data_parent = piece.find("CellData")
    if points_parent is None or polys_parent is None or cell_data_parent is None:
        raise ValueError(f"{path} is missing Points, Polys, or CellData")

    points = decode_data_array(find_data_array(points_parent, "Points"))
    connectivity = decode_data_array(find_data_array(polys_parent, "connectivity")).astype(int)
    offsets = decode_data_array(find_data_array(polys_parent, "offsets")).astype(int)
    temperature = decode_data_array(find_data_array(cell_data_parent, "T")).astype(float)

    centroids = []
    start = 0
    for stop in offsets:
        point_ids = connectivity[start:stop]
        centroids.append(points[point_ids].mean(axis=0))
        start = int(stop)
    centers = np.array(centroids)

    if len(centers) != len(temperature):
        raise ValueError(f"{path}: {len(centers)} cells but {len(temperature)} T values")
    return centers, temperature


def binned_profile(centers: np.ndarray, temperature: np.ndarray, bins: int) -> list[dict[str, float]]:
    x = centers[:, 0]
    edges = np.linspace(float(x.min()), float(x.max()), bins + 1)
    indices = np.digitize(x, edges, right=False) - 1
    indices = np.clip(indices, 0, bins - 1)

    rows: list[dict[str, float]] = []
    for i in range(bins):
        mask = indices == i
        if not np.any(mask):
            continue
        t = temperature[mask]
        rows.append(
            {
                "bin": float(i),
                "x": float((edges[i] + edges[i + 1]) * 0.5),
                "count": float(mask.sum()),
                "T_min": float(t.min()),
                "T_mean": float(t.mean()),
                "T_max": float(t.max()),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]], inlet_temperature: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        fieldnames = ["bin", "x", "count", "T_min", "T_mean", "T_max", "dT_min", "dT_mean", "dT_max"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["dT_min"] = row["T_min"] - inlet_temperature
            out["dT_mean"] = row["T_mean"] - inlet_temperature
            out["dT_max"] = row["T_max"] - inlet_temperature
            writer.writerow(out)


def plot_profile(path: Path, rows: list[dict[str, float]], inlet_temperature: float, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    x = [row["x"] for row in rows]
    mean_dt = [row["T_mean"] - inlet_temperature for row in rows]
    max_dt = [row["T_max"] - inlet_temperature for row in rows]

    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    ax.plot(x, mean_dt, color="#0f766e", linewidth=2.0, label="wall mean")
    ax.plot(x, max_dt, color="#b3261e", linewidth=1.8, linestyle="--", label="wall max")
    ax.set_title(title)
    ax.set_xlabel("x along heated wall [m]")
    ax.set_ylabel("Wall temperature rise T - Tin [K]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtp", type=Path)
    parser.add_argument("--name", required=True)
    parser.add_argument("--tin", type=float, default=800.0)
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--outdir", type=Path, default=Path("results/wall_profiles"))
    args = parser.parse_args()

    centers, temperature = read_vtp_wall_cells(args.vtp)
    rows = binned_profile(centers, temperature, args.bins)

    csv_path = args.outdir / f"{args.name}_heated_wall_surface_profile.csv"
    png_path = args.outdir / f"{args.name}_heated_wall_surface_profile.png"
    write_csv(csv_path, rows, args.tin)
    plot_profile(png_path, rows, args.tin, args.name)
    print(csv_path)
    print(png_path)


if __name__ == "__main__":
    main()

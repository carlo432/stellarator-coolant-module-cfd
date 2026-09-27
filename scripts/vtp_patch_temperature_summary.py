#!/usr/bin/env python3
"""Summarize temperature values on OpenFOAM VTP patch exports."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from sample_wall_temperature_profile import read_vtp_wall_cells


def summarize(path: Path, inlet_temperature: float) -> dict[str, str | float]:
    centers, temperature = read_vtp_wall_cells(path)
    return {
        "patch_file": str(path),
        "faces": float(len(temperature)),
        "x_min": float(centers[:, 0].min()),
        "x_max": float(centers[:, 0].max()),
        "T_min": float(temperature.min()),
        "T_mean": float(temperature.mean()),
        "T_max": float(temperature.max()),
        "dT_min": float(temperature.min() - inlet_temperature),
        "dT_mean": float(temperature.mean() - inlet_temperature),
        "dT_max": float(temperature.max() - inlet_temperature),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtp", nargs="+", type=Path)
    parser.add_argument("--tin", type=float, default=800.0)
    parser.add_argument("--out", type=Path, default=Path("results/visualization/patch_temperature_summary.csv"))
    args = parser.parse_args()

    rows = [summarize(path, args.tin) for path in args.vtp]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="ascii") as handle:
        fieldnames = ["patch_file", "faces", "x_min", "x_max", "T_min", "T_mean", "T_max", "dT_min", "dT_mean", "dT_max"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(args.out)
    for row in rows:
        print(
            f"{Path(str(row['patch_file'])).name}: "
            f"T_mean={row['T_mean']:.6g}, T_min={row['T_min']:.6g}, T_max={row['T_max']:.6g}"
        )


if __name__ == "__main__":
    main()

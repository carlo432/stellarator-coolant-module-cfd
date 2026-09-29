#!/usr/bin/env python3
"""B6 audit for promoted open-flow stellarator-pipeline scalar fields.

Unlike the closed periodic Kawamura channel, these ducts equilibrate by
advection. The applicable checks are temporal flatness, inlet/outlet mass
closure, and source-to-outflow balances for temperature and tritium.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import numpy as np


DEFAULT_CASES = [
    "geom_pipeline_plasma",
    "geom_nowave_pipeline_plasma",
    "geom_straight_pipeline_plasma",
    "geom_throat_pipeline_plasma",
    "geom_phase_000",
    "geom_phase_075",
    "geom_phase_146",
    "geom_phase_222",
]
RHO_CP = 4.65e6
T_IN = 900.0


def data_files(case: Path, object_name: str) -> list[Path]:
    return sorted((case / "postProcessing" / object_name).glob("*/*.dat"))


def series(case: Path, object_name: str) -> np.ndarray:
    rows: dict[float, list[float]] = {}
    for path in data_files(case, object_name):
        for line in path.read_text(encoding="ascii").splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = [float(value) for value in line.split()]
            rows[values[0]] = values[1:]
    if not rows:
        raise FileNotFoundError(f"no {object_name} rows under {case}")
    return np.asarray([[time, *values] for time, values in sorted(rows.items())])


def header_value(path: Path, label: str) -> float:
    match = re.search(rf"^# {re.escape(label)}\s+:\s+(\S+)", path.read_text(), re.M)
    if not match:
        raise ValueError(f"missing {label} in {path}")
    return float(match.group(1))


def source_integral(case: Path, field: str) -> float:
    control = (case / "system" / "controlDict").read_text(encoding="ascii")
    pattern = rf"qzone(\d+)_{field}.*?sources\s*\{{\s*{field}\s*\((\S+)"
    sources = {int(zone): float(value) for zone, value in re.findall(pattern, control)}
    total = 0.0
    for zone, source in sources.items():
        files = data_files(case, f"qzone{zone}_stats")
        if not files:
            raise FileNotFoundError(f"missing qzone{zone}_stats for {case}")
        total += source * header_value(files[0], "Volume")
    return total


def tail(values: np.ndarray, window_s: float) -> np.ndarray:
    return values[values[:, 0] >= values[-1, 0] - window_s - 1e-12]


def half_shift(values: np.ndarray) -> float:
    split = len(values) // 2
    return float(values[split:].mean() - values[:split].mean())


def audit_case(case: Path, window_s: float) -> dict[str, float | str | bool]:
    wall = tail(series(case, "fwT"), window_s)
    inlet_phi = tail(series(case, "inletPhi"), window_s)
    outlet_phi = tail(series(case, "outletPhi"), window_s)
    outlet_scalar = tail(series(case, "outletScal"), window_s)

    wall_t = wall[:, -1]
    wall_superheat = float(wall_t.mean() - T_IN)
    duration = float(wall[-1, 0] - wall[0, 0])
    wall_trend = float(np.polyfit(wall[:, 0], wall_t, 1)[0] * duration)
    wall_half = half_shift(wall_t)

    volume_flow_in = float(abs(inlet_phi[:, -1].mean()))
    volume_flow_out = float(outlet_phi[:, -1].mean())
    mass_error = abs(volume_flow_out - volume_flow_in) / volume_flow_in

    wall_area = header_value(data_files(case, "fwT")[0], "Area")
    # The phase sweep reports 0.5025 MW/m2; completed legacy-map cases carry
    # the measured 0.5026 MW/m2 endpoint-audit value.
    q_wall = 502500.0 if case.name.startswith("geom_phase_") else 502600.0
    expected_t_rise = (
        q_wall * wall_area / RHO_CP + source_integral(case, "T")
    ) / volume_flow_in
    measured_t_rise = float(outlet_scalar[:, 1].mean() - T_IN)
    energy_error = (measured_t_rise - expected_t_rise) / expected_t_rise

    expected_tritium = source_integral(case, "Ctrit") / volume_flow_in
    measured_tritium = float(outlet_scalar[:, 2].mean())
    tritium_error = (measured_tritium - expected_tritium) / expected_tritium
    tritium_half = half_shift(outlet_scalar[:, 2]) / measured_tritium

    zone_t_half_shifts = []
    zone_c_half_shifts = []
    for zone in range(6):
        zone_values = tail(series(case, f"qzone{zone}_stats"), window_s)
        zone_t_half_shifts.append(abs(half_shift(zone_values[:, 1])))
        zone_c_half_shifts.append(
            abs(half_shift(zone_values[:, 2])) / abs(zone_values[:, 2].mean())
        )
    max_zone_t_half_shift = float(max(zone_t_half_shifts))
    max_zone_c_half_shift = float(max(zone_c_half_shifts))

    temporal_pass = bool(
        abs(wall_trend) / wall_superheat <= 0.01
        and abs(wall_half) / wall_superheat <= 0.01
        and max_zone_t_half_shift <= 0.5
        and max_zone_c_half_shift <= 0.02
    )
    conservation_pass = bool(
        mass_error <= 0.001
        and abs(energy_error) <= 0.05
        and abs(tritium_error) <= 0.01
    )
    borderline = bool(abs(energy_error) <= 0.055 and not conservation_pass)
    return {
        "case": case.name,
        "window_start_s": float(wall[0, 0]),
        "window_end_s": float(wall[-1, 0]),
        "samples": len(wall),
        "wall_mean_K": float(wall_t.mean()),
        "wall_half_shift_percent_superheat": 100 * wall_half / wall_superheat,
        "wall_trend_percent_superheat": 100 * wall_trend / wall_superheat,
        "mass_balance_error_percent": 100 * mass_error,
        "expected_outlet_rise_K": expected_t_rise,
        "measured_outlet_rise_K": measured_t_rise,
        "thermal_energy_error_percent": 100 * energy_error,
        "tritium_balance_error_percent": 100 * tritium_error,
        "tritium_half_shift_percent": 100 * tritium_half,
        "max_zone_temperature_half_shift_K": max_zone_t_half_shift,
        "max_zone_tritium_half_shift_percent": 100 * max_zone_c_half_shift,
        "temporal_pass": temporal_pass,
        "conservation_pass": conservation_pass,
        "borderline": borderline,
        "status": "PASS" if temporal_pass and conservation_pass else (
            "BORDERLINE" if temporal_pass and borderline else "FAIL"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", default=DEFAULT_CASES)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--window", type=float, default=0.7)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rows = [audit_case(args.root / name, args.window) for name in args.cases]

    for row in rows:
        print(
            f"[{row['status']:10s}] {row['case']:32s} "
            f"wall trend={row['wall_trend_percent_superheat']:+.3f}%  "
            f"energy={row['thermal_energy_error_percent']:+.2f}%  "
            f"tritium={row['tritium_balance_error_percent']:+.2f}%"
        )

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="ascii") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps({"window_s": args.window, "cases": rows}, indent=2) + "\n")
    return 1 if any(row["status"] == "FAIL" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Audit C9 mean-profile stability across walls and late-window halves."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from independent_c9_gate_review import (  # noqa: E402
    accumulation_time,
    available_times,
    load_global,
    read_scalar,
)


def grouped_profile(coordinate: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rounded = np.round(coordinate, 9)
    levels = np.unique(rounded)
    profile = np.asarray([values[rounded == level].mean() for level in levels])
    return levels, profile


def window_mean(
    case: Path,
    field: str,
    base_field: str,
    start: tuple[float, str],
    end: tuple[float, str],
    size: int,
) -> tuple[np.ndarray, list[float]]:
    start_total = accumulation_time(case, start[1], base_field)
    end_total = accumulation_time(case, end[1], base_field)
    if end_total <= start_total:
        raise ValueError(f"non-increasing accumulator for {field}: {start_total} -> {end_total}")
    start_values = load_global(case, start[1], field, size)
    end_values = load_global(case, end[1], field, size)
    values = (
        end_values * end_total - start_values * start_total
    ) / (end_total - start_total)
    return values, [start_total, end_total]


def pressure_samples(log: Path) -> dict[float, float]:
    current: float | None = None
    values: dict[float, float] = {}
    for line in log.read_text(encoding="ascii", errors="replace").splitlines():
        time_match = re.match(r"\s*Time\s*=\s*([-+0-9.eE]+)\s*$", line)
        if time_match:
            current = float(time_match.group(1))
            continue
        gradient = re.search(r"pressure gradient =\s*([-+0-9.eE]+)", line)
        if gradient and current is not None:
            values[current] = float(gradient.group(1))
    return values


def load_dns(path: Path) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            xs.append(float(row["y_plus"]))
            ys.append(float(row["Theta_plus"]))
    return xs, ys


def interp(xs: list[float] | np.ndarray, ys: list[float] | np.ndarray, xq: float) -> float:
    return float(np.interp(xq, np.asarray(xs), np.asarray(ys)))


def wall_metrics(
    y: np.ndarray,
    profile: np.ndarray,
    side: str,
    dpdx: float,
    dns_y: list[float],
    dns_t: list[float],
    nu: float,
    pr: float,
    delta: float,
) -> dict[str, Any]:
    utau = math.sqrt(abs(dpdx) * delta)
    alpha = nu / pr
    if side == "lower":
        mask = y <= delta
        distance = y[mask]
        scalar = profile[mask]
        gradient = float(np.polyfit(np.r_[0.0, distance[:3]], np.r_[0.0, scalar[:3]], 1)[0])
    else:
        mask = y >= delta
        distance = 2.0 * delta - y[mask]
        scalar = profile[mask]
        order = np.argsort(distance)
        distance = distance[order]
        scalar = scalar[order]
        gradient = float(np.polyfit(np.r_[0.0, distance[:3]], np.r_[0.0, scalar[:3]], 1)[0])
    order = np.argsort(distance)
    distance = distance[order]
    scalar = scalar[order]
    yplus = distance * utau / nu
    ttau = alpha * gradient / utau
    tplus = scalar / ttau
    points = [30.0, 50.0, 75.0, 100.0, 150.0]
    errors = {
        str(point): 100.0
        * abs(interp(yplus, tplus, point) - interp(dns_y, dns_t, point))
        / interp(dns_y, dns_t, point)
        for point in points
    }
    center_dns = dns_t[-1]
    center = float(tplus[-1])
    return {
        "dpdx": dpdx,
        "u_tau": utau,
        "Re_tau": utau * delta / nu,
        "wall_gradient": gradient,
        "T_tau": ttau,
        "centerline_T_plus": center,
        "centerline_error_percent": 100.0 * (center - center_dns) / center_dns,
        "log_layer_errors_percent": errors,
        "log_layer_max_error_percent": max(errors.values()),
        "sublayer_slope_ratio_at_1.79": interp(yplus, tplus, 1.79) / (pr * 1.79),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--solver-log", required=True, type=Path)
    parser.add_argument("--dns", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--nu", type=float, default=3.09278e-6)
    parser.add_argument("--pr", type=float, default=5.0)
    parser.add_argument("--delta", type=float, default=0.01)
    args = parser.parse_args()

    case = args.case.resolve()
    contract = json.loads((case / "constant/sourceDrivenScalarGate.json").read_text())
    start = float(contract["balance_window_start"])
    times = [row for row in available_times(case, "HMean") if row[0] + 1e-12 >= start]
    if len(times) < 5:
        raise ValueError(f"need at least five window writes, found {len(times)}")
    opening = times[0]
    endpoint = times[-1]
    midpoint_target = 0.5 * (opening[0] + endpoint[0])
    midpoint = min(times[1:-1], key=lambda row: abs(row[0] - midpoint_target))

    endpoint_values = read_scalar(case / endpoint[1] / "HMean")
    size = len(endpoint_values)
    coordinate = None
    coordinate_path = None
    for _, name in reversed(available_times(case, "Cy")):
        candidate = case / name / "Cy"
        values = read_scalar(candidate)
        if len(values) == size:
            coordinate = values
            coordinate_path = candidate
            break
    if coordinate is None:
        raise FileNotFoundError("no matching root Cy field")

    dns_y, dns_t = load_dns(args.dns.resolve())
    gradients = pressure_samples(args.solver_log.resolve())
    windows = {
        "full": (opening, endpoint),
        "early": (opening, midpoint),
        "late": (midpoint, endpoint),
    }
    result: dict[str, Any] = {
        "case": str(case),
        "solver_log": str(args.solver_log.resolve()),
        "coordinate": str(coordinate_path),
        "dns": str(args.dns.resolve()),
        "windows": {},
    }
    for label, (first, last) in windows.items():
        values, totals = window_mean(case, "HMean", "H", first, last, size)
        y, profile = grouped_profile(coordinate, values)
        selected_dpdx = [
            value for time, value in gradients.items() if first[0] <= time <= last[0]
        ]
        if not selected_dpdx:
            raise ValueError(f"no pressure-gradient samples in {first[0]}..{last[0]}")
        dpdx = float(np.mean(selected_dpdx))
        result["windows"][label] = {
            "write_window": [first[0], last[0]],
            "accumulation_times": totals,
            "pressure_sample_count": len(selected_dpdx),
            "pressure_gradient_mean": dpdx,
            "lower": wall_metrics(
                y, profile, "lower", dpdx, dns_y, dns_t, args.nu, args.pr, args.delta
            ),
            "upper": wall_metrics(
                y, profile, "upper", dpdx, dns_y, dns_t, args.nu, args.pr, args.delta
            ),
        }

    full = result["windows"]["full"]
    early = result["windows"]["early"]
    late = result["windows"]["late"]
    result["stability"] = {
        "full_wall_centerline_spread_percent": abs(
            full["upper"]["centerline_T_plus"] - full["lower"]["centerline_T_plus"]
        )
        / (0.5 * (full["upper"]["centerline_T_plus"] + full["lower"]["centerline_T_plus"]))
        * 100.0,
        "lower_early_late_centerline_change_percent": 100.0
        * (late["lower"]["centerline_T_plus"] - early["lower"]["centerline_T_plus"])
        / early["lower"]["centerline_T_plus"],
        "upper_early_late_centerline_change_percent": 100.0
        * (late["upper"]["centerline_T_plus"] - early["upper"]["centerline_T_plus"])
        / early["upper"]["centerline_T_plus"],
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps(result["stability"], indent=2))
    for label, window in result["windows"].items():
        print(
            f"{label}: lower center err={window['lower']['centerline_error_percent']:+.3f}% "
            f"P1={window['lower']['log_layer_max_error_percent']:.3f}%; "
            f"upper center err={window['upper']['centerline_error_percent']:+.3f}% "
            f"P1={window['upper']['log_layer_max_error_percent']:.3f}%"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Reference-free convergence gate for source-driven scalar statistics.

The gate is configured by ``constant/sourceDrivenScalarGate.json``. It checks
the exact source-to-sink balance and verifies that the scalar bulk is flat over
the averaging window. A failed gate returns exit code 2 so callers cannot
silently continue to profile or DNS comparisons.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

import numpy as np


FIELD_RE = {
    "scalar": re.compile(
        r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;",
        re.S,
    ),
    "vector": re.compile(
        r"internalField\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;",
        re.S,
    ),
}


@dataclass
class Snapshot:
    values: np.ndarray
    distributed: bool
    part_lengths: list[int]


def numeric_time_dirs(base: Path) -> list[tuple[float, str]]:
    rows: list[tuple[float, str]] = []
    if not base.exists():
        return rows
    for path in base.iterdir():
        if not path.is_dir():
            continue
        try:
            rows.append((float(path.name), path.name))
        except ValueError:
            continue
    return sorted(rows)


def processor_dirs(case: Path) -> list[Path]:
    def index(path: Path) -> int:
        match = re.fullmatch(r"processor(\d+)", path.name)
        return int(match.group(1)) if match else 10**9

    return sorted(case.glob("processor[0-9]*"), key=index)


def read_field(path: Path, kind: str = "scalar", component: int = 0) -> np.ndarray:
    text = path.read_text(encoding="ascii", errors="replace")
    match = FIELD_RE[kind].search(text)
    if not match:
        uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
        if uniform:
            raise ValueError(f"uniform field unsupported for convergence statistics: {path}")
        raise ValueError(f"could not parse {kind} internalField: {path}")
    if kind == "scalar":
        return np.fromstring(match.group(1), sep=" ")
    vectors = re.findall(r"\(([^()]*)\)", match.group(1))
    return np.asarray([float(row.split()[component]) for row in vectors])


def available_times(case: Path, field: str) -> list[tuple[float, str]]:
    choices: dict[float, str] = {}
    procs = processor_dirs(case)
    if procs:
        for value, name in numeric_time_dirs(procs[0]):
            if (procs[0] / name / field).exists():
                choices[value] = name
    for value, name in numeric_time_dirs(case):
        if (case / name / field).exists():
            choices[value] = name
    return sorted(choices.items())


def load_snapshot(
    case: Path,
    time_name: str,
    field: str,
    kind: str = "scalar",
    component: int = 0,
) -> Snapshot:
    root_path = case / time_name / field
    if root_path.exists():
        values = read_field(root_path, kind=kind, component=component)
        return Snapshot(values=values, distributed=False, part_lengths=[len(values)])

    arrays = []
    for proc in processor_dirs(case):
        path = proc / time_name / field
        if not path.exists():
            raise FileNotFoundError(path)
        arrays.append(read_field(path, kind=kind, component=component))
    if not arrays:
        raise FileNotFoundError(root_path)
    return Snapshot(
        values=np.concatenate(arrays),
        distributed=True,
        part_lengths=[len(array) for array in arrays],
    )


def find_matching_field(
    base: Path,
    preferred_time: str,
    field: str,
    expected_length: int,
) -> np.ndarray:
    candidates = [(float("inf"), preferred_time)]
    candidates.extend(reversed(numeric_time_dirs(base)))
    seen = set()
    for _, name in candidates:
        if name in seen:
            continue
        seen.add(name)
        path = base / name / field
        if not path.exists():
            continue
        values = read_field(path)
        if len(values) == expected_length:
            return values
    raise FileNotFoundError(
        f"no {field} field with {expected_length} values under {base}"
    )


def load_coordinate(
    case: Path,
    time_name: str,
    coordinate_field: str,
    snapshot: Snapshot,
) -> np.ndarray:
    if not snapshot.distributed:
        return find_matching_field(
            case,
            time_name,
            coordinate_field,
            snapshot.part_lengths[0],
        )
    arrays = []
    for proc, length in zip(processor_dirs(case), snapshot.part_lengths):
        arrays.append(
            find_matching_field(proc, time_name, coordinate_field, length)
        )
    return np.concatenate(arrays)


def plane_profile(
    coordinate: np.ndarray,
    values: np.ndarray,
    decimals: int,
) -> tuple[np.ndarray, np.ndarray]:
    rounded = np.round(coordinate, decimals)
    levels = np.unique(rounded)
    profile = np.asarray([values[rounded == level].mean() for level in levels])
    return levels, profile


def wall_flux(
    coordinate: np.ndarray,
    values: np.ndarray,
    scalar: dict[str, Any],
    contract: dict[str, Any],
) -> dict[str, float]:
    walls = contract["walls_m"]
    wall_values = contract.get("wall_values", [0.0, 0.0])
    fit_levels = int(contract.get("near_wall_fit_levels", 3))
    decimals = int(contract.get("coordinate_round_decimals", 10))
    levels, profile = plane_profile(coordinate, values, decimals)
    if len(levels) < 2 * fit_levels:
        raise ValueError("not enough wall-normal levels for two wall-gradient fits")

    lower_x = np.r_[float(walls[0]), levels[:fit_levels]]
    lower_y = np.r_[float(wall_values[0]), profile[:fit_levels]]
    upper_x = np.r_[levels[-fit_levels:], float(walls[1])]
    upper_y = np.r_[profile[-fit_levels:], float(wall_values[1])]
    lower_gradient = float(np.polyfit(lower_x, lower_y, 1)[0])
    upper_gradient = float(np.polyfit(upper_x, upper_y, 1)[0])
    diffusivity = float(scalar["molecular_diffusivity_m2_s"])
    lower_flux = diffusivity * abs(lower_gradient)
    upper_flux = diffusivity * abs(upper_gradient)
    return {
        "lower_gradient": lower_gradient,
        "upper_gradient": upper_gradient,
        "lower_flux": lower_flux,
        "upper_flux": upper_flux,
        "mean_wall_flux": 0.5 * (lower_flux + upper_flux),
    }


def source_mean(
    case: Path,
    time_name: str,
    scalar: dict[str, Any],
    snapshot: Snapshot,
) -> float:
    source = scalar["source"]
    kind = source["type"]
    if kind == "constant_mean":
        return float(source["mean_value_per_s"])
    if kind == "velocity_scaled_component":
        velocity = load_snapshot(
            case,
            time_name,
            source.get("field", "U"),
            kind="vector",
            component=int(source.get("component", 0)),
        )
        if velocity.distributed != snapshot.distributed:
            raise ValueError("source velocity/scalar layouts differ")
        coefficient = float(source["coefficient_per_s"])
        reference = float(source["reference_velocity_m_s"])
        # A volumeMode-specific source integrates against cell volume, so the
        # bulk velocity must be volume-weighted; a plain cell mean is biased
        # low on wall-graded meshes.
        volumes = load_coordinate(
            case,
            time_name,
            str(source.get("volume_field", "V")),
            velocity,
        )
        bulk = float(np.average(velocity.values, weights=volumes))
        return coefficient * bulk / reference
    raise ValueError(f"unsupported source type: {kind}")


def accumulated_time(case: Path, time_name: str, source_field: str) -> float:
    """Exact fieldAverage accumulation time for source_field at a write.

    Read from uniform/functionObjects/functionObjectProperties so the window
    mean needs no assumption about when accumulation actually started.
    """
    rel = Path("uniform") / "functionObjects" / "functionObjectProperties"
    candidates = [case / time_name / rel]
    procs = processor_dirs(case)
    if procs:
        candidates.append(procs[0] / time_name / rel)
    pattern = re.compile(
        rf"(?<![\w]){re.escape(source_field)}(?![\w])\s*\{{[^{{}}]*?"
        rf"totalTime\s+([\d.eE+-]+);",
        re.S,
    )
    for path in candidates:
        if not path.exists():
            continue
        match = pattern.search(path.read_text(encoding="ascii", errors="replace"))
        if match:
            return float(match.group(1))
    raise FileNotFoundError(
        f"no accumulated totalTime for {source_field} at t={time_name}"
    )


def flatness_history(
    case: Path,
    field: str,
    time_start: float,
    volume_field: str | None = None,
) -> list[dict[str, float]]:
    rows = []
    for value, name in available_times(case, field):
        if value + 1e-12 < time_start:
            continue
        snapshot = load_snapshot(case, name, field)
        row = {"time": value, "cell_mean": float(snapshot.values.mean())}
        if volume_field is not None:
            # Wall-graded meshes over-weight small near-wall cells in a plain
            # cell mean; the physical scalar inventory is volume-weighted.
            weights = load_coordinate(case, name, volume_field, snapshot)
            row["volume_weighted_mean"] = float(
                np.average(snapshot.values, weights=weights)
            )
        rows.append(row)
    return rows


def evaluate_contract(case: Path, contract_path: Path | None = None) -> dict[str, Any]:
    case = Path(case)
    contract_path = contract_path or case / "constant" / "sourceDrivenScalarGate.json"
    contract = json.loads(contract_path.read_text(encoding="ascii"))
    balance_tolerance = float(contract.get("balance_tolerance_fraction", 0.05))
    flatness_tolerance = float(contract.get("flatness_tolerance_fraction", 0.01))
    flatness_start = float(contract["flatness_time_start"])
    minimum_snapshots = int(contract.get("minimum_flatness_snapshots", 3))
    flatness_mode = str(contract.get("flatness_mean", "arithmetic"))
    if flatness_mode not in ("arithmetic", "volume"):
        raise ValueError(f"unsupported flatness_mean: {flatness_mode}")
    volume_field = str(contract.get("volume_field", "V"))
    results = []

    for scalar in contract["scalars"]:
        field = scalar["field"]
        times = available_times(case, field)
        if not times:
            raise FileNotFoundError(f"no time directories contain {field}")
        latest_time, latest_name = times[-1]
        snapshot = load_snapshot(case, latest_name, field)
        coordinate = load_coordinate(
            case,
            latest_name,
            contract.get("coordinate_field", "Cy"),
            snapshot,
        )
        instantaneous_flux = wall_flux(coordinate, snapshot.values, scalar, contract)
        # The steady balance holds for the time-mean flux; an instantaneous
        # endpoint gradient carries turbulent fluctuations of several percent,
        # so a balance_field (e.g. HMean over the averaging window) governs
        # when the contract provides one.
        balance_field = scalar.get("balance_field")
        if balance_field:
            balance_times = available_times(case, balance_field)
            if not balance_times:
                raise FileNotFoundError(
                    f"no time directories contain balance field {balance_field}"
                )
            balance_time, balance_name = balance_times[-1]
            if abs(balance_time - latest_time) > 1e-9:
                raise ValueError(
                    f"latest {field} time {latest_time:.12g} does not have a "
                    f"matching {balance_field}; newest mean field is "
                    f"{balance_time:.12g}"
                )
            balance_snapshot = load_snapshot(case, balance_name, balance_field)
            balance_coordinate = load_coordinate(
                case,
                balance_name,
                contract.get("coordinate_field", "Cy"),
                balance_snapshot,
            )
            cumulative_flux = wall_flux(
                balance_coordinate, balance_snapshot.values, scalar, contract
            )
            window_start = contract.get("balance_window_start")
            if window_start is not None:
                # A cumulative accumulator drags any early filling transient
                # into the mean forever. The late-window mean over [w1, w2] is
                # exact from two cumulative writes and their recorded
                # accumulation times: (M2*T2 - M1*T1) / (T2 - T1).
                window_start = float(window_start)
                early = [
                    (value, name)
                    for value, name in balance_times
                    if value + 1e-12 >= window_start
                    and value < balance_time - 1e-9
                ]
                if not early:
                    raise ValueError(
                        f"no {balance_field} write in [{window_start:.12g}, "
                        f"{balance_time:.12g}) to open the balance window"
                    )
                window_time, window_name = early[0]
                window_snapshot = load_snapshot(case, window_name, balance_field)
                if window_snapshot.part_lengths != balance_snapshot.part_lengths:
                    raise ValueError(
                        f"{balance_field} layouts differ between window writes"
                    )
                t_early = accumulated_time(case, window_name, field)
                t_late = accumulated_time(case, balance_name, field)
                if not t_late > t_early > 0.0:
                    raise ValueError(
                        f"non-increasing accumulation times for {field}: "
                        f"{t_early:.12g} -> {t_late:.12g}"
                    )
                window_values = (
                    balance_snapshot.values * t_late
                    - window_snapshot.values * t_early
                ) / (t_late - t_early)
                flux = wall_flux(
                    balance_coordinate, window_values, scalar, contract
                )
                balance_mode = "window_mean"
                balance_window = [window_time, balance_time]
            else:
                flux = cumulative_flux
                balance_mode = "cumulative_mean"
                balance_window = None
        else:
            balance_time = latest_time
            flux = instantaneous_flux
            cumulative_flux = instantaneous_flux
            balance_mode = "instantaneous"
            balance_window = None
        mean_source = source_mean(case, latest_name, scalar, snapshot)
        expected_flux = mean_source * float(contract["source_to_sink_area_length_m"])
        balance_error = abs(flux["mean_wall_flux"] - expected_flux) / max(
            abs(expected_flux), 1e-30
        )

        history = flatness_history(
            case,
            field,
            flatness_start,
            volume_field if flatness_mode == "volume" else None,
        )
        mean_key = "volume_weighted_mean" if flatness_mode == "volume" else "cell_mean"
        if len(history) >= minimum_snapshots:
            means = np.asarray([row[mean_key] for row in history])
            flatness_error = float(
                (means.max() - means.min()) / max(abs(means.mean()), 1e-30)
            )
            arithmetic_means = np.asarray([row["cell_mean"] for row in history])
            arithmetic_flatness_error = float(
                (arithmetic_means.max() - arithmetic_means.min())
                / max(abs(arithmetic_means.mean()), 1e-30)
            )
        else:
            flatness_error = float("inf")
            arithmetic_flatness_error = float("inf")
        passed = balance_error <= balance_tolerance and flatness_error <= flatness_tolerance
        results.append(
            {
                "field": field,
                "latest_time": latest_time,
                "balance_field": balance_field or field,
                "balance_time": balance_time,
                "balance_mode": balance_mode,
                "balance_window": balance_window,
                "source_mean_per_s": mean_source,
                "expected_mean_wall_flux": expected_flux,
                "measured_mean_wall_flux": flux["mean_wall_flux"],
                "measured_mean_wall_flux_cumulative": cumulative_flux[
                    "mean_wall_flux"
                ],
                "measured_mean_wall_flux_instantaneous": instantaneous_flux[
                    "mean_wall_flux"
                ],
                "balance_relative_error": balance_error,
                "balance_fraction_of_target": flux["mean_wall_flux"] / expected_flux,
                "flatness_history": history,
                "flatness_mean_type": flatness_mode,
                "flatness_relative_range": flatness_error,
                "flatness_relative_range_arithmetic": arithmetic_flatness_error,
                "minimum_flatness_snapshots": minimum_snapshots,
                "balance_tolerance_fraction": balance_tolerance,
                "flatness_tolerance_fraction": flatness_tolerance,
                "wall_flux_detail": flux,
                "instantaneous_wall_flux_detail": instantaneous_flux,
                "passed": passed,
            }
        )

    return {
        "case": str(case),
        "contract": str(contract_path),
        "status": "PASS" if all(row["passed"] for row in results) else "FAIL",
        "statistics_allowed": all(row["passed"] for row in results),
        "scalars": results,
    }


def print_report(result: dict[str, Any]) -> None:
    print(f"=== source-driven scalar convergence gate: {result['case']} ===")
    for row in result["scalars"]:
        mark = "PASS" if row["passed"] else "FAIL"
        window = row.get("balance_window")
        window_tag = (
            f" window [{window[0]:.6g}, {window[1]:.6g}]" if window else ""
        )
        print(
            f"[{mark}] {row['field']} t={row['latest_time']:.9g}: "
            f"q_wall={row['measured_mean_wall_flux']:.6g} "
            f"({row['balance_field']} {row['balance_mode']}{window_tag} "
            f"@ t={row['balance_time']:.9g}), "
            f"q_exact={row['expected_mean_wall_flux']:.6g}, "
            f"balance error={100*row['balance_relative_error']:.2f}% "
            f"(limit {100*row['balance_tolerance_fraction']:.1f}%; "
            f"instantaneous q_wall={row['measured_mean_wall_flux_instantaneous']:.6g})"
        )
        print(
            f"       bulk range ({row['flatness_mean_type']})="
            f"{100*row['flatness_relative_range']:.3f}% over "
            f"{len(row['flatness_history'])} snapshots "
            f"(limit {100*row['flatness_tolerance_fraction']:.1f}%; "
            f"arithmetic={100*row['flatness_relative_range_arithmetic']:.3f}%)"
        )
    if result["statistics_allowed"]:
        print("STATISTICS ALLOWED: source balance and bulk flatness pass.")
    else:
        print("STATISTICS REFUSED: convergence gate failed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = evaluate_contract(args.case, args.contract)
    print_report(result)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    return 0 if result["statistics_allowed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Independent C9 scalar-balance review without using gate_verdict.json."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import numpy as np


SCALAR_RE = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;",
    re.S,
)
LABEL_RE = re.compile(r"\n\s*(\d+)\s*\n\s*\((.*?)\)\s*", re.S)
TOTAL_TIME_RE_TEMPLATE = (
    r"(?m)^\s*{field}\s*\{{\s*"
    r"(?:[^{{}}]|\{{[^{{}}]*\}})*?"
    r"^\s*totalTime\s+([\d.eE+-]+);"
)


def read_scalar(path: Path) -> np.ndarray:
    text = path.read_text(encoding="ascii", errors="replace")
    match = SCALAR_RE.search(text)
    if not match:
        raise ValueError(f"cannot parse scalar field: {path}")
    return np.fromstring(match.group(1), sep=" ")


def read_labels(path: Path) -> np.ndarray:
    text = path.read_text(encoding="ascii", errors="replace")
    match = LABEL_RE.search(text)
    if not match:
        raise ValueError(f"cannot parse label list: {path}")
    expected = int(match.group(1))
    values = np.fromstring(match.group(2), sep=" ", dtype=np.int64)
    if len(values) != expected:
        raise ValueError(f"{path}: parsed {len(values)} labels, expected {expected}")
    return values


def numeric_dirs(base: Path) -> list[tuple[float, str]]:
    rows: list[tuple[float, str]] = []
    if not base.exists():
        return rows
    for path in base.iterdir():
        if not path.is_dir():
            continue
        try:
            rows.append((float(path.name), path.name))
        except ValueError:
            pass
    return sorted(rows)


def processors(case: Path) -> list[Path]:
    paths = list(case.glob("processor[0-9]*"))
    return sorted(paths, key=lambda path: int(path.name.removeprefix("processor")))


def available_times(case: Path, field: str) -> list[tuple[float, str]]:
    choices: dict[float, str] = {}
    for value, name in numeric_dirs(case):
        if (case / name / field).exists():
            choices[value] = name
    procs = processors(case)
    if procs:
        for value, name in numeric_dirs(procs[0]):
            if (procs[0] / name / field).exists():
                choices[value] = name
    return sorted(choices.items())


def load_global(case: Path, time_name: str, field: str, size: int) -> np.ndarray:
    root = case / time_name / field
    if root.exists():
        values = read_scalar(root)
        if len(values) != size:
            raise ValueError(f"{root}: {len(values)} values, expected {size}")
        return values

    result = np.empty(size, dtype=float)
    covered = np.zeros(size, dtype=bool)
    for proc in processors(case):
        local_path = proc / time_name / field
        address_path = proc / "constant/polyMesh/cellProcAddressing"
        if not local_path.exists():
            raise FileNotFoundError(local_path)
        local = read_scalar(local_path)
        address = read_labels(address_path)
        if len(local) != len(address):
            raise ValueError(
                f"{proc.name}: {field} has {len(local)} cells but addressing has {len(address)}"
            )
        if np.any(address < 0) or np.any(address >= size):
            raise ValueError(f"{address_path}: global cell index outside 0..{size - 1}")
        if covered[address].any():
            raise ValueError(f"{address_path}: duplicate global cell ownership")
        result[address] = local
        covered[address] = True
    if not covered.all():
        raise ValueError(f"decomposition covers {covered.sum()} of {size} cells")
    return result


def find_root_scalar(case: Path, field: str, size: int) -> np.ndarray:
    for _, name in reversed(numeric_dirs(case)):
        path = case / name / field
        if not path.exists():
            continue
        values = read_scalar(path)
        if len(values) == size:
            return values
    raise FileNotFoundError(f"no root {field} field with {size} values under {case}")


def plane_profile(coordinate: np.ndarray, values: np.ndarray, decimals: int) -> tuple[np.ndarray, np.ndarray]:
    rounded = np.round(coordinate, decimals)
    levels = np.unique(rounded)
    profile = np.asarray([values[rounded == level].mean() for level in levels])
    return levels, profile


def wall_flux(
    coordinate: np.ndarray,
    values: np.ndarray,
    scalar: dict[str, Any],
    contract: dict[str, Any],
    fit_levels: int,
) -> dict[str, float]:
    levels, profile = plane_profile(
        coordinate, values, int(contract.get("coordinate_round_decimals", 9))
    )
    walls = [float(value) for value in contract["walls_m"]]
    wall_values = [float(value) for value in contract.get("wall_values", [0.0, 0.0])]
    lower_gradient = float(
        np.polyfit(
            np.r_[walls[0], levels[:fit_levels]],
            np.r_[wall_values[0], profile[:fit_levels]],
            1,
        )[0]
    )
    upper_gradient = float(
        np.polyfit(
            np.r_[levels[-fit_levels:], walls[1]],
            np.r_[profile[-fit_levels:], wall_values[1]],
            1,
        )[0]
    )
    diffusivity = float(scalar["molecular_diffusivity_m2_s"])
    lower = diffusivity * abs(lower_gradient)
    upper = diffusivity * abs(upper_gradient)
    return {
        "lower": lower,
        "upper": upper,
        "mean": 0.5 * (lower + upper),
        "wall_asymmetry_fraction": abs(upper - lower) / max(0.5 * (lower + upper), 1e-30),
    }


def relative_range(values: list[float]) -> float:
    array = np.asarray(values)
    return float((array.max() - array.min()) / max(abs(array.mean()), 1e-30))


def source_mean(scalar: dict[str, Any]) -> float:
    source = scalar["source"]
    if source["type"] != "constant_mean":
        raise ValueError(f"independent checker expects constant_mean, got {source['type']}")
    return float(source["mean_value_per_s"])


def accumulation_time(case: Path, time_name: str, field: str) -> float:
    """Read fieldAverage totalTime without sharing the canonical gate parser."""
    relative = Path("uniform/functionObjects/functionObjectProperties")
    candidates = [case / time_name / relative]
    candidates.extend(proc / time_name / relative for proc in processors(case))
    pattern = re.compile(
        TOTAL_TIME_RE_TEMPLATE.format(field=re.escape(field)), re.S
    )
    values: list[float] = []
    for path in candidates:
        if not path.exists():
            continue
        match = pattern.search(path.read_text(encoding="ascii", errors="replace"))
        if match:
            values.append(float(match.group(1)))
    if not values:
        raise ValueError(f"no fieldAverage totalTime for {field} at {time_name}")
    if max(values) - min(values) > 1e-10:
        raise ValueError(
            f"inconsistent fieldAverage totalTime for {field} at {time_name}: {values}"
        )
    return values[0]


def review(args: argparse.Namespace) -> dict[str, Any]:
    case = args.case.resolve()
    history_case = args.history_case.resolve() if args.history_case else None
    contract = json.loads((case / "constant/sourceDrivenScalarGate.json").read_text())
    volumes = read_scalar(args.volume_field.resolve())
    size = len(volumes)
    if np.any(volumes <= 0):
        raise ValueError("cell volumes must all be positive")
    coordinate = find_root_scalar(case, contract.get("coordinate_field", "Cy"), size)
    registered_start = float(contract["flatness_time_start"])
    balance_window_start = contract.get("balance_window_start")
    if balance_window_start is not None:
        balance_window_start = float(balance_window_start)
    result: dict[str, Any] = {
        "case": str(case),
        "history_case": str(history_case) if history_case else None,
        "volume_field": str(args.volume_field.resolve()),
        "cell_count": size,
        "volume_max_min_ratio": float(volumes.max() / volumes.min()),
        "registered_start": registered_start,
        "balance_window_start": balance_window_start,
        "strict_start": args.strict_start,
        "scalars": [],
    }

    for scalar in contract["scalars"]:
        field = scalar["field"]
        balance_field = scalar.get("balance_field", field)
        time_sources: dict[float, tuple[Path, str]] = {}
        if history_case:
            for value, name in available_times(history_case, field):
                if value <= args.history_max + 1e-12:
                    time_sources[value] = (history_case, name)
        for value, name in available_times(case, field):
            time_sources[value] = (case, name)
        times = sorted((value, source, name) for value, (source, name) in time_sources.items())
        if args.endpoint is not None:
            times = [row for row in times if row[0] <= args.endpoint + 1e-12]
        if not times:
            raise FileNotFoundError(f"no completed {field} writes found")
        exact = source_mean(scalar) * float(contract["source_to_sink_area_length_m"])
        rows = []
        for time_value, source_case, time_name in times:
            values = load_global(source_case, time_name, field, size)
            flux = wall_flux(
                coordinate,
                values,
                scalar,
                contract,
                int(contract.get("near_wall_fit_levels", 3)),
            )
            fit_sensitivity = {
                str(level): wall_flux(coordinate, values, scalar, contract, level)["mean"]
                for level in (1, 2, 3)
            }
            rows.append(
                {
                    "time": time_value,
                    "source_case": str(source_case),
                    "cell_mean": float(values.mean()),
                    "volume_mean": float(np.average(values, weights=volumes)),
                    "mean_wall_flux": flux["mean"],
                    "lower_wall_flux": flux["lower"],
                    "upper_wall_flux": flux["upper"],
                    "wall_asymmetry_fraction": flux["wall_asymmetry_fraction"],
                    "balance_error_fraction": (flux["mean"] - exact) / exact,
                    "fit_sensitivity": fit_sensitivity,
                }
            )

        balance_times = available_times(case, balance_field)
        if args.endpoint is not None:
            balance_times = [
                row for row in balance_times if row[0] <= args.endpoint + 1e-12
            ]
        if not balance_times:
            raise FileNotFoundError(
                f"no completed governing balance field {balance_field} found"
            )
        balance_time, balance_name = balance_times[-1]
        if abs(balance_time - rows[-1]["time"]) > 1e-9:
            raise ValueError(
                f"latest {field} time {rows[-1]['time']:.12g} does not have a "
                f"matching {balance_field}; newest mean field is "
                f"{balance_time:.12g}"
            )
        cumulative_values = load_global(case, balance_name, balance_field, size)
        cumulative_flux = wall_flux(
            coordinate,
            cumulative_values,
            scalar,
            contract,
            int(contract.get("near_wall_fit_levels", 3)),
        )
        balance_values = cumulative_values
        balance_mode = "cumulative_mean"
        balance_window = None
        accumulation_times = None
        if balance_window_start is not None:
            opening_writes = [
                (value, name)
                for value, name in balance_times
                if value + 1e-12 >= balance_window_start
                and value < balance_time - 1e-12
            ]
            if not opening_writes:
                raise ValueError(
                    f"no {balance_field} write in [{balance_window_start}, "
                    f"{balance_time}) to open independent balance window"
                )
            opening_time, opening_name = opening_writes[0]
            opening_values = load_global(case, opening_name, balance_field, size)
            early_total = accumulation_time(case, opening_name, field)
            late_total = accumulation_time(case, balance_name, field)
            if late_total <= early_total:
                raise ValueError(
                    f"non-increasing accumulation times for {field}: "
                    f"{early_total} -> {late_total}"
                )
            balance_values = (
                cumulative_values * late_total - opening_values * early_total
            ) / (late_total - early_total)
            balance_mode = "window_mean"
            balance_window = [opening_time, balance_time]
            accumulation_times = [early_total, late_total]
        governing_flux = wall_flux(
            coordinate,
            balance_values,
            scalar,
            contract,
            int(contract.get("near_wall_fit_levels", 3)),
        )
        governing_sensitivity = {
            str(level): wall_flux(
                coordinate, balance_values, scalar, contract, level
            )["mean"]
            for level in (1, 2, 3)
        }

        def summarize_window(start: float) -> dict[str, Any]:
            selected = [row for row in rows if row["time"] + 1e-12 >= start]
            if not selected:
                return {"start": start, "snapshots": 0}
            fluxes = [row["mean_wall_flux"] for row in selected]
            return {
                "start": start,
                "snapshots": len(selected),
                "cell_mean_relative_range": relative_range([row["cell_mean"] for row in selected]),
                "volume_mean_relative_range": relative_range([row["volume_mean"] for row in selected]),
                "mean_wall_flux": float(np.mean(fluxes)),
                "mean_balance_error_fraction": float(np.mean(fluxes) / exact - 1.0),
                "wall_flux_relative_range": relative_range(fluxes),
            }

        result["scalars"].append(
            {
                "field": field,
                "balance_field": balance_field,
                "expected_wall_flux": exact,
                "governing_balance": {
                    "time": balance_time,
                    "mode": balance_mode,
                    "window": balance_window,
                    "accumulation_times": accumulation_times,
                    "mean_wall_flux": governing_flux["mean"],
                    "lower_wall_flux": governing_flux["lower"],
                    "upper_wall_flux": governing_flux["upper"],
                    "wall_asymmetry_fraction": governing_flux[
                        "wall_asymmetry_fraction"
                    ],
                    "balance_error_fraction": governing_flux["mean"] / exact - 1.0,
                    "fit_sensitivity": governing_sensitivity,
                },
                "cumulative_balance": {
                    "time": balance_time,
                    "mean_wall_flux": cumulative_flux["mean"],
                    "lower_wall_flux": cumulative_flux["lower"],
                    "upper_wall_flux": cumulative_flux["upper"],
                    "wall_asymmetry_fraction": cumulative_flux[
                        "wall_asymmetry_fraction"
                    ],
                    "balance_error_fraction": cumulative_flux["mean"] / exact - 1.0,
                },
                "endpoint": rows[-1],
                "strict_window": summarize_window(args.strict_start),
                "registered_window": summarize_window(registered_start),
                "history": rows,
            }
        )
    return result


def print_report(result: dict[str, Any]) -> None:
    print(f"Independent C9 review: {result['case']}")
    print(
        f"cells={result['cell_count']}, Vmax/Vmin={result['volume_max_min_ratio']:.3f}, "
        f"strict_start={result['strict_start']}, registered_start={result['registered_start']}"
    )
    for scalar in result["scalars"]:
        endpoint = scalar["endpoint"]
        governing = scalar["governing_balance"]
        mode = governing.get("mode", "cumulative_mean")
        window = governing.get("window")
        window_text = f" window={window}" if window else ""
        print(
            f"{scalar['field']} governing {scalar['balance_field']} {mode}{window_text} "
            f"t={governing['time']:.9g}: q={governing['mean_wall_flux']:.8g}, "
            f"balance={100 * governing['balance_error_fraction']:+.3f}%, "
            f"wall asym={100 * governing['wall_asymmetry_fraction']:.3f}%"
        )
        cumulative = scalar.get("cumulative_balance")
        if cumulative and mode == "window_mean":
            print(
                f"  cumulative endpoint mean: q={cumulative['mean_wall_flux']:.8g}, "
                f"balance={100 * cumulative['balance_error_fraction']:+.3f}%"
            )
        print(
            f"  instantaneous endpoint t={endpoint['time']:.9g}: "
            f"q={endpoint['mean_wall_flux']:.8g}, "
            f"balance={100 * endpoint['balance_error_fraction']:+.3f}%, "
            f"wall asym={100 * endpoint['wall_asymmetry_fraction']:.3f}%"
        )
        for label in ("strict_window", "registered_window"):
            window = scalar[label]
            if not window["snapshots"]:
                print(f"  {label}: no snapshots")
                continue
            print(
                f"  {label}: n={window['snapshots']}, "
                f"flat(cell)={100 * window['cell_mean_relative_range']:.3f}%, "
                f"flat(volume)={100 * window['volume_mean_relative_range']:.3f}%, "
                f"mean-q balance={100 * window['mean_balance_error_fraction']:+.3f}%, "
                f"q range={100 * window['wall_flux_relative_range']:.3f}%"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--volume-field", required=True, type=Path)
    parser.add_argument("--strict-start", type=float, default=2.4)
    parser.add_argument("--endpoint", type=float)
    parser.add_argument("--history-case", type=Path)
    parser.add_argument("--history-max", type=float, default=4.3999142004455196)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = review(args)
    print_report(result)
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Minimal regression test for the source-driven scalar convergence gate."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile

import numpy as np

from scalar_convergence_gate import evaluate_contract


def write_scalar(path: Path, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"{value:.12g}" for value in values)
    path.write_text(
        "FoamFile { format ascii; class volScalarField; object test; }\n"
        f"internalField nonuniform List<scalar> {len(values)}\n(\n{body}\n)\n;\n"
        "boundaryField {}\n",
        encoding="ascii",
    )


def build_case(root: Path, drifting: bool) -> Path:
    case = root / ("drifting" if drifting else "balanced")
    coordinates = np.asarray([0.001, 0.002, 0.003, 0.017, 0.018, 0.019])
    base = np.minimum(coordinates, 0.02 - coordinates)
    for index, time in enumerate((1.0, 2.0, 3.0)):
        scale = 1.0 + (0.02 * index if drifting else 0.0)
        write_scalar(case / f"{time:g}" / "H", scale * base)
        write_scalar(case / f"{time:g}" / "Cy", coordinates)

    contract = {
        "version": 1,
        "balance_model": "closed_periodic_channel_with_two_fixed-value_wall_sinks",
        "coordinate_field": "Cy",
        "walls_m": [0.0, 0.02],
        "wall_values": [0.0, 0.0],
        "source_to_sink_area_length_m": 0.01,
        "near_wall_fit_levels": 3,
        "balance_tolerance_fraction": 0.05,
        "flatness_time_start": 1.0,
        "flatness_tolerance_fraction": 0.01,
        "minimum_flatness_snapshots": 3,
        "scalars": [
            {
                "field": "H",
                "molecular_diffusivity_m2_s": 1.0,
                "source": {"type": "constant_mean", "mean_value_per_s": 100.0},
            }
        ],
    }
    path = case / "constant" / "sourceDrivenScalarGate.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(contract), encoding="ascii")
    return case


def build_volume_case(root: Path, flatness_mean: str) -> Path:
    """Arithmetic bulk mean drifts, volume-weighted bulk mean is exactly flat.

    Cell 0 (small volume) is perturbed up while cell 5 (9x volume) compensates
    down by 1/9 of the amount, preserving sum(V*H) at every snapshot.
    """
    case = root / f"volume_{flatness_mean}"
    coordinates = np.asarray([0.001, 0.002, 0.003, 0.017, 0.018, 0.019])
    base = np.minimum(coordinates, 0.02 - coordinates)
    volumes = np.asarray([1.0, 1.0, 1.0, 1.0, 1.0, 9.0])
    for index, time in enumerate((1.0, 2.0, 3.0)):
        delta = 2.0e-4 * index
        values = base.copy()
        values[0] += delta
        values[5] -= delta / 9.0
        write_scalar(case / f"{time:g}" / "H", values)
        write_scalar(case / f"{time:g}" / "Cy", coordinates)
        write_scalar(case / f"{time:g}" / "V", volumes)

    contract = {
        "version": 2,
        "balance_model": "closed_periodic_channel_with_two_fixed-value_wall_sinks",
        "coordinate_field": "Cy",
        "walls_m": [0.0, 0.02],
        "wall_values": [0.0, 0.0],
        "source_to_sink_area_length_m": 0.01,
        "near_wall_fit_levels": 3,
        "balance_tolerance_fraction": 0.05,
        "flatness_time_start": 1.0,
        "flatness_tolerance_fraction": 0.01,
        "minimum_flatness_snapshots": 3,
        "flatness_mean": flatness_mean,
        "volume_field": "V",
        "scalars": [
            {
                "field": "H",
                "molecular_diffusivity_m2_s": 1.0,
                "source": {"type": "constant_mean", "mean_value_per_s": 100.0},
            }
        ],
    }
    path = case / "constant" / "sourceDrivenScalarGate.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(contract), encoding="ascii")
    return case


def build_balance_field_case(root: Path) -> Path:
    """Instantaneous endpoint flux is noisy; the mean field carries the balance.

    The final H snapshot moves value from cell 3 to cell 1 (equal volume), so
    the bulk mean is conserved (flatness passes) while the two-wall average
    gradient shifts ~10%, outside +/-5%. HMean stays exact, so a balance_field
    gate passes where an instantaneous gate would fail.
    """
    case = root / "balance_field"
    coordinates = np.asarray([0.001, 0.002, 0.003, 0.017, 0.018, 0.019])
    base = np.minimum(coordinates, 0.02 - coordinates)
    for index, time in enumerate((1.0, 2.0, 3.0)):
        values = base.copy()
        if index == 2:
            values[1] += 0.001
            values[3] -= 0.001
        write_scalar(case / f"{time:g}" / "H", values)
        write_scalar(case / f"{time:g}" / "HMean", base)
        write_scalar(case / f"{time:g}" / "Cy", coordinates)

    contract = {
        "version": 3,
        "balance_model": "closed_periodic_channel_with_two_fixed-value_wall_sinks",
        "coordinate_field": "Cy",
        "walls_m": [0.0, 0.02],
        "wall_values": [0.0, 0.0],
        "source_to_sink_area_length_m": 0.01,
        "near_wall_fit_levels": 3,
        "balance_tolerance_fraction": 0.05,
        "flatness_time_start": 1.0,
        "flatness_tolerance_fraction": 0.01,
        "minimum_flatness_snapshots": 3,
        "scalars": [
            {
                "field": "H",
                "balance_field": "HMean",
                "molecular_diffusivity_m2_s": 1.0,
                "source": {"type": "constant_mean", "mean_value_per_s": 100.0},
            }
        ],
    }
    path = case / "constant" / "sourceDrivenScalarGate.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(contract), encoding="ascii")
    return case


def write_accumulation_time(case: Path, time: float, field: str, total: float) -> None:
    path = case / f"{time:g}" / "uniform" / "functionObjects" / "functionObjectProperties"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "fieldAverage1\n{\n"
        f"    {field}\n    {{\n        totalIter 100;\n"
        f"        totalTime {total:.12g};\n    }}\n}}\n",
        encoding="ascii",
    )


def build_window_mean_case(root: Path, window_start: float) -> Path:
    """Cumulative HMean drags a filling transient; the late-window mean is exact.

    Accumulation runs [1, 3]. Cumulative means: M(2) = 0.8*base over T=1,
    M(3) = 0.9*base over T=2. The [2, 3] window mean is (0.9*2 - 0.8*1)*base
    = 1.0*base exactly, so window mode passes while the cumulative endpoint
    mean alone is 10% low.
    """
    case = root / f"window_{window_start:g}"
    coordinates = np.asarray([0.001, 0.002, 0.003, 0.017, 0.018, 0.019])
    base = np.minimum(coordinates, 0.02 - coordinates)
    for time, mean_scale, total in ((2.0, 0.8, 1.0), (3.0, 0.9, 2.0)):
        write_scalar(case / f"{time:g}" / "H", base)
        write_scalar(case / f"{time:g}" / "HMean", mean_scale * base)
        write_scalar(case / f"{time:g}" / "Cy", coordinates)
        write_accumulation_time(case, time, "H", total)
    write_scalar(case / "1" / "H", base)
    write_scalar(case / "1" / "Cy", coordinates)

    contract = {
        "version": 6,
        "balance_model": "closed_periodic_channel_with_two_fixed-value_wall_sinks",
        "coordinate_field": "Cy",
        "walls_m": [0.0, 0.02],
        "wall_values": [0.0, 0.0],
        "source_to_sink_area_length_m": 0.01,
        "near_wall_fit_levels": 3,
        "balance_tolerance_fraction": 0.05,
        "balance_window_start": window_start,
        "flatness_time_start": 1.0,
        "flatness_tolerance_fraction": 0.01,
        "minimum_flatness_snapshots": 3,
        "scalars": [
            {
                "field": "H",
                "balance_field": "HMean",
                "molecular_diffusivity_m2_s": 1.0,
                "source": {"type": "constant_mean", "mean_value_per_s": 100.0},
            }
        ],
    }
    path = case / "constant" / "sourceDrivenScalarGate.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(contract), encoding="ascii")
    return case


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="scalar-gate-") as directory:
        root = Path(directory)
        passing = evaluate_contract(build_case(root, drifting=False))
        drifting = evaluate_contract(build_case(root, drifting=True))
        assert passing["statistics_allowed"], passing
        assert not drifting["statistics_allowed"], drifting
        assert drifting["scalars"][0]["flatness_relative_range"] > 0.01

        volume = evaluate_contract(build_volume_case(root, "volume"))
        arithmetic = evaluate_contract(build_volume_case(root, "arithmetic"))
        assert volume["statistics_allowed"], volume
        assert volume["scalars"][0]["flatness_relative_range"] < 1e-12, volume
        assert volume["scalars"][0]["flatness_relative_range_arithmetic"] > 0.01
        assert not arithmetic["statistics_allowed"], arithmetic

        balanced_mean = evaluate_contract(build_balance_field_case(root))
        row = balanced_mean["scalars"][0]
        assert balanced_mean["statistics_allowed"], balanced_mean
        assert row["balance_field"] == "HMean", row
        assert row["balance_relative_error"] <= 0.05, row
        instantaneous_error = abs(
            row["measured_mean_wall_flux_instantaneous"]
            / row["expected_mean_wall_flux"]
            - 1.0
        )
        assert instantaneous_error > 0.05, row

        stale_mean = build_balance_field_case(root / "stale")
        (stale_mean / "3" / "HMean").unlink()
        try:
            evaluate_contract(stale_mean)
        except ValueError as exc:
            assert "does not have a matching HMean" in str(exc), exc
        else:
            raise AssertionError("gate accepted a stale balance_field")

        window = evaluate_contract(build_window_mean_case(root, 1.5))
        row = window["scalars"][0]
        assert window["statistics_allowed"], window
        assert row["balance_mode"] == "window_mean", row
        assert row["balance_window"] == [2.0, 3.0], row
        assert row["balance_relative_error"] <= 1e-9, row
        cumulative_error = abs(
            row["measured_mean_wall_flux_cumulative"]
            / row["expected_mean_wall_flux"]
            - 1.0
        )
        assert abs(cumulative_error - 0.10) < 1e-9, row

        try:
            evaluate_contract(build_window_mean_case(root / "late", 2.5))
        except ValueError as exc:
            assert "to open the balance window" in str(exc), exc
        else:
            raise AssertionError("gate accepted an unopenable balance window")
    print("scalar convergence gate regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

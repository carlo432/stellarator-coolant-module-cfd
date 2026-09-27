#!/usr/bin/env python3
"""Compute a three-grid Grid Convergence Index summary from a CSV table.

Input CSV columns:
    case,mesh_level,cells,h,metric_name,value

Rows may be in any order. The script sorts by grid spacing h, with the
smallest h treated as the finest grid. If h is blank, h is estimated as
cells^(-1/3), which is only a representative spacing.

The script refuses formal GCI for non-monotonic three-grid sequences.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value:
        return None
    return float(value)


def read_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"case", "mesh_level", "cells", "h", "metric_name", "value"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")
        for raw in reader:
            cells = parse_float(raw["cells"])
            h = parse_float(raw["h"])
            value = parse_float(raw["value"])
            if cells is None or value is None:
                raise ValueError(f"Row must include cells and value: {raw}")
            if h is None:
                h = cells ** (-1.0 / 3.0)
            rows.append(
                {
                    "case": raw["case"].strip(),
                    "mesh_level": raw["mesh_level"].strip(),
                    "cells": cells,
                    "h": h,
                    "metric_name": raw["metric_name"].strip(),
                    "value": value,
                }
            )
    if len(rows) != 3:
        raise ValueError(f"Expected exactly 3 grid rows, found {len(rows)}")
    return sorted(rows, key=lambda row: float(row["h"]))


def observed_order(phi1: float, phi2: float, phi3: float, r21: float, r32: float) -> float:
    e21 = phi2 - phi1
    e32 = phi3 - phi2
    ratio = abs(e32 / e21)
    if abs(r21 - r32) / max(r21, r32) < 0.05:
        return abs(math.log(ratio) / math.log((r21 + r32) / 2.0))

    sign = 1.0 if (e32 / e21) > 0 else -1.0
    p = abs(math.log(ratio) / math.log(math.sqrt(r21 * r32)))
    for _ in range(100):
        numerator = (r21**p) - sign
        denominator = (r32**p) - sign
        if numerator <= 0 or denominator <= 0:
            break
        q = math.log(numerator / denominator)
        next_p = abs((math.log(ratio) + q) / math.log(r21))
        if abs(next_p - p) < 1e-8:
            return next_p
        p = next_p
    return p


def summarize(rows: list[dict[str, object]], safety_factor: float) -> list[tuple[str, str]]:
    fine, medium, coarse = rows
    phi1 = float(fine["value"])
    phi2 = float(medium["value"])
    phi3 = float(coarse["value"])
    h1 = float(fine["h"])
    h2 = float(medium["h"])
    h3 = float(coarse["h"])
    r21 = h2 / h1
    r32 = h3 / h2
    e21 = phi2 - phi1
    e32 = phi3 - phi2

    out: list[tuple[str, str]] = [
        ("case", str(fine["case"])),
        ("metric_name", str(fine["metric_name"])),
        ("finest_mesh_level", str(fine["mesh_level"])),
        ("medium_mesh_level", str(medium["mesh_level"])),
        ("coarsest_mesh_level", str(coarse["mesh_level"])),
        ("finest_value", f"{phi1:.10g}"),
        ("medium_value", f"{phi2:.10g}"),
        ("coarse_value", f"{phi3:.10g}"),
        ("r21", f"{r21:.10g}"),
        ("r32", f"{r32:.10g}"),
    ]

    if e21 == 0 or e32 == 0:
        out.extend(
            [
                ("status", "zero_difference"),
                ("formal_gci_valid", "false"),
                ("reason", "At least one grid-pair difference is zero; observed order is undefined."),
            ]
        )
        return out

    if e21 * e32 <= 0:
        out.extend(
            [
                ("status", "non_monotonic"),
                ("formal_gci_valid", "false"),
                ("reason", "Three-grid sequence is non-monotonic; do not compute formal GCI."),
            ]
        )
        return out

    p = observed_order(phi1, phi2, phi3, r21, r32)
    denominator = (r21**p) - 1.0
    if denominator <= 0:
        out.extend(
            [
                ("status", "invalid_order"),
                ("formal_gci_valid", "false"),
                ("reason", "Observed order produced a non-positive GCI denominator."),
            ]
        )
        return out

    extrapolated = ((r21**p) * phi1 - phi2) / denominator
    approx_error_21 = abs((phi1 - phi2) / phi1)
    extrapolated_error = abs((extrapolated - phi1) / extrapolated) if extrapolated else math.nan
    gci_21 = safety_factor * approx_error_21 / denominator

    out.extend(
        [
            ("status", "monotonic"),
            ("formal_gci_valid", "true"),
            ("observed_order_p", f"{p:.10g}"),
            ("extrapolated_value", f"{extrapolated:.10g}"),
            ("approx_relative_error_21_percent", f"{100.0 * approx_error_21:.10g}"),
            ("extrapolated_relative_error_percent", f"{100.0 * extrapolated_error:.10g}"),
            ("gci_21_percent", f"{100.0 * gci_21:.10g}"),
            ("safety_factor", f"{safety_factor:.10g}"),
        ]
    )
    return out


def write_summary(path: Path, summary: list[tuple[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value"])
        writer.writerows(summary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--safety-factor", type=float, default=1.25)
    args = parser.parse_args()

    rows = read_rows(args.input_csv)
    metric_names = {str(row["metric_name"]) for row in rows}
    case_names = {str(row["case"]) for row in rows}
    if len(metric_names) != 1:
        raise ValueError(f"Input rows must describe one metric, found {sorted(metric_names)}")
    if len(case_names) != 1:
        raise ValueError(f"Input rows must describe one case, found {sorted(case_names)}")

    summary = summarize(rows, args.safety_factor)
    write_summary(args.output, summary)


if __name__ == "__main__":
    main()

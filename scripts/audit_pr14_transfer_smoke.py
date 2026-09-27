#!/usr/bin/env python3
"""Audit the Pr=14.4 transformed-scalar smoke without promoting physics."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re


REQUIRED_FIELDS = (
    "H", "HMean", "HPrime2Mean", "U", "UMean", "UPrime2Mean",
    "yPlus", "Cy", "V",
)


def numeric_times(case: Path) -> list[tuple[float, Path]]:
    rows = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            rows.append((float(path.name), path))
        except ValueError:
            pass
    return sorted(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    case = args.case.resolve()
    transfer = json.loads(
        (case / "constant/highPrTransferContract.json").read_text(encoding="ascii")
    )
    gate_contract = json.loads(
        (case / "constant/sourceDrivenScalarGate.json").read_text(encoding="ascii")
    )
    gate = json.loads((case / "gate_verdict_smoke.json").read_text(encoding="ascii"))
    latest_value, latest = numeric_times(case)[-1]
    control = (case / "system/controlDict").read_text(encoding="ascii")
    solver_log = (case / "log.pimpleFoam").read_text(encoding="ascii", errors="replace")
    mesh_log = (case / "log.checkMesh").read_text(encoding="ascii", errors="replace")

    missing = [field for field in REQUIRED_FIELDS if not (latest / field).is_file()]
    nonfinite = []
    token = re.compile(r"(?<![A-Za-z])(?:nan|[+-]?inf)(?![A-Za-z])", re.I)
    for field in REQUIRED_FIELDS:
        path = latest / field
        if path.is_file() and token.search(path.read_text(encoding="ascii", errors="replace")):
            nonfinite.append(field)

    expected_alpha_d = 1.0 / 14.4
    expected_alpha_dt = 1.0 / 0.9
    alpha_d_match = re.search(r"\balphaD\s+([0-9.eE+-]+)\s*;", control)
    alpha_dt_match = re.search(r"\balphaDt\s+([0-9.eE+-]+)\s*;", control)
    alpha_d = float(alpha_d_match.group(1)) if alpha_d_match else math.nan
    alpha_dt = float(alpha_dt_match.group(1)) if alpha_dt_match else math.nan

    checks = {
        "mesh_ok": "Mesh OK." in mesh_log,
        "solver_clean_end": solver_log.rstrip().endswith("End") and "FOAM FATAL" not in solver_log,
        "endpoint_reached": latest_value >= 0.01,
        "required_fields_present": not missing,
        "required_fields_finite": not nonfinite,
        "alphaD_is_1_over_14p4": math.isclose(alpha_d, expected_alpha_d, rel_tol=0, abs_tol=1e-12),
        "alphaDt_is_1_over_0p9": math.isclose(alpha_dt, expected_alpha_dt, rel_tol=0, abs_tol=1e-11),
        "v7_gate_contract": (
            gate_contract.get("version") == 7
            and gate_contract.get("flatness_mean") == "volume"
            and gate_contract["scalars"][0].get("balance_field") == "HMean"
        ),
        "smoke_statistics_refused_as_expected": (
            gate.get("status") == "FAIL" and gate.get("statistics_allowed") is False
        ),
        "no_pr14_kawamura_claim": (
            transfer["claim_limits"]["kawamura_p1_p4_certification_allowed"] is False
            and transfer["claim_limits"]["stellarator_geometry_validation"] is False
        ),
        "production_seed_registered": (
            transfer["production_requirements"]["default_initialization"]
            == "Kader-profile seeded near equilibrium"
        ),
    }
    result = {
        "case": str(case),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "interpretation": "Pr=14.4 transfer machinery only; no physics certification",
        "latest_time": latest_value,
        "cells": 27648,
        "alphaD": alpha_d,
        "alphaDt": alpha_dt,
        "physical_gate_status": gate["status"],
        "physical_gate_balance_error_percent": 100.0 * gate["scalars"][0]["balance_relative_error"],
        "physical_gate_flatness_percent": 100.0 * gate["scalars"][0]["flatness_relative_range"],
        "missing_fields": missing,
        "nonfinite_fields": nonfinite,
        "checks": checks,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(f"Pr=14.4 transfer smoke: {result['status']}")
    for name, passed in checks.items():
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
    print(f"wrote {args.json}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

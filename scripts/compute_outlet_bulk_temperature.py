#!/usr/bin/env python3
"""Compute area- and flux-weighted outlet temperatures from saved OpenFOAM fields."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(".")
TIN_K = 800.0
RHO_KG_M3 = 1940.0
CP_J_KG_K = 2400.0
HEATED_AREA_M2 = 0.0080


@dataclass(frozen=True)
class CaseSpec:
    label: str
    case_dir: Path
    heat_flux_w_m2: float
    time_name: str | None = None


DEFAULT_CASES = [
    CaseSpec(
        label="v14 100 kW/m^2",
        case_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025",
        heat_flux_w_m2=100000.0,
        time_name="2",
    ),
    CaseSpec(
        label="v15 100 kW/m^2 long-run",
        case_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025",
        heat_flux_w_m2=100000.0,
        time_name="4.8",
    ),
    CaseSpec(
        label="v16 500 kW/m^2",
        case_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025",
        heat_flux_w_m2=500000.0,
        time_name="2",
    ),
    CaseSpec(
        label="v23 offset 100 kW/m^2",
        case_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v23_scalar_temperature_flibe_re10000_heatflux100kw_bent_recirculation_dt00025",
        heat_flux_w_m2=100000.0,
        time_name="2",
    ),
    CaseSpec(
        label="v24 offset 500 kW/m^2",
        case_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v24_scalar_temperature_flibe_re10000_heatflux500kw_bent_recirculation_dt00025",
        heat_flux_w_m2=500000.0,
        time_name="2",
    ),
]


def parse_time_name(name: str) -> float | None:
    try:
        return float(name)
    except ValueError:
        return None


def latest_time_name(case_dir: Path) -> str:
    times: list[tuple[float, str]] = []
    for child in case_dir.iterdir():
        if not child.is_dir():
            continue
        value = parse_time_name(child.name)
        if value is not None:
            times.append((value, child.name))
    if not times:
        raise ValueError(f"No numeric time directories found in {case_dir}")
    return max(times)[1]


def read_patch_scalar_values(field_path: Path, patch_name: str) -> list[float]:
    text = field_path.read_text(encoding="ascii")
    patch_match = re.search(
        rf"\b{re.escape(patch_name)}\s*\{{(?P<body>.*?)\n\s*\}}",
        text,
        re.DOTALL,
    )
    if not patch_match:
        raise ValueError(f"Patch {patch_name!r} not found in {field_path}")

    body = patch_match.group("body")
    value_match = re.search(
        r"value\s+nonuniform\s+List<scalar>\s+(?P<count>\d+)\s*\((?P<values>.*?)\)\s*;",
        body,
        re.DOTALL,
    )
    if not value_match:
        raise ValueError(f"Patch {patch_name!r} in {field_path} does not have nonuniform scalar values")

    expected_count = int(value_match.group("count"))
    values = [float(item) for item in value_match.group("values").split()]
    if len(values) != expected_count:
        raise ValueError(
            f"Patch {patch_name!r} in {field_path} declared {expected_count} values but parsed {len(values)}"
        )
    return values


def compute_case(case: CaseSpec, patch_name: str) -> dict[str, str | float | int]:
    time_name = case.time_name or latest_time_name(case.case_dir)
    time_dir = case.case_dir / time_name
    temperatures = read_patch_scalar_values(time_dir / "T", patch_name)
    fluxes = read_patch_scalar_values(time_dir / "phi", patch_name)

    if len(temperatures) != len(fluxes):
        raise ValueError(
            f"Outlet T count {len(temperatures)} does not match phi count {len(fluxes)} for {case.case_dir}"
        )

    flux_sum = sum(fluxes)
    if abs(flux_sum) <= 1e-18:
        raise ValueError(f"Outlet flux is too close to zero for {case.case_dir}")

    area_average_t = sum(temperatures) / len(temperatures)
    flux_weighted_t = sum(t * phi for t, phi in zip(temperatures, fluxes)) / flux_sum
    mass_flow = RHO_KG_M3 * abs(flux_sum)
    heat_input = case.heat_flux_w_m2 * HEATED_AREA_M2
    energy_balance_delta_t = heat_input / (mass_flow * CP_J_KG_K)

    return {
        "case": case.label,
        "case_dir": str(case.case_dir),
        "time": time_name,
        "patch": patch_name,
        "heat_flux_w_m2": case.heat_flux_w_m2,
        "heat_input_w": heat_input,
        "outlet_face_count": len(temperatures),
        "outlet_volumetric_flow_m3_s": flux_sum,
        "outlet_mass_flow_kg_s": mass_flow,
        "area_avg_T_K": area_average_t,
        "area_avg_deltaT_K": area_average_t - TIN_K,
        "mass_flow_weighted_T_K": flux_weighted_t,
        "mass_flow_weighted_deltaT_K": flux_weighted_t - TIN_K,
        "energy_balance_deltaT_K": energy_balance_delta_t,
        "area_minus_mass_deltaT_K": area_average_t - flux_weighted_t,
        "mass_minus_energy_deltaT_K": (flux_weighted_t - TIN_K) - energy_balance_delta_t,
    }


def write_csv(path: Path, rows: list[dict[str, str | float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case",
        "case_dir",
        "time",
        "patch",
        "heat_flux_w_m2",
        "heat_input_w",
        "outlet_face_count",
        "outlet_volumetric_flow_m3_s",
        "outlet_mass_flow_kg_s",
        "area_avg_T_K",
        "area_avg_deltaT_K",
        "mass_flow_weighted_T_K",
        "mass_flow_weighted_deltaT_K",
        "energy_balance_deltaT_K",
        "area_minus_mass_deltaT_K",
        "mass_minus_energy_deltaT_K",
    ]
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, help="OpenFOAM case directory for a single-case calculation")
    parser.add_argument("--label", default="custom case", help="Label for --case-dir output")
    parser.add_argument("--time", dest="time_name", help="Time directory to read; defaults to latest numeric time")
    parser.add_argument("--heat-flux", type=float, help="Heat flux in W/m^2 for --case-dir")
    parser.add_argument("--patch", default="outlet", help="Patch name to evaluate")
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results/outlet_temperature_metrics/mass_weighted_outlet_temperature.csv",
        help="Output CSV path",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.case_dir:
        if args.heat_flux is None:
            raise SystemExit("--heat-flux is required when --case-dir is used")
        cases = [CaseSpec(args.label, args.case_dir, args.heat_flux, args.time_name)]
    else:
        cases = DEFAULT_CASES

    rows = [compute_case(case, args.patch) for case in cases]
    write_csv(args.out, rows)
    print(args.out)
    for row in rows:
        print(
            f"{row['case']}: area dT={float(row['area_avg_deltaT_K']):.6f} K, "
            f"mass-weighted dT={float(row['mass_flow_weighted_deltaT_K']):.6f} K, "
            f"energy-balance dT={float(row['energy_balance_deltaT_K']):.6f} K"
        )


if __name__ == "__main__":
    main()

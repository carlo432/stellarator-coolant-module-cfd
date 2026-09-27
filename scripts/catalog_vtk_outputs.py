#!/usr/bin/env python3
"""Catalog top-level VTK exports from the OpenFOAM cases."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASES_DIR = ROOT / "cases" / "baseline_module" / "openfoam_cases"
OUT_DIR = ROOT / "results" / "visualization"
OUT_CSV = OUT_DIR / "vtk_manifest.csv"


NOTES = {
    "v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025": (
        1,
        "Best first thermal view: refined FLiBe-like T field and heated wall.",
    ),
    "v13_isothermal_rans_flibe_re10000_thermal_refined": (
        2,
        "Best refined flow view: U, p, and recirculation structure.",
    ),
    "v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025": (
        3,
        "High heat-flux thermal view: 500 kW/m^2 bracket on the refined mesh.",
    ),
    "v17_unsteady_rans_flibe_re10000_thermal_refined_pimple": (
        4,
        "First transient RANS smoke-test export: U, p, and turbulence fields at the final short-run time.",
    ),
    "v18_unsteady_rans_flibe_re10000_thermal_refined_pimple_0p1s": (
        5,
        "Transient RANS continuation to about 0.1 s for early unsteady stability and visualization.",
    ),
    "v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005": (
        6,
        "First refined thermal comparison point.",
    ),
    "v9_scalar_temperature_flibe_re10000_heatflux100kw": (
        7,
        "Coarse thermal comparison point; useful for showing mesh sensitivity.",
    ),
    "v8_isothermal_rans_flibe_re10000_clean": (
        8,
        "Coarse FLiBe-like RANS baseline.",
    ),
    "v7_scalar_temperature_water_re10000_heatflux100kw": (
        9,
        "Water 100 kW/m^2 thermal debug case.",
    ),
    "v5_isothermal_rans_water_re10000_clean": (
        10,
        "Water Re=10000 RANS debug baseline.",
    ),
}


def classify(path: Path) -> str:
    if path.name.endswith(".vtm.series"):
        return "time_series"
    return "single_time"


def time_label(path: Path) -> str:
    if path.name.endswith(".vtm.series"):
        return "series"
    stem = path.stem
    maybe_time = stem.rsplit("_", 1)[-1]
    return maybe_time if maybe_time.isdigit() else "unknown"


def main() -> None:
    rows = []

    for case_dir in sorted(CASES_DIR.iterdir()):
        vtk_dir = case_dir / "VTK"
        if not vtk_dir.is_dir():
            continue

        priority, note = NOTES.get(case_dir.name, (99, "OpenFOAM VTK export."))
        for path in sorted(vtk_dir.iterdir()):
            if not path.is_file():
                continue
            if not (path.name.endswith(".vtm") or path.name.endswith(".vtm.series")):
                continue
            rows.append(
                {
                    "priority": priority,
                    "case": case_dir.name,
                    "kind": classify(path),
                    "time": time_label(path),
                    "path": str(path.relative_to(ROOT)),
                    "note": note,
                }
            )

    rows.sort(key=lambda row: (int(row["priority"]), row["case"], row["kind"], row["time"]))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["priority", "case", "kind", "time", "path", "note"])
        writer.writeheader()
        writer.writerows(rows)

    print(OUT_CSV.relative_to(ROOT))
    print(f"{len(rows)} top-level VTK entries")


if __name__ == "__main__":
    main()

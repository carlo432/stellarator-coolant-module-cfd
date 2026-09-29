#!/usr/bin/env python3
"""Create a C3 CHT case with material-split volumetric heating.

This is the run-ready follow-on to the C3 all-fluid volumetric smoke case. It
reuses the C2 curved CHT geometry, disables the external base flux, and injects
same-total-power volumetric sources into both the FLiBe-like fluid and the solid
heater region. The split is read from results/c3_material_split_targets.csv.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path


def polygon_area(points: list[list[float]]) -> float:
    if len(points) < 3:
        return 0.0
    p0 = points[0]
    area_vec = [0.0, 0.0, 0.0]
    for i in range(1, len(points) - 1):
        a = [points[i][j] - p0[j] for j in range(3)]
        b = [points[i + 1][j] - p0[j] for j in range(3)]
        cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
        area_vec = [area_vec[j] + cross[j] for j in range(3)]
    return 0.5 * math.sqrt(sum(x * x for x in area_vec))


def patch_area(geometry: dict, patch_name: str) -> float:
    vertices = geometry["vertices"]
    return sum(polygon_area([vertices[idx] for idx in face]) for face in geometry["patches"][patch_name])


def load_target(path: Path, target_id: str) -> dict[str, str]:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["target_id"] == target_id:
                return row
    raise ValueError(f"target_id {target_id!r} not found in {path}")


def write_region_fv_options(path: Path, *, region: str, option_name: str, power_w: float) -> None:
    path.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2512/v2506                           |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system/{region}";
    object      fvOptions;
}}

{option_name}
{{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   all;
    volumeMode      absolute;

    sources
    {{
        h           ( {power_w:.12g} 0 );
    }}
}}
"""
    )


def main() -> None:
    work = Path(__file__).resolve().parent
    repo = Path(__file__).resolve().parents[2]
    default_targets = repo / "results/c3_material_split_targets.csv"

    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="c3_curved_material_split_55_45_cht")
    parser.add_argument("--target-id", default="c3_ms_total_blanket_55_45")
    parser.add_argument("--target-csv", default=str(default_targets))
    parser.add_argument("--radial-cells", type=int, default=50)
    parser.add_argument("--solid-cells", type=int, default=8)
    parser.add_argument("--wall-adjacent-cell", type=float, default=5e-5)
    parser.add_argument("--end-time", type=int, default=400)
    parser.add_argument("--write-interval", type=int, default=400)
    parser.add_argument("--target-heat-flux", type=float, default=100000.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    target_path = Path(args.target_csv)
    target = load_target(target_path, args.target_id)
    fluid_fraction = float(target["fluid_fraction"])
    solid_fraction = float(target["solid_fraction"])
    omitted_fraction = float(target["omitted_fraction"])
    if abs(fluid_fraction + solid_fraction + omitted_fraction - 1.0) > 1e-8:
        raise ValueError(f"target fractions do not sum to one: {target}")
    if omitted_fraction:
        raise ValueError(
            f"{args.target_id} omits {omitted_fraction:.6g} of power; choose a two-region target for this generator"
        )

    generator = work / "gen_c2_curved_cht.py"
    cmd = [
        sys.executable,
        str(generator),
        "--case",
        args.case,
        "--radial-cells",
        str(args.radial_cells),
        "--solid-cells",
        str(args.solid_cells),
        "--wall-adjacent-cell",
        str(args.wall_adjacent_cell),
        "--end-time",
        str(args.end_time),
        "--write-interval",
        str(args.write_interval),
        "--heat-flux",
        "0",
    ]
    if args.overwrite:
        cmd.append("--overwrite")
    subprocess.run(cmd, cwd=work, check=True)

    case = work / args.case
    geometry = json.loads((case / "geometry_c2cht.json").read_text())
    heated_area = patch_area(geometry, "heatedBase")
    total_power = args.target_heat_flux * heated_area
    fluid_power = total_power * fluid_fraction
    solid_power = total_power * solid_fraction

    write_region_fv_options(
        case / "system/bottomWater/fvOptions",
        region="bottomWater",
        option_name="c3MaterialSplitFluidHeating",
        power_w=fluid_power,
    )
    write_region_fv_options(
        case / "system/heater/fvOptions",
        region="heater",
        option_name="c3MaterialSplitSolidHeating",
        power_w=solid_power,
    )

    metadata = {
        "case": str(case),
        "basis": "C2 curved CHT geometry with external base flux disabled",
        "target_id": args.target_id,
        "target_csv": str(target_path),
        "literature_basis": target.get("literature_basis"),
        "solver_use": target.get("solver_use"),
        "target_heat_flux_W_m2": args.target_heat_flux,
        "heated_base_area_m2": heated_area,
        "total_power_W": total_power,
        "target_card_total_reference_power_W": float(target["total_reference_power_W"]),
        "fluid_fraction": fluid_fraction,
        "solid_fraction": solid_fraction,
        "fluid_power_W": fluid_power,
        "solid_power_W": solid_power,
        "target_card_fluid_power_W": float(target["fluid_power_W"]),
        "target_card_solid_power_W": float(target["solid_power_W"]),
        "source_method": "scalarSemiImplicitSource on h in bottomWater and heater regions; volumeMode absolute",
        "radial_cells": args.radial_cells,
        "solid_cells": args.solid_cells,
        "end_time": args.end_time,
        "write_interval": args.write_interval,
    }
    (case / "c3_material_split_cht_manifest.json").write_text(json.dumps(metadata, indent=2))
    (case / "README.c3_material_split.md").write_text(
        "\n".join(
            [
                "# C3 Curved Material-Split Volumetric CHT",
                "",
                "This case reuses the C2 curved CHT geometry, disables the external",
                "heated-base gradient, and injects same-total-power volumetric sources",
                "into the FLiBe-like fluid and solid heater regions.",
                "",
                f"- target id: `{args.target_id}`",
                f"- target equivalent heat flux: `{args.target_heat_flux:g} W/m^2`",
                f"- heated-base reference area: `{heated_area:.8g} m^2`",
                f"- total power: `{total_power:.8g} W`",
                f"- fluid power: `{fluid_power:.8g} W`",
                f"- solid power: `{solid_power:.8g} W`",
                f"- target source table: `{target_path}`",
                "",
                "Run:",
                "",
                "```bash",
                "./Allpre.c2cht",
                "./Allrun.c2cht",
                "```",
                "",
                "This is a literature-guided source placement case, not a validated",
                "stellarator neutronics prediction.",
            ]
        )
        + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create a C3 CHT case with same-total-power volumetric FLiBe heating.

The case reuses the C2 curved CHT geometry but turns off the external base
flux and injects heat into the fluid region via fvOptions.  This is a C3
plumbing check for one-way neutronics-to-thermal coupling, not a final blanket
calculation.
"""
from __future__ import annotations

import argparse
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


def write_fv_options(path: Path, total_power_w: float) -> None:
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
    location    "system/bottomWater";
    object      fvOptions;
}}

flibeVolumetricHeating
{{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   all;
    volumeMode      absolute;

    sources
    {{
        h           ( {total_power_w:.12g} 0 );
    }}
}}
"""
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="c3_curved_openmc_volumetric_cht_smoke")
    parser.add_argument("--radial-cells", type=int, default=50)
    parser.add_argument("--solid-cells", type=int, default=8)
    parser.add_argument("--wall-adjacent-cell", type=float, default=5e-5)
    parser.add_argument("--end-time", type=int, default=400)
    parser.add_argument("--write-interval", type=int, default=400)
    parser.add_argument("--target-heat-flux", type=float, default=100000.0)
    parser.add_argument(
        "--openmc-summary",
        default="figs/c3_openmc_flibe_slab_heating_profile.json",
        help="Optional OpenMC heating-profile summary to link in metadata.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    work = Path(__file__).resolve().parent
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
    geom_path = case / "geometry_c2cht.json"
    geometry = json.loads(geom_path.read_text())
    heated_area = patch_area(geometry, "heatedBase")
    total_power = args.target_heat_flux * heated_area
    write_fv_options(case / "system/bottomWater/fvOptions", total_power)

    openmc_summary_path = work / args.openmc_summary
    openmc_summary = None
    if openmc_summary_path.exists():
        openmc_summary = json.loads(openmc_summary_path.read_text())

    metadata = {
        "case": str(case),
        "basis": "C2 curved CHT geometry with external base flux disabled",
        "volumetric_region": "bottomWater",
        "source_method": "scalarSemiImplicitSource on h, volumeMode absolute",
        "target_heat_flux_W_m2": args.target_heat_flux,
        "heated_base_area_m2": heated_area,
        "total_power_W": total_power,
        "radial_cells": args.radial_cells,
        "solid_cells": args.solid_cells,
        "end_time": args.end_time,
        "write_interval": args.write_interval,
        "openmc_summary": str(openmc_summary_path) if openmc_summary_path.exists() else None,
        "openmc_source_kind": openmc_summary.get("source_kind") if openmc_summary else None,
        "openmc_profile_peak_to_mean": openmc_summary.get("normalized_qvol_peak_to_mean") if openmc_summary else None,
    }
    (case / "c3_volumetric_cht_manifest.json").write_text(json.dumps(metadata, indent=2))
    readme = case / "README.c3cht.md"
    readme.write_text(
        "\n".join(
            [
                "# C3 Curved OpenMC-Normalized Volumetric CHT Smoke",
                "",
                "This case reuses the C2 curved CHT geometry, sets the external",
                "heated-base gradient to zero, and injects the same total power into",
                "the FLiBe-like fluid region with `fvOptions`.",
                "",
                f"- target equivalent heat flux: `{args.target_heat_flux:g} W/m^2`",
                f"- heated-base reference area: `{heated_area:.8g} m^2`",
                f"- total power injected into fluid: `{total_power:.8g} W`",
                f"- OpenMC profile linked: `{metadata['openmc_summary']}`",
                "",
                "Run:",
                "",
                "```bash",
                "./Allpre.c2cht",
                "./Allrun.c2cht",
                "```",
                "",
                "This is a one-way coupling smoke case. The OpenMC depth profile is",
                "not yet mapped into per-bin OpenFOAM cell zones; only the same total",
                "power normalization is injected here.",
            ]
        )
        + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()


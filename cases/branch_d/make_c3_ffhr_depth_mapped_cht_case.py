#!/usr/bin/env python3
"""Create a C3 CHT case with FFHR-like depth-mapped fluid/solid heating.

This is the mapped-source follow-on to the aggregated FFHR-like Be/FLiBe C3
case. It still uses the simplified C2 curved CHT mesh, so the OpenMC stack
depth is resampled onto streamwise stations as a plumbing/sensitivity exercise,
not as a physical wall-normal source map.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import stat
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


def load_region_depth_profile(path: Path) -> tuple[float, list[dict[str, float | str]]]:
    rows: list[dict[str, float | str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            shape = max(float(row["normalized_shape"]), 0.0)
            rows.append(
                {
                    "region": row["region"],
                    "material_class": row["material_class"],
                    "x0_m": float(row["x0_m"]),
                    "x1_m": float(row["x1_m"]),
                    "normalized_shape": shape,
                }
            )
    if not rows:
        raise ValueError(f"no profile rows found in {path}")
    total_shape = sum(float(row["normalized_shape"]) for row in rows)
    if total_shape <= 0.0:
        raise ValueError(f"profile shape is non-positive in {path}")
    thickness = max(float(row["x1_m"]) for row in rows)
    for row in rows:
        row["normalized_shape"] = float(row["normalized_shape"]) / total_shape
        row["x0_hat"] = float(row["x0_m"]) / thickness
        row["x1_hat"] = float(row["x1_m"]) / thickness
    return thickness, rows


def map_class_to_segments(
    rows: list[dict[str, float | str]],
    *,
    material_class: str,
    segments: int,
) -> list[float]:
    fractions: list[float] = []
    selected = [row for row in rows if row["material_class"] == material_class]
    for seg in range(segments):
        z0 = seg / segments
        z1 = (seg + 1) / segments
        frac = 0.0
        for row in selected:
            x0_hat = float(row["x0_hat"])
            x1_hat = float(row["x1_hat"])
            overlap = max(0.0, min(z1, x1_hat) - max(z0, x0_hat))
            width = x1_hat - x0_hat
            if overlap > 0.0 and width > 0.0:
                frac += float(row["normalized_shape"]) * overlap / width
        fractions.append(frac)
    total = sum(fractions)
    if total <= 0.0:
        raise ValueError(f"{material_class} mapped fractions sum to zero")
    return fractions


def zone_block(name: str, start: int, stop: int) -> list[str]:
    labels = "\n".join(str(i) for i in range(start, stop))
    return [
        name,
        "{",
        "    type cellZone;",
        "    cellLabels List<label>",
        str(stop - start),
        "(",
        labels,
        ");",
        "}",
    ]


def write_region_zone_helper(
    path: Path,
    *,
    fluid_cells: int,
    solid_cells: int,
    fluid_block_cells: int,
    solid_block_cells: int,
    fluid_zone_names: list[str],
    solid_zone_names: list[str],
) -> None:
    fluid_zones = [
        {"name": name, "start": i * fluid_block_cells, "stop": (i + 1) * fluid_block_cells}
        for i, name in enumerate(fluid_zone_names)
    ]
    solid_zones = [
        {"name": name, "start": i * solid_block_cells, "stop": (i + 1) * solid_block_cells}
        for i, name in enumerate(solid_zone_names)
    ]
    path.write_text(
        f"""#!/usr/bin/env python3
from pathlib import Path

fluid_cells = {fluid_cells}
solid_cells = {solid_cells}
fluid_zones = {fluid_zones!r}
solid_zones = {solid_zones!r}

def labels(start, stop):
    return "\\n".join(str(i) for i in range(start, stop))

def zone_block(name, start, stop):
    return [
        name,
        "{{",
        "    type cellZone;",
        "    cellLabels List<label>",
        str(stop - start),
        "(",
        labels(start, stop),
        ");",
        "}}",
    ]

def write_cell_zones(region, region_total_cells, zones):
    out = Path(f"constant/{{region}}/polyMesh/cellZones")
    names = [region] + [zone["name"] for zone in zones]
    parts = [
        "/*--------------------------------*- C++ -*----------------------------------*/",
        "FoamFile",
        "{{",
        "    version 2.0;",
        "    format ascii;",
        "    class regIOobject;",
        f"    location \\"constant/{{region}}/polyMesh\\";",
        "    object cellZones;",
        "    meta {{ names ( " + " ".join(names) + " ); }}",
        "}}",
        str(len(names)),
        "(",
    ]
    parts.extend(zone_block(region, 0, region_total_cells))
    for zone in zones:
        parts.extend(zone_block(zone["name"], zone["start"], zone["stop"]))
    parts.extend([")", ""])
    out.write_text("\\n".join(parts))
    print(f"wrote {{out}} with {{len(zones)}} mapped source zones")

write_cell_zones("bottomWater", fluid_cells, fluid_zones)
write_cell_zones("heater", solid_cells, solid_zones)
"""
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_fv_options(path: Path, *, region: str, prefix: str, zone_powers: list[tuple[str, float]]) -> None:
    entries = []
    for i, (zone_name, power_w) in enumerate(zone_powers):
        entries.append(
            f"""{prefix}_{i:02d}
{{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   cellZone;
    cellZone        {zone_name};
    volumeMode      absolute;

    sources
    {{
        h           ( {power_w:.12g} 0 );
    }}
}}
"""
        )
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

"""
        + "\n".join(entries)
    )


def patch_allpre(case: Path) -> None:
    allpre = case / "Allpre.c2cht"
    text = allpre.read_text()
    old = "splitMeshRegions -cellZones -overwrite > log.split 2>&1\ncheckMesh -allRegions > log.checkMesh 2>&1\n"
    new = (
        "splitMeshRegions -cellZones -overwrite > log.split 2>&1\n"
        "./write_c3_ffhr_depth_mapped_cell_zones.py > log.c3DepthMappedZones 2>&1\n"
        "checkMesh -allRegions > log.checkMesh 2>&1\n"
    )
    if old not in text:
        raise RuntimeError(f"could not patch {allpre}")
    allpre.write_text(text.replace(old, new))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="c3_curved_openmc_ffhr_be_stack_depth_mapped_cht")
    parser.add_argument("--radial-cells", type=int, default=50)
    parser.add_argument("--solid-cells", type=int, default=8)
    parser.add_argument("--wall-adjacent-cell", type=float, default=5e-5)
    parser.add_argument("--end-time", type=int, default=400)
    parser.add_argument("--write-interval", type=int, default=400)
    parser.add_argument("--target-heat-flux", type=float, default=100000.0)
    parser.add_argument(
        "--region-depth-profile",
        default="figs/c3_openmc_ffhr_be_stack_heating_split_region_depth.csv",
        help="OpenMC FFHR-like region-depth CSV with material_class and normalized_shape columns.",
    )
    parser.add_argument(
        "--openmc-summary",
        default="figs/c3_openmc_ffhr_be_stack_heating_split.json",
        help="OpenMC FFHR-like stack summary JSON to link in metadata.",
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
    geometry = json.loads((case / "geometry_c2cht.json").read_text())
    params = geometry["parameters"]
    segments = int(params["segments"])
    stream_cells = int(params["stream_cells"])
    height_cells = int(params["height_cells"])
    radial_cells = int(params["radial_cells"])
    solid_cells_per_block = int(params["solid_cells"])
    fluid_block_cells = stream_cells * height_cells * radial_cells
    solid_block_cells = stream_cells * height_cells * solid_cells_per_block
    fluid_cells = int(geometry["fluid_cells"])
    solid_cells = int(geometry["solid_cells"])
    if fluid_cells != segments * fluid_block_cells:
        raise ValueError("unexpected fluid-cell ordering; refusing to write mapped zones")
    if solid_cells != segments * solid_block_cells:
        raise ValueError("unexpected solid-cell ordering; refusing to write mapped zones")

    profile_path = work / args.region_depth_profile
    stack_thickness, profile = load_region_depth_profile(profile_path)
    fluid_fractions = map_class_to_segments(profile, material_class="fluid", segments=segments)
    solid_fractions = map_class_to_segments(profile, material_class="solid", segments=segments)

    heated_area = patch_area(geometry, "heatedBase")
    total_power = args.target_heat_flux * heated_area
    fluid_power_total = total_power * sum(fluid_fractions)
    solid_power_total = total_power * sum(solid_fractions)

    fluid_zone_names = [f"c3FFHRFluidSource_{i:02d}" for i in range(segments)]
    solid_zone_names = [f"c3FFHRSolidSource_{i:02d}" for i in range(segments)]
    fluid_zone_powers = [(name, total_power * fluid_fractions[i]) for i, name in enumerate(fluid_zone_names)]
    solid_zone_powers = [(name, total_power * solid_fractions[i]) for i, name in enumerate(solid_zone_names)]

    write_region_zone_helper(
        case / "write_c3_ffhr_depth_mapped_cell_zones.py",
        fluid_cells=fluid_cells,
        solid_cells=solid_cells,
        fluid_block_cells=fluid_block_cells,
        solid_block_cells=solid_block_cells,
        fluid_zone_names=fluid_zone_names,
        solid_zone_names=solid_zone_names,
    )
    write_fv_options(
        case / "system/bottomWater/fvOptions",
        region="bottomWater",
        prefix="c3FFHRFluidMappedHeating",
        zone_powers=fluid_zone_powers,
    )
    write_fv_options(
        case / "system/heater/fvOptions",
        region="heater",
        prefix="c3FFHRSolidMappedHeating",
        zone_powers=solid_zone_powers,
    )
    patch_allpre(case)

    summary_path = work / args.openmc_summary
    openmc_summary = json.loads(summary_path.read_text()) if summary_path.exists() else None
    mapping = []
    for i in range(segments):
        mapping.append(
            {
                "station": i,
                "fluid_zone": fluid_zone_names[i],
                "solid_zone": solid_zone_names[i],
                "fluid_fraction_of_total": fluid_fractions[i],
                "solid_fraction_of_total": solid_fractions[i],
                "fluid_power_W": fluid_zone_powers[i][1],
                "solid_power_W": solid_zone_powers[i][1],
            }
        )

    metadata = {
        "case": str(case),
        "basis": "C2 curved CHT geometry with external base flux disabled",
        "source_method": (
            "FFHR-like OpenMC region-depth heating resampled to streamwise "
            "fluid and solid cellZones; scalarSemiImplicitSource on h"
        ),
        "mapping_caveat": (
            "OpenMC stack depth is mapped onto streamwise stations as a local "
            "plumbing/sensitivity exercise, not as a physical wall-normal source map."
        ),
        "target_heat_flux_W_m2": args.target_heat_flux,
        "heated_base_area_m2": heated_area,
        "total_power_W": total_power,
        "mapped_power_sum_W": fluid_power_total + solid_power_total,
        "fluid_power_W": fluid_power_total,
        "solid_power_W": solid_power_total,
        "fluid_fraction": sum(fluid_fractions),
        "solid_fraction": sum(solid_fractions),
        "segments": segments,
        "fluid_block_cells_per_source_zone": fluid_block_cells,
        "solid_block_cells_per_source_zone": solid_block_cells,
        "stack_thickness_m": stack_thickness,
        "region_depth_profile": str(profile_path),
        "openmc_summary": str(summary_path) if summary_path.exists() else None,
        "openmc_max_depth_relative_uncertainty": (
            openmc_summary.get("max_depth_relative_uncertainty") if openmc_summary else None
        ),
        "fluid_source_power_min_W": min(power for _, power in fluid_zone_powers),
        "fluid_source_power_max_W": max(power for _, power in fluid_zone_powers),
        "solid_source_power_min_W": min(power for _, power in solid_zone_powers),
        "solid_source_power_max_W": max(power for _, power in solid_zone_powers),
        "mapping": mapping,
    }
    manifest_path = case / "c3_ffhr_depth_mapped_cht_manifest.json"
    manifest_path.write_text(json.dumps(metadata, indent=2))
    (case / "README.c3_ffhr_depth_mapped.md").write_text(
        "\n".join(
            [
                "# C3 Curved FFHR-Like Depth-Mapped CHT",
                "",
                "This case reuses the C2 curved CHT geometry, sets the external",
                "heated-base gradient to zero, and injects the same total power using",
                "mapped fluid and solid cell-zone source distributions from the",
                "FFHR-like OpenMC region-depth profile.",
                "",
                f"- total power: `{total_power:.8g} W`",
                f"- fluid power: `{fluid_power_total:.8g} W`",
                f"- solid power: `{solid_power_total:.8g} W`",
                f"- source stations per region: `{segments}`",
                f"- OpenMC profile: `{profile_path}`",
                f"- fluid source power range: `{metadata['fluid_source_power_min_W']:.6g}` to `{metadata['fluid_source_power_max_W']:.6g} W`",
                f"- solid source power range: `{metadata['solid_source_power_min_W']:.6g}` to `{metadata['solid_source_power_max_W']:.6g} W`",
                "",
                "Run:",
                "",
                "```bash",
                "./Allpre.c2cht",
                "./Allrun.c2cht",
                "```",
                "",
                "Caveat: this maps OpenMC stack depth onto streamwise mesh stations.",
                "It exercises profile-to-CHT plumbing and source-placement sensitivity;",
                "it is not a physical wall-normal source map or reactor neutronics.",
            ]
        )
        + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

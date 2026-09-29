#!/usr/bin/env python3
"""Create a C3 CHT case with OpenMC-profile-mapped volumetric heating.

This is the next step after the uniform C3 volumetric smoke case.  It reuses
the C2 curved CHT geometry, turns off external base heat flux, then maps the
OpenMC slab heating profile onto streamwise fluid cell zones.  The mapping is
still one-way and approximate, but it exercises the profile-to-CHT plumbing.
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


def load_openmc_profile(csv_path: Path) -> tuple[float, list[dict[str, float]]]:
    rows: list[dict[str, float]] = []
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "x0_m": float(row["x0_m"]),
                    "x1_m": float(row["x1_m"]),
                    "normalized_shape": float(row["normalized_shape"]),
                    "normalized_qvol_W_m3": float(row["normalized_qvol_W_m3"]),
                }
            )
    if not rows:
        raise ValueError(f"no OpenMC profile rows found in {csv_path}")
    thickness = max(row["x1_m"] for row in rows)
    total_shape = sum(row["normalized_shape"] for row in rows)
    if total_shape <= 0.0:
        raise ValueError(f"OpenMC profile shape is non-positive in {csv_path}")
    for row in rows:
        row["normalized_shape"] /= total_shape
        row["x0_hat"] = row["x0_m"] / thickness
        row["x1_hat"] = row["x1_m"] / thickness
    return thickness, rows


def map_profile_to_segments(profile: list[dict[str, float]], segments: int) -> list[float]:
    """Integrate equal-depth OpenMC bins onto curved-mesh stream segments."""
    fractions: list[float] = []
    for seg in range(segments):
        z0 = seg / segments
        z1 = (seg + 1) / segments
        frac = 0.0
        for row in profile:
            overlap = max(0.0, min(z1, row["x1_hat"]) - max(z0, row["x0_hat"]))
            bin_width = row["x1_hat"] - row["x0_hat"]
            if overlap > 0.0 and bin_width > 0.0:
                frac += row["normalized_shape"] * overlap / bin_width
        fractions.append(frac)

    total = sum(fractions)
    if total <= 0.0:
        raise ValueError("mapped source fractions sum to zero")
    return [frac / total for frac in fractions]


def labels_expr(start: int, stop: int) -> str:
    return f"labels({start}, {stop})"


def write_mapped_zone_helper(
    path: Path,
    *,
    fluid_cells: int,
    block_cells: int,
    zone_names: list[str],
) -> None:
    zones_literal = [
        {"name": name, "start": i * block_cells, "stop": (i + 1) * block_cells}
        for i, name in enumerate(zone_names)
    ]
    path.write_text(
        f"""#!/usr/bin/env python3
from pathlib import Path

fluid_cells = {fluid_cells}
zones = {zones_literal!r}
out = Path("constant/bottomWater/polyMesh/cellZones")

def labels(start, stop):
    return "\\n".join(str(i) for i in range(start, stop))

def zone_block(name, start, stop):
    count = stop - start
    return [
        name,
        "{{",
        "    type cellZone;",
        "    cellLabels List<label>",
        str(count),
        "(",
        labels(start, stop),
        ");",
        "}}",
    ]

names = ["bottomWater"] + [zone["name"] for zone in zones]
parts = [
    "/*--------------------------------*- C++ -*----------------------------------*/",
    "FoamFile",
    "{{",
    "    version 2.0;",
    "    format ascii;",
    "    class regIOobject;",
    "    location \\"constant/bottomWater/polyMesh\\";",
    "    object cellZones;",
    "    meta {{ names ( " + " ".join(names) + " ); }}",
    "}}",
    str(len(names)),
    "(",
]
parts.extend(zone_block("bottomWater", 0, fluid_cells))
for zone in zones:
    parts.extend(zone_block(zone["name"], zone["start"], zone["stop"]))
parts.extend([")", ""])
out.write_text("\\n".join(parts))
print(f"wrote {{out}} with {{len(zones)}} mapped C3 source zones")
"""
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_fv_options(path: Path, zone_powers: list[tuple[str, float]]) -> None:
    entries = []
    for i, (zone_name, power_w) in enumerate(zone_powers):
        entries.append(
            f"""c3OpenMCMappedHeating_{i:02d}
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
        """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2512/v2506                           |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system/bottomWater";
    object      fvOptions;
}

"""
        + "\n".join(entries)
    )


def patch_allpre(case: Path) -> None:
    allpre = case / "Allpre.c2cht"
    text = allpre.read_text()
    old = "splitMeshRegions -cellZones -overwrite > log.split 2>&1\ncheckMesh -allRegions > log.checkMesh 2>&1\n"
    new = (
        "splitMeshRegions -cellZones -overwrite > log.split 2>&1\n"
        "./write_c3_mapped_cell_zones.py > log.c3MappedZones 2>&1\n"
        "checkMesh -allRegions > log.checkMesh 2>&1\n"
    )
    if old not in text:
        raise RuntimeError(f"could not patch {allpre}")
    allpre.write_text(text.replace(old, new))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="c3_curved_openmc_mapped_volumetric_cht_smoke")
    parser.add_argument("--radial-cells", type=int, default=50)
    parser.add_argument("--solid-cells", type=int, default=8)
    parser.add_argument("--wall-adjacent-cell", type=float, default=5e-5)
    parser.add_argument("--end-time", type=int, default=400)
    parser.add_argument("--write-interval", type=int, default=400)
    parser.add_argument("--target-heat-flux", type=float, default=100000.0)
    parser.add_argument(
        "--openmc-profile",
        default="figs/c3_openmc_flibe_slab_heating_profile.csv",
        help="OpenMC depth-profile CSV with normalized_shape and normalized_qvol_W_m3 columns.",
    )
    parser.add_argument(
        "--openmc-summary",
        default="figs/c3_openmc_flibe_slab_heating_profile.json",
        help="OpenMC profile summary JSON to link in metadata.",
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
    block_cells = stream_cells * height_cells * radial_cells
    fluid_cells = int(geometry["fluid_cells"])
    if fluid_cells != segments * block_cells:
        raise ValueError("unexpected C2/C3 fluid-cell ordering; refusing to write mapped zones")

    profile_path = work / args.openmc_profile
    slab_thickness, profile = load_openmc_profile(profile_path)
    fractions = map_profile_to_segments(profile, segments)
    heated_area = patch_area(geometry, "heatedBase")
    total_power = args.target_heat_flux * heated_area

    zone_names = [f"c3SourceStation_{i:02d}" for i in range(segments)]
    zone_powers = [(name, total_power * fractions[i]) for i, name in enumerate(zone_names)]
    write_mapped_zone_helper(
        case / "write_c3_mapped_cell_zones.py",
        fluid_cells=fluid_cells,
        block_cells=block_cells,
        zone_names=zone_names,
    )
    write_fv_options(case / "system/bottomWater/fvOptions", zone_powers)
    patch_allpre(case)

    summary_path = work / args.openmc_summary
    openmc_summary = json.loads(summary_path.read_text()) if summary_path.exists() else None
    mapping = [
        {
            "zone": zone_names[i],
            "cell_start": i * block_cells,
            "cell_stop": (i + 1) * block_cells,
            "source_fraction": fractions[i],
            "source_power_W": zone_powers[i][1],
        }
        for i in range(segments)
    ]
    metadata = {
        "case": str(case),
        "basis": "C2 curved CHT geometry with external base flux disabled",
        "volumetric_region": "bottomWater",
        "source_method": "OpenMC profile resampled to streamwise cellZones; scalarSemiImplicitSource on h",
        "target_heat_flux_W_m2": args.target_heat_flux,
        "heated_base_area_m2": heated_area,
        "total_power_W": total_power,
        "mapped_power_sum_W": sum(power for _, power in zone_powers),
        "segments": segments,
        "stream_cells": stream_cells,
        "height_cells": height_cells,
        "radial_cells": radial_cells,
        "block_cells_per_source_zone": block_cells,
        "slab_thickness_m": slab_thickness,
        "openmc_profile": str(profile_path),
        "openmc_summary": str(summary_path) if summary_path.exists() else None,
        "openmc_source_kind": openmc_summary.get("source_kind") if openmc_summary else None,
        "openmc_profile_peak_to_mean": openmc_summary.get("normalized_qvol_peak_to_mean") if openmc_summary else None,
        "source_fraction_min": min(fractions),
        "source_fraction_max": max(fractions),
        "source_power_min_W": min(power for _, power in zone_powers),
        "source_power_max_W": max(power for _, power in zone_powers),
        "mapping": mapping,
    }
    (case / "c3_mapped_volumetric_cht_manifest.json").write_text(json.dumps(metadata, indent=2))
    (case / "README.c3mapped.md").write_text(
        "\n".join(
            [
                "# C3 Curved OpenMC-Mapped Volumetric CHT Smoke",
                "",
                "This case reuses the C2 curved CHT geometry, sets the external",
                "heated-base gradient to zero, and injects the same total power into",
                "streamwise FLiBe-like fluid cell zones using the OpenMC slab profile.",
                "",
                f"- target equivalent heat flux: `{args.target_heat_flux:g} W/m^2`",
                f"- heated-base reference area: `{heated_area:.8g} m^2`",
                f"- total power injected into fluid: `{total_power:.8g} W`",
                f"- source zones: `{segments}`",
                f"- OpenMC profile: `{profile_path}`",
                f"- mapped source power range: `{metadata['source_power_min_W']:.6g}` to `{metadata['source_power_max_W']:.6g} W`",
                "",
                "Run:",
                "",
                "```bash",
                "./Allpre.c2cht",
                "./Allrun.c2cht",
                "```",
                "",
                "This is still a one-way coupling smoke case. The OpenMC slab depth",
                "coordinate is resampled onto streamwise mesh stations as a local",
                "mapping exercise, not as a validated stellarator blanket source.",
            ]
        )
        + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

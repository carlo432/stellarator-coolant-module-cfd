#!/usr/bin/env python3
"""Postprocess a C2 curved-slice CHT smoke case."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


def latest_time_dir(case: Path) -> Path:
    time_dirs = []
    for path in case.iterdir():
        if path.is_dir():
            try:
                time_dirs.append((float(path.name), path))
            except ValueError:
                continue
    if not time_dirs:
        raise ValueError(f"{case}: no numeric time directories found")
    return max(time_dirs, key=lambda item: item[0])[1]


def parse_solver_log(case: Path) -> dict[str, float | str]:
    text = (case / "log.cht").read_text()
    result: dict[str, float | str] = {"log": "log.cht"}
    times = re.findall(r"^Time =\s+(\S+)", text, re.M)
    if times:
        result["final_time_or_iteration"] = float(times[-1])
    execs = re.findall(r"ExecutionTime =\s+([-+0-9.eE]+)\s+s", text)
    if execs:
        result["execution_time_s"] = float(execs[-1])
    cont = re.findall(
        r"time step continuity errors : sum local =\s+([-+0-9.eE]+), "
        r"global =\s+([-+0-9.eE]+), cumulative =\s+([-+0-9.eE]+)",
        text,
    )
    if cont:
        local, global_, cumulative = cont[-1]
        result["continuity_local_final"] = float(local)
        result["continuity_global_final"] = float(global_)
        result["continuity_cumulative_final"] = float(cumulative)
    h_matches = re.findall(
        r"Solving for h, Initial residual =\s+([-+0-9.eE]+), "
        r"Final residual =\s+([-+0-9.eE]+), No Iterations\s+(\d+)",
        text,
    )
    if len(h_matches) >= 2:
        fluid_h = h_matches[-2]
        solid_h = h_matches[-1]
        result["fluid_h_initial_residual_final_step"] = float(fluid_h[0])
        result["fluid_h_final_residual_final_step"] = float(fluid_h[1])
        result["solid_h_initial_residual_final_step"] = float(solid_h[0])
        result["solid_h_final_residual_final_step"] = float(solid_h[1])
    minmax = re.findall(r"Min/max T:([-+0-9.eE]+)\s+([-+0-9.eE]+)", text)
    if len(minmax) >= 2:
        fluid = minmax[-2]
        solid = minmax[-1]
        result["fluid_min_T_K"] = float(fluid[0])
        result["fluid_max_T_K"] = float(fluid[1])
        result["solid_min_T_K"] = float(solid[0])
        result["solid_max_T_K"] = float(solid[1])
    return result


def patch_average(case: Path, region: str, patch: str, field: str = "T") -> dict[str, float]:
    cmd = (
        "source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 "
        "|| source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1; "
        f"chtMultiRegionSimpleFoam -postProcess -region {region} "
        f"-func 'patchAverage(name={patch},fields=({field}))' -latestTime"
    )
    completed = subprocess.run(
        ["bash", "-lc", cmd],
        cwd=case,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    text = completed.stdout
    area = re.search(r"total area\s+=\s+([-+0-9.eE]+)", text)
    avg = re.search(rf"areaAverage\({re.escape(patch)}\) of {field} =\s+([-+0-9.eE]+)", text)
    if not area or not avg:
        raise ValueError(f"could not parse patchAverage output for {region}/{patch}/{field}")
    return {
        "area_m2": float(area.group(1)),
        "area_average": float(avg.group(1)),
    }


def parse_checkmesh(case: Path) -> dict[str, float | int]:
    text = (case / "log.checkMesh").read_text()
    cells = re.findall(r"cells:\s+(\d+)", text)
    aspect = re.findall(r"Max aspect ratio =\s+([-+0-9.eE]+)", text)
    nonorth = re.findall(r"Mesh non-orthogonality Max:\s+([-+0-9.eE]+)", text)
    skew = re.findall(r"Max skewness =\s+([-+0-9.eE]+)", text)
    return {
        "bottomWater_cells": int(cells[0]) if len(cells) > 0 else None,
        "heater_cells": int(cells[1]) if len(cells) > 1 else None,
        "bottomWater_max_aspect_ratio": float(aspect[0]) if len(aspect) > 0 else None,
        "heater_max_aspect_ratio": float(aspect[1]) if len(aspect) > 1 else None,
        "bottomWater_max_nonorthogonality": float(nonorth[0]) if len(nonorth) > 0 else None,
        "heater_max_nonorthogonality": float(nonorth[1]) if len(nonorth) > 1 else None,
        "bottomWater_max_skewness": float(skew[0]) if len(skew) > 0 else None,
        "heater_max_skewness": float(skew[1]) if len(skew) > 1 else None,
    }


def write_markdown(summary: dict, path: Path) -> None:
    tin = summary["tin_K"]
    patches = summary["patch_averages"]
    solver = summary["solver"]
    mesh = summary["mesh"]
    setup = summary.get("setup", {})
    case_name = str(summary["case"])
    if "ffhr_be_stack_depth_mapped" in case_name:
        title = "# C3 Curved OpenMC FFHR-Like Depth-Mapped CHT Summary"
        interpretation = (
            "This is a C3 one-way OpenMC-tallied FFHR-like Be/FLiBe stack CHT "
            "result with the region-depth profile resampled into separate "
            "streamwise fluid and solid source zones. It is a plumbing and "
            "source-placement sensitivity exercise: the OpenMC stack depth is "
            "mapped onto streamwise stations, not a validated physical wall-normal "
            "blanket source map."
        )
    elif "ffhr_be_stack" in case_name:
        title = "# C3 Curved OpenMC FFHR-Like Be-Stack CHT Summary"
        interpretation = (
            "This is a C3 one-way OpenMC-tallied FFHR-like Be/FLiBe stack CHT "
            "result. It reuses the simplified C2 curved geometry and applies "
            "aggregated fluid/solid source powers from a steel-proxy/front-FLiBe/"
            "Be-proxy/rear-FLiBe/steel-proxy OpenMC heating-local tally. It is "
            "a source-placement demonstration, not validated reactor neutronics "
            "or full blanket CHT."
        )
    elif "material_stack" in case_name:
        title = "# C3 Curved OpenMC Material-Stack CHT Summary"
        interpretation = (
            "This is a C3 one-way OpenMC-tallied material-stack CHT result. "
            "It reuses the simplified C2 curved geometry and applies separate "
            "fluid/solid source powers from a tiny steel-proxy/FLiBe/steel-proxy "
            "OpenMC heating-local tally. It is a source-placement demonstration, "
            "not validated reactor neutronics or full blanket CHT."
        )
    elif "material_split" in case_name:
        title = "# C3 Curved Material-Split CHT Summary"
        interpretation = (
            "This is a C3 one-way material-split volumetric-heating CHT result. "
            "It reuses the simplified C2 curved geometry and applies separate "
            "fluid/solid source powers from the Yamanishi-guided target card. It "
            "is a source-placement demonstration, not validated reactor neutronics "
            "or full blanket CHT."
        )
    elif "c3_curved_openmc_mapped" in case_name:
        title = "# C3 Curved OpenMC-Mapped Volumetric CHT Summary"
        interpretation = (
            "This is a short C3 one-way OpenMC-to-CHT mapping smoke result. It "
            "proves that an OpenMC depth profile can be resampled onto streamwise "
            "OpenFOAM fluid cell zones and used as a volumetric enthalpy source. "
            "It is not a validated neutronics model, reactor-geometry source map, "
            "or grid-converged thermal result."
        )
    elif "c3_curved_openmc" in case_name:
        title = "# C3 Curved Volumetric CHT Summary"
        interpretation = (
            "This is a short C3 one-way volumetric-heating CHT smoke result. It "
            "reuses the C2 curved CHT geometry, disables the external heated-base "
            "flux, and injects same-total power into the FLiBe-like fluid. It is "
            "not validated reactor neutronics or full blanket CHT."
        )
    elif "xfine" in case_name:
        title = "# C2 Curved CHT Extra-Fine-Grid Summary"
        interpretation = (
            "This is a C2 curved-slice wall-resolved CHT extra-fine-grid result. It "
            "supports a local three-level wall-normal sensitivity check for this "
            "simplified geometry, but it is not a formal GCI, validation, LES, or "
            "reactor-geometry CHT result."
        )
    elif "fine" in case_name:
        title = "# C2 Curved CHT Fine-Grid Summary"
        interpretation = (
            "This is a C2 curved-slice wall-resolved CHT fine-grid result. It "
            "supports a local two-grid sensitivity check for this simplified "
            "geometry, but it is not a formal GCI, validation, LES, or "
            "reactor-geometry CHT result."
        )
    elif "longrun" in case_name:
        title = "# C2 Curved CHT Longrun Summary"
        interpretation = (
            "This is a C2 curved-slice wall-resolved CHT continuation result. It "
            "gives a locally stable wall-temperature scale for this simplified "
            "geometry, but it is not a grid-converged, validated, LES, or "
            "reactor-geometry CHT result."
        )
    else:
        title = "# C2 Curved CHT Smoke Summary"
        interpretation = (
            "This is a short C2 curved-slice CHT smoke result. It proves the "
            "multi-region fluid-solid path and gives a first wall-temperature "
            "scale, but it is not a grid-converged, validated, or reactor-geometry "
            "CHT result."
        )
    lines = [
        title,
        "",
        f"- case: `{summary['case']}`",
        f"- latest time: `{summary['latest_time']}`",
        f"- fluid cells: `{mesh['bottomWater_cells']}`",
        f"- solid cells: `{mesh['heater_cells']}`",
        f"- radial grading: `{setup.get('radial_grading', 'n/a')}`",
        f"- requested wall-adjacent fluid cell: `{setup.get('parameters', {}).get('wall_adjacent_cell', 'n/a')}`",
        f"- solver log: `{solver.get('log', 'n/a')}`",
        f"- final time/iteration: `{solver.get('final_time_or_iteration', 'n/a')}`",
        f"- final execution time: `{solver.get('execution_time_s', 'n/a')}` s",
        f"- final continuity local: `{solver.get('continuity_local_final', 'n/a')}`",
        f"- final fluid h residual: `{solver.get('fluid_h_final_residual_final_step', 'n/a')}`",
        f"- final solid h residual: `{solver.get('solid_h_final_residual_final_step', 'n/a')}`",
        "",
        "## Patch Temperatures",
        "",
        "| Patch | Area, m^2 | Average T, K | Average dT over Tin, K |",
        "|---|---:|---:|---:|",
    ]
    for label, metric in patches.items():
        lines.append(
            f"| {label} | {metric['area_m2']:.8g} | {metric['area_average']:.5f} | "
            f"{metric['area_average'] - tin:.5f} |"
        )
    lines.extend(
        [
            "",
            "## Field Min/Max",
            "",
            f"- fluid T min/max: `{solver.get('fluid_min_T_K', 'n/a')}` / `{solver.get('fluid_max_T_K', 'n/a')}` K",
            f"- solid T min/max: `{solver.get('solid_min_T_K', 'n/a')}` / `{solver.get('solid_max_T_K', 'n/a')}` K",
            "",
            "## Mesh Checks",
            "",
            f"- bottomWater max aspect/nonorth/skew: `{mesh['bottomWater_max_aspect_ratio']}` / "
            f"`{mesh['bottomWater_max_nonorthogonality']}` / `{mesh['bottomWater_max_skewness']}`",
            f"- heater max aspect/nonorth/skew: `{mesh['heater_max_aspect_ratio']}` / "
            f"`{mesh['heater_max_nonorthogonality']}` / `{mesh['heater_max_skewness']}`",
            "",
            "## Interpretation",
            "",
            interpretation,
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="c2_curved_stellarator_slice_cht_wallres_smoke")
    ap.add_argument("--out-dir", default="figs")
    ap.add_argument("--tin", type=float, default=800.0)
    args = ap.parse_args()

    case = Path(args.case)
    latest = latest_time_dir(case)
    setup_path = case / "geometry_c2cht.json"
    setup = json.loads(setup_path.read_text()) if setup_path.exists() else {}
    summary = {
        "case": str(case),
        "latest_time": latest.name,
        "tin_K": args.tin,
        "setup": setup,
        "mesh": parse_checkmesh(case),
        "solver": parse_solver_log(case),
        "patch_averages": {
            "heater_to_bottomWater": patch_average(case, "heater", "heater_to_bottomWater"),
            "bottomWater_to_heater": patch_average(case, "bottomWater", "bottomWater_to_heater"),
            "heatedBase": patch_average(case, "heater", "heatedBase"),
        },
    }

    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = case.parent / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{case.name}_summary"
    json_path = out_dir / f"{stem}.json"
    md_path = out_dir / f"{stem}.md"
    json_path.write_text(json.dumps(summary, indent=2))
    write_markdown(summary, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

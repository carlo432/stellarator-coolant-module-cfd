#!/usr/bin/env python3
"""Compare baseline, fine, and extra-fine C2 CHT wall-resolution cases."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def patch_dt(data: dict, patch: str) -> float:
    return data["patch_averages"][patch]["area_average"] - data["tin_K"]


def metric_value(data: dict, name: str) -> float:
    if name == "interface_mean_dT_K":
        return patch_dt(data, "heater_to_bottomWater")
    if name == "heatedBase_mean_dT_K":
        return patch_dt(data, "heatedBase")
    if name in data["solver"]:
        return data["solver"][name]
    raise KeyError(name)


def pct(new: float, old: float) -> float:
    return 100.0 * (new - old) / old if old else 0.0


def span_pct(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return 100.0 * (max(values) - min(values)) / mean if mean else 0.0


def mesh_info(data: dict) -> dict:
    params = data.get("setup", {}).get("parameters", {})
    return {
        "case": data["case"],
        "latest_time": data["latest_time"],
        "fluid_cells": data["mesh"]["bottomWater_cells"],
        "solid_cells": data["mesh"]["heater_cells"],
        "radial_cells": params.get("radial_cells"),
        "solid_cells_through": params.get("solid_cells"),
        "wall_adjacent_cell": params.get("wall_adjacent_cell"),
        "radial_grading": data.get("setup", {}).get("radial_grading"),
    }


def write_markdown(summary: dict, path: Path) -> None:
    cases = summary["cases"]
    rows = summary["metrics"]
    lines = [
        "# C2 Curved CHT Three-Level Wall-Normal Sensitivity",
        "",
        "This compares the baseline, fine, and extra-fine wall-resolved C2 CHT",
        "cases. The refinement is concentrated in radial fluid resolution and solid",
        "thickness resolution, so this is not a formal all-direction GCI.",
        "",
        "| Quantity | Baseline | Fine | Extra-fine |",
        "|---|---:|---:|---:|",
        (
            f"| Fluid cells | {cases['baseline']['fluid_cells']} | "
            f"{cases['fine']['fluid_cells']} | {cases['xfine']['fluid_cells']} |"
        ),
        (
            f"| Solid cells | {cases['baseline']['solid_cells']} | "
            f"{cases['fine']['solid_cells']} | {cases['xfine']['solid_cells']} |"
        ),
        (
            f"| Radial cells | {cases['baseline']['radial_cells']} | "
            f"{cases['fine']['radial_cells']} | {cases['xfine']['radial_cells']} |"
        ),
        (
            f"| Solid-thickness cells | {cases['baseline']['solid_cells_through']} | "
            f"{cases['fine']['solid_cells_through']} | {cases['xfine']['solid_cells_through']} |"
        ),
    ]
    for row in rows.values():
        lines.append(
            f"| {row['label']} | {row['baseline']:.6g} | {row['fine']:.6g} | "
            f"{row['xfine']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Changes",
            "",
            "| Metric | Fine - baseline | Extra-fine - fine | Three-level span |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in rows.values():
        lines.append(
            f"| {row['label']} | {row['fine_minus_baseline']:+.6g} "
            f"({row['fine_minus_baseline_pct']:+.3f}%) | "
            f"{row['xfine_minus_fine']:+.6g} "
            f"({row['xfine_minus_fine_pct']:+.3f}%) | "
            f"{row['span_pct']:.3f}% |"
        )
    lines.extend(["", "## Interpretation", "", summary["interpretation"], ""])
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", default="figs/c2_curved_stellarator_slice_cht_wallres_longrun_summary.json")
    ap.add_argument("--fine", default="figs/c2_curved_stellarator_slice_cht_wallres_fine_summary.json")
    ap.add_argument("--xfine", default="figs/c2_curved_stellarator_slice_cht_wallres_xfine_summary.json")
    ap.add_argument("--out", default="figs/c2_curved_cht_three_grid_sensitivity")
    args = ap.parse_args()

    baseline = load(Path(args.baseline))
    fine = load(Path(args.fine))
    xfine = load(Path(args.xfine))
    metric_specs = {
        "interface_mean_dT_K": "Interface mean dT, K",
        "heatedBase_mean_dT_K": "Heated-base mean dT, K",
        "fluid_max_T_K": "OpenFOAM log fluid max T, K",
        "solid_max_T_K": "OpenFOAM log solid max T, K",
    }
    metrics = {}
    for key, label in metric_specs.items():
        values = [
            metric_value(baseline, key),
            metric_value(fine, key),
            metric_value(xfine, key),
        ]
        metrics[key] = {
            "label": label,
            "baseline": values[0],
            "fine": values[1],
            "xfine": values[2],
            "fine_minus_baseline": values[1] - values[0],
            "fine_minus_baseline_pct": pct(values[1], values[0]),
            "xfine_minus_fine": values[2] - values[1],
            "xfine_minus_fine_pct": pct(values[2], values[1]),
            "span_pct": span_pct(values),
        }

    interface_span = metrics["interface_mean_dT_K"]["span_pct"]
    base_span = metrics["heatedBase_mean_dT_K"]["span_pct"]
    if interface_span <= 1.0 and base_span <= 1.0:
        interpretation = (
            "All three wall-resolved C2 grids keep interface and heated-base mean "
            "dT within a 1% span. This is strong local wall-normal mesh-sensitivity "
            "evidence for the simplified C2 geometry. It remains not-a-formal-GCI "
            "because streamwise/height resolution and geometry are not refined as a "
            "matched all-direction sequence, and it is not validation or LES."
        )
    else:
        interpretation = (
            "The three-grid span exceeds 1% for at least one wall-temperature metric. "
            "Keep the mesh-sensitivity caveat and avoid quoting this as a stable "
            "wall-temperature scale without further refinement."
        )

    summary = {
        "cases": {
            "baseline": mesh_info(baseline),
            "fine": mesh_info(fine),
            "xfine": mesh_info(xfine),
        },
        "metrics": metrics,
        "interpretation": interpretation,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json_path = out.with_suffix(".json")
    md_path = out.with_suffix(".md")
    json_path.write_text(json.dumps(summary, indent=2))
    write_markdown(summary, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

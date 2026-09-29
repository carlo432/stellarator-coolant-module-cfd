#!/usr/bin/env python3
"""Compare baseline and fine wall-resolved C2 CHT longrun summaries."""
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


def percent_change(new: float, old: float) -> float:
    return 100.0 * (new - old) / old if old else 0.0


def write_markdown(summary: dict, path: Path) -> None:
    coarse = summary["coarse"]
    fine = summary["fine"]
    rows = summary["metrics"]
    lines = [
        "# C2 Curved CHT Grid Sensitivity",
        "",
        "This compares the locally stable wall-resolved C2 CHT longrun against a finer",
        "same-physics grid. It is a two-grid sensitivity check, not a formal GCI.",
        "",
        "| Quantity | Baseline wall-res | Fine wall-res | Fine - baseline |",
        "|---|---:|---:|---:|",
        (
            f"| Fluid cells | {coarse['fluid_cells']} | {fine['fluid_cells']} | "
            f"{fine['fluid_cells'] - coarse['fluid_cells']} |"
        ),
        (
            f"| Solid cells | {coarse['solid_cells']} | {fine['solid_cells']} | "
            f"{fine['solid_cells'] - coarse['solid_cells']} |"
        ),
        (
            f"| Radial cells | {coarse['radial_cells']} | {fine['radial_cells']} | "
            f"{fine['radial_cells'] - coarse['radial_cells']} |"
        ),
        (
            f"| Solid-thickness cells | {coarse['solid_cells_through']} | "
            f"{fine['solid_cells_through']} | "
            f"{fine['solid_cells_through'] - coarse['solid_cells_through']} |"
        ),
    ]
    for label, row in rows.items():
        lines.append(
            f"| {row['label']} | {row['coarse']:.6g} | {row['fine']:.6g} | "
            f"{row['delta']:+.6g} ({row['delta_percent']:+.3f}%) |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            summary["interpretation"],
            "",
        ]
    )
    path.write_text("\n".join(lines))


def case_mesh_info(data: dict) -> dict[str, int | float | str | None]:
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--coarse",
        default="figs/c2_curved_stellarator_slice_cht_wallres_longrun_summary.json",
        help="Baseline wall-resolved longrun summary JSON.",
    )
    ap.add_argument(
        "--fine",
        default="figs/c2_curved_stellarator_slice_cht_wallres_fine_summary.json",
        help="Fine-grid wall-resolved summary JSON.",
    )
    ap.add_argument("--out", default="figs/c2_curved_cht_grid_sensitivity")
    args = ap.parse_args()

    coarse = load(Path(args.coarse))
    fine = load(Path(args.fine))
    metric_specs = {
        "interface_mean_dT_K": "Interface mean dT, K",
        "heatedBase_mean_dT_K": "Heated-base mean dT, K",
        "fluid_max_T_K": "OpenFOAM log fluid max T, K",
        "solid_max_T_K": "OpenFOAM log solid max T, K",
        "continuity_local_final": "Final local continuity",
        "fluid_h_final_residual_final_step": "Final fluid h residual",
        "solid_h_final_residual_final_step": "Final solid h residual",
    }
    metrics = {}
    for name, label in metric_specs.items():
        coarse_value = metric_value(coarse, name)
        fine_value = metric_value(fine, name)
        metrics[name] = {
            "label": label,
            "coarse": coarse_value,
            "fine": fine_value,
            "delta": fine_value - coarse_value,
            "delta_percent": percent_change(fine_value, coarse_value),
        }

    interface_pct = abs(metrics["interface_mean_dT_K"]["delta_percent"])
    base_pct = abs(metrics["heatedBase_mean_dT_K"]["delta_percent"])
    if interface_pct <= 2.0 and base_pct <= 2.0:
        interpretation = (
            "The fine grid changes both the interface and heated-base mean dT by "
            "<= 2%, so the C2 wall-temperature scale is locally mesh-stable across "
            "this two-grid check. This is stronger than a smoke result, but still "
            "not a formal GCI, validation, LES, or reactor-geometry calculation."
        )
    else:
        interpretation = (
            "The fine grid changes at least one wall-temperature metric by > 2%, so "
            "C2 should keep a grid-sensitivity caveat and consider a third grid before "
            "quoting a report-grade wall-temperature scale."
        )

    summary = {
        "coarse": case_mesh_info(coarse),
        "fine": case_mesh_info(fine),
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

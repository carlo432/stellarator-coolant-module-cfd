#!/usr/bin/env python3
"""Compare coarse and wall-resolved C2 curved CHT smoke summaries."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def dT(data: dict, patch: str) -> float:
    return data["patch_averages"][patch]["area_average"] - data["tin_K"]


def pct(new: float, ref: float) -> float:
    return 100.0 * (new - ref) / ref


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coarse", default="figs/c2_curved_stellarator_slice_cht_preflight_summary.json")
    ap.add_argument("--wallres", default="figs/c2_curved_stellarator_slice_cht_wallres_smoke_summary.json")
    ap.add_argument("--out", default="figs/c2_curved_cht_resolution_comparison")
    args = ap.parse_args()

    coarse = load(Path(args.coarse))
    wallres = load(Path(args.wallres))
    metrics = {
        "interface_dT_K": {
            "coarse": dT(coarse, "heater_to_bottomWater"),
            "wallres": dT(wallres, "heater_to_bottomWater"),
        },
        "heatedBase_dT_K": {
            "coarse": dT(coarse, "heatedBase"),
            "wallres": dT(wallres, "heatedBase"),
        },
        "fluid_max_T_K": {
            "coarse": coarse["solver"]["fluid_max_T_K"],
            "wallres": wallres["solver"]["fluid_max_T_K"],
        },
        "solid_max_T_K": {
            "coarse": coarse["solver"]["solid_max_T_K"],
            "wallres": wallres["solver"]["solid_max_T_K"],
        },
    }
    for metric in metrics.values():
        metric["delta"] = metric["wallres"] - metric["coarse"]
        metric["delta_percent"] = pct(metric["wallres"], metric["coarse"])

    summary = {
        "coarse_case": coarse["case"],
        "wallres_case": wallres["case"],
        "coarse_cells": coarse["mesh"],
        "wallres_cells": wallres["mesh"],
        "metrics": metrics,
        "interpretation": (
            "The coarse CHT run is a plumbing proof only. The wall-resolved radial "
            "grading cuts the interface/base superheat scale by roughly 80%, so "
            "C2 CHT claims should use the wall-resolved smoke as the first credible "
            "temperature-scale artifact, still with short-run/non-converged caveats."
        ),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json_path = out.with_suffix(".json")
    md_path = out.with_suffix(".md")
    json_path.write_text(json.dumps(summary, indent=2))
    md_path.write_text(
        "# C2 Curved CHT Resolution Comparison\n\n"
        "This compares the first coarse CHT plumbing run against the radially graded "
        "wall-resolved CHT smoke run on the same curved C2 slice.\n\n"
        "| Metric | Coarse | Wall-resolved | Wall-resolved - Coarse |\n"
        "|---|---:|---:|---:|\n"
        f"| Fluid cells | {coarse['mesh']['bottomWater_cells']} | {wallres['mesh']['bottomWater_cells']} | {wallres['mesh']['bottomWater_cells'] - coarse['mesh']['bottomWater_cells']} |\n"
        f"| Solid cells | {coarse['mesh']['heater_cells']} | {wallres['mesh']['heater_cells']} | {wallres['mesh']['heater_cells'] - coarse['mesh']['heater_cells']} |\n"
        f"| Interface mean dT, K | {metrics['interface_dT_K']['coarse']:.6g} | {metrics['interface_dT_K']['wallres']:.6g} | {metrics['interface_dT_K']['delta']:.6g} ({metrics['interface_dT_K']['delta_percent']:+.2f}%) |\n"
        f"| Heated-base mean dT, K | {metrics['heatedBase_dT_K']['coarse']:.6g} | {metrics['heatedBase_dT_K']['wallres']:.6g} | {metrics['heatedBase_dT_K']['delta']:.6g} ({metrics['heatedBase_dT_K']['delta_percent']:+.2f}%) |\n"
        f"| Fluid max T, K | {metrics['fluid_max_T_K']['coarse']:.6g} | {metrics['fluid_max_T_K']['wallres']:.6g} | {metrics['fluid_max_T_K']['delta']:.6g} |\n"
        f"| Solid max T, K | {metrics['solid_max_T_K']['coarse']:.6g} | {metrics['solid_max_T_K']['wallres']:.6g} | {metrics['solid_max_T_K']['delta']:.6g} |\n\n"
        "Interpretation: the coarse CHT run is a plumbing proof only. The "
        "wall-resolved radial grading cuts the apparent wall-temperature scale by "
        "about 80%, so C2 CHT claims should use the wall-resolved smoke as the first "
        "credible temperature-scale artifact, still with short-run and non-converged "
        "caveats.\n"
    )
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

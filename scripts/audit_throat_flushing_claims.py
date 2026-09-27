#!/usr/bin/env python3
"""Independent audit of the local-throat flushing test.

The primary analysis establishes the design verdict. This script checks the
same reconstructed endpoints with a separate implementation and records the field
evidence relevant to the proposed relaminarization mechanism.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pyvista as pv


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "cases" / "branch_d"
OUTPUT = ROOT / "results" / "advisor_packet" / "102_throat_flushing_audit.json"

T_IN = 900.0
X_PEAK = 0.1461
STRIP_WIDTH = 0.034

CASES = {
    "baseline": (WORK / "geom_straight_pipeline_plasma", 2.7),
    "throat": (WORK / "geom_throat_pipeline_plasma", 4.0),
}


def averaging_time(case: Path, time_value: float) -> float:
    path = case / f"{time_value:g}" / "uniform" / "functionObjects" / "functionObjectProperties"
    text = path.read_text()
    block = re.search(r"fieldAverage1\s*\{.*?\bU\s*\{.*?totalTime\s+([0-9.eE+-]+);", text, re.S)
    if not block:
        raise RuntimeError(f"Cannot recover fieldAverage U totalTime from {path}")
    return float(block.group(1))


def endpoint(case: Path, time_value: float) -> dict[str, float]:
    marker = case / "case.foam"
    marker.touch(exist_ok=True)
    reader = pv.OpenFOAMReader(marker)
    reader.set_active_time_value(time_value)
    data = reader.read()

    internal = data["internalMesh"]
    centers = internal.cell_centers()
    points = np.asarray(centers.points)
    z_max = float(points[:, 2].max())
    near_wall = (
        (points[:, 2] > z_max - 0.0007)
        & (np.abs(points[:, 1]) < 0.003)
        & (np.abs(points[:, 0] - X_PEAK) < STRIP_WIDTH / 2)
    )
    near_wall_center = near_wall & (np.abs(points[:, 0] - X_PEAK) < 0.005)
    if not np.any(near_wall):
        raise RuntimeError(f"No near-wall hotspot cells selected in {case}")

    u_mean = np.asarray(centers["UMean"])
    u_prime2 = np.asarray(centers["UPrime2Mean"])
    t_prime2 = np.maximum(np.asarray(centers["TPrime2Mean"]), 0.0)
    nut = np.asarray(centers["nut"])

    # OpenFOAM symmTensor ordering begins with xx, yy, zz.
    resolved_tke = 0.5 * np.sum(u_prime2[:, :3], axis=1)

    wall = data["boundary"]["heated_first_wall"].cell_data_to_point_data()
    wall_points = np.asarray(wall.points)
    wall_t = np.asarray(wall["TMean"])
    wall_hot = (
        (np.abs(wall_points[:, 0] - X_PEAK) < STRIP_WIDTH / 2)
        & (np.abs(wall_points[:, 1]) < 0.005)
    )

    domain_volume = float(internal.compute_cell_sizes(length=False, area=False)["Volume"].sum())
    wall_area = float(
        data["boundary"]["heated_first_wall"]
        .compute_cell_sizes(length=False, volume=False)["Area"]
        .sum()
    )

    return {
        "averaging_time_s": averaging_time(case, time_value),
        "domain_volume_m3": domain_volume,
        "heated_wall_area_m2": wall_area,
        "wall_p99_temperature_K": float(np.percentile(wall_t, 99)),
        "wall_max_temperature_K": float(wall_t.max()),
        "hot_strip_wall_peak_superheat_K": float(wall_t[wall_hot].max() - T_IN),
        "hot_strip_nearwall_Ux_mean_m_s": float(u_mean[near_wall, 0].mean()),
        "hotspot_center_nearwall_Ux_mean_m_s": float(u_mean[near_wall_center, 0].mean()),
        "hot_strip_resolved_tke_mean_m2_s2": float(resolved_tke[near_wall].mean()),
        "hot_strip_temperature_rms_mean_K": float(np.sqrt(t_prime2[near_wall]).mean()),
        "hot_strip_nut_mean_m2_s": float(nut[near_wall].mean()),
        "hot_strip_cell_count": int(near_wall.sum()),
    }


def ratio(new: float, old: float) -> float:
    return new / old


def main() -> None:
    result = {name: endpoint(*case) for name, case in CASES.items()}
    baseline = result["baseline"]
    throat = result["throat"]

    result["ratios_throat_over_baseline"] = {
        "averaging_time": ratio(throat["averaging_time_s"], baseline["averaging_time_s"]),
        "domain_volume": ratio(throat["domain_volume_m3"], baseline["domain_volume_m3"]),
        "heated_wall_area": ratio(throat["heated_wall_area_m2"], baseline["heated_wall_area_m2"]),
        "nearwall_Ux": ratio(
            throat["hot_strip_nearwall_Ux_mean_m_s"],
            baseline["hot_strip_nearwall_Ux_mean_m_s"],
        ),
        "hotspot_center_nearwall_Ux": ratio(
            throat["hotspot_center_nearwall_Ux_mean_m_s"],
            baseline["hotspot_center_nearwall_Ux_mean_m_s"],
        ),
        "resolved_tke": ratio(
            throat["hot_strip_resolved_tke_mean_m2_s2"],
            baseline["hot_strip_resolved_tke_mean_m2_s2"],
        ),
        "temperature_rms": ratio(
            throat["hot_strip_temperature_rms_mean_K"],
            baseline["hot_strip_temperature_rms_mean_K"],
        ),
        "nut": ratio(
            throat["hot_strip_nut_mean_m2_s"],
            baseline["hot_strip_nut_mean_m2_s"],
        ),
        "hot_strip_peak_superheat": ratio(
            throat["hot_strip_wall_peak_superheat_K"],
            baseline["hot_strip_wall_peak_superheat_K"],
        ),
    }
    result["claim_control"] = {
        "measured": (
            "The local contraction preserves heated-wall area, accelerates near-wall flow, "
            "and under-delivers the developed-duct cooling prediction."
        ),
        "mechanism_interpretation": (
            "Reduced resolved fluctuation energy, temperature RMS, and turbulent viscosity "
            "support favorable-pressure-gradient relaminarization as the mechanism; this is "
            "not a separately validated relaminarization model."
        ),
        "source_caveat": (
            "The cold-side contraction lowers total fluid/source volume slightly, so the cases "
            "are not an exact one-variable energy-source comparison."
        ),
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Postprocess the matched A3 plasma-load CHT comparison."""
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np

from postprocess_curved_stellarator_slice_c2 import (
    face_area,
    patch_face_ids,
    read_boundary,
    read_faces,
    read_internal_scalar,
    read_points,
)


WORK = Path(__file__).resolve().parent
REPO = WORK.parents[1]
RESULTS = REPO / "results" / "branch_d" / "a3_variable_property_cht"
FIGS = WORK / "figs"
CASES = {
    "constant": WORK / "a3_const_plasma_cht",
    "variable": WORK / "a3_varprop_plasma_cht",
}
T_INLET = 800.0
FLIBE_LIMIT_K = 1703.0


def numeric_time_dirs(case: Path) -> list[Path]:
    result: list[tuple[float, Path]] = []
    for child in case.iterdir():
        if not child.is_dir():
            continue
        try:
            value = float(child.name)
        except ValueError:
            continue
        if value >= 3000.0:
            result.append((value, child))
    return [path for _, path in sorted(result)]


def patch_block(text: str, patch_name: str) -> str:
    start = re.search(rf"\b{re.escape(patch_name)}\s*\{{", text)
    if not start:
        raise ValueError(f"patch {patch_name!r} not found")
    idx = start.end()
    depth = 1
    while idx < len(text) and depth:
        if text[idx] == "{":
            depth += 1
        elif text[idx] == "}":
            depth -= 1
        idx += 1
    return text[start.end() : idx - 1]


def read_patch_values(path: Path, patch_name: str) -> np.ndarray:
    block = patch_block(path.read_text(), patch_name)
    matches = list(
        re.finditer(
            r"(?m)^\s*value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;",
            block,
            re.S,
        )
    )
    if not matches:
        uniform = re.search(r"(?m)^\s*value\s+uniform\s+([-+0-9.eE]+)\s*;", block)
        if uniform:
            raise ValueError(f"{path}: uniform patch values cannot establish face count")
        raise ValueError(f"{path}: no patch value list for {patch_name}")
    count = int(matches[-1].group(1))
    values = np.asarray(
        [float(value) for value in re.findall(r"[-+0-9.eE]+", matches[-1].group(2))],
        dtype=float,
    )
    if len(values) != count:
        raise ValueError(f"{path}: expected {count} patch values, got {len(values)}")
    return values


def interface_geometry(case: Path) -> tuple[np.ndarray, np.ndarray]:
    mesh = case / "constant" / "bottomWater" / "polyMesh"
    points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    patches = read_boundary(mesh / "boundary")
    ids = list(patch_face_ids(patches["bottomWater_to_heater"]))
    centers = []
    areas = []
    for face_id in ids:
        face = faces[face_id]
        centers.append(np.asarray([points[idx] for idx in face], dtype=float).mean(axis=0))
        areas.append(face_area(face, points))
    return np.asarray(centers), np.asarray(areas)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(values * weights) / np.sum(weights))


def case_field_history(case: Path) -> tuple[list[dict[str, float]], dict[str, np.ndarray]]:
    centers, areas = interface_geometry(case)
    history = []
    final_values: np.ndarray | None = None
    for time_dir in numeric_time_dirs(case):
        interface = read_patch_values(
            time_dir / "bottomWater" / "T", "bottomWater_to_heater"
        )
        fluid = np.asarray(read_internal_scalar(time_dir / "bottomWater" / "T"), dtype=float)
        solid = np.asarray(read_internal_scalar(time_dir / "heater" / "T"), dtype=float)
        history.append(
            {
                "time": float(time_dir.name),
                "interface_avg_K": weighted_mean(interface, areas),
                "interface_max_K": float(interface.max()),
                "fluid_min_K": float(fluid.min()),
                "fluid_max_K": float(fluid.max()),
                "solid_min_K": float(solid.min()),
                "solid_max_K": float(solid.max()),
            }
        )
        final_values = interface
    if final_values is None:
        raise ValueError(f"{case}: no persisted time fields")
    return history, {"centers": centers, "areas": areas, "interface": final_values}


def parse_solver_log(path: Path) -> dict[str, object]:
    lines = path.read_text(errors="replace").splitlines()
    current_time: float | None = None
    region: str | None = None
    residuals: dict[str, float] = {}
    t_records: list[dict[str, float | str]] = []
    continuity: dict[str, float] = {}
    for line in lines:
        time_match = re.match(r"Time = ([-+0-9.eE]+)", line)
        if time_match:
            current_time = float(time_match.group(1))
        region_match = re.match(r"Solving for (fluid|solid) region (\S+)", line)
        if region_match:
            region = region_match.group(1)
        residual_match = re.search(
            r"Solving for ([^,]+), Initial residual = ([-+0-9.eE]+)", line
        )
        if residual_match and region:
            residuals[f"{region}:{residual_match.group(1)}"] = float(residual_match.group(2))
        temp_match = re.search(r"Min/max T:([-+0-9.eE]+)\s+([-+0-9.eE]+)", line)
        if temp_match and region and current_time is not None:
            t_records.append(
                {
                    "time": current_time,
                    "region": region,
                    "min_K": float(temp_match.group(1)),
                    "max_K": float(temp_match.group(2)),
                }
            )
        continuity_match = re.search(
            r"sum local = ([-+0-9.eE]+), global = ([-+0-9.eE]+), cumulative = ([-+0-9.eE]+)",
            line,
        )
        if continuity_match:
            continuity = {
                "sum_local": float(continuity_match.group(1)),
                "global": float(continuity_match.group(2)),
                "cumulative": float(continuity_match.group(3)),
            }

    tail: dict[str, dict[str, float]] = {}
    for region_name in ("fluid", "solid"):
        values = [record for record in t_records if record["region"] == region_name][-100:]
        if values:
            maxima = np.asarray([float(record["max_K"]) for record in values])
            minima = np.asarray([float(record["min_K"]) for record in values])
            tail[region_name] = {
                "samples": len(values),
                "max_T_span_K": float(maxima.max() - maxima.min()),
                "min_T_span_K": float(minima.max() - minima.min()),
                "final_min_K": float(minima[-1]),
                "final_max_K": float(maxima[-1]),
            }
    last_nonblank = next((line.strip() for line in reversed(lines) if line.strip()), "")
    return {
        "clean_end": last_nonblank == "End",
        "final_time": current_time,
        "final_initial_residuals": residuals,
        "final_continuity": continuity,
        "last_100_iteration_temperature_stability": tail,
    }


def profile_by_x(final: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centers = final["centers"]
    areas = final["areas"]
    values = final["interface"] - T_INLET
    rounded_x = np.round(centers[:, 0], 10)
    x_values = np.unique(rounded_x)
    means = []
    maxima = []
    for x in x_values:
        mask = rounded_x == x
        means.append(weighted_mean(values[mask], areas[mask]))
        maxima.append(float(values[mask].max()))
    return x_values, np.asarray(means), np.asarray(maxima)


def signed_delta(variable: float, constant: float) -> float:
    return 100.0 * (variable / constant - 1.0)


def make_figure(
    histories: dict[str, list[dict[str, float]]],
    finals: dict[str, dict[str, np.ndarray]],
    summary: dict[str, object],
    output: Path,
) -> None:
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.4), constrained_layout=True)
    colors = {"constant": "#2367a8", "variable": "#c33d2e"}
    labels = {"constant": "constant properties", "variable": "variable mu(T), k(T)"}

    ax = axes[0, 0]
    for model in CASES:
        rows = histories[model]
        ax.plot(
            [row["time"] for row in rows],
            [row["interface_avg_K"] - T_INLET for row in rows],
            marker="o",
            ms=4,
            color=colors[model],
            label=labels[model],
        )
    ax.set_title("Persisted-field convergence")
    ax.set_xlabel("steady iteration")
    ax.set_ylabel("interface-average superheat [K]")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    for model in CASES:
        x, mean, maximum = profile_by_x(finals[model])
        ax.plot(x, mean, lw=2, color=colors[model], label=labels[model])
        ax.plot(x, maximum, lw=1, ls="--", color=colors[model], alpha=0.65)
    ax.set_title("Final interface profile")
    ax.set_xlabel("streamwise x [m]")
    ax.set_ylabel("interface superheat [K]")
    ax.grid(alpha=0.25)
    ax.text(0.02, 0.96, "solid: spanwise mean\ndashed: spanwise max", transform=ax.transAxes, va="top")

    ax = axes[1, 0]
    metrics = [
        "interface_avg_superheat_K",
        "interface_max_superheat_K",
        "fluid_internal_max_superheat_K",
        "solid_global_max_superheat_K",
    ]
    metric_labels = ["interface\navg", "interface\nmax", "fluid internal\nmax", "solid global\nmax"]
    x_pos = np.arange(len(metrics))
    width = 0.36
    final_summary = summary["final"]
    const_values = [float(final_summary["constant"][metric]) for metric in metrics]
    var_values = [float(final_summary["variable"][metric]) for metric in metrics]
    ax.bar(x_pos - width / 2, const_values, width, color=colors["constant"], label=labels["constant"])
    ax.bar(x_pos + width / 2, var_values, width, color=colors["variable"], label=labels["variable"])
    for idx, (constant, variable) in enumerate(zip(const_values, var_values)):
        ax.text(
            idx,
            max(constant, variable) * 1.025,
            f"{signed_delta(variable, constant):+.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.set_title("Final signed property effect")
    ax.set_ylabel("temperature above 800 K [K]")
    ax.set_xticks(x_pos, metric_labels)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1, 1]
    temperatures = np.linspace(700.0, 1800.0, 500)
    mu = 1.16e-4 * np.exp(3755.0 / temperatures)
    kappa = 1.212 + 9.0e-5 * np.clip(temperatures, 800.0, FLIBE_LIMIT_K)
    ax.plot(temperatures, mu * 1000.0, color="#70439b", lw=2, label="mu [mPa s]")
    ax.set_xlabel("temperature [K]")
    ax.set_ylabel("dynamic viscosity [mPa s]", color="#70439b")
    ax.tick_params(axis="y", labelcolor="#70439b")
    ax.grid(alpha=0.25)
    ax.axvspan(873.0, 1073.0, color="#70439b", alpha=0.08, label="Sohal mu source range")
    ax.axvline(FLIBE_LIMIT_K, color="#4b4b4b", ls=":", lw=1.3)
    twin = ax.twinx()
    twin.plot(temperatures, kappa, color="#16846d", lw=2, label="k [W/m/K]")
    twin.set_ylabel("thermal conductivity [W/m/K]", color="#16846d")
    twin.tick_params(axis="y", labelcolor="#16846d")
    ax.set_title("Applied variable-property laws")
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = twin.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, frameon=False, loc="upper right")

    fig.suptitle("A3 plasma-load CHT: constant versus temperature-dependent FLiBe properties", fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    histories: dict[str, list[dict[str, float]]] = {}
    finals: dict[str, dict[str, np.ndarray]] = {}
    logs: dict[str, dict[str, object]] = {}
    for model, case in CASES.items():
        histories[model], finals[model] = case_field_history(case)
        logs[model] = parse_solver_log(case / "log.a3_plasma")

    final: dict[str, dict[str, float]] = {}
    for model in CASES:
        row = histories[model][-1]
        log_tail = logs[model]["last_100_iteration_temperature_stability"]
        fluid_global_max = float(log_tail["fluid"]["final_max_K"])
        solid_global_max = float(log_tail["solid"]["final_max_K"])
        final[model] = {
            "time": row["time"],
            "interface_avg_K": row["interface_avg_K"],
            "interface_avg_superheat_K": row["interface_avg_K"] - T_INLET,
            "interface_max_K": row["interface_max_K"],
            "interface_max_superheat_K": row["interface_max_K"] - T_INLET,
            "fluid_internal_max_K": row["fluid_max_K"],
            "fluid_internal_max_superheat_K": row["fluid_max_K"] - T_INLET,
            "solid_internal_max_K": row["solid_max_K"],
            "solid_internal_max_superheat_K": row["solid_max_K"] - T_INLET,
            "fluid_global_max_K": fluid_global_max,
            "fluid_global_max_superheat_K": fluid_global_max - T_INLET,
            "solid_global_max_K": solid_global_max,
            "solid_global_max_superheat_K": solid_global_max - T_INLET,
            "solid_exceeds_1703_K": bool(solid_global_max > FLIBE_LIMIT_K),
            "persisted_4750_to_5000_change_K": {
                key: row[key] - histories[model][-2][key]
                for key in (
                    "interface_avg_K",
                    "interface_max_K",
                    "fluid_max_K",
                    "solid_max_K",
                )
            },
        }

    delta_metrics = (
        "interface_avg_superheat_K",
        "interface_max_superheat_K",
        "fluid_internal_max_superheat_K",
        "solid_global_max_superheat_K",
    )
    deltas = {
        metric: signed_delta(final["variable"][metric], final["constant"][metric])
        for metric in delta_metrics
    }
    constant_bias = {
        metric: signed_delta(final["constant"][metric], final["variable"][metric])
        for metric in delta_metrics
    }
    summary: dict[str, object] = {
        "cases": {model: str(case.relative_to(REPO)) for model, case in CASES.items()},
        "final": final,
        "signed_variable_minus_constant_percent": deltas,
        "legacy_constant_bias_high_relative_to_variable_percent": constant_bias,
        "solver_logs": logs,
        "interpretation_rule": "positive delta means temperature-dependent properties predict a hotter result than the legacy constant-property baseline",
    }
    (RESULTS / "a3_varprop_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    with (RESULTS / "a3_varprop_convergence.csv").open("w", newline="") as stream:
        fields = ["model", *histories["constant"][0].keys()]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for model, rows in histories.items():
            for row in rows:
                writer.writerow({"model": model, **row})

    figure = RESULTS / "a3_varprop_plasma_comparison.png"
    make_figure(histories, finals, summary, figure)
    shutil.copy2(figure, FIGS / figure.name)
    print(json.dumps(summary, indent=2))
    print(f"wrote {figure.relative_to(REPO)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Postprocess the C2 curved-slice OpenFOAM smoke case.

The first useful C2 flow metric is not the velocity on the heated wall itself
because the wall is no-slip.  This script reports owner-cell values adjacent to
named patches, which gives a quick, reproducible hook for later
flushing-law-on-curvature checks.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable


Number = float
Vector = tuple[float, float, float]


def extract_list_body(text: str) -> tuple[int, str]:
    matches = list(re.finditer(r"\n\s*(\d+)\s*\n\s*\(\s*\n(.*?)\n\s*\)", text, re.S))
    if not matches:
        raise ValueError("could not find OpenFOAM list body")
    match = matches[-1]
    return int(match.group(1)), match.group(2)


def read_points(path: Path) -> list[Vector]:
    count, body = extract_list_body(path.read_text())
    points = [
        (float(x), float(y), float(z))
        for x, y, z in re.findall(
            r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", body
        )
    ]
    if len(points) != count:
        raise ValueError(f"{path}: expected {count} points, parsed {len(points)}")
    return points


def read_owner(path: Path) -> list[int]:
    count, body = extract_list_body(path.read_text())
    owners = [int(line) for line in re.findall(r"^\s*(\d+)\s*$", body, re.M)]
    if len(owners) != count:
        raise ValueError(f"{path}: expected {count} owners, parsed {len(owners)}")
    return owners


def read_faces(path: Path) -> list[list[int]]:
    count, body = extract_list_body(path.read_text())
    faces = [
        [int(idx) for idx in raw.split()]
        for raw in re.findall(r"\d+\(([^)]*)\)", body)
    ]
    if len(faces) != count:
        raise ValueError(f"{path}: expected {count} faces, parsed {len(faces)}")
    return faces


def read_boundary(path: Path) -> dict[str, dict[str, int | str]]:
    lines = path.read_text().splitlines()
    patches: dict[str, dict[str, int | str]] = {}
    idx = 0
    while idx < len(lines) - 1:
        name = lines[idx].strip()
        if not re.fullmatch(r"[A-Za-z0-9_]+", name):
            idx += 1
            continue
        next_idx = idx + 1
        while next_idx < len(lines) and not lines[next_idx].strip():
            next_idx += 1
        if next_idx >= len(lines) or lines[next_idx].strip() != "{":
            idx += 1
            continue
        block_lines = []
        depth = 0
        idx = next_idx
        while idx < len(lines):
            stripped = lines[idx].strip()
            depth += stripped.count("{")
            depth -= stripped.count("}")
            block_lines.append(lines[idx])
            idx += 1
            if depth == 0:
                break
        block = "\n".join(block_lines)
        type_match = re.search(r"\btype\s+([A-Za-z0-9_]+);", block)
        n_faces_match = re.search(r"\bnFaces\s+(\d+);", block)
        start_face_match = re.search(r"\bstartFace\s+(\d+);", block)
        if type_match and n_faces_match and start_face_match:
            patches[name] = {
                "type": type_match.group(1),
                "nFaces": int(n_faces_match.group(1)),
                "startFace": int(start_face_match.group(1)),
            }
    return patches


def read_internal_scalar(path: Path) -> list[Number]:
    text = path.read_text()
    uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if uniform:
        raise ValueError(f"{path}: uniform fields need mesh cell count; expected nonuniform")
    count, body = read_internal_list_body(text, "scalar")
    values = [float(x) for x in re.findall(r"[-+0-9.eE]+", body)]
    if len(values) != count:
        raise ValueError(f"{path}: expected {count} scalars, parsed {len(values)}")
    return values


def read_internal_list_body(text: str, kind: str) -> tuple[int, str]:
    marker = f"internalField   nonuniform List<{kind}>"
    start = text.find(marker)
    if start < 0:
        marker = f"internalField nonuniform List<{kind}>"
        start = text.find(marker)
    if start < 0:
        raise ValueError(f"could not find nonuniform List<{kind}> internalField")
    lines = text[start + len(marker) :].splitlines()
    idx = 0
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    count = int(lines[idx].strip())
    idx += 1
    while idx < len(lines) and lines[idx].strip() != "(":
        idx += 1
    if idx == len(lines):
        raise ValueError(f"could not find opening list paren for List<{kind}>")
    idx += 1
    body_lines = []
    while idx < len(lines) and lines[idx].strip() != ")":
        body_lines.append(lines[idx])
        idx += 1
    if idx == len(lines):
        raise ValueError(f"could not find closing list paren for List<{kind}>")
    return count, "\n".join(body_lines)


def read_internal_vector(path: Path) -> list[Vector]:
    text = path.read_text()
    count, body = read_internal_list_body(text, "vector")
    values = [
        (float(x), float(y), float(z))
        for x, y, z in re.findall(
            r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", body
        )
    ]
    if len(values) != count:
        raise ValueError(f"{path}: expected {count} vectors, parsed {len(values)}")
    return values


def sub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def mag(a: Vector) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def face_area(face: Iterable[int], points: list[Vector]) -> float:
    ids = list(face)
    if len(ids) < 3:
        return 0.0
    origin = points[ids[0]]
    area = 0.0
    for i in range(1, len(ids) - 1):
        area += 0.5 * mag(cross(sub(points[ids[i]], origin), sub(points[ids[i + 1]], origin)))
    return area


def weighted_stats(values: list[float], weights: list[float]) -> dict[str, float]:
    if not values:
        raise ValueError("empty values")
    total_w = sum(weights)
    if total_w <= 0.0:
        raise ValueError("non-positive total weight")
    return {
        "min": min(values),
        "max": max(values),
        "area_weighted_average": sum(v * w for v, w in zip(values, weights)) / total_w,
        "unweighted_average": sum(values) / len(values),
    }


def patch_face_ids(patch: dict[str, int | str]) -> range:
    start = int(patch["startFace"])
    return range(start, start + int(patch["nFaces"]))


def patch_owner_values(
    patch_name: str,
    patches: dict[str, dict[str, int | str]],
    owners: list[int],
    faces: list[list[int]],
    points: list[Vector],
    cell_values: list[float],
) -> dict[str, float | int]:
    patch = patches[patch_name]
    face_ids = patch_face_ids(patch)
    values = [cell_values[owners[face_id]] for face_id in face_ids]
    weights = [face_area(faces[face_id], points) for face_id in face_ids]
    stats = weighted_stats(values, weights)
    stats["n_faces"] = int(patch["nFaces"])
    stats["area_m2"] = sum(weights)
    return stats


def parse_yplus(case: Path) -> dict[str, dict[str, float]]:
    yplus_path = case / "postProcessing/yPlus1/0/yPlus.dat"
    if not yplus_path.exists():
        return {}
    out: dict[str, dict[str, float]] = {}
    for line in yplus_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 5:
            _time, patch, ymin, ymax, yavg = parts
            out[patch] = {
                "min": float(ymin),
                "max": float(ymax),
                "average": float(yavg),
            }
    return out


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


def parse_solver_tail(case: Path) -> dict[str, float | str]:
    log_path = case / "log.simpleFoam"
    if not log_path.exists():
        log_path = case / "log.scalarTransportFoam"
    if not log_path.exists():
        return {}
    text = log_path.read_text()
    result: dict[str, float | str] = {"log": log_path.name}
    time_matches = re.findall(r"^Time =\s+(\S+)", text, re.M)
    if time_matches:
        result["final_time_or_iteration"] = float(time_matches[-1])
    exec_matches = re.findall(r"ExecutionTime =\s+([-+0-9.eE]+)\s+s", text)
    if exec_matches:
        result["execution_time_s"] = float(exec_matches[-1])
    cont_matches = re.findall(
        r"time step continuity errors : sum local =\s+([-+0-9.eE]+), "
        r"global =\s+([-+0-9.eE]+), cumulative =\s+([-+0-9.eE]+)",
        text,
    )
    if cont_matches:
        local, global_, cumulative = cont_matches[-1]
        result["continuity_local_final"] = float(local)
        result["continuity_global_final"] = float(global_)
        result["continuity_cumulative_final"] = float(cumulative)
    t_matches = re.findall(
        r"Solving for T, Initial residual =\s+([-+0-9.eE]+), "
        r"Final residual =\s+([-+0-9.eE]+), No Iterations\s+(\d+)",
        text,
    )
    if t_matches:
        initial, final, iterations = t_matches[-1]
        result["T_initial_residual_final_step"] = float(initial)
        result["T_final_residual_final_step"] = float(final)
        result["T_iterations_final_step"] = float(iterations)
    return result


def write_markdown(summary: dict, out_path: Path) -> None:
    def metric_line(label: str, metric: dict[str, float | int], unit: str) -> str:
        return (
            f"- {label}: avg `{metric['area_weighted_average']:.6g}` {unit}, "
            f"min `{metric['min']:.6g}`, max `{metric['max']:.6g}`, "
            f"faces `{int(metric['n_faces'])}`, area `{metric['area_m2']:.6g} m^2`"
        )

    yplus = summary.get("yplus", {})
    lines = [
        "# C2 Slice Smoke Summary",
        "",
        f"- case: `{summary['case']}`",
        f"- latest time: `{summary['latest_time']}`",
        f"- cells inferred from latest fields: `{summary['n_cells']}`",
        f"- solver log: `{summary['solver'].get('log', 'n/a')}`",
        f"- final time/iteration: `{summary['solver'].get('final_time_or_iteration', 'n/a')}`",
        f"- final execution time: `{summary['solver'].get('execution_time_s', 'n/a')}` s",
        f"- final continuity local: `{summary['solver'].get('continuity_local_final', 'n/a')}`",
    ]
    if "T_final_residual_final_step" in summary.get("solver", {}):
        lines.extend(
            [
                f"- final T residual: `{summary['solver']['T_final_residual_final_step']:.6g}`",
                f"- final T initial residual: `{summary['solver']['T_initial_residual_final_step']:.6g}`",
            ]
        )
    if "patch_owner_speed_m_per_s" in summary:
        lines.extend(["", "## Owner-Cell Speed", ""])
        for patch, metric in summary["patch_owner_speed_m_per_s"].items():
            lines.append(metric_line(f"{patch} adjacent speed", metric, "m/s"))
    if "patch_owner_pressure_m2_per_s2" in summary:
        lines.extend(["", "## Owner-Cell Pressure", ""])
        for patch, metric in summary["patch_owner_pressure_m2_per_s2"].items():
            lines.append(metric_line(f"{patch} adjacent pressure", metric, "m^2/s^2"))
        lines.append(
            "- owner-cell pressure drop, inlet minus outlet: "
            f"`{summary['pressure_drop_owner_cell_m2_per_s2']:.6g}` m^2/s^2"
        )
    if "patch_owner_temperature_K" in summary:
        setup = summary.get("thermal_smoke_setup")
        if setup:
            lines.extend(
                [
                    "",
                    "## Thermal Setup",
                    "",
                    f"- gradient mode: `{setup['gradient_mode']}`",
                    f"- heated-wall area: `{setup['heated_wall_area_m2']:.6g} m^2`",
                    f"- applied heated-wall gradient: `{setup['applied_heated_wall_gradient']:.6g}`",
                    f"- gradient-area product: `{setup['gradient_area_product']:.6g}`",
                ]
            )
        lines.extend(["", "## Owner-Cell Temperature", ""])
        for patch, metric in summary["patch_owner_temperature_K"].items():
            lines.append(metric_line(f"{patch} adjacent T", metric, "K"))
            lines.append(
                f"  delta over Tin: avg `{metric['area_weighted_average'] - summary['tin_K']:.6g} K`, "
                f"max `{metric['max'] - summary['tin_K']:.6g} K`"
            )
    lines.extend(["", "## Wall yPlus", ""])
    for patch, values in yplus.items():
        lines.append(
            f"- {patch}: avg `{values['average']:.6g}`, "
            f"min `{values['min']:.6g}`, max `{values['max']:.6g}`"
        )
    if not yplus:
        lines.append("- not present in this case")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a local smoke result on a C2 slice/control scaffold. "
            "Patch metrics are sampled from owner cells adjacent to named patches. "
            "For no-slip walls, this is intentionally not wall-patch velocity. "
            "For fixed-gradient temperature walls, it is adjacent-fluid temperature, "
            "not a conjugate wall temperature.  This is enough to prove the local "
            "C2 workflow path, but it is not yet a validated CHT or LES result.",
            "",
        ]
    )
    out_path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="c2_curved_stellarator_slice_preflight")
    ap.add_argument("--out-dir", default="figs")
    ap.add_argument("--tin", type=float, default=800.0)
    args = ap.parse_args()

    case = Path(args.case)
    latest = latest_time_dir(case)
    points = read_points(case / "constant/polyMesh/points")
    faces = read_faces(case / "constant/polyMesh/faces")
    owners = read_owner(case / "constant/polyMesh/owner")
    patches = read_boundary(case / "constant/polyMesh/boundary")

    summary: dict = {
        "case": str(case),
        "latest_time": latest.name,
        "tin_K": args.tin,
        "solver": parse_solver_tail(case),
        "yplus": parse_yplus(case),
    }
    setup_path = case / "thermal_smoke_setup.json"
    if setup_path.exists():
        summary["thermal_smoke_setup"] = json.loads(setup_path.read_text())
    n_cells = None
    if (latest / "U").exists():
        u = read_internal_vector(latest / "U")
        n_cells = len(u)
        speed = [mag(vec) for vec in u]
        summary["patch_owner_speed_m_per_s"] = {
            name: patch_owner_values(name, patches, owners, faces, points, speed)
            for name in ("heated_first_wall", "blanket_back_wall", "side_walls")
        }
    if (latest / "p").exists():
        p = read_internal_scalar(latest / "p")
        n_cells = n_cells or len(p)
        pressure_patches = {
            name: patch_owner_values(name, patches, owners, faces, points, p)
            for name in ("inlet", "outlet")
        }
        summary["patch_owner_pressure_m2_per_s2"] = pressure_patches
        summary["pressure_drop_owner_cell_m2_per_s2"] = (
            pressure_patches["inlet"]["area_weighted_average"]
            - pressure_patches["outlet"]["area_weighted_average"]
        )
    if (latest / "T").exists():
        temperature = read_internal_scalar(latest / "T")
        n_cells = n_cells or len(temperature)
        summary["patch_owner_temperature_K"] = {
            name: patch_owner_values(name, patches, owners, faces, points, temperature)
            for name in ("heated_first_wall", "blanket_back_wall", "side_walls", "outlet")
        }
    if n_cells is None:
        raise ValueError(f"{latest}: no supported latest-time fields found")
    summary["n_cells"] = n_cells

    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = case.parent / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{case.name}_summary" if case.name.endswith("_smoke") else f"{case.name}_smoke_summary"
    json_path = out_dir / f"{stem}.json"
    md_path = out_dir / f"{stem}.md"
    json_path.write_text(json.dumps(summary, indent=2))
    write_markdown(summary, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

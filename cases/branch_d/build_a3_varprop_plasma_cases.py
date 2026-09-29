#!/usr/bin/env python3
"""Build matched constant- and variable-property A3 plasma-load CHT cases."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np

from postprocess_curved_stellarator_slice_c2 import (
    face_area,
    patch_face_ids,
    read_boundary,
    read_faces,
    read_points,
)


WORK = Path(__file__).resolve().parent
REPO = WORK.parents[1]
SOURCE_CASE = WORK / "cht_wallres_500"
QMAP_PATH = WORK / "figs" / "plasma_qmap.json"
CASE_NAMES = {
    "constant": "a3_const_plasma_cht",
    "variable": "a3_varprop_plasma_cht",
}
TARGET_Q_MEAN = 500_000.0
T_INLET = 800.0
T_LIQUID_MIN = 732.0
T_LIQUID_MAX = 1703.0


def named_block_span(text: str, name: str) -> tuple[int, int]:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*\{{", text)
    if not match:
        raise ValueError(f"block {name!r} not found")
    open_brace = text.find("{", match.start())
    depth = 0
    for idx in range(open_brace, len(text)):
        if text[idx] == "{":
            depth += 1
        elif text[idx] == "}":
            depth -= 1
            if depth == 0:
                return match.start(), idx + 1
    raise ValueError(f"unterminated block {name!r}")


def set_control_entry(text: str, name: str, value: str) -> str:
    pattern = rf"(?m)^(\s*{re.escape(name)}\s+)[^;]+;"
    updated, count = re.subn(pattern, rf"\g<1>{value};", text, count=1)
    if count != 1:
        raise ValueError(f"controlDict entry {name!r} not found exactly once")
    return updated


def configure_control_dict(path: Path, end_time: int, write_interval: int) -> None:
    text = path.read_text()
    for name, value in (
        ("startFrom", "latestTime"),
        ("endTime", str(end_time)),
        ("deltaT", "1"),
        ("writeControl", "timeStep"),
        ("writeInterval", str(write_interval)),
        ("writePrecision", "10"),
        ("writeCompression", "off"),
    ):
        text = set_control_entry(text, name, value)
    start, end = named_block_span(text, "functions")
    text = text[:start] + "functions\n{\n}\n" + text[end:]
    path.write_text(text)


def replace_patch_gradient(path: Path, patch: str, gradients: np.ndarray) -> None:
    text = path.read_text()
    start, end = named_block_span(text, patch)
    block = text[start:end]
    values = "\n".join(f"{value:.12g}" for value in gradients)
    replacement = (
        "gradient        nonuniform List<scalar>\n"
        f"{len(gradients)}\n"
        "(\n"
        f"{values}\n"
        ")\n"
        "        ;"
    )
    pattern = r"\bgradient\s+(?:uniform\s+[-+0-9.eE]+|nonuniform\s+List<scalar>\s+\d+\s*\(.*?\)\s*);"
    updated, count = re.subn(pattern, replacement, block, count=1, flags=re.S)
    if count != 1:
        raise ValueError(f"{path}: expected one gradient entry in patch {patch}")
    path.write_text(text[:start] + updated + text[end:])


def copy_case(destination: Path, force: bool) -> None:
    if destination.exists():
        if not force:
            raise FileExistsError(f"{destination} exists; rerun with --force to rebuild it")
        shutil.rmtree(destination)

    def ignore(_directory: str, names: list[str]) -> set[str]:
        ignored = {
            name
            for name in names
            if name == "postProcessing"
            or name.startswith("processor")
            or name.startswith("log.")
        }
        return ignored

    shutil.copytree(SOURCE_CASE, destination, ignore=ignore)


def heat_flux_map() -> tuple[list[dict[str, float | int]], dict[str, float]]:
    mesh = SOURCE_CASE / "constant" / "heater" / "polyMesh"
    points = read_points(mesh / "points")
    faces = read_faces(mesh / "faces")
    patches = read_boundary(mesh / "boundary")
    face_ids = list(patch_face_ids(patches["heatedBase"]))

    qmap = json.loads(QMAP_PATH.read_text())
    source_s = np.asarray(qmap["s_m"], dtype=float)
    source_q = np.asarray(qmap["q_W_m2"], dtype=float)

    records: list[dict[str, float | int]] = []
    centers_x: list[float] = []
    areas: list[float] = []
    for local_idx, face_id in enumerate(face_ids):
        face = faces[face_id]
        xyz = np.asarray([points[point_id] for point_id in face], dtype=float)
        center = xyz.mean(axis=0)
        area = face_area(face, points)
        centers_x.append(float(center[0]))
        areas.append(float(area))
        records.append(
            {
                "patch_face_index": local_idx,
                "mesh_face_id": face_id,
                "x_m": float(center[0]),
                "y_m": float(center[1]),
                "z_m": float(center[2]),
                "area_m2": float(area),
            }
        )

    x = np.asarray(centers_x)
    area = np.asarray(areas)
    patch_vertices = np.asarray(
        [points[point_id][0] for face_id in face_ids for point_id in faces[face_id]],
        dtype=float,
    )
    x_min = float(patch_vertices.min())
    x_max = float(patch_vertices.max())
    source_coordinate = (x - x_min) / (x_max - x_min) * float(source_s[-1])
    q = np.interp(source_coordinate, source_s, source_q)
    raw_mean = float(np.sum(q * area) / np.sum(area))
    q *= TARGET_Q_MEAN / raw_mean

    solid_text = (SOURCE_CASE / "constant" / "heater" / "thermophysicalProperties").read_text()
    k_match = re.search(r"\bkappa\s+([-+0-9.eE]+)\s*;", solid_text)
    if not k_match:
        raise ValueError("solid kappa not found")
    solid_k = float(k_match.group(1))

    for record, source_s_value, q_value in zip(records, source_coordinate, q):
        record["source_s_m"] = float(source_s_value)
        record["q_W_m2"] = float(q_value)
        record["gradient_K_m"] = float(q_value / solid_k)

    q_mean = float(np.sum(q * area) / np.sum(area))
    stats = {
        "n_patch_faces": len(records),
        "patch_area_m2": float(np.sum(area)),
        "source_period_m": float(source_s[-1]),
        "channel_length_m": x_max - x_min,
        "raw_interpolated_mean_W_m2": raw_mean,
        "applied_mean_W_m2": q_mean,
        "applied_min_W_m2": float(q.min()),
        "applied_max_W_m2": float(q.max()),
        "applied_peak_over_mean": float(q.max() / q_mean),
        "normalization_factor": float(TARGET_Q_MEAN / raw_mean),
        "solid_k_W_mK": solid_k,
    }
    return records, stats


def mu_sohal(T: float) -> float:
    return 1.16e-4 * math.exp(3755.0 / T)


def k_attarian_fit(T: float) -> float:
    # Least-squares line through Attarian Appendix A.10 values at 800-1200 K.
    return 1.212 + 9.0e-5 * T


def property_tables() -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    mu_temperatures = [650.0, T_LIQUID_MIN]
    mu_temperatures.extend(float(T) for T in range(750, 1701, 25))
    mu_temperatures.extend([T_LIQUID_MAX, 1800.0, 2200.0])
    mu_at_min = mu_sohal(T_LIQUID_MIN)
    mu_at_max = mu_sohal(T_LIQUID_MAX)
    mu_table = []
    for T in mu_temperatures:
        if T <= T_LIQUID_MIN:
            value = mu_at_min
        elif T >= T_LIQUID_MAX:
            value = mu_at_max
        else:
            value = mu_sohal(T)
        mu_table.append((T, value))

    k_temperatures = [650.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1400.0, 1600.0, T_LIQUID_MAX, 1800.0, 2200.0]
    k_at_low = k_attarian_fit(800.0)
    k_at_high = k_attarian_fit(T_LIQUID_MAX)
    k_table = []
    for T in k_temperatures:
        if T <= 800.0:
            value = k_at_low
        elif T >= T_LIQUID_MAX:
            value = k_at_high
        else:
            value = k_attarian_fit(T)
        k_table.append((T, value))
    return mu_table, k_table


def table_text(values: list[tuple[float, float]]) -> str:
    return "\n".join(f"            ({T:.8g} {value:.12g})" for T, value in values)


def write_variable_thermo(path: Path, mu_table: list[tuple[float, float]], k_table: list[tuple[float, float]]) -> None:
    text = f"""FoamFile
{{
    version 2.0;
    format ascii;
    class dictionary;
    object thermophysicalProperties;
}}

thermoType
{{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       tabulated;
    thermo          hPolynomial;
    equationOfState icoPolynomial;
    specie          specie;
    energy          sensibleEnthalpy;
}}

mixture
{{
    specie
    {{
        molWeight 33;
    }}
    equationOfState
    {{
        rhoCoeffs<8> (1940 0 0 0 0 0 0 0);
    }}
    thermodynamics
    {{
        Hf 0;
        Sf 0;
        CpCoeffs<8> (2400 0 0 0 0 0 0 0);
    }}
    transport
    {{
        mu
        (
{table_text(mu_table)}
        );
        kappa
        (
{table_text(k_table)}
        );
    }}
}}
"""
    path.write_text(text)


def write_audit_csv(path: Path, records: list[dict[str, float | int]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_property_csv(path: Path, mu_table: list[tuple[float, float]], k_table: list[tuple[float, float]]) -> None:
    temperatures = sorted({T for T, _ in mu_table} | {T for T, _ in k_table})
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["T_K", "mu_Pa_s", "k_W_mK"])
        for T in temperatures:
            mu = next((value for temp, value in mu_table if temp == T), "")
            kappa = next((value for temp, value in k_table if temp == T), "")
            writer.writerow([f"{T:.8g}", f"{mu:.12g}" if mu != "" else "", f"{kappa:.12g}" if kappa != "" else ""])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="replace only the two generated A3 case directories")
    parser.add_argument("--end-time", type=int, default=5000)
    parser.add_argument("--write-interval", type=int, default=250)
    args = parser.parse_args()
    if args.end_time <= 3000 or (args.end_time - 3000) % args.write_interval:
        raise ValueError("end-time must exceed 3000 and land on the write interval")

    records, heat_stats = heat_flux_map()
    gradients = np.asarray([float(row["gradient_K_m"]) for row in records])
    mu_table, k_table = property_tables()

    for model, case_name in CASE_NAMES.items():
        case = WORK / case_name
        copy_case(case, args.force)
        configure_control_dict(case / "system" / "controlDict", args.end_time, args.write_interval)
        for time_name in ("0", "3000"):
            replace_patch_gradient(case / time_name / "heater" / "T", "heatedBase", gradients)
        write_audit_csv(case / "A3_APPLIED_HEAT_FLUX.csv", records)

        if model == "variable":
            write_variable_thermo(
                case / "constant" / "bottomWater" / "thermophysicalProperties",
                mu_table,
                k_table,
            )
            write_property_csv(case / "A3_PROPERTY_TABLES.csv", mu_table, k_table)

        metadata = {
            "track": "A3 temperature-dependent-property CHT",
            "model": model,
            "source_case": str(SOURCE_CASE.relative_to(REPO)),
            "common_start_time": 3000,
            "end_time": args.end_time,
            "write_interval": args.write_interval,
            "inlet_temperature_K": T_INLET,
            "heat_flux_source": str(QMAP_PATH.relative_to(REPO)),
            "heat_flux_mapping": "one toroidal q-map period mapped linearly onto the 0.4 m streamwise heated base and area-normalized",
            "heat_flux": heat_stats,
            "property_model": (
                {
                    "mu": "1.16e-4 exp(3755/T) Pa s; Sohal/INL (Williams), source range 873-1073 K; explicitly extrapolated over 732-1703 K and held flat outside",
                    "kappa": "least-squares k=1.212+9e-5*T W/m/K through Attarian Appendix A.10 800-1200 K values; extrapolated to 1703 K and held flat outside",
                    "rho_kg_m3": 1940.0,
                    "Cp_J_kgK": 2400.0,
                    "table_interpolation": "OpenFOAM tabulatedTransport linear interpolation; flat first/last table segments prevent unbounded extrapolation",
                }
                if model == "variable"
                else {
                    "mu_Pa_s": 0.006,
                    "kappa_W_mK": 1.0,
                    "rho_kg_m3": 1940.0,
                    "Cp_J_kgK": 2400.0,
                }
            ),
        }
        (case / "A3_CASE_METADATA.json").write_text(json.dumps(metadata, indent=2) + "\n")
        (case / "README_A3.txt").write_text(
            "Generated by build_a3_varprop_plasma_cases.py.\n"
            "Run after sourcing OpenFOAM v2512: chtMultiRegionSimpleFoam -case <case>.\n"
            "The 3000 field is the common persisted warm start; writes are forced at 250-step intervals.\n"
        )
        print(f"built {case.relative_to(REPO)}")

    print(json.dumps(heat_stats, indent=2))


if __name__ == "__main__":
    main()

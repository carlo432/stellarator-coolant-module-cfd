#!/usr/bin/env python3
"""Compute Branch D dimensionless heat-transfer and loss diagnostics.

The current module is an idealized expansion/bend/step geometry with passive
scalar heat transport, so these values are reporting diagnostics rather than
validated channel correlations. They are still useful because they translate
the wall-temperature and near-wall-flushing result into Nu, St, Re, and loss
coefficient language.
"""

from __future__ import annotations

import csv
import base64
import os
import re
import struct
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"

TIN_K = 800.0
HEAT_FLUX_W_M2 = 100_000.0
RHO_KG_M3 = 1940.0
MU_PA_S = 0.006
K_W_M_K = 1.0
CP_J_KG_K = 2400.0
INLET_WIDTH_M = 0.02
INLET_DEPTH_M = 0.04
MODULE_LENGTH_M = 0.40
HYDRAULIC_DIAMETER_M = 2.0 * INLET_WIDTH_M * INLET_DEPTH_M / (INLET_WIDTH_M + INLET_DEPTH_M)
INLET_SPEED_M_S = 1.1598


@dataclass(frozen=True)
class CaseSpec:
    geometry: str
    case_pair: str
    rans_time_dir: Path
    thermal_time_dir: Path
    heated_wall_vtp: Path
    nearwall_speed_m_s: float
    pressure_drop_pa: float | None
    pressure_drop_source: str
    marker_color: str


CASES = [
    CaseSpec(
        geometry="90-deg bend",
        case_pair="v28/v29",
        rans_time_dir=ROOT / "cases/baseline_module/openfoam_cases/v28_isothermal_rans_flibe_re10000_bend90/1000",
        thermal_time_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v29_scalar_temperature_flibe_re10000_heatflux100kw_bend90_dt00025/2",
        heated_wall_vtp=ROOT
        / "cases/baseline_module/openfoam_cases/v29_scalar_temperature_flibe_re10000_heatflux100kw_bend90_dt00025/VTK/v29_scalar_temperature_flibe_re10000_heatflux100kw_bend90_dt00025_8000/boundary/heated_wall.vtp",
        nearwall_speed_m_s=0.86,
        pressure_drop_pa=None,
        pressure_drop_source="patch-average inlet pressure from v28 postProcessing",
        marker_color="#2f9e44",
    ),
    CaseSpec(
        geometry="straight",
        case_pair="v13/v14",
        rans_time_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v13_isothermal_rans_flibe_re10000_thermal_refined/1000",
        thermal_time_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025/2",
        heated_wall_vtp=ROOT
        / "cases/baseline_module/openfoam_cases/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025/VTK/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025_8000/boundary/heated_wall.vtp",
        nearwall_speed_m_s=0.47,
        pressure_drop_pa=None,
        pressure_drop_source="patch-average inlet minus outlet pressure from v13 postProcessing",
        marker_color="#1c7ed6",
    ),
    CaseSpec(
        geometry="outlet-offset",
        case_pair="v22/v23",
        rans_time_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v22_isothermal_rans_flibe_re10000_bent_recirculation/1000",
        thermal_time_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v23_scalar_temperature_flibe_re10000_heatflux100kw_bent_recirculation_dt00025/2",
        heated_wall_vtp=ROOT
        / "cases/baseline_module/openfoam_cases/v23_scalar_temperature_flibe_re10000_heatflux100kw_bent_recirculation_dt00025/VTK/v23_scalar_temperature_flibe_re10000_heatflux100kw_bent_recirculation_dt00025_8000/boundary/heated_wall.vtp",
        nearwall_speed_m_s=0.39,
        pressure_drop_pa=1983.0,
        pressure_drop_source="documented Branch D v22 result; patch-average files not retained",
        marker_color="#e03131",
    ),
    CaseSpec(
        geometry="backward-facing step",
        case_pair="v30/v31",
        rans_time_dir=ROOT / "cases/baseline_module/openfoam_cases/v30_isothermal_rans_flibe_re10000_bfs/1000",
        thermal_time_dir=ROOT
        / "cases/baseline_module/openfoam_cases/v31_scalar_temperature_flibe_re10000_heatflux100kw_bfs_dt00025/2",
        heated_wall_vtp=ROOT
        / "cases/baseline_module/openfoam_cases/v31_scalar_temperature_flibe_re10000_heatflux100kw_bfs_dt00025/VTK/v31_scalar_temperature_flibe_re10000_heatflux100kw_bfs_dt00025_8000/boundary/heated_wall.vtp",
        nearwall_speed_m_s=0.296,
        pressure_drop_pa=None,
        pressure_drop_source="not reported; no patch-average pressure-drop file retained",
        marker_color="#f08c00",
    ),
]

VTK_DTYPES = {
    "Float32": np.dtype("<f4"),
    "Float64": np.dtype("<f8"),
    "Int32": np.dtype("<i4"),
    "UInt64": np.dtype("<u8"),
}


def _internal_block(text: str) -> tuple[str, int, str]:
    match = re.search(
        r"internalField\s+nonuniform\s+List<(?P<kind>\w+)>\s+"
        r"(?P<count>\d+)\s*\((?P<body>.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError("Could not find nonuniform internalField block")
    return match.group("kind"), int(match.group("count")), match.group("body")


def read_scalar_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "scalar":
        raise ValueError(f"{path} is {kind}, not scalar")
    data = np.fromstring(body, sep=" ")
    if data.size != count:
        raise ValueError(f"{path}: expected {count}, got {data.size}")
    return data


def read_vector_field(path: Path) -> np.ndarray:
    kind, count, body = _internal_block(path.read_text(encoding="ascii"))
    if kind != "vector":
        raise ValueError(f"{path} is {kind}, not vector")
    rows = re.findall(r"\(([^()]+)\)", body)
    data = np.array([[float(x) for x in row.split()] for row in rows], dtype=float)
    if data.shape != (count, 3):
        raise ValueError(f"{path}: expected {(count, 3)}, got {data.shape}")
    return data


def read_patch_scalar_values(field_path: Path, patch_name: str) -> np.ndarray:
    text = field_path.read_text(encoding="ascii")
    patch_match = re.search(
        rf"\b{re.escape(patch_name)}\s*\{{(?P<body>.*?)\n\s*\}}",
        text,
        re.DOTALL,
    )
    if not patch_match:
        raise ValueError(f"Patch {patch_name!r} not found in {field_path}")
    body = patch_match.group("body")
    value_match = re.search(
        r"value\s+nonuniform\s+List<scalar>\s+(?P<count>\d+)\s*\((?P<values>.*?)\)\s*;",
        body,
        re.DOTALL,
    )
    if not value_match:
        raise ValueError(f"Patch {patch_name!r} in {field_path} does not have nonuniform scalar values")
    count = int(value_match.group("count"))
    data = np.fromstring(value_match.group("values"), sep=" ")
    if data.size != count:
        raise ValueError(f"{field_path}: patch {patch_name!r} expected {count}, got {data.size}")
    return data


def decode_vtk_data_array(element: ET.Element) -> np.ndarray:
    vtk_type = element.attrib["type"]
    dtype = VTK_DTYPES[vtk_type]
    if element.attrib.get("format", "ascii") == "ascii":
        values = np.fromstring(element.text or "", sep=" ", dtype=dtype)
        components = int(element.attrib.get("NumberOfComponents", "1"))
        if components > 1:
            values = values.reshape((-1, components))
        return values

    encoded = "".join((element.text or "").split())
    blob = base64.b64decode(encoded)
    if len(blob) < 8:
        raise ValueError(f"DataArray {element.attrib.get('Name')} is missing UInt64 payload header")
    payload_length = struct.unpack("<Q", blob[:8])[0]
    payload = blob[8 : 8 + payload_length]
    values = np.frombuffer(payload, dtype=dtype).copy()
    components = int(element.attrib.get("NumberOfComponents", "1"))
    if components > 1:
        values = values.reshape((-1, components))
    return values


def read_vtp_cell_temperature(path: Path) -> np.ndarray:
    root = ET.parse(path).getroot()
    piece = root.find("./PolyData/Piece")
    if piece is None:
        raise ValueError(f"{path} does not contain a PolyData Piece")
    cell_data_parent = piece.find("CellData")
    if cell_data_parent is None:
        raise ValueError(f"{path} is missing CellData")
    for element in cell_data_parent.iter("DataArray"):
        if element.attrib.get("Name") == "T":
            return decode_vtk_data_array(element).astype(float)
    raise ValueError(f"{path} does not contain CellData array T")


def read_surface_average(path: Path) -> float:
    last_value: float | None = None
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) >= 2:
            last_value = float(parts[1])
    if last_value is None:
        raise ValueError(f"No surface-average value found in {path}")
    return last_value


def postprocessed_pressure_drop(case: CaseSpec) -> float | None:
    if case.case_pair == "v13/v14":
        base = case.rans_time_dir.parents[0]
        inlet = read_surface_average(base / "postProcessing/inletPressure/1000/surfaceFieldValue.dat")
        outlet = read_surface_average(base / "postProcessing/outletPressure/1000/surfaceFieldValue.dat")
        return (inlet - outlet) * RHO_KG_M3
    if case.case_pair == "v28/v29":
        base = case.rans_time_dir.parents[0]
        inlet = read_surface_average(
            base / "postProcessing/patchAverage(name=inlet,fields=(p))/1000/surfaceFieldValue.dat"
        )
        return inlet * RHO_KG_M3
    return None


def duct_cell_pressure_proxy(case: CaseSpec) -> tuple[float | None, int, int]:
    c_path = case.rans_time_dir / "C"
    p_path = case.rans_time_dir / "p"
    if not c_path.exists() or not p_path.exists():
        return None, 0, 0
    centers = read_vector_field(c_path)
    pressure = read_scalar_field(p_path)
    x = centers[:, 0]
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    inlet_mask = x <= xmin + 0.03
    outlet_mask = x >= xmax - 0.03
    if np.count_nonzero(inlet_mask) == 0 or np.count_nonzero(outlet_mask) == 0:
        return None, int(np.count_nonzero(inlet_mask)), int(np.count_nonzero(outlet_mask))
    delta_p_kinematic = float(np.mean(pressure[inlet_mask]) - np.mean(pressure[outlet_mask]))
    return delta_p_kinematic * RHO_KG_M3, int(np.count_nonzero(inlet_mask)), int(np.count_nonzero(outlet_mask))


def module_loss_coefficient(delta_p_pa: float) -> float:
    return 2.0 * delta_p_pa / (RHO_KG_M3 * INLET_SPEED_M_S**2)


def global_darcy_f_diagnostic(loss_k: float) -> float:
    return loss_k * HYDRAULIC_DIAMETER_M / MODULE_LENGTH_M


def dimensionless_row(case: CaseSpec) -> dict[str, str | float | int]:
    wall_t = read_vtp_cell_temperature(case.heated_wall_vtp)
    wall_delta_t = wall_t - TIN_K
    wall_mean_delta_t = float(np.mean(wall_delta_t))
    wall_peak_delta_t = float(np.max(wall_delta_t))
    h_mean = HEAT_FLUX_W_M2 / wall_mean_delta_t
    h_peak = HEAT_FLUX_W_M2 / wall_peak_delta_t
    nu_mean = h_mean * HYDRAULIC_DIAMETER_M / K_W_M_K
    nu_peak = h_peak * HYDRAULIC_DIAMETER_M / K_W_M_K
    re_local = RHO_KG_M3 * case.nearwall_speed_m_s * HYDRAULIC_DIAMETER_M / MU_PA_S
    re_inlet = RHO_KG_M3 * INLET_SPEED_M_S * HYDRAULIC_DIAMETER_M / MU_PA_S
    st_peak_local = h_peak / (RHO_KG_M3 * CP_J_KG_K * case.nearwall_speed_m_s)
    st_mean_local = h_mean / (RHO_KG_M3 * CP_J_KG_K * case.nearwall_speed_m_s)
    st_peak_inlet = h_peak / (RHO_KG_M3 * CP_J_KG_K * INLET_SPEED_M_S)
    st_mean_inlet = h_mean / (RHO_KG_M3 * CP_J_KG_K * INLET_SPEED_M_S)

    post_dp = postprocessed_pressure_drop(case)
    proxy_dp, proxy_inlet_cells, proxy_outlet_cells = duct_cell_pressure_proxy(case)
    primary_dp = case.pressure_drop_pa if case.pressure_drop_pa is not None else post_dp
    primary_source = case.pressure_drop_source
    if primary_dp is None and proxy_dp is not None and proxy_dp > 0.0:
        primary_dp = proxy_dp
        primary_source = "duct-cell pressure proxy; no patch-average file retained"
    elif primary_dp is None and proxy_dp is not None and proxy_dp <= 0.0:
        primary_source = "not reported; duct-cell proxy failed sign check"
    loss_k = module_loss_coefficient(primary_dp) if primary_dp is not None else np.nan
    f_diag = global_darcy_f_diagnostic(loss_k) if primary_dp is not None else np.nan

    return {
        "geometry": case.geometry,
        "case_pair": case.case_pair,
        "thermal_time_dir": str(case.thermal_time_dir.relative_to(ROOT)),
        "heated_wall_vtp": str(case.heated_wall_vtp.relative_to(ROOT)),
        "wall_face_count": int(wall_t.size),
        "heat_flux_W_m2": HEAT_FLUX_W_M2,
        "wall_mean_deltaT_K": wall_mean_delta_t,
        "wall_peak_deltaT_K": wall_peak_delta_t,
        "wall_mean_metric_status": "preferred diagnostic; wall average is more grid-stable than single-cell peak",
        "wall_peak_metric_status": "mesh-sensitive same-mesh hotspot diagnostic; not grid-converged",
        "nearwall_mean_speed_m_s": case.nearwall_speed_m_s,
        "local_Re_Dh_nearwall": re_local,
        "inlet_Re_Dh_reference": re_inlet,
        "h_mean_W_m2K": h_mean,
        "h_peak_W_m2K": h_peak,
        "Nu_mean_Dh": nu_mean,
        "Nu_peak_Dh": nu_peak,
        "St_mean_local_nearwall": st_mean_local,
        "St_peak_local_nearwall": st_peak_local,
        "St_mean_inlet_reference": st_mean_inlet,
        "St_peak_inlet_reference": st_peak_inlet,
        "pressure_drop_primary_Pa": primary_dp if primary_dp is not None else "",
        "pressure_drop_primary_source": primary_source,
        "pressure_drop_patch_or_documented_Pa": case.pressure_drop_pa if case.pressure_drop_pa is not None else post_dp or "",
        "pressure_drop_cell_proxy_Pa": proxy_dp if proxy_dp is not None else "",
        "pressure_proxy_inlet_cell_count": proxy_inlet_cells,
        "pressure_proxy_outlet_cell_count": proxy_outlet_cells,
        "module_loss_coefficient_K": loss_k if np.isfinite(loss_k) else "",
        "global_darcy_f_diagnostic": f_diag if np.isfinite(f_diag) else "",
        "Pr_flibe_like": CP_J_KG_K * MU_PA_S / K_W_M_K,
        "Dh_m": HYDRAULIC_DIAMETER_M,
        "inlet_speed_m_s": INLET_SPEED_M_S,
    }


def write_csv(path: Path, rows: list[dict[str, str | float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_plot(path: Path, advisor_path: Path, rows: list[dict[str, str | float | int]]) -> None:
    labels = [str(row["geometry"]) for row in rows]
    local_re = np.array([float(row["local_Re_Dh_nearwall"]) for row in rows])
    mean_dt = np.array([float(row["wall_mean_deltaT_K"]) for row in rows])
    peak_dt = np.array([float(row["wall_peak_deltaT_K"]) for row in rows])
    nu_mean = np.array([float(row["Nu_mean_Dh"]) for row in rows])
    nu_peak = np.array([float(row["Nu_peak_Dh"]) for row in rows])
    st_mean = np.array([float(row["St_mean_local_nearwall"]) for row in rows])
    loss_k = np.array(
        [
            float(row["module_loss_coefficient_K"]) if row["module_loss_coefficient_K"] != "" else np.nan
            for row in rows
        ]
    )
    colors = [case.marker_color for case in CASES]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.0), dpi=180)

    ax = axes[0]
    for label, x_value, y_value, color in zip(labels, local_re, mean_dt, colors):
        display = "BFS" if label == "backward-facing step" else label
        ax.scatter(x_value, y_value, s=110, color=color, edgecolor="black", linewidth=0.8, zorder=3)
        ax.text(x_value + 90, y_value + 0.08, display, fontsize=8.5, va="bottom")
    ax.set_title("Wall-average dT vs local near-wall Re", fontsize=12, fontweight="bold")
    ax.set_xlabel("Local Re_Dh from near-wall speed")
    ax.set_ylabel("Mean heated-wall dT [K]")
    ax.set_ylim(42.5, 47.4)
    ax.grid(True, alpha=0.28)

    ax2 = axes[1]
    x = np.arange(len(rows))
    bars = ax2.bar(x - 0.18, nu_mean, width=0.36, color=colors, edgecolor="black", linewidth=0.6, label="Nu_mean")
    ax2.scatter(x + 0.18, nu_peak, color="#212529", marker="x", s=58, linewidths=1.6, label="Nu_peak diagnostic")
    ax2.set_xticks(x)
    ax2.set_xticklabels(["bend", "straight", "offset", "BFS"], rotation=0)
    ax2.set_ylabel("Nu_Dh")
    ax2.set_title("Wall-average Nu, local St, and module loss", fontsize=12, fontweight="bold")
    ax2.grid(True, axis="y", alpha=0.25)
    for bar, value in zip(bars, nu_mean):
        ax2.text(bar.get_x() + bar.get_width() / 2.0, value + 0.45, f"{value:.1f}", ha="center", fontsize=8)

    ax3 = ax2.twinx()
    ax3.plot(x + 0.18, st_mean * 1e3, color="#343a40", marker="o", linewidth=1.4, label="St_mean x1e3")
    ax3.plot(x, loss_k, color="#6741d9", marker="s", linewidth=1.2, linestyle="--", label="loss K")
    ax3.set_ylabel("St_mean x 1e3 / loss coefficient K")
    ax3.set_ylim(bottom=0.0)

    lines, line_labels = ax2.get_legend_handles_labels()
    lines2, line_labels2 = ax3.get_legend_handles_labels()
    ax2.legend(lines + lines2, line_labels + line_labels2, loc="lower left", fontsize=8, frameon=True)

    fig.suptitle("Branch D Dimensionless Heat-Transfer Diagnostics", fontsize=14, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "Wall-average Nu/St are preferred after peak-dT GCI warning; peak values are same-mesh hotspot diagnostics only.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    advisor_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(advisor_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rows = [dimensionless_row(case) for case in CASES]
    csv_path = OUT_DIR / "branch_d_dimensionless_metrics.csv"
    fig_path = OUT_DIR / "branch_d_dimensionless_metrics.png"
    advisor_path = ADVISOR_DIR / "25_branch_d_dimensionless_metrics.png"
    write_csv(csv_path, rows)
    make_plot(fig_path, advisor_path, rows)
    print(csv_path)
    print(fig_path)
    print(advisor_path)
    for row in rows:
        print(
            f"{row['case_pair']} {row['geometry']}: "
            f"Nu_mean={float(row['Nu_mean_Dh']):.2f}, "
            f"Nu_peak_diag={float(row['Nu_peak_Dh']):.2f}, "
            f"St_mean_local={float(row['St_mean_local_nearwall']):.6f}, "
            f"K={row['module_loss_coefficient_K']}"
        )


if __name__ == "__main__":
    main()

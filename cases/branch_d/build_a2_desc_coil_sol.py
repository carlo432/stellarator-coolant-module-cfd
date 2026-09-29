#!/usr/bin/env python3
"""Build the A2 W7-X-equilibrium-derived coil/SOL surrogate.

This is deliberately not labelled as the W7-X engineering coil set.  DESC 0.17.2
ships a fixed-boundary W7-X equilibrium but no W7-X coil geometry.  The script
therefore fits a REGCOIL current-potential surface to that equilibrium, optionally
discretizes the sheet current into filament coils, and includes the equilibrium's
plasma-current contribution through a virtual-casing surface integral.

Stages are restartable because the REGCOIL fit and converged wall-field audit are
the costly parts::

    python build_a2_desc_coil_sol.py --stage fit
    python build_a2_desc_coil_sol.py --stage spline
    python build_a2_desc_coil_sol.py --stage trace
    python build_a2_desc_coil_sol.py --stage all

The SOL-width calculation is a reduced transport closure, not a first-principles
edge-plasma solution.  Connection length supplies the parallel residence time;
perpendicular thermal diffusivity and LCFS temperature remain explicit inputs:

    lambda_far = sqrt(chi_perp * L_connection / c_s).

The 1--3 m2/s and 100--200 eV band is the reactor-study band used by Lion et al.
(Stellaris, Fusion Engineering and Design 214 (2025) 114868).  Mau et al.
(Fusion Science and Technology 54 (2008) 771--786, doi:10.13182/FST08-27)
independently used 1 m2/s for ARIES-CS diffusive field-line tracing.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import time
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-a2-desc")
warnings.filterwarnings("ignore", message=".*Matplotlib.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from desc.coils import CoilSet
from desc.examples import get
from desc.geometry import FourierRZToroidalSurface
from desc.grid import Grid, LinearGrid
from desc.magnetic_fields import (
    FourierCurrentPotentialField,
    field_line_integrate,
    solve_regularized_surface_current,
)
from desc.utils import rpz2xyz, rpz2xyz_vec, xyz2rpz_vec


ROOT = Path(__file__).resolve().parents[2]
CASE_DIR = ROOT / "cases" / "branch_d" / "a2_desc_coil_sol"
RESULT_DIR = ROOT / "results" / "branch_d" / "a2_desc_coil_sol"
FIG_DIR = ROOT / "cases" / "branch_d" / "figs"

FIELD_FILE = CASE_DIR / "w7x_regcoil_current_potential.h5"
COIL_FILE = CASE_DIR / "w7x_regcoil_filaments.h5"
FIT_METRICS_FILE = RESULT_DIR / "a2_fit_metrics.json"
SPLINE_METRICS_FILE = RESULT_DIR / "a2_spline_metrics.json"

NFP = 5
WINDING_OFFSET_M = 0.45
WINDING_FIT_M = 20
WINDING_FIT_N = 20
CURRENT_M = 7
CURRENT_N = 7
REG_LAMBDAS = np.array([1e-20, 1e-18, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12])

REG_GRID_M = 20
REG_GRID_N = 20
VALID_GRID_M = 30
VALID_GRID_N = 30
FIELD_SOURCE_M = 26
FIELD_SOURCE_N = 26
TRACE_SOURCE_M = 16
TRACE_SOURCE_N = 16
VC_SOURCE_M = 24
VC_SOURCE_N = 24
VC_WALL_SOURCE_LEVELS = (128, 160, 192)

NW = 192
CLEARANCE_M = 0.03
QWALL_W_M2 = 0.5e6
F_RAD = 0.5
S_DUCT_M = np.radians(95.0) * 0.18

M_DT_KG = 2.5 * 1.67262192595e-27
EV_J = 1.602176634e-19


def ensure_dirs() -> None:
    for path in (CASE_DIR, RESULT_DIR, FIG_DIR):
        path.mkdir(parents=True, exist_ok=True)


def json_ready(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    return value


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(json_ready(data), indent=2) + "\n", encoding="ascii")


def boundary_grid(m: int, n: int) -> LinearGrid:
    return LinearGrid(
        rho=np.array([1.0]), M=m, N=n, NFP=NFP, sym=False, endpoint=False
    )


def make_winding_surface(eq):
    """Fit a manufacturable Fourier surface outside the W7-X LCFS.

    FourierRZToroidalSurface cannot fit a true 3-D normal offset when zeta != phi.
    The offset therefore follows the outward normal projected into each constant-phi
    R-Z plane.  This is recorded as a surrogate geometry limitation.
    """
    sample_grid = boundary_grid(40, 40)
    data = eq.surface.compute(["R", "phi", "Z", "n_rho"], grid=sample_grid)
    normal_rz = np.asarray(data["n_rho"])[:, [0, 2]].copy()
    normal_rz /= np.linalg.norm(normal_rz, axis=1)[:, None]
    target = np.column_stack(
        [
            np.asarray(data["R"]) + WINDING_OFFSET_M * normal_rz[:, 0],
            np.asarray(data["phi"]),
            np.asarray(data["Z"]) + WINDING_OFFSET_M * normal_rz[:, 1],
        ]
    )
    surface = FourierRZToroidalSurface.from_values(
        target,
        sample_grid.nodes[:, 1],
        M=WINDING_FIT_M,
        N=WINDING_FIT_N,
        NFP=NFP,
        sym=True,
    )
    fitted = surface.compute(["R", "phi", "Z"], grid=sample_grid)
    residual_rz = np.column_stack(
        [np.asarray(fitted["R"]) - target[:, 0], np.asarray(fitted["Z"]) - target[:, 2]]
    )
    fit_norm = np.linalg.norm(residual_rz, axis=1)
    metrics = {
        "offset_definition": "R-Z-projected outward normal at fixed cylindrical phi",
        "requested_offset_m": WINDING_OFFSET_M,
        "fit_rms_m": float(np.sqrt(np.mean(fit_norm**2))),
        "fit_max_m": float(np.max(fit_norm)),
        "target_R_range_m": [float(target[:, 0].min()), float(target[:, 0].max())],
        "target_Z_range_m": [float(target[:, 2].min()), float(target[:, 2].max())],
        "orientation": float(surface._compute_orientation()),
    }
    return surface, metrics


def normalized_bn_stats(bn: np.ndarray, bmag: np.ndarray) -> dict:
    rel = np.abs(np.asarray(bn)) / np.maximum(np.asarray(bmag), 1e-12)
    return {
        "bn_abs_max_T": float(np.max(np.abs(bn))),
        "bn_over_b_rms": float(np.sqrt(np.mean(rel**2))),
        "bn_over_b_mean_abs": float(np.mean(rel)),
        "bn_over_b_max_abs": float(np.max(rel)),
        "bn_over_b_p95_abs": float(np.quantile(rel, 0.95)),
    }


def select_regularization(rows: list[dict]) -> int:
    """Select the largest regularization meeting a 1% max-|Bn|/|B| gate."""
    acceptable = [i for i, row in enumerate(rows) if row["bn_over_b_max_abs"] <= 0.01]
    if acceptable:
        return acceptable[-1]
    relaxed = [i for i, row in enumerate(rows) if row["bn_over_b_max_abs"] <= 0.02]
    if relaxed:
        return relaxed[-1]
    return int(np.argmin([row["bn_over_b_rms"] for row in rows]))


def direct_lcfs_standoff(eq):
    phi = np.linspace(0.0, 2 * np.pi / NFP, NW, endpoint=False)
    theta = np.linspace(0.0, 2 * np.pi, 513, endpoint=False)
    th, ph = np.meshgrid(theta, phi, indexing="ij")
    nodes = np.column_stack([np.ones(th.size), th.ravel(), ph.ravel()])
    data = eq.compute(["R", "Z"], grid=Grid(nodes, sort=False))
    rr = np.asarray(data["R"]).reshape(theta.size, NW)
    zz = np.asarray(data["Z"]).reshape(theta.size, NW)
    imax = np.argmax(rr, axis=0)
    rout = rr[imax, np.arange(NW)]
    zout = zz[imax, np.arange(NW)]
    rwall = float(rout.max() + CLEARANCE_M)
    return phi, rout, zout, rwall, rwall - rout


def fit_stage() -> dict:
    ensure_dirs()
    eq = get("W7-X")
    if eq.NFP != NFP:
        raise RuntimeError(f"Expected W7-X NFP={NFP}, got {eq.NFP}")

    print("[fit] constructing projected-normal winding surface", flush=True)
    winding, winding_metrics = make_winding_surface(eq)
    winding.save(CASE_DIR / "w7x_winding_surface.h5")

    current_field = FourierCurrentPotentialField.from_surface(
        winding, M_Phi=CURRENT_M, N_Phi=CURRENT_N, sym_Phi="sin"
    )
    source_grid = boundary_grid(REG_GRID_M, REG_GRID_N)
    eval_grid = boundary_grid(REG_GRID_M, REG_GRID_N)

    print(
        f"[fit] REGCOIL sweep: M_Phi={CURRENT_M}, N_Phi={CURRENT_N}, "
        f"grid=({REG_GRID_M},{REG_GRID_N})",
        flush=True,
    )
    started = time.time()
    fields, data = solve_regularized_surface_current(
        current_field,
        eq,
        lambda_regularization=REG_LAMBDAS,
        current_helicity=(1, 0),
        vacuum=False,
        regularization_type="regcoil",
        source_grid=source_grid,
        eval_grid=eval_grid,
        vc_source_grid=source_grid,
        verbose=1,
        chunk_size=32,
        B_plasma_chunk_size=32,
    )
    elapsed = time.time() - started
    bmag = np.asarray(eq.compute("|B|", grid=eval_grid)["|B|"])
    sweep_rows = []
    for i, lam in enumerate(REG_LAMBDAS):
        row = {
            "lambda_regularization": float(lam),
            "chi2_B": float(data["chi^2_B"][i]),
            "chi2_K": float(data["chi^2_K"][i]),
            "K_max_A_m": float(np.max(np.asarray(data["|K|"][i]))),
            "K_rms_A_m": float(np.sqrt(np.mean(np.asarray(data["|K|"][i]) ** 2))),
        }
        row.update(normalized_bn_stats(np.asarray(data["Bn_total"][i]), bmag))
        sweep_rows.append(row)

    coarse_selected_idx = select_regularization(sweep_rows)
    print(
        f"[fit] coarse-grid candidate lambda={REG_LAMBDAS[coarse_selected_idx]:.1e}; "
        "auditing every candidate independently",
        flush=True,
    )

    print(
        f"[fit] independent {VALID_GRID_M}x{VALID_GRID_N} boundary-normal sweep",
        flush=True,
    )
    valid_grid = boundary_grid(VALID_GRID_M, VALID_GRID_N)
    valid_source = boundary_grid(FIELD_SOURCE_M, FIELD_SOURCE_N)
    bmag_valid = np.asarray(eq.compute("|B|", grid=valid_grid)["|B|"])
    independent_rows = []
    for i, candidate in enumerate(fields):
        bn_valid, _ = candidate.compute_Bnormal(
            eq,
            eval_grid=valid_grid,
            source_grid=valid_source,
            vc_source_grid=valid_source,
            chunk_size=32,
            B_plasma_chunk_size=32,
        )
        stats = normalized_bn_stats(np.asarray(bn_valid), bmag_valid)
        stats["lambda_regularization"] = float(REG_LAMBDAS[i])
        stats["gate"] = "p95 <= 1% and max <= 2%"
        stats["passed"] = bool(
            stats["bn_over_b_p95_abs"] <= 0.01
            and stats["bn_over_b_max_abs"] <= 0.02
        )
        independent_rows.append(stats)
        sweep_rows[i]["independent"] = stats
        print(
            f"      lambda={REG_LAMBDAS[i]:.1e}: "
            f"p95={stats['bn_over_b_p95_abs']:.3%}, "
            f"max={stats['bn_over_b_max_abs']:.3%}",
            flush=True,
        )
    passing = [i for i, stats in enumerate(independent_rows) if stats["passed"]]
    selected_idx = passing[-1] if passing else int(
        np.argmin([stats["bn_over_b_p95_abs"] for stats in independent_rows])
    )
    field = fields[selected_idx]
    validation = independent_rows[selected_idx]
    if not passing:
        raise RuntimeError(
            "No REGCOIL candidate passed the independent boundary-normal gate; "
            "best candidate: "
            f"p95={validation['bn_over_b_p95_abs']:.3%}, "
            f"max={validation['bn_over_b_max_abs']:.3%}"
        )
    field.save(FIELD_FILE)
    print(
        f"[fit] selected lambda={REG_LAMBDAS[selected_idx]:.1e}; "
        f"independent p95={validation['bn_over_b_p95_abs']:.3%}, "
        f"max={validation['bn_over_b_max_abs']:.3%}",
        flush=True,
    )

    phi_wall, rout, zout, rwall, standoff = direct_lcfs_standoff(eq)
    wall_coords = np.column_stack(
        [np.full(NW, rwall), phi_wall, np.zeros(NW)]
    )
    wall_source = boundary_grid(FIELD_SOURCE_M, FIELD_SOURCE_N)
    b_wall_sheet = np.asarray(
        field.compute_magnetic_field(
            wall_coords, source_grid=wall_source, chunk_size=32
        )
    )

    coil_rows = []
    selected_coils = None
    print("[fit] filament discretization sweep (5/10/20 contours per half-period)", flush=True)
    for count in (5, 10, 20):
        started_coils = time.time()
        coils = field.to_CoilSet(
            num_coils=count,
            npts=192,
            stell_sym=True,
            check_intersection=False,
        )
        b_wall_coils = np.asarray(
            coils.compute_magnetic_field(
                wall_coords, source_grid=192, chunk_size=32
            )
        )
        denom = np.maximum(np.linalg.norm(b_wall_sheet, axis=1), 1e-12)
        rel = np.linalg.norm(b_wall_coils - b_wall_sheet, axis=1) / denom
        row = {
            "contours_per_half_period": count,
            "represented_full_torus_coils": int(2 * NFP * count),
            "filament_current_A": float(abs(list(coils)[0].current)),
            "wall_vector_error_mean": float(np.mean(rel)),
            "wall_vector_error_p95": float(np.quantile(rel, 0.95)),
            "wall_vector_error_max": float(np.max(rel)),
            "build_seconds": float(time.time() - started_coils),
        }
        coil_rows.append(row)
        print(
            f"      {count:2d}: mean wall-vector error {row['wall_vector_error_mean']:.2%}, "
            f"p95 {row['wall_vector_error_p95']:.2%}",
            flush=True,
        )
        if selected_coils is None and row["wall_vector_error_p95"] <= 0.10:
            selected_coils = coils
            selected_count = count
        if count == 20 and selected_coils is None:
            selected_coils = coils
            selected_count = count

    selected_coils.save(COIL_FILE)
    metrics = {
        "status": "W7-X-equilibrium-derived current-potential and filament surrogate",
        "not_claimed": "published W7-X engineering coil geometry or a new free-boundary equilibrium",
        "desc_version": "0.17.2",
        "equilibrium": "DESC bundled W7-X fixed-boundary example",
        "winding_surface": winding_metrics,
        "regcoil": {
            "regularization_type": "regcoil",
            "includes_virtual_casing_plasma_field": True,
            "current_helicity": [1, 0],
            "M_Phi": CURRENT_M,
            "N_Phi": CURRENT_N,
            "G_A": float(field.G),
            "I_A": float(field.I),
            "elapsed_seconds": elapsed,
            "selected_index": selected_idx,
            "selected_lambda": float(REG_LAMBDAS[selected_idx]),
            "coarse_selected_index": coarse_selected_idx,
            "sweep": sweep_rows,
            "independent_validation": validation,
        },
        "filament_discretization": {
            "selected_contours_per_half_period": selected_count,
            "selected_full_torus_coils": int(2 * NFP * selected_count),
            "use_for_physics": "continuous current-potential field; filaments are a geometry/fidelity check",
            "sweep": coil_rows,
        },
        "baseline_wall": {
            "R_wall_m": rwall,
            "standoff_min_m": float(standoff.min()),
            "standoff_max_m": float(standoff.max()),
            "standoff_mean_m": float(standoff.mean()),
        },
        "artifacts": {
            "current_potential": FIELD_FILE,
            "filament_coils": COIL_FILE,
            "winding_surface": CASE_DIR / "w7x_winding_surface.h5",
        },
    }
    write_json(FIT_METRICS_FILE, metrics)
    np.savez_compressed(
        CASE_DIR / "a2_baseline_wall_geometry.npz",
        phi=phi_wall,
        R_out=rout,
        Z_out=zout,
        R_wall=np.array([rwall]),
        standoff=standoff,
    )
    print(f"[fit] wrote {FIT_METRICS_FILE}", flush=True)
    return metrics


def virtual_casing_sources(eq, m: int = VC_SOURCE_M, n: int = VC_SOURCE_N):
    """Discretize the LCFS virtual-casing sheet current over the full torus."""
    grid = boundary_grid(m, n)
    data = eq.compute(
        ["R", "phi", "Z", "K_vc", "|e_theta x e_zeta|"], grid=grid
    )
    ht = 2 * np.pi / grid.num_theta
    hz = 2 * np.pi / (grid.num_zeta * NFP)
    xyz_parts = []
    kxyz_parts = []
    da_parts = []
    for period in range(NFP):
        phi = np.asarray(data["phi"]) + period * 2 * np.pi / NFP
        rpz = np.column_stack([data["R"], phi, data["Z"]])
        xyz_parts.append(np.asarray(rpz2xyz(rpz)))
        kxyz_parts.append(
            np.asarray(rpz2xyz_vec(np.asarray(data["K_vc"]), phi=phi))
        )
        da_parts.append(np.asarray(data["|e_theta x e_zeta|"]) * ht * hz)
    return {
        "xyz": np.concatenate(xyz_parts),
        "K_xyz": np.concatenate(kxyz_parts),
        "dA": np.concatenate(da_parts),
        "num_theta": grid.num_theta,
        "num_zeta_per_period": grid.num_zeta,
    }


def plasma_field_direct(coords_rpz: np.ndarray, sources: dict, chunk: int = 64):
    """External plasma field from the virtual-casing sheet (ordinary quadrature).

    Evaluation points are in the SOL/wall and therefore do not require DESC's
    singular on-surface quadrature or the B/2 jump term.
    """
    coords_rpz = np.atleast_2d(np.asarray(coords_rpz, dtype=float))
    eval_xyz = np.asarray(rpz2xyz(coords_rpz))
    src_xyz = sources["xyz"]
    k_da = sources["K_xyz"] * sources["dA"][:, None]
    out_xyz = np.empty_like(eval_xyz)
    prefactor = 1e-7  # mu0/(4*pi)
    for i0 in range(0, len(eval_xyz), chunk):
        i1 = min(i0 + chunk, len(eval_xyz))
        dx = eval_xyz[i0:i1, None, :] - src_xyz[None, :, :]
        r2 = np.einsum("ijk,ijk->ij", dx, dx)
        inv_r3 = np.maximum(r2, 1e-16) ** -1.5
        cross = np.cross(k_da[None, :, :], dx, axis=-1)
        out_xyz[i0:i1] = prefactor * np.sum(cross * inv_r3[:, :, None], axis=1)
    return np.asarray(xyz2rpz_vec(out_xyz, phi=coords_rpz[:, 1]))


def total_field_direct(
    coords_rpz: np.ndarray,
    current_field,
    current_source_grid,
    vc_sources: dict,
):
    coil = np.asarray(
        current_field.compute_magnetic_field(
            coords_rpz,
            source_grid=current_source_grid,
            chunk_size=64,
        )
    )
    plasma = plasma_field_direct(coords_rpz, vc_sources, chunk=48)
    return coil + plasma, coil, plasma


def spline_stage() -> dict:
    """Audit the direct trace field and certify vacuum+plasma values at the wall.

    The stage name is retained for restart compatibility.  A rectangular field spline
    was tested and rejected because its box crossed the REGCOIL sheet-current surface.
    All trajectories therefore use the continuous Biot-Savart field directly.
    """
    ensure_dirs()
    if not FIELD_FILE.exists():
        raise FileNotFoundError(f"Run --stage fit first; missing {FIELD_FILE}")
    eq = get("W7-X")
    field = FourierCurrentPotentialField.load(FIELD_FILE)
    wall = np.load(CASE_DIR / "a2_baseline_wall_geometry.npz")
    rwall = float(wall["R_wall"][0])

    started = time.time()
    source_grid = boundary_grid(FIELD_SOURCE_M, FIELD_SOURCE_N)

    phi_wall = wall["phi"]
    wall_coords = np.column_stack(
        [np.full_like(phi_wall, rwall), phi_wall, np.zeros_like(phi_wall)]
    )
    bcoil_wall = np.asarray(
        field.compute_magnetic_field(
            wall_coords, source_grid=source_grid, chunk_size=32
        )
    )
    # Near the 3 cm closest approach, ordinary virtual-casing quadrature needs a
    # much finer source grid than the smooth external field.  Certify the wall
    # correction separately instead of contaminating a rectangular volume spline.
    vc_convergence = []
    previous_plasma = None
    bplasma_wall = None
    for level in VC_WALL_SOURCE_LEVELS:
        print(f"[spline] wall plasma-field quadrature M=N={level}", flush=True)
        vc_sources = virtual_casing_sources(eq, level, level)
        candidate = plasma_field_direct(wall_coords, vc_sources, chunk=2)
        row = {
            "M_N": level,
            "full_torus_sources": int(len(vc_sources["xyz"])),
            "magnitude_min_T": float(np.linalg.norm(candidate, axis=1).min()),
            "magnitude_max_T": float(np.linalg.norm(candidate, axis=1).max()),
        }
        if previous_plasma is not None:
            delta = np.linalg.norm(candidate - previous_plasma, axis=1)
            row["delta_from_previous_p95_T"] = float(np.quantile(delta, 0.95))
            row["delta_from_previous_max_T"] = float(np.max(delta))
        vc_convergence.append(row)
        previous_plasma = candidate
        bplasma_wall = candidate
        del vc_sources
    btot_wall = bcoil_wall + bplasma_wall
    btot_mag = np.linalg.norm(btot_wall, axis=1)
    bplasma_mag = np.linalg.norm(bplasma_wall, axis=1)
    incidence = np.abs(btot_wall[:, 0]) / np.maximum(btot_mag, 1e-12)

    # Bound the plasma term on the LCFS from the exact equilibrium field minus
    # the fitted external field.  This is the sensitivity carried by direct
    # external-field trajectories.
    audit_grid = boundary_grid(VALID_GRID_M, VALID_GRID_N)
    audit_data = eq.compute(["B", "|B|", "R", "phi", "Z"], grid=audit_grid)
    audit_coords = np.column_stack(
        [audit_data["R"], audit_data["phi"], audit_data["Z"]]
    )
    boundary_external = np.asarray(
        field.compute_magnetic_field(
            audit_coords, source_grid=source_grid, chunk_size=32
        )
    )
    boundary_plasma = np.asarray(audit_data["B"]) - boundary_external
    boundary_ratio = np.linalg.norm(boundary_plasma, axis=1) / np.asarray(
        audit_data["|B|"]
    )
    trace_source_grid = boundary_grid(TRACE_SOURCE_M, TRACE_SOURCE_N)
    boundary_trace = np.asarray(
        field.compute_magnetic_field(
            audit_coords, source_grid=trace_source_grid, chunk_size=32
        )
    )
    wall_trace = np.asarray(
        field.compute_magnetic_field(
            wall_coords, source_grid=trace_source_grid, chunk_size=32
        )
    )
    boundary_trace_error = np.linalg.norm(
        boundary_trace - boundary_external, axis=1
    ) / np.maximum(np.linalg.norm(boundary_external, axis=1), 1e-12)
    wall_trace_error = np.linalg.norm(wall_trace - bcoil_wall, axis=1) / np.maximum(
        np.linalg.norm(bcoil_wall, axis=1), 1e-12
    )
    np.savetxt(
        RESULT_DIR / "a2_wall_magnetic_field.csv",
        np.column_stack(
            [
                phi_wall,
                wall["standoff"],
                btot_wall,
                bcoil_wall,
                bplasma_wall,
                incidence,
            ]
        ),
        delimiter=",",
        header=(
            "phi_rad,baseline_standoff_m,Btotal_R_T,Btotal_phi_T,Btotal_Z_T,"
            "Bcoil_R_T,Bcoil_phi_T,Bcoil_Z_T,Bplasma_R_T,Bplasma_phi_T,"
            "Bplasma_Z_T,sin_incidence"
        ),
        comments="",
    )
    metrics = {
        "stage_name_compatibility": (
            "called spline stage, but the tested rectangular spline was rejected"
        ),
        "trace_field": "continuous external REGCOIL current-potential Biot-Savart field",
        "rejected_method": (
            "rectangular R-phi-Z spline crossed the winding sheet and produced an "
            "81% p95 interpolation error"
        ),
        "plasma_field_treatment": (
            "included at wall with converged virtual-casing quadrature; omitted from "
            "trajectories and bounded by the LCFS equilibrium-minus-external audit"
        ),
        "wall_virtual_casing_convergence": vc_convergence,
        "elapsed_seconds": float(time.time() - started),
        "wall_field": {
            "Btotal_min_T": float(btot_mag.min()),
            "Btotal_max_T": float(btot_mag.max()),
            "plasma_fraction_mean": float(np.mean(bplasma_mag / btot_mag)),
            "plasma_fraction_max": float(np.max(bplasma_mag / btot_mag)),
            "sin_incidence_mean": float(incidence.mean()),
            "sin_incidence_max": float(incidence.max()),
        },
        "trajectory_plasma_sensitivity_bound": {
            "LCFS_full_vector_over_B_rms": float(
                np.sqrt(np.mean(boundary_ratio**2))
            ),
            "LCFS_full_vector_over_B_max": float(np.max(boundary_ratio)),
            "LCFS_full_vector_over_B_p95": float(np.quantile(boundary_ratio, 0.95)),
            "wall_plasma_over_total_max": float(np.max(bplasma_mag / btot_mag)),
        },
        "trace_source_quadrature": {
            "M_N": TRACE_SOURCE_M,
            "reference_M_N": FIELD_SOURCE_M,
            "LCFS_relative_vector_error_p95": float(
                np.quantile(boundary_trace_error, 0.95)
            ),
            "LCFS_relative_vector_error_max": float(np.max(boundary_trace_error)),
            "wall_relative_vector_error_p95": float(
                np.quantile(wall_trace_error, 0.95)
            ),
            "wall_relative_vector_error_max": float(np.max(wall_trace_error)),
        },
        "artifact": FIELD_FILE,
    }
    write_json(SPLINE_METRICS_FILE, metrics)
    print(
        f"[spline] direct trace field retained; wall |B|={btot_mag.min():.2f}-"
        f"{btot_mag.max():.2f} T; LCFS plasma correction p95="
        f"{metrics['trajectory_plasma_sensitivity_bound']['LCFS_full_vector_over_B_p95']:.2%}",
        flush=True,
    )
    return metrics


def _surface_plane(eq, phi: float, ntheta: int = 96):
    theta = np.linspace(0.0, 2 * np.pi, ntheta, endpoint=False)
    nodes = np.column_stack(
        [np.ones_like(theta), theta, np.full_like(theta, phi)]
    )
    data = eq.surface.compute(["R", "Z", "n_rho"], grid=Grid(nodes, sort=False))
    return theta, np.asarray(data["R"]), np.asarray(data["Z"]), np.asarray(data["n_rho"])


def trace_boundary(eq, field, source_grid, rwall: float):
    """Reconstruct the outboard boundary by transporting one phi=0 surface ring."""
    theta, r0, z0, _ = _surface_plane(eq, 0.0, ntheta=64)
    phis = np.linspace(0.0, 2 * np.pi / NFP, NW, endpoint=False)
    r, z = field_line_integrate(
        r0,
        z0,
        phis,
        field,
        source_grid=source_grid,
        rtol=2e-6,
        atol=2e-6,
        max_steps=12000,
        bounds_R=(3.70, rwall + 0.02),
        bounds_Z=(-1.40, 1.40),
        chunk_size=16,
        bs_chunk_size=32,
    )
    r = np.asarray(r)
    z = np.asarray(z)
    if r.shape[0] == len(phis):
        # DESC returns (nphi,nline); use (nline,nphi) below.
        r = r.T
        z = z.T
    if r.shape != (len(theta), len(phis)):
        raise RuntimeError(f"Unexpected field-line shape {r.shape}")
    valid = np.isfinite(r)
    if not np.all(valid):
        raise RuntimeError("A nominal LCFS seed left the direct-field audit bounds")
    idx = np.argmax(r, axis=0)
    rout = r[idx, np.arange(len(phis))]
    zout = z[idx, np.arange(len(phis))]
    return phis, r, z, rout, zout


def path_length(r: np.ndarray, z: np.ndarray, phi: np.ndarray) -> float:
    finite = np.isfinite(r) & np.isfinite(z)
    if np.count_nonzero(finite) < 2:
        return 0.0
    last = np.flatnonzero(finite)[-1]
    r = r[: last + 1]
    z = z[: last + 1]
    phi = phi[: last + 1]
    dr = np.diff(r)
    dz = np.diff(z)
    dphi = np.diff(phi)
    rmid = 0.5 * (r[:-1] + r[1:])
    return float(np.sum(np.sqrt(dr**2 + dz**2 + (rmid * dphi) ** 2)))


def trace_seed_group(
    field,
    source_grid,
    r0: np.ndarray,
    z0: np.ndarray,
    phi0: float,
    rwall: float,
    direction: int,
):
    max_turns = 10
    dphi = direction * 0.10
    phis = phi0 + np.arange(int(max_turns * 2 * np.pi / abs(dphi)) + 1) * dphi
    r, z = field_line_integrate(
        np.asarray(r0),
        np.asarray(z0),
        phis,
        field,
        source_grid=source_grid,
        rtol=2e-6,
        atol=2e-6,
        max_steps=250000,
        bounds_R=(0.1, rwall),
        bounds_Z=(-np.inf, np.inf),
        chunk_size=8,
        bs_chunk_size=32,
    )
    r = np.asarray(r)
    z = np.asarray(z)
    if r.shape[0] == len(phis):
        r = r.T
        z = z.T
    rows = []
    for ri, zi in zip(r, z):
        finite = np.isfinite(ri) & np.isfinite(zi)
        nfinite = int(np.count_nonzero(finite))
        hit = nfinite < len(phis)
        rows.append(
            {
                "hit_wall": bool(hit),
                "length_m": path_length(ri, zi, phis),
                "turns_traced": float(
                    abs(phis[max(nfinite - 1, 0)] - phi0) / (2 * np.pi)
                ),
                "samples": nfinite,
            }
        )
    return rows


def derive_lambda(connection_lengths: np.ndarray):
    te_values = np.array([100.0, 150.0, 200.0])
    chi_values = np.array([1.0, 2.0, 3.0])
    rows = []
    for length in connection_lengths:
        for te in te_values:
            cs = math.sqrt(2 * te * EV_J / M_DT_KG)
            for chi in chi_values:
                lam = math.sqrt(chi * length / cs)
                rows.append((length, te, chi, cs, lam))
    arr = np.asarray(rows)
    nominal_length = float(np.median(connection_lengths))
    nominal_cs = math.sqrt(2 * 150.0 * EV_J / M_DT_KG)
    nominal_lambda = math.sqrt(2.0 * nominal_length / nominal_cs)
    return arr, {
        "nominal_m": nominal_lambda,
        "band_min_m": float(arr[:, 4].min()),
        "band_max_m": float(arr[:, 4].max()),
        "median_m": float(np.median(arr[:, 4])),
        "nominal_inputs": {
            "connection_length_m": nominal_length,
            "Te_eV": 150.0,
            "chi_perp_m2_s": 2.0,
            "sound_speed_m_s": nominal_cs,
        },
    }


def periodic_interp(x_old, y_old, x_new, period):
    x = np.r_[x_old, x_old[0] + period]
    y = np.r_[y_old, y_old[0]]
    return interp1d(x, y, kind="cubic")(np.mod(x_new, period))


def feed_load_model(phi, standoff, incidence, lambda_m):
    nearfield_path = FIG_DIR / "plasma_qmap_v2_nearfield.json"
    if not nearfield_path.exists():
        nearfield_path = ROOT / "proposal" / "05_DATA" / "plasma_qmap_v2_nearfield.json"
    baseline = json.loads(nearfield_path.read_text(encoding="utf-8"))
    qrad = np.asarray(baseline["q_rad_rel"], dtype=float)
    if len(qrad) != len(phi):
        source_phi = np.linspace(0.0, 2 * np.pi / NFP, len(qrad), endpoint=False)
        qrad = periodic_interp(source_phi, qrad, phi, 2 * np.pi / NFP)
    qrad /= qrad.mean()
    gperp = np.exp(-(standoff - standoff.min()) / lambda_m)
    qpar = incidence * gperp
    nrm = lambda a: np.asarray(a) / np.mean(a)
    qtransport = nrm(gperp) * (1 - incidence.mean()) + nrm(qpar) * incidence.mean()
    qtotal = QWALL_W_M2 * (F_RAD * qrad + (1 - F_RAD) * nrm(qtransport))
    return qrad, nrm(qtransport), qtotal


def trace_stage() -> dict:
    ensure_dirs()
    if not FIELD_FILE.exists() or not SPLINE_METRICS_FILE.exists():
        raise FileNotFoundError("Run --stage fit and --stage spline before tracing")
    eq = get("W7-X")
    field = FourierCurrentPotentialField.load(FIELD_FILE)
    source_grid = boundary_grid(TRACE_SOURCE_M, TRACE_SOURCE_N)
    wall = np.load(CASE_DIR / "a2_baseline_wall_geometry.npz")
    rwall = float(wall["R_wall"][0])

    print(
        "[trace] reconstructing the LCFS envelope with direct Biot-Savart integration",
        flush=True,
    )
    phi, ring_r, ring_z, traced_rout, traced_zout = trace_boundary(
        eq, field, source_grid, rwall
    )
    baseline_standoff = np.asarray(wall["standoff"])
    traced_standoff = rwall - traced_rout
    standoff_delta = traced_standoff - baseline_standoff

    wall_field = np.loadtxt(
        RESULT_DIR / "a2_wall_magnetic_field.csv", delimiter=",", skiprows=1
    )
    incidence = wall_field[:, -1]

    connection_rows = []
    print(
        "[trace] tracing two 16-point conformal SOL launch rings in both directions",
        flush=True,
    )
    seed_number = 0
    for plane_index, phi0 in enumerate((0.0, np.pi / NFP)):
        theta, rs, zs, normals = _surface_plane(eq, phi0, ntheta=16)
        normals_rz = normals[:, [0, 2]].copy()
        normals_rz /= np.linalg.norm(normals_rz, axis=1)[:, None]
        r0 = rs + 0.010 * normals_rz[:, 0]
        z0 = zs + 0.010 * normals_rz[:, 1]
        forward_rows = trace_seed_group(
            field, source_grid, r0, z0, phi0, rwall, +1
        )
        backward_rows = trace_seed_group(
            field, source_grid, r0, z0, phi0, rwall, -1
        )
        for local, (forward, backward) in enumerate(
            zip(forward_rows, backward_rows)
        ):
            complete = forward["hit_wall"] and backward["hit_wall"]
            total_length = forward["length_m"] + backward["length_m"]
            row = {
                "seed": seed_number,
                "launch_plane": plane_index,
                "theta_rad": float(theta[local]),
                "phi_rad": float(phi0),
                "R0_m": float(r0[local]),
                "Z0_m": float(z0[local]),
                "forward_hit": forward["hit_wall"],
                "forward_length_m": forward["length_m"],
                "forward_turns": forward["turns_traced"],
                "backward_hit": backward["hit_wall"],
                "backward_length_m": backward["length_m"],
                "backward_turns": backward["turns_traced"],
                "complete_connection": complete,
                "connection_length_m": total_length,
            }
            connection_rows.append(row)
            print(
                f"      seed {seed_number:02d}: L={total_length:8.1f} m "
                f"({'two-wall' if complete else 'censored'})",
                flush=True,
            )
            seed_number += 1

    complete_lengths = np.asarray(
        [r["connection_length_m"] for r in connection_rows if r["complete_connection"]]
    )
    all_lengths = np.asarray([r["connection_length_m"] for r in connection_rows])
    censored_lengths = np.asarray(
        [r["connection_length_m"] for r in connection_rows if not r["complete_connection"]]
    )
    complete_fraction = float(np.mean([r["complete_connection"] for r in connection_rows]))
    if len(complete_lengths) >= 4:
        lambda_lengths = complete_lengths
        length_status = (
            "conditional finite-connection estimate; censored majority reported separately"
        )
    else:
        lambda_lengths = all_lengths
        length_status = "censored lower-bound connections; lambda is a lower-bound band"
    lambda_grid, lambda_stats = derive_lambda(lambda_lengths)
    if len(censored_lengths):
        censored_lambda_grid, censored_lambda_stats = derive_lambda(censored_lengths)
    else:
        censored_lambda_grid = np.empty((0, 5))
        censored_lambda_stats = None

    qrad, qtransport, qtotal = feed_load_model(
        phi, traced_standoff, incidence, lambda_stats["nominal_m"]
    )
    if censored_lambda_stats is not None:
        _, qtransport_censored, qtotal_censored = feed_load_model(
            phi,
            traced_standoff,
            incidence,
            censored_lambda_stats["nominal_m"],
        )
    else:
        qtransport_censored = np.full_like(qtransport, np.nan)
        qtotal_censored = np.full_like(qtotal, np.nan)
    baseline_lambda = 0.03
    # Keep the envelope method fixed across the lambda sweep. The original A2
    # artifact used direct geometry only for the 3 cm row, which confounded a
    # width change with an envelope-method change.
    _, baseline_transport, baseline_total = feed_load_model(
        phi, traced_standoff, incidence, baseline_lambda
    )

    np.savetxt(
        RESULT_DIR / "a2_connection_lengths.csv",
        np.asarray(
            [
                [
                    r["seed"],
                    r["launch_plane"],
                    r["theta_rad"],
                    r["phi_rad"],
                    r["R0_m"],
                    r["Z0_m"],
                    int(r["forward_hit"]),
                    r["forward_length_m"],
                    r["forward_turns"],
                    int(r["backward_hit"]),
                    r["backward_length_m"],
                    r["backward_turns"],
                    int(r["complete_connection"]),
                    r["connection_length_m"],
                ]
                for r in connection_rows
            ]
        ),
        delimiter=",",
        header=(
            "seed,launch_plane,theta_rad,phi_rad,R0_m,Z0_m,forward_hit,forward_length_m,forward_turns,"
            "backward_hit,backward_length_m,backward_turns,complete_connection,"
            "connection_length_m"
        ),
        comments="",
    )
    np.savetxt(
        RESULT_DIR / "a2_lambda_transport_band.csv",
        np.vstack(
            [
                np.column_stack([np.zeros(len(lambda_grid)), lambda_grid]),
                np.column_stack(
                    [np.ones(len(censored_lambda_grid)), censored_lambda_grid]
                ),
            ]
        ),
        delimiter=",",
        header=(
            "population_code,connection_length_m,Te_eV,chi_perp_m2_s,"
            "sound_speed_m_s,lambda_far_m"
        ),
        comments="",
    )
    np.savetxt(
        RESULT_DIR / "a2_standoff_load_feedback.csv",
        np.column_stack(
            [
                phi,
                baseline_standoff,
                traced_standoff,
                incidence,
                qrad,
                baseline_transport,
                qtransport,
                qtransport_censored,
                baseline_total,
                qtotal,
                qtotal_censored,
            ]
        ),
        delimiter=",",
        header=(
            "phi_rad,baseline_standoff_m,coil_traced_standoff_m,sin_incidence,"
            "qrad_rel,traced_3cm_transport_rel,conditional_transport_rel,"
            "censored_lower_transport_rel,traced_3cm_qtotal_W_m2,"
            "conditional_qtotal_W_m2,censored_lower_qtotal_W_m2"
        ),
        comments="",
    )

    metrics = {
        "method_status": (
            "direct external current-potential field-line trace; converged virtual-casing "
            "plasma correction included at the wall and bounded for trajectories; reduced "
            "diffusive SOL closure"
        ),
        "standoff": {
            "baseline_min_max_m": [float(baseline_standoff.min()), float(baseline_standoff.max())],
            "coil_traced_min_max_m": [float(traced_standoff.min()), float(traced_standoff.max())],
            "trace_minus_boundary_rms_m": float(np.sqrt(np.mean(standoff_delta**2))),
            "trace_minus_boundary_max_abs_m": float(np.max(np.abs(standoff_delta))),
        },
        "connection_length": {
            "seed_offset_m": 0.010,
            "seeds": len(connection_rows),
            "complete_fraction": complete_fraction,
            "status": length_status,
            "complete_count": int(len(complete_lengths)),
            "censored_count": int(len(censored_lengths)),
            "used_min_m": float(lambda_lengths.min()),
            "used_median_m": float(np.median(lambda_lengths)),
            "used_max_m": float(lambda_lengths.max()),
            "trace_limit_turns_each_direction": 10,
            "launch_planes_rad": [0.0, float(np.pi / NFP)],
        },
        "lambda_far": {
            **lambda_stats,
            "interpretation": "conditional on the finite two-wall subset",
            "censored_majority_lower_bound": (
                {
                    **censored_lambda_stats,
                    "interpretation": (
                        "nominal and band are lower bounds because these trajectories "
                        "did not hit within the trace limit"
                    ),
                }
                if censored_lambda_stats is not None
                else None
            ),
            "closure": "sqrt(chi_perp * L_connection / c_s)",
            "assumption_band": {"Te_eV": [100, 200], "chi_perp_m2_s": [1, 3]},
            "old_assumed_nominal_m": 0.03,
            "old_assumed_band_m": [0.01, 0.05],
            "references": [
                "Lion et al., Stellaris, Fusion Engineering and Design 214 (2025) 114868",
                "Mau et al., Fusion Science and Technology 54 (2008) 771-786, doi:10.13182/FST08-27",
            ],
        },
        "load_feedback": {
            "comparison_method": (
                "method-consistent direct external-field envelope at 3.0, 5.9, "
                "and >=11.6 cm"
            ),
            "baseline_peak_over_mean": float(baseline_total.max() / baseline_total.mean()),
            "conditional_peak_over_mean": float(qtotal.max() / qtotal.mean()),
            "conditional_peak_change_percent": float(
                100
                * (
                    (qtotal.max() / qtotal.mean())
                    / (baseline_total.max() / baseline_total.mean())
                    - 1
                )
            ),
            "baseline_peak_MW_m2": float(baseline_total.max() / 1e6),
            "conditional_peak_MW_m2": float(qtotal.max() / 1e6),
            "censored_lower_peak_over_mean": (
                float(qtotal_censored.max() / qtotal_censored.mean())
                if censored_lambda_stats is not None
                else None
            ),
            "censored_lower_peak_MW_m2": (
                float(qtotal_censored.max() / 1e6)
                if censored_lambda_stats is not None
                else None
            ),
        },
        "limitations": [
            "No published W7-X engineering coil file was available locally; coils are REGCOIL-derived.",
            "The starting equilibrium is fixed-boundary; this is a supported-boundary field reconstruction, not a new equilibrium optimization.",
            "Trajectories use the external current-potential field; the omitted plasma contribution is <=2.1% of |B| on the LCFS and <0.5% at the wall.",
            "Connection length alone does not determine lambda_far; Te and chi_perp remain explicit transport inputs.",
            "The first wall remains the project's cylindrical comparison wall, not stellarator CAD.",
        ],
    }
    write_json(RESULT_DIR / "a2_metrics.json", metrics)
    make_summary_figure(
        eq,
        phi,
        ring_r,
        ring_z,
        rwall,
        baseline_standoff,
        traced_standoff,
        connection_rows,
        lambda_grid,
        censored_lambda_grid,
        baseline_total,
        qtotal,
        qtotal_censored,
        metrics,
    )
    print(
        f"[trace] lambda_far={lambda_stats['nominal_m']*100:.1f} cm "
        f"({lambda_stats['band_min_m']*100:.1f}-{lambda_stats['band_max_m']*100:.1f} cm band); "
        f"conditional load pk/mean {metrics['load_feedback']['baseline_peak_over_mean']:.2f} -> "
        f"{metrics['load_feedback']['conditional_peak_over_mean']:.2f}; "
        f"censored-majority lower-width bound -> "
        f"{metrics['load_feedback']['censored_lower_peak_over_mean']:.2f}",
        flush=True,
    )
    return metrics


def make_summary_figure(
    eq,
    phi,
    ring_r,
    ring_z,
    rwall,
    baseline_standoff,
    traced_standoff,
    connection_rows,
    lambda_grid,
    censored_lambda_grid,
    baseline_total,
    qtotal,
    qtotal_censored,
    metrics,
):
    coils = CoilSet.load(COIL_FILE)
    coil_xyz = np.asarray(coils._compute_position(grid=120, basis="xyz"))
    surface_grid = boundary_grid(32, 32)
    surface_data = eq.surface.compute(["R", "phi", "Z"], grid=surface_grid)
    surface_xyz = np.asarray(
        rpz2xyz(
            np.column_stack(
                [surface_data["R"], surface_data["phi"], surface_data["Z"]]
            )
        )
    )

    fig = plt.figure(figsize=(14.2, 9.4), facecolor="white")
    ax3 = fig.add_subplot(2, 2, 1, projection="3d")
    for curve in coil_xyz:
        ax3.plot(curve[:, 0], curve[:, 1], curve[:, 2], color="#b24a3a", lw=0.55, alpha=0.45)
    ax3.scatter(
        surface_xyz[::4, 0],
        surface_xyz[::4, 1],
        surface_xyz[::4, 2],
        s=0.8,
        color="#247a74",
        alpha=0.35,
    )
    ax3.set_title("W7-X equilibrium-derived coil surrogate", loc="left", fontsize=11)
    ax3.set_box_aspect((1, 1, 0.34))
    ax3.view_init(elev=24, azim=36)
    ax3.set_axis_off()

    ax = fig.add_subplot(2, 2, 2)
    # ``phi`` is endpoint-free over one field period, so its mapped plotting
    # coordinate must use the same periodic sampling convention.
    s = np.linspace(0.0, S_DUCT_M, len(phi), endpoint=False) * 1e3
    ax.plot(s, baseline_standoff * 100, color="0.45", ls="--", lw=1.7, label="fixed-boundary geometry")
    ax.plot(
        s,
        traced_standoff * 100,
        color="#247a74",
        lw=2.2,
        label="direct external-field line envelope",
    )
    ax.set_xlabel("mapped duct coordinate s [mm]")
    ax.set_ylabel("outboard standoff [cm]")
    ax.set_title("Field reconstruction preserves the standoff phase", loc="left", fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(2, 2, 3)
    lengths = np.asarray([r["connection_length_m"] for r in connection_rows])
    complete = np.asarray([r["complete_connection"] for r in connection_rows])
    colors = np.where(complete, "#2d6a9f", "#b24a3a")
    ax.scatter(np.arange(len(lengths)), lengths, c=colors, s=38, zorder=3)
    ax.axhline(np.median(lengths), color="0.35", ls="--", lw=1)
    ax.set_xlabel("SOL launch seed")
    ax.set_ylabel("two-direction traced length [m]")
    ax.set_title("Connection length (red = censored at trace limit)", loc="left", fontsize=11)
    ax.grid(alpha=0.25)
    inset = ax.inset_axes([0.56, 0.52, 0.40, 0.41])
    inset.hist(
        lambda_grid[:, 4] * 100,
        bins=10,
        color="#8b6f3d",
        alpha=0.78,
        label="finite",
    )
    if len(censored_lambda_grid):
        inset.hist(
            censored_lambda_grid[:, 4] * 100,
            bins=10,
            histtype="step",
            color="#b24a3a",
            lw=1.3,
            label="censored LB",
        )
    inset.axvline(3.0, color="0.2", ls=":", lw=1.3)
    inset.set_xlabel(r"$\lambda_{far}$ [cm]", fontsize=7)
    inset.set_ylabel("cases", fontsize=7)
    inset.tick_params(labelsize=7)
    inset.legend(fontsize=6, frameon=False)

    ax = fig.add_subplot(2, 2, 4)
    ax.plot(
        s,
        baseline_total / 1e6,
        color="0.45",
        ls="--",
        lw=1.8,
        label=(
            f"traced envelope, 3 cm: pk/mean "
            f"{metrics['load_feedback']['baseline_peak_over_mean']:.2f}"
        ),
    )
    ax.plot(
        s,
        qtotal / 1e6,
        color="#9b2d30",
        lw=2.3,
        label=(
            f"finite subset, 5.9 cm: pk/mean "
            f"{metrics['load_feedback']['conditional_peak_over_mean']:.2f}"
        ),
    )
    if np.all(np.isfinite(qtotal_censored)):
        ax.plot(
            s,
            qtotal_censored / 1e6,
            color="#247a74",
            lw=2.0,
            label=(
                f"censored majority, >=11.6 cm: pk/mean "
                f"{metrics['load_feedback']['censored_lower_peak_over_mean']:.2f}"
            ),
        )
    ax.set_xlabel("mapped duct coordinate s [mm]")
    ax.set_ylabel(r"surface load $q''$ [MW m$^{-2}$]")
    ax.set_title("A2 feedback into the near-field-corrected load", loc="left", fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    fig.suptitle(
        "A2: magnet-coupled W7-X surrogate and reduced SOL closure",
        x=0.055,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.015,
        "REGCOIL surrogate, not published engineering coils | direct external-field traces; converged plasma correction retained at wall | "
        r"$\lambda_{far}=\sqrt{\chi_\perp L_\parallel/c_s}$ keeps $T_e$ and $\chi_\perp$ explicit",
        fontsize=8.2,
        color="0.35",
    )
    fig.tight_layout(rect=(0.03, 0.045, 0.99, 0.95))
    out = RESULT_DIR / "a2_desc_coil_sol_summary.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[trace] wrote {out}", flush=True)


def refresh_load_feedback_stage() -> dict:
    """Refresh only A2 load feedback after a comparison-method correction."""
    ensure_dirs()
    feedback_path = RESULT_DIR / "a2_standoff_load_feedback.csv"
    metrics_path = RESULT_DIR / "a2_metrics.json"
    connection_path = RESULT_DIR / "a2_connection_lengths.csv"
    lambda_path = RESULT_DIR / "a2_lambda_transport_band.csv"
    required = (feedback_path, metrics_path, connection_path, lambda_path)
    if any(not path.exists() for path in required):
        raise FileNotFoundError("Run --stage trace before --stage refresh")

    old = np.loadtxt(feedback_path, delimiter=",", skiprows=1)
    phi = old[:, 0]
    baseline_standoff = old[:, 1]
    traced_standoff = old[:, 2]
    incidence = old[:, 3]
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    lambda_conditional = float(metrics["lambda_far"]["nominal_m"])
    lambda_censored = float(
        metrics["lambda_far"]["censored_majority_lower_bound"]["nominal_m"]
    )

    qrad, baseline_transport, baseline_total = feed_load_model(
        phi, traced_standoff, incidence, 0.03
    )
    _, _, direct_geometry_3cm_total = feed_load_model(
        phi, baseline_standoff, incidence, 0.03
    )
    _, qtransport, qtotal = feed_load_model(
        phi, traced_standoff, incidence, lambda_conditional
    )
    _, qtransport_censored, qtotal_censored = feed_load_model(
        phi, traced_standoff, incidence, lambda_censored
    )
    np.savetxt(
        feedback_path,
        np.column_stack(
            [
                phi,
                baseline_standoff,
                traced_standoff,
                incidence,
                qrad,
                baseline_transport,
                qtransport,
                qtransport_censored,
                baseline_total,
                qtotal,
                qtotal_censored,
            ]
        ),
        delimiter=",",
        header=(
            "phi_rad,baseline_standoff_m,coil_traced_standoff_m,sin_incidence,"
            "qrad_rel,traced_3cm_transport_rel,conditional_transport_rel,"
            "censored_lower_transport_rel,traced_3cm_qtotal_W_m2,"
            "conditional_qtotal_W_m2,censored_lower_qtotal_W_m2"
        ),
        comments="",
    )

    baseline_peak = float(baseline_total.max() / baseline_total.mean())
    conditional_peak = float(qtotal.max() / qtotal.mean())
    censored_peak = float(qtotal_censored.max() / qtotal_censored.mean())
    direct_geometry_3cm_peak = float(
        direct_geometry_3cm_total.max() / direct_geometry_3cm_total.mean()
    )
    metrics["load_feedback"] = {
        "comparison_method": (
            "method-consistent direct external-field envelope at 3.0, 5.9, "
            "and >=11.6 cm"
        ),
        "baseline_peak_over_mean": baseline_peak,
        "conditional_peak_over_mean": conditional_peak,
        "conditional_peak_change_percent": float(
            100 * (conditional_peak / baseline_peak - 1)
        ),
        "baseline_peak_MW_m2": float(baseline_total.max() / 1e6),
        "conditional_peak_MW_m2": float(qtotal.max() / 1e6),
        "censored_lower_peak_over_mean": censored_peak,
        "censored_lower_peak_MW_m2": float(qtotal_censored.max() / 1e6),
        "censored_lower_peak_change_percent": float(
            100 * (censored_peak / baseline_peak - 1)
        ),
        "direct_geometry_3cm_diagnostic_peak_over_mean": direct_geometry_3cm_peak,
        "envelope_method_spread_3cm_percent": float(
            100 * (baseline_peak / direct_geometry_3cm_peak - 1)
        ),
    }
    write_json(metrics_path, metrics)

    connection_data = np.genfromtxt(connection_path, delimiter=",", names=True)
    connection_rows = [
        {
            "connection_length_m": float(row["connection_length_m"]),
            "complete_connection": bool(row["complete_connection"]),
        }
        for row in np.atleast_1d(connection_data)
    ]
    lambda_data = np.loadtxt(lambda_path, delimiter=",", skiprows=1)
    lambda_grid = lambda_data[lambda_data[:, 0] == 0, 1:]
    censored_lambda_grid = lambda_data[lambda_data[:, 0] == 1, 1:]
    eq = get("W7-X")
    wall = np.load(CASE_DIR / "a2_baseline_wall_geometry.npz")
    rwall = float(wall["R_wall"][0])
    make_summary_figure(
        eq,
        phi,
        None,
        None,
        rwall,
        baseline_standoff,
        traced_standoff,
        connection_rows,
        lambda_grid,
        censored_lambda_grid,
        baseline_total,
        qtotal,
        qtotal_censored,
        metrics,
    )
    advisor = ROOT / "results" / "advisor_packet"
    advisor.mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        RESULT_DIR / "a2_desc_coil_sol_summary.png",
        advisor / "105_a2_desc_coil_sol_summary.png",
    )
    print(
        "[refresh] method-consistent peak/mean: "
        f"{baseline_peak:.3f}, {conditional_peak:.3f}, {censored_peak:.3f}",
        flush=True,
    )
    return metrics


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage", choices=("fit", "spline", "trace", "refresh", "all"), default="all"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.stage in ("fit", "all"):
        fit_stage()
    if args.stage in ("spline", "all"):
        spline_stage()
    if args.stage in ("trace", "all"):
        trace_stage()
    if args.stage == "refresh":
        refresh_load_feedback_stage()


if __name__ == "__main__":
    main()

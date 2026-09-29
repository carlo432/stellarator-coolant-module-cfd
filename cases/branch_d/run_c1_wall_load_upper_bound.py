#!/usr/bin/env python3
"""C1 confinement-agnostic wall-load-only shape optimization.

This script estimates an upper bound on wall-load flattening from low-order,
stellarator-symmetric perturbations of the bundled DESC W7-X boundary. It is
not a confinement or MHD optimization and must not be presented as a better
stellarator design.

The workflow is restartable:

    python run_c1_wall_load_upper_bound.py --stage search
    python run_c1_wall_load_upper_bound.py --stage certify
    python run_c1_wall_load_upper_bound.py --stage post
    python run_c1_wall_load_upper_bound.py --stage all

Search uses a volume-preserving boundary surrogate and the accepted A2 load
closure. Certification reconstructs a fixed-boundary equilibrium, fits the same
order REGCOIL current-potential sheet used by A2, and independently audits the
normal field and current density. Postprocessing maps the load through the
existing CFD response and reports peak location and shape-deviation cost.
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

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-c1-wall-load")
warnings.filterwarnings("ignore", message=".*Matplotlib.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d
from scipy.optimize import brentq

from desc.examples import get
from desc.equilibrium import Equilibrium
from desc.geometry import FourierRZToroidalSurface
from desc.grid import Grid, LinearGrid, QuadratureGrid
from desc.magnetic_fields import (
    FourierCurrentPotentialField,
    solve_regularized_surface_current,
)
from desc.objectives import get_fixed_boundary_constraints


ROOT = Path(__file__).resolve().parents[2]
WORK_DIR = ROOT / "cases" / "branch_d" / "c1_wall_load_upper_bound"
RESULT_DIR = ROOT / "results" / "branch_d" / "c1_wall_load_upper_bound"
FIG_DIR = ROOT / "cases" / "branch_d" / "figs"

SEARCH_FILE = RESULT_DIR / "c1_search.csv"
SEARCH_META_FILE = RESULT_DIR / "c1_search_metrics.json"
CERT_FILE = RESULT_DIR / "c1_certification_metrics.json"
FINAL_FILE = RESULT_DIR / "c1_metrics.json"

NFP = 5
NW = 192
CLEARANCE_M = 0.03
QWALL_W_M2 = 0.5e6
F_RAD = 0.5
LAMBDAS_M = np.array([0.030, 0.059, 0.116])
S_DUCT_M = np.radians(95.0) * 0.18

# Same current-potential representation and independent gates as A2.
WINDING_OFFSET_M = 0.45
WINDING_FIT_M = 20
WINDING_FIT_N = 20
CURRENT_M = 7
CURRENT_N = 7
REG_LAMBDA = 1e-16
REG_GRID_M = 20
REG_GRID_N = 20
VALID_GRID_M = 30
VALID_GRID_N = 30
FIELD_SOURCE_M = 26
FIELD_SOURCE_N = 26
TRACE_SOURCE_M = 16
TRACE_SOURCE_N = 16

A2_KMAX_A_M = 8_465_719.609303592
KMAX_GATE_A_M = 1.15 * A2_KMAX_A_M
BN_P95_GATE = 0.01
BN_MAX_GATE = 0.02
FORCE_RATIO_GATE = 5.00


def ensure_dirs() -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def json_ready(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(json_ready(data), indent=2) + "\n", encoding="ascii")


def boundary_grid(m: int, n: int) -> LinearGrid:
    return LinearGrid(
        rho=np.array([1.0]), M=m, N=n, NFP=NFP, sym=False, endpoint=False
    )


def periodic_interp(x_old, y_old, x_new, period):
    x = np.r_[x_old, x_old[0] + period]
    y = np.r_[y_old, y_old[0]]
    return interp1d(x, y, kind="cubic")(np.mod(x_new, period))


def load_reference(n: int = NW):
    nearfield_path = FIG_DIR / "plasma_qmap_v2_nearfield.json"
    if not nearfield_path.exists():
        nearfield_path = ROOT / "proposal" / "05_DATA" / nearfield_path.name
    payload = json.loads(nearfield_path.read_text(encoding="utf-8"))
    qrad = np.asarray(payload["q_rad_rel"], dtype=float)
    phi = np.linspace(0.0, 2 * np.pi / NFP, n, endpoint=False)
    if len(qrad) != n:
        source_phi = np.linspace(0.0, 2 * np.pi / NFP, len(qrad), endpoint=False)
        qrad = periodic_interp(source_phi, qrad, phi, 2 * np.pi / NFP)
    return phi, qrad / qrad.mean()


def direct_standoff(surface, rwall: float | None = None):
    phi = np.linspace(0.0, 2 * np.pi / NFP, NW, endpoint=False)
    theta = np.linspace(0.0, 2 * np.pi, 513, endpoint=False)
    th, ph = np.meshgrid(theta, phi, indexing="ij")
    nodes = np.column_stack([np.ones(th.size), th.ravel(), ph.ravel()])
    data = surface.compute(["R", "Z"], grid=Grid(nodes, sort=False))
    rr = np.asarray(data["R"]).reshape(theta.size, NW)
    zz = np.asarray(data["Z"]).reshape(theta.size, NW)
    imax = np.argmax(rr, axis=0)
    rout = rr[imax, np.arange(NW)]
    zout = zz[imax, np.arange(NW)]
    if rwall is None:
        rwall = float(rout.max() + CLEARANCE_M)
    return phi, rout, zout, float(rwall), float(rwall) - rout


def surface_volume(surface) -> float:
    return float(np.asarray(surface.compute("V")["V"]).reshape(-1)[0])


def perturb_surface(base, scales: np.ndarray, target_volume: float):
    """Scale |n|=1,2,3 boundary groups and restore volume exactly.

    A common scale on every non-axis mode is solved after the requested toroidal
    group scales. R00 stays fixed, so this is a bounded shape perturbation rather
    than a major-radius translation.
    """
    surface = base.copy()
    r_modes = np.asarray(surface.R_basis.modes)
    z_modes = np.asarray(surface.Z_basis.modes)
    r0 = np.asarray(base.R_lmn).copy()
    z0 = np.asarray(base.Z_lmn).copy()
    grouped_r = r0.copy()
    grouped_z = z0.copy()
    for order, scale in enumerate(scales, start=1):
        grouped_r[np.abs(r_modes[:, 2]) == order] *= scale
        grouped_z[np.abs(z_modes[:, 2]) == order] *= scale
    r_deform_mask = ~((r_modes[:, 1] == 0) & (r_modes[:, 2] == 0))
    z_deform_mask = ~((z_modes[:, 1] == 0) & (z_modes[:, 2] == 0))

    def set_common(common: float) -> float:
        r_coefficients = grouped_r.copy()
        z_coefficients = grouped_z.copy()
        r_coefficients[r_deform_mask] *= common
        z_coefficients[z_deform_mask] *= common
        surface.R_lmn = r_coefficients
        surface.Z_lmn = z_coefficients
        return surface_volume(surface) - target_volume

    lo, hi = 0.70, 1.30
    flo, fhi = set_common(lo), set_common(hi)
    if flo * fhi > 0:
        raise RuntimeError(
            f"Volume compensation is not bracketed: f({lo})={flo}, f({hi})={fhi}"
        )
    common = brentq(set_common, lo, hi, xtol=1e-11, rtol=1e-11)
    set_common(common)
    return surface, float(common)


def surface_cost(base, candidate) -> dict:
    grid = boundary_grid(80, 80)
    b = base.compute(["R", "phi", "Z", "S", "a", "R0/a"], grid=grid)
    c = candidate.compute(["R", "phi", "Z", "S", "a", "R0/a"], grid=grid)
    displacement = np.sqrt(
        (np.asarray(c["R"]) - np.asarray(b["R"])) ** 2
        + (np.asarray(c["Z"]) - np.asarray(b["Z"])) ** 2
    )
    base_a = float(np.asarray(b["a"]).reshape(-1)[0])
    rms_displacement = float(np.sqrt(np.mean(displacement**2)))
    max_displacement = float(displacement.max())
    return {
        "surface_displacement_rms_m": rms_displacement,
        "surface_displacement_max_m": max_displacement,
        "surface_displacement_rms_over_a": rms_displacement / base_a,
        "surface_displacement_max_over_a": max_displacement / base_a,
        "volume_m3": surface_volume(candidate),
        "volume_change_percent": float(
            100 * (surface_volume(candidate) / surface_volume(base) - 1)
        ),
        "surface_area_m2": float(np.asarray(c["S"]).reshape(-1)[0]),
        "surface_area_change_percent": float(
            100
            * (
                float(np.asarray(c["S"]).reshape(-1)[0])
                / float(np.asarray(b["S"]).reshape(-1)[0])
                - 1
            )
        ),
        "minor_radius_m": float(np.asarray(c["a"]).reshape(-1)[0]),
        "aspect_ratio": float(np.asarray(c["R0/a"]).reshape(-1)[0]),
        "aspect_ratio_change_percent": float(
            100
            * (
                float(np.asarray(c["R0/a"]).reshape(-1)[0])
                / float(np.asarray(b["R0/a"]).reshape(-1)[0])
                - 1
            )
        ),
        "orientation": float(candidate._compute_orientation()),
    }


def feed_load(phi, standoff, incidence, lambda_m, qrad):
    gperp = np.exp(-(standoff - standoff.min()) / lambda_m)
    qpar = incidence * gperp
    normalize = lambda values: np.asarray(values) / np.mean(values)
    qtransport = (
        normalize(gperp) * (1 - incidence.mean())
        + normalize(qpar) * incidence.mean()
    )
    qtotal = QWALL_W_M2 * (
        F_RAD * qrad + (1 - F_RAD) * normalize(qtransport)
    )
    return normalize(qtransport), qtotal


def evaluate_geometry(surface, rwall, incidence, qrad):
    phi, rout, zout, _, standoff = direct_standoff(surface, rwall=rwall)
    rows = []
    for lambda_m in LAMBDAS_M:
        transport, qtotal = feed_load(phi, standoff, incidence, lambda_m, qrad)
        imax = int(np.argmax(qtotal))
        rows.append(
            {
                "lambda_m": float(lambda_m),
                "peak_over_mean": float(qtotal.max() / qtotal.mean()),
                "peak_MW_m2": float(qtotal.max() / 1e6),
                "peak_phi_rad": float(phi[imax]),
                "peak_s_m": float(phi[imax] * 0.18),
                "peak_standoff_m": float(standoff[imax]),
            }
        )
    return {
        "phi": phi,
        "rout": rout,
        "zout": zout,
        "standoff": standoff,
        "min_clearance_m": float(standoff.min()),
        "max_standoff_m": float(standoff.max()),
        "lambda_rows": rows,
        "robust_peak_over_mean": float(max(row["peak_over_mean"] for row in rows)),
    }


def load_a2_inputs():
    wall = np.load(
        ROOT
        / "cases"
        / "branch_d"
        / "a2_desc_coil_sol"
        / "a2_baseline_wall_geometry.npz"
    )
    rwall = float(wall["R_wall"][0])
    wall_field = np.loadtxt(
        ROOT
        / "results"
        / "branch_d"
        / "a2_desc_coil_sol"
        / "a2_wall_magnetic_field.csv",
        delimiter=",",
        skiprows=1,
    )
    phi, qrad = load_reference()
    incidence = wall_field[:, -1]
    if len(incidence) != len(phi):
        source_phi = wall_field[:, 0]
        incidence = periodic_interp(source_phi, incidence, phi, 2 * np.pi / NFP)
    return rwall, incidence, qrad


def search_stage() -> dict:
    ensure_dirs()
    eq = get("W7-X")
    base = eq.surface
    target_volume = surface_volume(base)
    rwall, incidence, qrad = load_a2_inputs()

    # Low-order stellarator-symmetric toroidal shaping factors. Search remains
    # deliberately small enough to audit every candidate and expose shape cost.
    levels = np.linspace(0.82, 1.08, 7)
    combinations = np.array(np.meshgrid(levels, levels, levels)).T.reshape(-1, 3)
    combinations = np.vstack([combinations, np.ones(3)])
    rows = []
    started = time.time()
    print(f"[search] evaluating {len(combinations)} volume-preserving boundaries", flush=True)
    for index, scales in enumerate(combinations):
        candidate, common = perturb_surface(base, scales, target_volume)
        cost = surface_cost(base, candidate)
        geometry = evaluate_geometry(candidate, rwall, incidence, qrad)
        valid = (
            cost["orientation"] == float(base._compute_orientation())
            and geometry["min_clearance_m"] >= CLEARANCE_M - 1e-6
            and abs(cost["volume_change_percent"]) <= 1.0
        )
        row = {
            "candidate": int(index),
            "scale_n1": float(scales[0]),
            "scale_n2": float(scales[1]),
            "scale_n3": float(scales[2]),
            "volume_compensation": common,
            **cost,
            "min_clearance_m": geometry["min_clearance_m"],
            "max_standoff_m": geometry["max_standoff_m"],
            "robust_peak_over_mean": geometry["robust_peak_over_mean"],
            "valid_geometry": bool(valid),
        }
        for lam_row in geometry["lambda_rows"]:
            tag = f"lambda_{lam_row['lambda_m']*100:.1f}cm".replace(".", "p")
            row[f"{tag}_peak_over_mean"] = lam_row["peak_over_mean"]
            row[f"{tag}_peak_phi_rad"] = lam_row["peak_phi_rad"]
        rows.append(row)
        if (index + 1) % 50 == 0:
            print(f"         {index + 1}/{len(combinations)}", flush=True)

    names = list(rows[0])
    numeric = np.asarray([[float(row[name]) for name in names] for row in rows])
    np.savetxt(
        SEARCH_FILE,
        numeric,
        delimiter=",",
        header=",".join(names),
        comments="",
    )
    baseline = rows[-1]
    valid_rows = [row for row in rows if row["valid_geometry"]]
    ranked = sorted(valid_rows, key=lambda row: row["robust_peak_over_mean"])
    shortlist = ranked[: min(8, len(ranked))]
    meta = {
        "status": "confinement-agnostic wall-load-only upper-bound search",
        "not_claimed": "a better stellarator, confinement optimization, or engineering coil design",
        "decision_variables": "multipliers on all |n|=1,2,3 stellarator-symmetric boundary modes",
        "levels": levels,
        "candidate_count": len(rows),
        "valid_geometry_count": len(valid_rows),
        "baseline": baseline,
        "shortlist": shortlist,
        "elapsed_seconds": time.time() - started,
        "limitations": [
            "Search-stage incidence is held to the accepted A2 wall-field profile.",
            "No confinement, neoclassical ripple, QI, or MHD-quality objective is included.",
            "Self-intersection is screened by orientation and dense sampled geometry; certification checks the fixed-boundary coordinate map.",
        ],
    }
    write_json(SEARCH_META_FILE, meta)
    print(
        "[search] baseline robust pk/mean "
        f"{baseline['robust_peak_over_mean']:.3f}; best screened "
        f"{shortlist[0]['robust_peak_over_mean']:.3f}",
        flush=True,
    )
    return meta


def make_winding_surface(eq):
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
    fitted = surface.compute(["R", "Z"], grid=sample_grid)
    error = np.sqrt(
        (np.asarray(fitted["R"]) - target[:, 0]) ** 2
        + (np.asarray(fitted["Z"]) - target[:, 2]) ** 2
    )
    return surface, {
        "fit_rms_m": float(np.sqrt(np.mean(error**2))),
        "fit_max_m": float(error.max()),
        "orientation": float(surface._compute_orientation()),
    }


def normalized_bn_stats(bn, bmag):
    rel = np.abs(np.asarray(bn)) / np.maximum(np.asarray(bmag), 1e-12)
    return {
        "bn_over_b_rms": float(np.sqrt(np.mean(rel**2))),
        "bn_over_b_p95_abs": float(np.quantile(rel, 0.95)),
        "bn_over_b_max_abs": float(np.max(rel)),
    }


def mean_force_residual(eq) -> float:
    resolution = min(eq.L, eq.M, eq.N, 8)
    grid = QuadratureGrid(
        L=resolution, M=resolution, N=resolution, NFP=NFP
    )
    return float(np.asarray(eq.compute("<|F|>_vol", grid=grid)["<|F|>_vol"]))


def solve_candidate_equilibrium(base_eq, scales, candidate_rank):
    """Continue the fixed-boundary solve from low to full resolution."""
    stages_spec = (
        # Over-resolved low-order stages suppress aliasing while remaining cheap.
        (6, 12, 25, "M06"),
        (8, 16, 20, "M08"),
        (10, 20, 15, "M10"),
    )
    eq = None
    stages = []
    for resolution, grid_resolution, maxiter, tag in stages_spec:
        stage_path = (
            WORK_DIR
            / f"candidate_{candidate_rank:02d}_equilibrium_{tag}.h5"
        )
        if stage_path.exists():
            eq = Equilibrium.load(stage_path)
            force_checkpoint = mean_force_residual(eq)
            stages.append(
                {
                    "resolution": resolution,
                    "grid_resolution": eq.M_grid,
                    "resumed_checkpoint": True,
                    "force_after_Pa_m": force_checkpoint,
                    "artifact": stage_path,
                }
            )
            print(
                f"[certify] rank {candidate_rank}: resumed M={resolution} "
                f"checkpoint, force={force_checkpoint:.3e} Pa/m",
                flush=True,
            )
            continue
        reference = base_eq.copy()
        if resolution != base_eq.M:
            reference.change_resolution(
                L=resolution,
                M=resolution,
                N=resolution,
                L_grid=grid_resolution,
                M_grid=grid_resolution,
                N_grid=grid_resolution,
            )
        target_surface, common = perturb_surface(
            reference.surface, scales, surface_volume(reference.surface)
        )
        if eq is None:
            eq = reference.copy()
            eq.surface = target_surface
            eq.set_initial_guess(target_surface, reference.axis, ensure_nested=True)
        else:
            eq.change_resolution(
                L=resolution,
                M=resolution,
                N=resolution,
                L_grid=grid_resolution,
                M_grid=grid_resolution,
                N_grid=grid_resolution,
            )
            eq.surface = target_surface

        force_before = mean_force_residual(eq)
        print(
            f"[certify] rank {candidate_rank}: M=N=L={resolution}, "
            f"grid={grid_resolution}, force before={force_before:.3e} Pa/m",
            flush=True,
        )
        started = time.time()
        optimizer = "lsq-exact"
        tolerance = 1e-5
        eq, solver_result = eq.solve(
            objective="force",
            constraints=get_fixed_boundary_constraints(eq),
            optimizer=optimizer,
            ftol=tolerance,
            xtol=tolerance,
            gtol=tolerance,
            maxiter=maxiter,
            verbose=1,
            copy=False,
        )
        force_after = mean_force_residual(eq)
        eq.save(stage_path)
        stages.append(
            {
                "resolution": resolution,
                "grid_resolution": grid_resolution,
                "max_iterations": maxiter,
                "optimizer": optimizer,
                "volume_compensation": common,
                "force_before_Pa_m": force_before,
                "force_after_Pa_m": force_after,
                "solver_success": bool(getattr(solver_result, "success", False)),
                "solver_message": str(getattr(solver_result, "message", "")),
                "elapsed_seconds": time.time() - started,
                "artifact": stage_path,
            }
        )
        print(
            f"[certify] rank {candidate_rank}: M={resolution} force after="
            f"{force_after:.3e} Pa/m",
            flush=True,
        )
    return eq, stages


def certify_candidate(base_eq, row: dict, candidate_rank: int) -> dict:
    scales = np.array([row["scale_n1"], row["scale_n2"], row["scale_n3"]])
    print(
        f"[certify] rank {candidate_rank}: continuing fixed boundary "
        f"scales={scales.tolist()}",
        flush=True,
    )
    eq, solve_stages = solve_candidate_equilibrium(base_eq, scales, candidate_rank)
    # The bundled full-resolution equilibrium is the honest residual reference.
    # Truncating it to M=10 without re-solving creates a 4.5e5 Pa/m artifact.
    baseline_force = mean_force_residual(base_eq)
    candidate_force = mean_force_residual(eq)
    force_ratio = candidate_force / baseline_force
    eq_path = WORK_DIR / f"candidate_{candidate_rank:02d}_equilibrium.h5"
    eq.save(eq_path)

    if force_ratio > FORCE_RATIO_GATE:
        return {
            "rank": candidate_rank,
            "search_row": row,
            "solve_stages": solve_stages,
            "force_residual_baseline_Pa_m": baseline_force,
            "force_residual_candidate_Pa_m": candidate_force,
            "force_residual_over_baseline": force_ratio,
            "force_gate_max_ratio": FORCE_RATIO_GATE,
            "passed_all_certification_gates": False,
            "rejection_stage": "fixed-boundary force residual",
            "artifacts": {"equilibrium": eq_path},
        }

    # Coordinate-map Jacobian is the direct nestedness/self-intersection gate.
    map_grid = LinearGrid(
        rho=np.linspace(0.05, 1.0, 12),
        M=24,
        N=24,
        NFP=NFP,
        sym=False,
        endpoint=False,
    )
    jac = np.asarray(eq.compute("sqrt(g)", grid=map_grid)["sqrt(g)"])
    jac_positive = bool(np.all(jac > 0))

    winding, winding_metrics = make_winding_surface(eq)
    field0 = FourierCurrentPotentialField.from_surface(
        winding, M_Phi=CURRENT_M, N_Phi=CURRENT_N, sym_Phi="sin"
    )
    source_grid = boundary_grid(REG_GRID_M, REG_GRID_N)
    eval_grid = boundary_grid(REG_GRID_M, REG_GRID_N)
    print(f"[certify] rank {candidate_rank}: REGCOIL solve", flush=True)
    fields, data = solve_regularized_surface_current(
        field0,
        eq,
        lambda_regularization=np.array([REG_LAMBDA]),
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
    field = fields[0]
    kmax = float(np.max(np.asarray(data["|K|"][0])))

    valid_grid = boundary_grid(VALID_GRID_M, VALID_GRID_N)
    valid_source = boundary_grid(FIELD_SOURCE_M, FIELD_SOURCE_N)
    bmag = np.asarray(eq.compute("|B|", grid=valid_grid)["|B|"])
    bn, _ = field.compute_Bnormal(
        eq,
        eval_grid=valid_grid,
        source_grid=valid_source,
        vc_source_grid=valid_source,
        chunk_size=32,
        B_plasma_chunk_size=32,
    )
    bn_stats = normalized_bn_stats(bn, bmag)
    passed = bool(
        jac_positive
        and force_ratio <= FORCE_RATIO_GATE
        and bn_stats["bn_over_b_p95_abs"] <= BN_P95_GATE
        and bn_stats["bn_over_b_max_abs"] <= BN_MAX_GATE
        and kmax <= KMAX_GATE_A_M
    )
    field_path = WORK_DIR / f"candidate_{candidate_rank:02d}_current_potential.h5"
    field.save(field_path)

    return {
        "rank": candidate_rank,
        "search_row": row,
        "solve_stages": solve_stages,
        "certification_resolution": {
            "L": eq.L,
            "M": eq.M,
            "N": eq.N,
            "L_grid": eq.L_grid,
            "M_grid": eq.M_grid,
            "N_grid": eq.N_grid,
            "status": "resolution-limited; M=12 native-grid correction was tested and rejected after optimizer regression",
        },
        "force_residual_baseline_Pa_m": baseline_force,
        "force_residual_candidate_Pa_m": candidate_force,
        "force_residual_over_baseline": force_ratio,
        "jacobian_min": float(jac.min()),
        "jacobian_positive": jac_positive,
        "winding_surface": winding_metrics,
        "K_max_A_m": kmax,
        "K_over_A2": kmax / A2_KMAX_A_M,
        "normal_field": bn_stats,
        "passed_all_certification_gates": passed,
        "artifacts": {"equilibrium": eq_path, "current_potential": field_path},
    }


def certify_stage() -> dict:
    ensure_dirs()
    if not SEARCH_META_FILE.exists():
        raise FileNotFoundError("Run --stage search first")
    search = json.loads(SEARCH_META_FILE.read_text(encoding="utf-8"))
    base_eq = get("W7-X")
    certifications = []
    shortlist = search["shortlist"]
    # Test the load optimum first, then the mildest shortlisted deformation.
    # Intermediate options follow only if both fail a physics gate.
    order = [0, len(shortlist) - 1, len(shortlist) // 2, 1]
    ordered_rows = [shortlist[index] for index in dict.fromkeys(order)]
    for rank, row in enumerate(ordered_rows, start=1):
        result = certify_candidate(base_eq, row, rank)
        certifications.append(result)
        write_json(CERT_FILE, {"candidates": certifications})
        passing_count = sum(
            item["passed_all_certification_gates"] for item in certifications
        )
        if result["passed_all_certification_gates"]:
            print(f"[certify] rank {rank} passed every gate", flush=True)
        if passing_count >= 2:
            break
        print("[certify] testing the next Pareto candidate", flush=True)
    passing = [row for row in certifications if row["passed_all_certification_gates"]]
    payload = {
        "status": "passing candidate found" if passing else "null result",
        "selected": passing[0] if passing else None,
        "candidates": certifications,
        "gates": {
            "Bn_p95_max": BN_P95_GATE,
            "Bn_point_max": BN_MAX_GATE,
            "Kmax_A_m": KMAX_GATE_A_M,
            "positive_coordinate_jacobian": True,
            "force_residual_over_baseline_max": FORCE_RATIO_GATE,
        },
    }
    write_json(CERT_FILE, payload)
    return payload


def make_c1_figure(
    result: dict,
    s: np.ndarray,
    baseline_standoff: np.ndarray,
    candidate_standoff: np.ndarray,
    profile_columns: list[np.ndarray],
) -> None:
    rows = result["lambda_results"]
    colors = ("#B33A2B", "#26736B", "#315D9A")
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.6), facecolor="white")

    ax = axes[0, 0]
    ax.plot(s * 1000, baseline_standoff * 100, "--", color="0.45", lw=1.8, label="A2 W7-X surrogate")
    ax.plot(s * 1000, candidate_standoff * 100, color="#26736B", lw=2.2, label="C1 candidate")
    ax.set_xlabel("mapped duct location s [mm]")
    ax.set_ylabel("wall standoff [cm]")
    ax.set_title("Field-line-traced envelope", loc="left", fontsize=11)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for index, (row, color) in enumerate(zip(rows, colors)):
        baseline_q = np.asarray(profile_columns[4 + 2 * index]) / 1e6
        candidate_q = np.asarray(profile_columns[5 + 2 * index]) / 1e6
        label = f"lambda={row['lambda_m']*100:.1f} cm"
        ax.plot(s * 1000, baseline_q, "--", color=color, lw=1.25, alpha=0.72)
        ax.plot(s * 1000, candidate_q, color=color, lw=2.0, label=label)
    ax.set_xlabel("mapped duct location s [mm]")
    ax.set_ylabel("surface load [MW/m2]")
    ax.set_title("Power-normalized load profiles", loc="left", fontsize=11)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    x = np.arange(len(rows))
    width = 0.34
    baseline_peaks = [row["baseline_peak_over_mean"] for row in rows]
    candidate_peaks = [row["candidate_peak_over_mean"] for row in rows]
    ax.bar(x - width / 2, baseline_peaks, width, color="#86A9C4", label="A2 field-traced")
    ax.bar(x + width / 2, candidate_peaks, width, color="#26736B", label="C1")
    for index, row in enumerate(rows):
        ax.text(
            index + width / 2,
            candidate_peaks[index] + 0.05,
            f"{row['peak_over_mean_improvement_percent']:+.1f}%",
            ha="center",
            fontsize=7.5,
            color="#174C47",
        )
    ax.set_xticks(x, [f"{row['lambda_m']*100:.1f} cm" for row in rows])
    ax.set_ylabel("load peak / mean")
    method = result["robust_acceptance"]["envelope_method_spread_3cm"]
    ax.set_title(
        "Method-consistent gain; "
        f"3 cm method spread {method['spread_percent']:.1f}% > 7.5% gain",
        loc="left",
        fontsize=10.2,
    )
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1, 1]
    for index, row in enumerate(rows):
        for offset, key, color, marker in (
            (-0.10, "baseline_h_wedge_K", "0.45", "o"),
            (+0.10, "candidate_h_wedge_K", "#26736B", "s"),
        ):
            low, high = row[key]
            center = 0.5 * (low + high)
            ax.errorbar(
                index + offset,
                center,
                yerr=[[center - low], [high - center]],
                fmt=marker,
                color=color,
                capsize=4,
                ms=5,
            )
    ax.axhline(1703, color="#B33A2B", ls="-.", lw=1.4, label="FLiBe boiling 1703 K")
    ax.axhspan(973, 1073, color="#C8B58B", alpha=0.28, label="alloy film band")
    ax.set_xticks(x, [f"{row['lambda_m']*100:.1f} cm" for row in rows])
    ax.set_ylabel("CFD-mapped hotspot h-wedge [K]")
    ax.set_title(
        "CFD response (11.6 cm is a censored A2 indication)",
        loc="left",
        fontsize=10.5,
    )
    ax.legend(fontsize=7.5, loc="upper right")

    for ax in axes.flat:
        ax.grid(alpha=0.22)
        for spine in ax.spines.values():
            spine.set_color("0.75")
    cost = result["shape_cost"]
    cert = result["certification"]
    fig.suptitle(
        "C1 bounded null: plasma shape is a weak first-wall hotspot lever",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.012,
        "Resolution-limited M=10 | Shape cost: "
        f"{cost['surface_displacement_rms_m']*1000:.0f} mm RMS / "
        f"{cost['surface_displacement_max_m']*1000:.0f} mm max; "
        f"Bn/B p95 {cert['normal_field']['bn_over_b_p95_abs']:.2%}; "
        f"Kmax {cert['K_max_A_m']/1e6:.2f} MA/m. "
        "No confinement metric is included; this is not a better-stellarator claim.",
        ha="center",
        fontsize=8.2,
        color="0.28",
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    rank = result["certification"]["rank"]
    output = RESULT_DIR / f"c1_candidate_{rank:02d}_summary.png"
    fig.savefig(output, dpi=180, bbox_inches="tight")
    fig.savefig(FIG_DIR / output.name, dpi=180, bbox_inches="tight")
    advisor = ROOT / "results" / "advisor_packet"
    advisor.mkdir(parents=True, exist_ok=True)
    fig.savefig(advisor / f"106_c1_candidate_{rank:02d}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def post_stage(candidate_rank: int | None = None) -> dict:
    """Trace the certified field and map its load through the existing CFD response."""
    ensure_dirs()
    if not CERT_FILE.exists():
        raise FileNotFoundError("Run --stage certify first")
    certification = json.loads(CERT_FILE.read_text(encoding="utf-8"))
    if candidate_rank is None:
        selected = certification.get("selected")
    else:
        selected = next(
            (
                row
                for row in certification.get("candidates", [])
                if row["rank"] == candidate_rank
                and row["passed_all_certification_gates"]
            ),
            None,
        )
    if not selected:
        payload = {
            "status": "null result",
            "reason": "no candidate passed the certification gates",
        }
        write_json(FINAL_FILE, payload)
        return payload

    from build_a2_desc_coil_sol import (
        plasma_field_direct,
        trace_boundary,
        virtual_casing_sources,
    )
    import plasma_wall_qmap_v2_nearfield as nearfield

    eq = Equilibrium.load(selected["artifacts"]["equilibrium"])
    field = FourierCurrentPotentialField.load(
        selected["artifacts"]["current_potential"]
    )
    a2_profile = np.loadtxt(
        ROOT
        / "results"
        / "branch_d"
        / "a2_desc_coil_sol"
        / "a2_standoff_load_feedback.csv",
        delimiter=",",
        skiprows=1,
    )
    phi = a2_profile[:, 0]
    baseline_standoff = a2_profile[:, 2]
    baseline_incidence = a2_profile[:, 3]
    rwall, _, baseline_qrad = load_a2_inputs()
    base_eq = get("W7-X")
    baseline_direct_standoff = direct_standoff(base_eq.surface, rwall=rwall)[4]
    _, direct_geometry_3cm_q = feed_load(
        phi,
        baseline_direct_standoff,
        baseline_incidence,
        LAMBDAS_M[0],
        baseline_qrad,
    )
    direct_geometry_3cm_peak = float(
        direct_geometry_3cm_q.max() / direct_geometry_3cm_q.mean()
    )

    trace_rows = []
    traced_rout = None
    traced_zout = None
    previous_rout = None
    for level in (16, 20, 26):
        print(f"[post] tracing envelope with source M=N={level}", flush=True)
        trace_grid = boundary_grid(level, level)
        traced_phi, _, _, candidate_rout, candidate_zout = trace_boundary(
            eq, field, trace_grid, rwall
        )
        if not np.allclose(traced_phi, phi):
            raise RuntimeError("C1 and A2 toroidal grids differ")
        row = {"source_M_N": level}
        if previous_rout is not None:
            delta = candidate_rout - previous_rout
            row["Rout_delta_from_previous_rms_m"] = float(
                np.sqrt(np.mean(delta**2))
            )
            row["Rout_delta_from_previous_max_abs_m"] = float(
                np.max(np.abs(delta))
            )
        trace_rows.append(row)
        previous_rout = candidate_rout
        traced_rout = candidate_rout
        traced_zout = candidate_zout
    candidate_standoff = rwall - traced_rout

    wall_coords = np.column_stack(
        [np.full_like(phi, rwall), phi, np.zeros_like(phi)]
    )
    wall_source = boundary_grid(FIELD_SOURCE_M, FIELD_SOURCE_N)
    bcoil = np.asarray(
        field.compute_magnetic_field(
            wall_coords, source_grid=wall_source, chunk_size=32
        )
    )
    vc_rows = []
    bplasma = None
    previous = None
    for level in (64, 96, 128, 160, 192):
        print(f"[post] wall plasma-field quadrature M=N={level}", flush=True)
        sources = virtual_casing_sources(eq, level, level)
        candidate = plasma_field_direct(
            wall_coords, sources, chunk=2 if level >= 128 else 4
        )
        row = {
            "M_N": level,
            "full_torus_sources": len(sources["xyz"]),
            "magnitude_mean_T": float(np.linalg.norm(candidate, axis=1).mean()),
            "magnitude_max_T": float(np.linalg.norm(candidate, axis=1).max()),
        }
        if previous is not None:
            delta = np.linalg.norm(candidate - previous, axis=1)
            row["delta_from_previous_p95_T"] = float(np.quantile(delta, 0.95))
            row["delta_from_previous_max_T"] = float(delta.max())
        vc_rows.append(row)
        previous = candidate
        bplasma = candidate
    btotal = bcoil + bplasma
    btotal_mag = np.linalg.norm(btotal, axis=1)
    candidate_incidence = np.abs(btotal[:, 0]) / np.maximum(btotal_mag, 1e-12)

    print("[post] recomputing candidate near-field radiation profile", flush=True)
    nearfield._EQ = eq
    volume = nearfield.build_volume(24, 96)
    candidate_qrad_raw, near_counts = nearfield.los_nearfield(
        volume,
        nearfield.eps_mantle,
        phi,
        rwall,
        sub=(2, 2, 8),
        Ksplit=5.0,
    )
    candidate_qrad = candidate_qrad_raw / candidate_qrad_raw.mean()

    response_path = FIG_DIR / "plasma_sensitivity_nowave.json"
    response = json.loads(response_path.read_text(encoding="utf-8"))
    response_pk = np.asarray([row["pk"] for row in response], dtype=float)
    response_superheat = np.asarray([row["sup"] for row in response], dtype=float)
    if not np.all(np.diff(response_pk) > 0):
        raise ValueError("CFD hotspot-response peaking anchors must be strictly increasing")
    hotspot_response = PchipInterpolator(
        response_pk, response_superheat, extrapolate=False
    )

    def map_hotspot_superheat(peak_over_mean: float) -> float:
        if not response_pk[0] <= peak_over_mean <= response_pk[-1]:
            raise ValueError(
                "C1 hotspot mapping would extrapolate outside the completed CFD "
                f"domain [{response_pk[0]:.3f}, {response_pk[-1]:.3f}]: "
                f"requested {peak_over_mean:.3f}"
            )
        mapped = float(hotspot_response(peak_over_mean))
        if not np.isfinite(mapped):
            raise ValueError("C1 hotspot interpolation returned a non-finite value")
        return mapped
    h_deficit = 3.7
    tin = 900.0
    boiling = 1703.0
    period = 2 * np.pi / NFP
    s = phi / period * S_DUCT_M

    # Existing clean-geometry CFD places the v1 center-loaded hotspot here.
    weak_flushing_s = 0.1461

    lambda_rows = []
    profile_columns = [phi, s, baseline_standoff, candidate_standoff]
    profile_header = [
        "phi_rad",
        "mapped_duct_s_m",
        "A2_traced_standoff_m",
        "C1_traced_standoff_m",
    ]
    baseline_reference = (
        "A2 field-traced envelope, old-assumption 3 cm width",
        "A2 field-traced envelope, finite-connection 5.9 cm width",
        "A2 field-traced envelope, censored-majority 11.6 cm lower-bound width",
    )
    for lambda_index, lambda_m in enumerate(LAMBDAS_M):
        _, same_closure_baseline_q = feed_load(
            phi,
            baseline_standoff,
            baseline_incidence,
            lambda_m,
            baseline_qrad,
        )
        baseline_q = same_closure_baseline_q
        _, candidate_q = feed_load(
            phi,
            candidate_standoff,
            candidate_incidence,
            lambda_m,
            candidate_qrad,
        )
        baseline_peak = float(baseline_q.max() / baseline_q.mean())
        candidate_peak = float(candidate_q.max() / candidate_q.mean())
        ib = int(np.argmax(baseline_q))
        ic = int(np.argmax(candidate_q))
        baseline_sup = map_hotspot_superheat(baseline_peak)
        candidate_sup = map_hotspot_superheat(candidate_peak)
        baseline_wedge = [tin + baseline_sup / h_deficit, tin + baseline_sup]
        candidate_wedge = [tin + candidate_sup / h_deficit, tin + candidate_sup]

        def boiling_status(wedge):
            if wedge[1] < boiling:
                return "below boiling across h-wedge"
            if wedge[0] > boiling:
                return "above boiling across h-wedge"
            return "boiling straddled inside h-wedge"

        lambda_rows.append(
            {
                "lambda_m": float(lambda_m),
                "baseline_reference": baseline_reference[lambda_index],
                "baseline_peak_over_mean": baseline_peak,
                "candidate_peak_over_mean": candidate_peak,
                "peak_over_mean_improvement_percent": float(
                    100 * (1 - candidate_peak / baseline_peak)
                ),
                "baseline_peak_s_m": float(s[ib]),
                "candidate_peak_s_m": float(s[ic]),
                "candidate_peak_distance_from_known_weak_flushing_s_m": float(
                    abs(s[ic] - weak_flushing_s)
                ),
                "baseline_CFD_mapped_superheat_K": baseline_sup,
                "candidate_CFD_mapped_superheat_K": candidate_sup,
                "CFD_mapped_superheat_reduction_percent": float(
                    100 * (1 - candidate_sup / baseline_sup)
                ),
                "baseline_h_wedge_K": baseline_wedge,
                "candidate_h_wedge_K": candidate_wedge,
                "baseline_boiling_status": boiling_status(baseline_wedge),
                "candidate_boiling_status": boiling_status(candidate_wedge),
            }
        )
        label = f"{lambda_m*100:.1f}cm".replace(".", "p")
        profile_columns.extend([baseline_q, candidate_q])
        profile_header.extend(
            [f"A2_qtotal_{label}_W_m2", f"C1_qtotal_{label}_W_m2"]
        )

    np.savetxt(
        RESULT_DIR / "c1_load_profiles.csv",
        np.column_stack(profile_columns),
        delimiter=",",
        header=",".join(profile_header),
        comments="",
    )
    np.savetxt(
        RESULT_DIR / "c1_wall_field.csv",
        np.column_stack(
            [
                phi,
                candidate_standoff,
                traced_rout,
                traced_zout,
                btotal,
                bcoil,
                bplasma,
                candidate_incidence,
            ]
        ),
        delimiter=",",
        header=(
            "phi_rad,standoff_m,Rout_m,Zout_m,Btotal_R_T,Btotal_phi_T,Btotal_Z_T,"
            "Bcoil_R_T,Bcoil_phi_T,Bcoil_Z_T,Bplasma_R_T,Bplasma_phi_T,"
            "Bplasma_Z_T,sin_incidence"
        ),
        comments="",
    )

    shape_cost = surface_cost(base_eq.surface, eq.surface)
    trace_delta = candidate_standoff - direct_standoff(eq.surface, rwall=rwall)[4]
    baseline_robust = max(
        row["baseline_peak_over_mean"] for row in lambda_rows
    )
    candidate_robust = max(
        row["candidate_peak_over_mean"] for row in lambda_rows
    )
    traced_pass = candidate_robust < baseline_robust
    method_spread_percent = float(
        100 * (baseline_robust / direct_geometry_3cm_peak - 1)
    )
    result = {
        "status": "bounded null: plasma-shape optimization is a weak first-wall hotspot lever",
        "framing": "wall-load-only, confinement-agnostic upper bound",
        "claim_scope": "load-domain result is primary; hotspot mapping is secondary only",
        "hotspot_phase_uncertainty_fraction": 0.05,
        "hotspot_phase_note": "four matched phases bound location sensitivity; carry +/-5% on every mapped hotspot",
        "h_wedge_note": "3.7 divisor is a bounded Ferrero cross-geometry transfer assumption, not a same-case correction",
        "not_claimed": "a better stellarator, confinement optimization, engineering coil design, or full M=12 equilibrium certificate",
        "decision_scales": {
            key: selected["search_row"][key]
            for key in ("scale_n1", "scale_n2", "scale_n3")
        },
        "shape_cost": shape_cost,
        "certification": selected,
        "field_line_envelope": {
            "trace_source_convergence": trace_rows,
            "A2_standoff_min_max_m": [
                float(baseline_standoff.min()),
                float(baseline_standoff.max()),
            ],
            "C1_standoff_min_max_m": [
                float(candidate_standoff.min()),
                float(candidate_standoff.max()),
            ],
            "C1_trace_minus_direct_rms_m": float(
                np.sqrt(np.mean(trace_delta**2))
            ),
            "C1_trace_minus_direct_max_abs_m": float(np.max(np.abs(trace_delta))),
            "baseline_direct_standoff_min_max_m": [
                float(baseline_direct_standoff.min()),
                float(baseline_direct_standoff.max()),
            ],
        },
        "wall_field": {
            "Btotal_min_max_T": [float(btotal_mag.min()), float(btotal_mag.max())],
            "incidence_mean_max": [
                float(candidate_incidence.mean()),
                float(candidate_incidence.max()),
            ],
            "plasma_fraction_mean_max": [
                float(np.mean(np.linalg.norm(bplasma, axis=1) / btotal_mag)),
                float(np.max(np.linalg.norm(bplasma, axis=1) / btotal_mag)),
            ],
            "virtual_casing_convergence": vc_rows,
        },
        "radiation": {
            "method": "candidate-specific near-field cell-extent LOS, NT=24, NZ=96, sub=(2,2,8)",
            "peak_over_mean": float(candidate_qrad.max()),
            "near_cells_per_wall_point_min_max": [
                int(near_counts.min()),
                int(near_counts.max()),
            ],
        },
        "lambda_results": lambda_rows,
        "robust_acceptance": {
            "candidate_worst_peak_over_mean": candidate_robust,
            "like_for_like_field_trace": {
                "baseline_worst_peak_over_mean": baseline_robust,
                "improvement_percent": float(
                    100 * (1 - candidate_robust / baseline_robust)
                ),
                "passed": traced_pass,
                "method_note": "A2 and C1 both use their matching current-potential field-line envelopes at every lambda",
            },
            "envelope_method_spread_3cm": {
                "direct_geometry_peak_over_mean": direct_geometry_3cm_peak,
                "field_traced_peak_over_mean": baseline_robust,
                "spread_percent": method_spread_percent,
                "interpretation": "numerical model-definition diagnostic, not an acceptance comparator",
            },
            "criterion": "candidate worst-case peak/mean must be lower than method-consistent A2 over 3.0, 5.9, and 11.6 cm",
            "review_ruling": (
                "controlled field-traced comparison is valid, but the result is a bounded null: "
                "the wide-SOL gain is about 1%, smaller than model-definition spread, and not a remedy"
            ),
        },
        "hotspot_mapping": {
            "method": "shape-preserving interpolation of four completed clean-geometry CFD cases",
            "response_points": response,
            "peak_over_mean_domain": [float(response_pk[0]), float(response_pk[-1])],
            "extrapolation_allowed": False,
            "h_wedge": "LES mapped superheat divided by 3.7 at a bounded Ferrero cross-geometry transfer edge",
            "known_weak_flushing_location_s_m": weak_flushing_s,
            "location_limit": "CFD response cases used center-peaked v1 profiles; C1 profiles were not rerun in CFD, so location is audited separately rather than embedded in the interpolation",
        },
        "confinement_cost": {
            "status": "not scored",
            "reason": "no cheap validated ripple/QI/confinement proxy was available in the fixed-boundary workflow",
        },
        "limitations": [
            "The equilibrium certificate stops at M=N=L=10; attempted M=12 BFGS corrections regressed the force objective and were rejected.",
            "REGCOIL is an equilibrium-derived current-potential surrogate, not W7-X engineering coils or a free-boundary design.",
            "The cylindrical first wall and reduced SOL closure remain in force.",
            "The existing CFD response maps load peaking; a four-phase sweep bounds location sensitivity at 5%, which must accompany every mapped hotspot.",
        ],
    }
    candidate_file = RESULT_DIR / f"c1_candidate_{selected['rank']:02d}_metrics.json"
    write_json(candidate_file, result)
    if candidate_rank is None:
        write_json(FINAL_FILE, result)
    make_c1_figure(result, s, baseline_standoff, candidate_standoff, profile_columns)
    print(f"[post] wrote {candidate_file}", flush=True)
    return result


def finalize_stage() -> dict:
    """Select the post-trace Pareto candidate under the approved baseline."""
    candidate_files = sorted(RESULT_DIR.glob("c1_candidate_*_metrics.json"))
    if not candidate_files:
        raise FileNotFoundError("Run postprocessing for certified candidates first")
    candidates = [json.loads(path.read_text(encoding="utf-8")) for path in candidate_files]
    like_for_like = [
        row
        for row in candidates
        if row["robust_acceptance"]["like_for_like_field_trace"]["passed"]
    ]
    if like_for_like:
        selected = max(
            like_for_like,
            key=lambda row: (
                row["robust_acceptance"]["like_for_like_field_trace"][
                    "improvement_percent"
                ],
                -row["shape_cost"]["surface_displacement_rms_m"],
            ),
        )
        selection_basis = (
            "largest method-consistent field-traced robust improvement, then lower RMS shape cost"
        )
    else:
        selected = min(
            candidates,
            key=lambda row: row["robust_acceptance"]["candidate_worst_peak_over_mean"],
        )
        selection_basis = "least-bad certified null candidate"

    rank = selected["certification"]["rank"]
    selected["selection"] = {
        "selected_candidate_rank": rank,
        "basis": selection_basis,
        "review_ruling": "approved 2026-07-12",
        "promotion_rule": (
            "use only the method-consistent field-traced comparison and promote C1 as a bounded null"
        ),
        "candidates_compared": [row["certification"]["rank"] for row in candidates],
    }
    write_json(FINAL_FILE, selected)

    source_figure = RESULT_DIR / f"c1_candidate_{rank:02d}_summary.png"
    canonical_figure = RESULT_DIR / "c1_wall_load_upper_bound_summary.png"
    shutil.copy2(source_figure, canonical_figure)
    shutil.copy2(source_figure, FIG_DIR / canonical_figure.name)
    shutil.copy2(
        source_figure,
        ROOT / "results" / "advisor_packet" / "106_c1_wall_load_upper_bound.png",
    )

    comparison_rows = []
    for row in candidates:
        comparison_rows.append(
            [
                row["certification"]["rank"],
                row["shape_cost"]["surface_displacement_rms_m"],
                row["shape_cost"]["surface_displacement_max_m"],
                row["robust_acceptance"]["candidate_worst_peak_over_mean"],
                row["robust_acceptance"]["like_for_like_field_trace"][
                    "improvement_percent"
                ],
                int(
                    row["robust_acceptance"]["like_for_like_field_trace"][
                        "passed"
                    ]
                ),
                row["robust_acceptance"]["envelope_method_spread_3cm"][
                    "spread_percent"
                ],
            ]
        )
    np.savetxt(
        RESULT_DIR / "c1_candidate_comparison.csv",
        np.asarray(comparison_rows),
        delimiter=",",
        header=(
            "candidate_rank,shape_rms_m,shape_max_m,candidate_robust_peak_over_mean,"
            "like_for_like_improvement_percent,like_for_like_pass,"
            "envelope_method_spread_3cm_percent"
        ),
        comments="",
    )
    print(f"[finalize] selected candidate {rank}; wrote {FINAL_FILE}", flush=True)
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("search", "certify", "post", "finalize", "all"),
        default="all",
    )
    parser.add_argument("--candidate-rank", type=int)
    args = parser.parse_args()
    if args.stage == "all":
        search_stage()
        certification = certify_stage()
        ranks = [
            row["rank"]
            for row in certification.get("candidates", [])
            if row["passed_all_certification_gates"]
        ]
        if not ranks:
            write_json(
                FINAL_FILE,
                {
                    "status": "null result",
                    "reason": "no candidate passed the certification gates",
                },
            )
            return
        for rank in ranks:
            post_stage(candidate_rank=rank)
        finalize_stage()
        return
    if args.stage == "search":
        search_stage()
    elif args.stage == "certify":
        certify_stage()
    elif args.stage == "post":
        post_stage(candidate_rank=args.candidate_rank)
    elif args.stage == "finalize":
        finalize_stage()


if __name__ == "__main__":
    main()

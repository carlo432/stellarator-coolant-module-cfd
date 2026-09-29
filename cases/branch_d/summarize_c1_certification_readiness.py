#!/usr/bin/env python3
"""Summarize C1 thermal LES readiness against Kawamura-style targets.

This does not certify C1.  It turns the existing smoke/mid profile CSVs into a
clear readiness table and writes the benchmark gates needed before the project
can claim thermal-channel LES certification.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
FIGS = WORK / "figs"
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
DIGITIZATION = ROOT / "references" / "digitization" / "kawamura_1998"
DIGITIZATION_MANIFEST = DIGITIZATION / "target_manifest.csv"

PR5_PREFLIGHT = {
    "case": "les_channel_c1_thermal_pr5_preflight",
    "path": WORK / "les_channel_c1_thermal_pr5_preflight",
    "cells": 18432,
    "pr": 5.0,
    "nu": 3.09278e-6,
    "alpha": 3.09278e-6 / 5.0,
    "mesh_status": "checkMesh OK",
    "run_status": "short smoke run completed; not certification",
}

TRANSFORMED_SMOKE = {
    "case": "les_channel_c1_kawamura_transformed_smoke",
    "path": WORK / "les_channel_c1_kawamura_transformed_smoke",
    "cells": 2304,
    "pr": 5.0,
    "field": "H",
    "latest_time": "0.0010184",
    "mesh_status": "checkMesh OK",
    "run_status": "tiny transformed-source smoke completed; not certification",
    "source": "5.567004*(U.x()/1.1598)",
}

TRANSFORMED_PR5_PREFLIGHT = {
    "case": "les_channel_c1_kawamura_transformed_pr5_preflight",
    "path": WORK / "les_channel_c1_kawamura_transformed_pr5_preflight",
    "cells": 18432,
    "pr": 5.0,
    "field": "H",
    "latest_time": "0.009892992",
    "mesh_status": "checkMesh OK",
    "run_status": "18k transformed Pr=5 preflight completed; not certification",
    "source": "5.567004*(U.x()/1.1598)",
}

RE180_CALIBRATION_ROWS = [
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_preflight",
        "u_bar": 1.1598,
        "latest_time": "0.009892992",
        "re_tau": 97.29004322396099,
        "status": "completed",
        "use": "original transformed preflight; undershot target",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_cal_ubar136",
        "u_bar": 1.36,
        "latest_time": "0.009892992",
        "re_tau": 178.1935728541673,
        "status": "completed",
        "use": "smoke-scale seed before larger shakedown",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_cal_ubar140",
        "u_bar": 1.40,
        "latest_time": "0.009892992",
        "re_tau": 190.27579190073666,
        "status": "completed",
        "use": "near-target upper bracket",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_cal_ubar160",
        "u_bar": 1.60,
        "latest_time": "0.009892992",
        "re_tau": 241.79150901718208,
        "status": "completed",
        "use": "strong overshoot",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_cal_ubar220",
        "u_bar": 2.20,
        "latest_time": "n/a",
        "re_tau": math.nan,
        "status": "generated_not_run",
        "use": "skipped after u_bar 1.60 overshot",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_cal_ubar280",
        "u_bar": 2.80,
        "latest_time": "n/a",
        "re_tau": math.nan,
        "status": "generated_not_run",
        "use": "skipped after u_bar 1.60 overshot",
    },
]

STAGE1_SHAKEDOWN_ROWS = [
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_re180_stage1_shakedown_ubar136",
        "u_bar": 1.36,
        "latest_time": "0.019889956",
        "re_tau": 165.75073525432018,
        "status": "completed",
        "use": "larger-mesh shakedown undershot target",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_re180_stage1_shakedown_ubar140",
        "u_bar": 1.40,
        "latest_time": "0.019924859",
        "re_tau": 174.85722962961077,
        "status": "completed",
        "use": "near-target lower bracket",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_re180_stage1_shakedown_ubar142",
        "u_bar": 1.42,
        "latest_time": "0.020038972",
        "re_tau": 179.15010850890164,
        "status": "completed",
        "use": "current larger-mesh Stage 1 seed",
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_re180_stage1_ubar142",
        "u_bar": 1.42,
        "latest_time": "n/a",
        "re_tau": math.nan,
        "status": "preflight_only",
        "use": "current full 0.2 s candidate; checkMesh OK, solver not launched",
    },
]


PROFILE_INPUTS = [
    {
        "case": "les_channel_c1_thermal_smoke",
        "role": "smoke/plumbing",
        "path": FIGS / "c1_thermal_profile_les_channel_c1_thermal_smoke.csv",
        "certification_use": (
            "Verifies scalarTransport, TMean/TPrime2Mean profile plumbing, and "
            "postprocessor I/O only."
        ),
    },
    {
        "case": "les_channel_c1_thermal_mid",
        "role": "short covariance route",
        "path": FIGS / "c1_thermal_profile_les_channel_c1_thermal_mid.csv",
        "certification_use": (
            "Verifies sparse write control and resolved u'v'/v'T'/Pr_t extraction "
            "route only."
        ),
    },
    {
        "case": "les_channel_c1_thermal_pr5_preflight",
        "role": "Pr=5 benchmark-path smoke",
        "path": FIGS / "c1_thermal_profile_les_channel_c1_thermal_pr5_preflight.csv",
        "certification_use": (
            "Verifies the direct Kawamura-range Pr=5 scalar setting, solver "
            "smoke path, and harness route only."
        ),
    },
    {
        "case": "les_channel_c1_kawamura_transformed_smoke",
        "role": "Kawamura transformed-source smoke",
        "path": FIGS / "c1_thermal_profile_les_channel_c1_kawamura_transformed_smoke.csv",
        "certification_use": (
            "Verifies the transformed H-field fixed-wall setup, inline U.x() "
            "source expression, and harness --field route only."
        ),
    },
    {
        "case": "les_channel_c1_kawamura_transformed_pr5_preflight",
        "role": "Kawamura transformed Pr=5 preflight",
        "path": FIGS / "c1_thermal_profile_les_channel_c1_kawamura_transformed_pr5_preflight.csv",
        "certification_use": (
            "Verifies the transformed H-field setup at the normal 18k smoke "
            "scale with sparse writes, HMean/HPrime2Mean, covariance, and Pr_t "
            "plotting.  Still short-run machinery evidence only."
        ),
    },
]


TARGET_ROWS = [
    {
        "gate": "reference_dns",
        "benchmark_target": (
            "Kawamura et al. 1998 turbulent channel heat-transfer DNS; "
            "Re_tau = 180; Pr = 0.025 to 5; uniform heating from both walls."
        ),
        "current_state": (
            "Source note and extracted text exist locally; project Pr = 14.4 is "
            "above the direct DNS Pr range."
        ),
        "needed_next": (
            "Digitize published profiles and keep Pr = 14.4 as an extrapolated "
            "project case unless a direct Pr = 5 benchmark run is added."
        ),
        "status": "partial",
    },
    {
        "gate": "thermal_boundary_condition",
        "benchmark_target": (
            "Fully developed channel flow with both walls heated in the DNS "
            "temperature transformation."
        ),
        "current_state": (
            "Generated C1 scaffold uses opposed fixed scalar gradients to prevent "
            "bulk scalar drift during local smoke tests.  A tiny transformed H-field "
            "source-expression smoke and an 18k transformed Pr=5 preflight now "
            "exist and run locally, but neither is a benchmark run."
        ),
        "needed_next": (
            "Scale the transformed H-field setup beyond expression smoke and run a "
            "statistically meaningful matched Pr = 5 benchmark before claiming "
            "comparison."
        ),
        "status": "route_exists_smoke_only",
    },
    {
        "gate": "friction_reynolds_number",
        "benchmark_target": "Re_tau = 180 for direct Kawamura comparison.",
        "current_state": (
            "Smoke case is about Re_tau = 97.3; mid covariance case is about "
            "Re_tau = 120; transformed smoke is about Re_tau = 67.3; transformed "
            "Pr=5 preflight is about Re_tau = 97.3; a short transformed calibration "
            "probe at u_bar = 1.36 measures Re_tau = 178.2; larger Stage 1 "
            "shakedown at u_bar = 1.42 measures Re_tau = 179.15; O-like case is "
            "meshed but not run."
        ),
        "needed_next": (
            "Run a statistically meaningful case near Re_tau = 180 or document "
            "why the comparison is only qualitative."
        ),
        "status": "open",
    },
    {
        "gate": "prandtl_number",
        "benchmark_target": (
            "Direct DNS profile comparison is available up to Pr = 5 in the "
            "local Kawamura source."
        ),
        "current_state": (
            "Current C1 scalar target is Pr = 14.4 for FLiBe relevance.  A "
            "Pr = 5 direct-benchmark-path smoke case now exists and ran locally, "
            "and Pr = 5 transformed-source smoke/preflight cases now exist, but "
            "none are statistically meaningful."
        ),
        "needed_next": (
            "Best certification path: run a Pr = 5 matched benchmark case, then "
            "run Pr = 14.4 as the high-Pr project extrapolation."
        ),
        "status": "open",
    },
    {
        "gate": "mean_temperature_profile",
        "benchmark_target": "Digitized DNS T+ versus y+ profile.",
        "current_state": "Local smoke/mid T+ CSVs and figures exist but are short-run scaffold data.",
        "needed_next": "Digitize the Kawamura mean-temperature profile and overlay only after a long run.",
        "status": "open",
    },
    {
        "gate": "temperature_variance_profile",
        "benchmark_target": "Digitized DNS theta_rms+ or equivalent temperature-variance profile.",
        "current_state": "TPrime2Mean/T_rms+ extraction exists; current values are not statistically converged.",
        "needed_next": "Digitize the variance profile and compare after stable thermal averaging.",
        "status": "open",
    },
    {
        "gate": "turbulent_heat_flux_and_prt",
        "benchmark_target": (
            "DNS wall-normal turbulent heat flux and turbulent Prandtl-number "
            "behavior, especially near-wall and log-layer trends."
        ),
        "current_state": (
            "Snapshot covariance route exists, but current Pr_t samples are sparse "
            "and noisy with many invalid/negative-diffusivity points."
        ),
        "needed_next": (
            "Use many saved statistically independent snapshots after the flow and "
            "thermal scalar are developed; compare only finite positive eddy "
            "viscosity/diffusivity regions."
        ),
        "status": "route_exists_data_open",
    },
    {
        "gate": "statistical_averaging",
        "benchmark_target": "Stationary, fully developed thermal statistics over a meaningful sampling window.",
        "current_state": (
            "Pr=14.4 smoke latest time 0.009892992 s; Pr=5 smoke latest time "
            "0.009892992 s; transformed smoke latest time 0.0010184 s; mid latest "
            "time 0.040121929 s; transformed Pr=5 preflight/latest calibration "
            "probes latest time 0.009892992 s; Stage 1 shakedowns latest time "
            "about 0.02004 s."
        ),
        "needed_next": "Run a long O-like/HPC-staged case with stable TMean, TPrime2Mean, and covariance windows.",
        "status": "open",
    },
]


def as_float(value: str | None) -> float:
    if value is None or value == "":
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def finite(values: list[float]) -> list[float]:
    return [v for v in values if math.isfinite(v)]


def median_or_nan(values: list[float]) -> float:
    vals = finite(values)
    return statistics.median(vals) if vals else math.nan


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{value:.6g}"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def summary_path_for_profile(path: Path) -> Path:
    name = path.name
    if name.startswith("c1_thermal_profile_"):
        return path.with_name(name.replace("c1_thermal_profile_", "c1_thermal_summary_", 1)).with_suffix(".json")
    return path.with_suffix(".json")


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def count_data_rows(path: Path) -> int:
    if not path.exists():
        return -1
    with path.open(newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return 0
    return max(len(rows) - 1, 0)


def digitization_status() -> tuple[list[dict[str, object]], int]:
    targets = read_rows(DIGITIZATION_MANIFEST)
    rows: list[dict[str, object]] = []
    missing_required = 0
    for target in targets:
        data_rows = count_data_rows(DIGITIZATION / target["digitized_csv"])
        if data_rows < 0:
            status = "missing_file"
        elif data_rows == 0:
            status = "empty"
        else:
            status = "has_data"
        if target["priority"] == "required" and status != "has_data":
            missing_required += 1
        rows.append(
            {
                "figure_id": target["figure_id"],
                "priority": target["priority"],
                "profile": target["profile"],
                "digitized_csv": target["digitized_csv"],
                "data_rows": data_rows,
                "status": status,
            }
        )
    return rows, missing_required


def summarize_profile(spec: dict[str, object]) -> dict[str, object]:
    path = Path(spec["path"])
    summary_path = summary_path_for_profile(path)
    summary = read_json(summary_path)
    rows = read_rows(path)
    y_plus = [as_float(row.get("y_plus")) for row in rows]
    t_plus = [as_float(row.get("T_plus")) for row in rows]
    t_rms_plus = [as_float(row.get("T_rms_plus")) for row in rows]
    pr_t = [as_float(row.get("Pr_t_resolved")) for row in rows]

    finite_y = finite(y_plus)
    finite_t = finite(t_plus)
    finite_rms = finite(t_rms_plus)
    finite_prt = finite(pr_t)
    positive_prt = [v for v in finite_prt if v > 0.0]
    positive_midband = [
        p
        for yp, p in zip(y_plus, pr_t)
        if math.isfinite(yp) and 10.0 <= yp <= 80.0 and math.isfinite(p) and p > 0.0
    ]

    max_t_plus = max(finite_t) if finite_t else math.nan
    y_at_max_t_plus = math.nan
    if math.isfinite(max_t_plus):
        for yp, tp in zip(y_plus, t_plus):
            if math.isfinite(tp) and tp == max_t_plus:
                y_at_max_t_plus = yp
                break

    total = len(rows)
    return {
        "case": spec["case"],
        "role": spec["role"],
        "profile_csv": str(path.relative_to(ROOT)),
        "summary_json": str(summary_path.relative_to(ROOT)) if summary_path.exists() else "",
        "summary_Pr": summary.get("Pr", math.nan),
        "summary_Re_tau": summary.get("Re_tau", math.nan),
        "summary_T_tau": summary.get("T_tau", math.nan),
        "summary_covariance_snapshots": summary.get("covariance", {}).get("snapshot_count", 0) if summary else 0,
        "profile_points": total,
        "first_cell_y_plus": min(finite_y) if finite_y else math.nan,
        "max_y_plus": max(finite_y) if finite_y else math.nan,
        "max_T_plus": max_t_plus,
        "y_plus_at_max_T_plus": y_at_max_t_plus,
        "max_T_rms_plus": max(finite_rms) if finite_rms else math.nan,
        "finite_Pr_t_points": len(finite_prt),
        "positive_Pr_t_points": len(positive_prt),
        "positive_Pr_t_fraction": len(positive_prt) / total if total else math.nan,
        "positive_Pr_t_median": median_or_nan(positive_prt),
        "positive_Pr_t_midband_median": median_or_nan(positive_midband),
        "certification_use": spec["certification_use"],
    }


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def readiness_markdown(readiness: list[dict[str, object]]) -> str:
    pr5_exists = PR5_PREFLIGHT["path"].exists()
    pr5_checkmesh = (PR5_PREFLIGHT["path"] / "log.checkMesh").exists()
    digitized_rows, required_missing = digitization_status()
    lines = [
        "# C1 Kawamura Certification Targets And Readiness",
        "",
        "Date: 2026-06-22",
        "",
        "## Verdict",
        "",
        "C1 is ready at the scaffold/target-definition level, not at the physics-certification level. "
        "The local cases prove that the high-Pr scalar path, profile extraction, sparse writes, and "
        "resolved covariance machinery work.  They do not yet prove agreement with thermal-channel DNS.",
        "",
        "The honest certification route is: first match the Kawamura DNS conditions as closely as the "
        "OpenFOAM scaffold permits, then label `Pr = 14.4` as a FLiBe-relevant extrapolation above the "
        "paper's direct `Pr <= 5` range.",
        "",
        "## Kawamura Anchor",
        "",
        "| Quantity | Value |",
        "|---|---|",
        "| Citation | Kawamura, Ohsaka, Abe, and Yamamoto, IJHFF 19(5), 482-491, 1998 |",
        "| DOI | `10.1016/S0142-727X(98)10026-7` |",
        "| Flow | fully developed turbulent channel |",
        "| Thermal condition | uniform heating from both walls |",
        "| Re_tau | `180` |",
        "| Pr range in local source | `0.025` to `5` |",
        "| Domain | `6.4d x 2d x 3.2d` |",
        "| Grid | `128 x 66 x 128` for `Pr <= 1.5`; `256 x 128 x 256` for `Pr = 5` |",
        "| Local project caveat | project `Pr = 14.4` is extrapolated beyond the direct DNS range |",
        "",
        "## Digitization Workspace",
        "",
        f"- Workspace: `{DIGITIZATION.relative_to(ROOT)}`",
        "- Rendered pages: `pages/kawamura_1998_page-03.png` through `pages/kawamura_1998_page-08.png`",
        "- Checker: `scripts/check_kawamura_digitization_status.py`",
        "- Boundary/source run-card: `docs/c1_kawamura_boundary_condition_matching_plan.md`",
        "- Source implementation note: `docs/c1_kawamura_source_term_implementation_note.md`",
        "- Transformed smoke status: `docs/c1_kawamura_transformed_smoke_status.md`",
        "- Transformed Pr=5 preflight status: `docs/c1_kawamura_transformed_pr5_preflight_status.md`",
        "- Re_tau 180 calibration status: `docs/c1_kawamura_re180_calibration_status.md`",
        "- Stage 1 shakedown status: `docs/c1_kawamura_stage1_shakedown_status.md`",
        "- Re_tau 180 next-run card: `docs/c1_kawamura_re180_next_run_card.md`",
        "- Re_tau 180 calibration CSV/plot: `results/c1_kawamura_re180_calibration.csv`, `results/c1_kawamura_re180_calibration.png`",
        "- Stage 1 shakedown CSV/plot: `results/c1_kawamura_stage1_shakedown.csv`, `results/c1_kawamura_stage1_shakedown.png`",
        f"- Required targets missing data rows: `{required_missing}`",
        "",
        "| Figure | Priority | Profile | Rows | Status |",
        "|---|---|---|---:|---|",
    ]
    for row in digitized_rows:
        lines.append(
            f"| `{row['figure_id']}` | {row['priority']} | {row['profile']} | "
            f"{row['data_rows']} | `{row['status']}` |"
        )
    lines.extend(
        [
            "",
            "## Certification Gates",
            "",
            "| Gate | Benchmark/Target | Current State | Needed Next | Status |",
            "|---|---|---|---|---|",
        ]
    )
    for row in TARGET_ROWS:
        lines.append(
            f"| `{row['gate']}` | {row['benchmark_target']} | {row['current_state']} | "
            f"{row['needed_next']} | `{row['status']}` |"
        )
    lines.extend(
        [
            "",
            "## Direct Pr=5 Benchmark-Path Smoke",
            "",
            "| Quantity | Value |",
            "|---|---|",
            f"| Case | `{PR5_PREFLIGHT['case']}` |",
            f"| Path | `{PR5_PREFLIGHT['path'].relative_to(ROOT)}` |",
            f"| Exists | `{pr5_exists}` |",
            f"| Cells | `{PR5_PREFLIGHT['cells']}` |",
            f"| Pr | `{PR5_PREFLIGHT['pr']:.6g}` |",
            f"| Scalar diffusivity | `{PR5_PREFLIGHT['alpha']:.6g} m^2/s` |",
            f"| Mesh/preflight | `{PR5_PREFLIGHT['mesh_status'] if pr5_checkmesh else 'not checked'}` |",
            f"| Solver status | `{PR5_PREFLIGHT['run_status']}` |",
            "",
            "This case is useful because it puts the next direct Kawamura-range run "
            "one command away.  It is still the opposed-gradient scaffold, so it "
            "does not remove the boundary-condition matching task.",
            "",
            "## Kawamura Transformed-Mode Source Smoke",
            "",
            "| Quantity | Value |",
            "|---|---|",
            f"| Case | `{TRANSFORMED_SMOKE['case']}` |",
            f"| Path | `{TRANSFORMED_SMOKE['path'].relative_to(ROOT)}` |",
            f"| Exists | `{TRANSFORMED_SMOKE['path'].exists()}` |",
            f"| Cells | `{TRANSFORMED_SMOKE['cells']}` |",
            f"| Field | `{TRANSFORMED_SMOKE['field']}` |",
            f"| Pr | `{TRANSFORMED_SMOKE['pr']:.6g}` |",
            f"| Latest time | `{TRANSFORMED_SMOKE['latest_time']} s` |",
            f"| Source expression | `{TRANSFORMED_SMOKE['source']}` |",
            f"| Mesh/preflight | `{TRANSFORMED_SMOKE['mesh_status']}` |",
            f"| Solver status | `{TRANSFORMED_SMOKE['run_status']}` |",
            "",
            "This case confirms that the inline transformed-source expression and "
            "`H`-field harness route work locally under OpenFOAM 2512.  It is only "
            "machinery evidence.",
            "",
            "## Kawamura Transformed Pr=5 Preflight",
            "",
            "| Quantity | Value |",
            "|---|---|",
            f"| Case | `{TRANSFORMED_PR5_PREFLIGHT['case']}` |",
            f"| Path | `{TRANSFORMED_PR5_PREFLIGHT['path'].relative_to(ROOT)}` |",
            f"| Exists | `{TRANSFORMED_PR5_PREFLIGHT['path'].exists()}` |",
            f"| Cells | `{TRANSFORMED_PR5_PREFLIGHT['cells']}` |",
            f"| Field | `{TRANSFORMED_PR5_PREFLIGHT['field']}` |",
            f"| Pr | `{TRANSFORMED_PR5_PREFLIGHT['pr']:.6g}` |",
            f"| Latest time | `{TRANSFORMED_PR5_PREFLIGHT['latest_time']} s` |",
            f"| Source expression | `{TRANSFORMED_PR5_PREFLIGHT['source']}` |",
            f"| Mesh/preflight | `{TRANSFORMED_PR5_PREFLIGHT['mesh_status']}` |",
            f"| Solver status | `{TRANSFORMED_PR5_PREFLIGHT['run_status']}` |",
            "",
            "This case scales the transformed `H` route to the normal local C1 smoke "
            "mesh size and verifies sparse `HMean`/`HPrime2Mean`, snapshot "
            "covariance, and resolved-`Pr_t` plotting.  It is still a short, "
            "low-`Re_tau` machinery check, not thermal-channel DNS validation.",
            "",
            "## Short Re_tau Calibration Sweep",
            "",
            "| Case | u_bar | Latest time | Measured Re_tau | Status | Use |",
            "|---|---:|---:|---:|---|---|",
        ]
    )
    for row in RE180_CALIBRATION_ROWS:
        lines.append(
            f"| `{row['case']}` | {row['u_bar']:.4g} | {row['latest_time']} | "
            f"{fmt(float(row['re_tau']))} | `{row['status']}` | {row['use']} |"
        )
    lines.extend(
        [
            "",
            "The short sweep identifies `u_bar = 1.36 m/s` as the smoke-scale "
            "seed because it measured `Re_tau = 178.19` on the smoke-scale "
            "transformed `Pr = 5` case.  This is calibration evidence only; it "
            "does not make the case statistically converged.",
            "",
            "## Larger-Mesh Stage 1 Shakedown",
            "",
            "| Case | u_bar | Latest time | Measured Re_tau | Status | Use |",
            "|---|---:|---:|---:|---|---|",
        ]
    )
    for row in STAGE1_SHAKEDOWN_ROWS:
        lines.append(
            f"| `{row['case']}` | {row['u_bar']:.4g} | {row['latest_time']} | "
            f"{fmt(float(row['re_tau']))} | `{row['status']}` | {row['use']} |"
        )
    lines.extend(
        [
            "",
            "The larger `147,456`-cell shakedown supersedes the smoke-scale seed: "
            "`u_bar = 1.42 m/s` measured `Re_tau = 179.15` and is the current "
            "full Stage 1 candidate.  The full `0.2 s` case is generated and "
            "`checkMesh OK`, but the solver has not been launched.",
            "",
            "## Current Local Profile Readiness",
            "",
            "| Case | Role | Pr | Re_tau | Snapshots | Points | y+ min | y+ max | max T+/H+ | max RMS+ | finite Pr_t | positive Pr_t | Use |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in readiness:
        lines.append(
            f"| `{row['case']}` | {row['role']} | "
            f"{fmt(float(row['summary_Pr']))} | {fmt(float(row['summary_Re_tau']))} | "
            f"{row['summary_covariance_snapshots']} | {row['profile_points']} | "
            f"{fmt(float(row['first_cell_y_plus']))} | {fmt(float(row['max_y_plus']))} | "
            f"{fmt(float(row['max_T_plus']))} | {fmt(float(row['max_T_rms_plus']))} | "
            f"{row['finite_Pr_t_points']} | {row['positive_Pr_t_points']} | "
            f"{row['certification_use']} |"
        )
    lines.extend(
        [
            "",
            "## Next Run Card",
            "",
            "1. Digitize Kawamura thermal profiles before making DNS pass/fail plots.",
            "2. Use the completed transformed-mode smoke as the local source-expression proof, not as validation.",
            "3. Use the completed `18,432`-cell transformed `Pr = 5` preflight as the local normal-smoke-scale proof, not as validation.",
            "4. Use the completed larger Stage 1 shakedown seed `u_bar = 1.42`; do not treat it as validation.",
            "5. Keep `Pr = 14.4` as the FLiBe-relevant high-Pr project case, with the extrapolation caveat visible in captions.",
            "6. Run the O-like or HPC-staged case long enough that transformed scalar means/variances and saved-snapshot covariance are stationary.",
            "7. Report current smoke, mid, Pr=5 smoke, transformed smoke, and transformed preflight outputs only as machinery checks.",
            "",
            "## Do Not Claim Yet",
            "",
            "- Do not claim C1 thermal LES certification.",
            "- Do not claim direct DNS validation at `Pr = 14.4` from Kawamura.",
            "- Do not use current `Pr_t` values as physics evidence; they are extraction-route evidence.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    readiness = [summarize_profile(spec) for spec in PROFILE_INPUTS]

    readiness_fields = [
        "case",
        "role",
        "profile_csv",
        "summary_json",
        "summary_Pr",
        "summary_Re_tau",
        "summary_T_tau",
        "summary_covariance_snapshots",
        "profile_points",
        "first_cell_y_plus",
        "max_y_plus",
        "max_T_plus",
        "y_plus_at_max_T_plus",
        "max_T_rms_plus",
        "finite_Pr_t_points",
        "positive_Pr_t_points",
        "positive_Pr_t_fraction",
        "positive_Pr_t_median",
        "positive_Pr_t_midband_median",
        "certification_use",
    ]
    target_fields = ["gate", "benchmark_target", "current_state", "needed_next", "status"]

    write_csv(RESULTS / "c1_certification_readiness.csv", readiness, readiness_fields)
    write_csv(RESULTS / "c1_kawamura_certification_targets.csv", TARGET_ROWS, target_fields)
    md_path = DOCS / "c1_kawamura_certification_targets.md"
    md_path.write_text(readiness_markdown(readiness))

    print(f"wrote {RESULTS / 'c1_certification_readiness.csv'}")
    print(f"wrote {RESULTS / 'c1_kawamura_certification_targets.csv'}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

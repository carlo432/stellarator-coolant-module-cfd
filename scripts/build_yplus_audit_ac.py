#!/usr/bin/env python3
"""Build the AC heated-wall y+ audit table and figure."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path
from statistics import mean

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/branch_d"
ADVISOR_DIR = ROOT / "results/report_figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
ADVISOR_DIR.mkdir(parents=True, exist_ok=True)


SOURCE_CASES = {
    "v13": {
        "case": "cases/baseline_module/openfoam_cases/v13_isothermal_rans_flibe_re10000_thermal_refined",
        "time": "1000",
        "solver_family": "RANS kOmegaSST",
        "role": "straight source flow for v14",
        "reported_deltaT": "source flow only",
    },
    "v22": {
        "case": "cases/branch_d/v22_isothermal_rans_flibe_re10000_bent_recirculation",
        "time": "1000",
        "solver_family": "RANS kOmegaSST",
        "role": "outlet-offset source flow for v23/v24",
        "reported_deltaT": "source flow only",
    },
    "v25": {
        "case": "cases/branch_d/v25_unsteady_rans_flibe_re10000_bent_recirculation_pimple",
        "time": "0.400106",
        "solver_family": "URANS kOmegaSST",
        "role": "outlet-offset unsteady RANS check",
        "reported_deltaT": "flow-only unsteady check",
    },
    "v26": {
        "case": "cases/branch_d/v26_les_smoke_flibe_re10000_bent_recirculation_wale",
        "time": "0.800034",
        "solver_family": "LES WALE",
        "role": "outlet-offset LES smoke test",
        "reported_deltaT": "isothermal LES",
    },
    "v27": {
        "case": "cases/branch_d/v27_thermal_les_flibe_re10000_bent_recirculation_wale_100kw",
        "time": "0.500073",
        "solver_family": "LES WALE + passive T",
        "role": "outlet-offset thermal LES extension",
        "reported_deltaT": "near-wall fluid T fluctuation only",
    },
    "v28": {
        "case": "cases/baseline_module/openfoam_cases/v28_isothermal_rans_flibe_re10000_bend90",
        "time": "1000",
        "solver_family": "RANS kOmegaSST",
        "role": "90-degree-bend source flow for v29",
        "reported_deltaT": "source flow only",
    },
    "v30": {
        "case": "cases/baseline_module/openfoam_cases/v30_isothermal_rans_flibe_re10000_bfs",
        "time": "1000",
        "solver_family": "RANS kOmegaSST",
        "role": "BFS source flow for v31",
        "reported_deltaT": "source flow only",
    },
}

INHERITED_CASES = {
    "v14": ("v13", "passive scalar wall dT inherits v13 flow y+"),
    "v23": ("v22", "passive scalar wall dT inherits v22 flow y+"),
    "v24": ("v22", "passive scalar 500 kW/m2 wall dT inherits v22 flow y+"),
    "v29": ("v28", "passive scalar wall dT inherits v28 flow y+"),
    "v31": ("v30", "passive scalar wall dT inherits v30 flow y+"),
}


def read_patch_values(path: Path, patch: str = "heated_wall") -> list[float]:
    text = path.read_text()
    patch_match = re.search(rf"\b{re.escape(patch)}\s*\{{(?P<body>.*?)\n\s*\}}", text, re.S)
    if not patch_match:
        raise ValueError(f"patch {patch!r} not found in {path}")
    body = patch_match.group("body")
    if "value           uniform" in body or "value uniform" in body:
        uniform_match = re.search(r"value\s+uniform\s+([-+0-9.eE]+)", body)
        if not uniform_match:
            raise ValueError(f"uniform value missing for {patch!r} in {path}")
        return [float(uniform_match.group(1))]

    list_match = re.search(r"value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", body, re.S)
    if not list_match:
        raise ValueError(f"nonuniform values missing for {patch!r} in {path}")
    expected = int(list_match.group(1))
    values = [float(item) for item in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", list_match.group(2))]
    if len(values) != expected:
        raise ValueError(f"expected {expected} values in {path}, got {len(values)}")
    return values


def percentile(values: list[float], pct: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    idx = (len(ordered) - 1) * pct / 100.0
    lower = int(idx)
    upper = min(lower + 1, len(ordered) - 1)
    frac = idx - lower
    return ordered[lower] * (1.0 - frac) + ordered[upper] * frac


def band_fractions(values: list[float]) -> dict[str, float]:
    n = len(values)
    return {
        "frac_lt_1": sum(value < 1.0 for value in values) / n,
        "frac_1_to_5": sum(1.0 <= value < 5.0 for value in values) / n,
        "frac_5_to_30": sum(5.0 <= value < 30.0 for value in values) / n,
        "frac_30_to_300": sum(30.0 <= value <= 300.0 for value in values) / n,
        "frac_gt_300": sum(value > 300.0 for value in values) / n,
    }


def classify(row: dict[str, float | str]) -> str:
    if float(row["frac_30_to_300"]) >= 0.8:
        return "mostly high-Re wall-function band"
    if float(row["yplus_max"]) < 30.0:
        return "outside high-Re wall-function band; heated wall below y+=30"
    if float(row["frac_5_to_30"]) > 0.5:
        return "mostly buffer-layer"
    return "mixed near-wall regime"


def source_row(case_id: str, meta: dict[str, str]) -> dict[str, str]:
    path = ROOT / meta["case"] / meta["time"] / "yPlus"
    values = read_patch_values(path)
    bands = band_fractions(values)
    raw: dict[str, float | str] = {
        "case": case_id,
        "source_flow_case": case_id,
        "solver_family": meta["solver_family"],
        "role": meta["role"],
        "time": meta["time"],
        "n_heated_wall_faces": len(values),
        "yplus_min": min(values),
        "yplus_mean": mean(values),
        "yplus_p95": percentile(values, 95.0),
        "yplus_max": max(values),
        **bands,
        "wall_function_validity": "",
        "reported_deltaT_status": meta["reported_deltaT"],
        "note": "direct yPlus post-processing",
    }
    raw["wall_function_validity"] = classify(raw)
    return stringify(raw)


def inherited_row(case_id: str, source_id: str, note: str, source: dict[str, str]) -> dict[str, str]:
    row = dict(source)
    row["case"] = case_id
    row["source_flow_case"] = source_id
    row["role"] = note
    row["reported_deltaT_status"] = "inherits source-flow y+; passive scalar fixedGradient T"
    row["note"] = "inherited from source flow; thermal case lacks turbulence fields"
    return row


def stringify(row: dict[str, float | str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, value in row.items():
        if isinstance(value, float):
            if key.startswith("frac_"):
                out[key] = f"{100.0 * value:.1f}"
            else:
                out[key] = f"{value:.3f}"
        else:
            out[key] = str(value)
    return out


def build_rows() -> list[dict[str, str]]:
    source_rows = {case_id: source_row(case_id, meta) for case_id, meta in SOURCE_CASES.items()}
    rows: list[dict[str, str]] = []
    order = ["v13", "v14", "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31"]
    for case_id in order:
        if case_id in source_rows:
            rows.append(source_rows[case_id])
        else:
            source_id, note = INHERITED_CASES[case_id]
            rows.append(inherited_row(case_id, source_id, note, source_rows[source_id]))
    return rows


def write_csv(rows: list[dict[str, str]]) -> Path:
    out = OUT_DIR / "yplus_audit_ac.csv"
    fieldnames = list(rows[0])
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out


def plot(rows: list[dict[str, str]]) -> tuple[Path, Path]:
    unique = [row for row in rows if row["case"] == row["source_flow_case"]]
    labels = [row["case"] for row in unique]
    bands = [
        ("frac_lt_1", "<1", "#2f9e44"),
        ("frac_1_to_5", "1-5", "#74b816"),
        ("frac_5_to_30", "5-30 buffer", "#f08c00"),
        ("frac_30_to_300", "30-300 wall-fn", "#1971c2"),
        ("frac_gt_300", ">300", "#c92a2a"),
    ]

    fig, ax = plt.subplots(figsize=(11.2, 5.8))
    fig.patch.set_facecolor("white")
    left = [0.0] * len(unique)
    for key, label, color in bands:
        values = [float(row[key]) for row in unique]
        ax.barh(labels, values, left=left, label=label, color=color, alpha=0.85)
        left = [base + value for base, value in zip(left, values)]

    means = [float(row["yplus_mean"]) for row in unique]
    maxes = [float(row["yplus_max"]) for row in unique]
    for idx, (mean_value, max_value) in enumerate(zip(means, maxes)):
        ax.text(101.0, idx, f"avg {mean_value:.1f}, max {max_value:.1f}", va="center", fontsize=8)

    ax.set_xlim(0, 118)
    ax.set_xlabel("heated-wall faces by y+ band [%]")
    ax.set_title("AC: Heated-Wall y+ Audit for Branch D Source Flow Cases")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=5, fontsize=8)
    ax.grid(axis="x", alpha=0.25)
    ax.text(
        0.0,
        -0.18,
        "Rows show unique source-flow cases; passive thermal cases v14/v23/v24/v29/v31 inherit these y+ fields. "
        "High-Re wall-function target shown as 30-300.",
        transform=ax.transAxes,
        fontsize=9,
        color="#495057",
    )
    fig.tight_layout(rect=[0, 0.08, 1, 1])

    out = OUT_DIR / "yplus_audit_ac.png"
    advisor = ADVISOR_DIR / "33_yplus_audit.png"
    fig.savefig(out, dpi=180)
    fig.savefig(advisor, dpi=180)
    plt.close(fig)
    return out, advisor


def main() -> None:
    rows = build_rows()
    csv_out = write_csv(rows)
    fig_out, advisor_out = plot(rows)
    print(csv_out)
    print(fig_out)
    print(advisor_out)


if __name__ == "__main__":
    main()

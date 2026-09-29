#!/usr/bin/env python3
"""Overlay C1 profiles against digitized Kawamura DNS profiles when available."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
FIGS = WORK / "figs"
DIGITIZATION = ROOT / "references" / "digitization" / "kawamura_1998"
DIGITIZED = DIGITIZATION / "digitized"


TARGETS = [
    {
        "id": "fig2_mean_temperature",
        "dns_csv": DIGITIZED / "kawamura_fig2_mean_temperature.csv",
        "dns_y": "Theta_plus",
        "c1_y": "T_plus",
        "ylabel": "Theta+",
        "out_png": "c1_vs_kawamura_fig2_mean_temperature.png",
        "normalization": "direct C1 T_plus vs digitized Theta_plus",
    },
    {
        "id": "fig3_near_wall_temperature",
        "dns_csv": DIGITIZED / "kawamura_fig3_near_wall_temperature.csv",
        "dns_y": "Theta_plus_over_Pr",
        "c1_y": "T_plus_over_Pr",
        "ylabel": "Theta+ / Pr",
        "out_png": "c1_vs_kawamura_fig3_near_wall_temperature.png",
        "normalization": "C1 T_plus divided by profile Pr",
    },
    {
        "id": "fig6_temperature_rms",
        "dns_csv": DIGITIZED / "kawamura_fig6_temperature_rms.csv",
        "dns_y": "theta_rms_plus_over_Pr",
        "c1_y": "T_rms_plus_over_Pr",
        "ylabel": "theta_rms+ / Pr",
        "out_png": "c1_vs_kawamura_fig6_temperature_rms.png",
        "normalization": "C1 T_rms_plus divided by profile Pr",
    },
    {
        "id": "fig8_wall_normal_heat_flux",
        "dns_csv": DIGITIZED / "kawamura_fig8_wall_normal_heat_flux.csv",
        "dns_y": "minus_v_theta_plus_over_Pr",
        "c1_y": "minus_v_theta_plus_over_Pr",
        "ylabel": "-v+ theta+ / Pr",
        "out_png": "c1_vs_kawamura_fig8_wall_normal_heat_flux.png",
        "normalization": "C1 -vT_resolved / (u_tau * T_tau * Pr), using harness summary sidecar when available",
    },
    {
        "id": "fig9_turbulent_prandtl",
        "dns_csv": DIGITIZED / "kawamura_fig9_turbulent_prandtl.csv",
        "dns_y": "Pr_t",
        "c1_y": "Pr_t_resolved",
        "ylabel": "Pr_t",
        "out_png": "c1_vs_kawamura_fig9_turbulent_prandtl.png",
        "normalization": "resolved positive eddy-diffusivity Pr_t points only",
    },
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def fval(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, ""))
    except ValueError:
        return math.nan


def finite_xy(rows: list[dict[str, str]], xkey: str, ykey: str) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        x = fval(row, xkey)
        y = fval(row, ykey)
        if math.isfinite(x) and math.isfinite(y):
            xs.append(x)
            ys.append(y)
    return xs, ys


def infer_utau(rows: list[dict[str, str]], nu: float) -> float:
    values: list[float] = []
    for row in rows:
        y = fval(row, "y")
        y_plus = fval(row, "y_plus")
        if math.isfinite(y) and y > 0 and math.isfinite(y_plus):
            values.append(y_plus * nu / y)
    if not values:
        return math.nan
    values.sort()
    return values[len(values) // 2]


def infer_wall_gradient(rows: list[dict[str, str]]) -> float:
    pairs = [
        (fval(row, "y"), fval(row, "T_mean"))
        for row in rows
        if math.isfinite(fval(row, "y")) and math.isfinite(fval(row, "T_mean"))
    ]
    pairs.sort()
    if len(pairs) < 2:
        return math.nan
    dy = pairs[1][0] - pairs[0][0]
    if dy == 0:
        return math.nan
    return (pairs[1][1] - pairs[0][1]) / dy


def add_derived_c1_columns(rows: list[dict[str, str]], pr: float, nu: float) -> dict[str, float]:
    alpha = nu / pr if pr else math.nan
    utau = infer_utau(rows, nu)
    wall_grad = infer_wall_gradient(rows)
    ttau = abs(alpha * wall_grad) / utau if math.isfinite(alpha) and math.isfinite(wall_grad) and math.isfinite(utau) and utau > 0 else math.nan
    for row in rows:
        t_plus = fval(row, "T_plus")
        t_rms_plus = fval(row, "T_rms_plus")
        vt = fval(row, "vT_resolved")
        row["T_plus_over_Pr"] = str(t_plus / pr) if math.isfinite(t_plus) and pr else ""
        row["T_rms_plus_over_Pr"] = str(t_rms_plus / pr) if math.isfinite(t_rms_plus) and pr else ""
        if math.isfinite(vt) and math.isfinite(utau) and math.isfinite(ttau) and utau > 0 and ttau > 0 and pr:
            row["minus_v_theta_plus_over_Pr"] = str(-vt / (utau * ttau * pr))
        else:
            row["minus_v_theta_plus_over_Pr"] = ""
    return {"nu": nu, "alpha": alpha, "u_tau": utau, "wall_dTdy": wall_grad, "T_tau": ttau}


def default_summary_path(profile: Path) -> Path:
    name = profile.name
    if name.startswith("c1_thermal_profile_"):
        return profile.with_name(name.replace("c1_thermal_profile_", "c1_thermal_summary_", 1)).with_suffix(".json")
    return profile.with_suffix(".json")


def load_harness_summary(path: Path | None) -> dict[str, float | str] | None:
    if path is None or not path.exists():
        return None
    return json.loads(path.read_text())


def add_derived_from_summary_or_profile(
    rows: list[dict[str, str]],
    pr: float,
    nu: float,
    summary: dict[str, float | str] | None,
) -> dict[str, float | str | None]:
    profile_derived = add_derived_c1_columns(rows, pr, nu)
    if not summary:
        profile_derived["source"] = "profile_inference"
        return profile_derived

    utau = float(summary.get("u_tau", math.nan))
    ttau = float(summary.get("T_tau", math.nan))
    summary_pr = float(summary.get("Pr", pr))
    for row in rows:
        vt = fval(row, "vT_resolved")
        if math.isfinite(vt) and math.isfinite(utau) and math.isfinite(ttau) and utau > 0 and ttau > 0 and summary_pr:
            row["minus_v_theta_plus_over_Pr"] = str(-vt / (utau * ttau * summary_pr))
    return {
        "source": str(path) if (path := summary.get("summary_path")) else "harness_summary",
        "nu": float(summary.get("nu", nu)),
        "alpha": float(summary.get("alpha", nu / pr if pr else math.nan)),
        "u_tau": utau,
        "wall_dTdy": float(summary.get("wall_dTdy", math.nan)),
        "T_tau": ttau,
        "Re_tau": float(summary.get("Re_tau", math.nan)),
    }


def plot_target(c1_rows: list[dict[str, str]], dns_rows: list[dict[str, str]], target: dict[str, str], outdir: Path) -> str:
    import matplotlib.pyplot as plt

    c1_x, c1_y = finite_xy(c1_rows, "y_plus", target["c1_y"])
    dns_x, dns_y = finite_xy(dns_rows, "y_plus", target["dns_y"])
    if not c1_x or not dns_x:
        return "skipped_no_data"

    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    ax.semilogx(dns_x, dns_y, "o", ms=4, label="Kawamura digitized")
    ax.semilogx(c1_x, c1_y, "-", lw=1.8, label="C1 local profile")
    ax.set_xlabel("y+")
    ax.set_ylabel(target["ylabel"])
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    ax.set_title(target["id"])
    out = outdir / target["out_png"]
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return str(out.relative_to(ROOT))


def write_markdown(summary: dict, path: Path) -> None:
    lines = [
        "# C1 Kawamura Overlay Readiness",
        "",
        "This file is generated by `cases/branch_d/compare_c1_kawamura_profiles.py`.",
        "",
        f"- C1 profile: `{summary['c1_profile']}`",
        f"- C1 summary: `{summary['c1_summary']}`",
        f"- Assumed C1 Pr: `{summary['pr']}`",
        f"- Normalization source: `{summary['derived_normalization']['source']}`",
        f"- `u_tau`: `{summary['derived_normalization']['u_tau']:.6g}`",
        f"- `T_tau`: `{summary['derived_normalization']['T_tau']:.6g}`",
        f"- near-wall `dT/dy`: `{summary['derived_normalization']['wall_dTdy']:.6g}`",
        "",
        "| Target | DNS rows | C1 rows | Plot status | Output |",
        "|---|---:|---:|---|---|",
    ]
    for row in summary["targets"]:
        lines.append(
            f"| `{row['id']}` | {row['dns_rows']} | {row['c1_rows']} | "
            f"`{row['plot_status']}` | {row['output']} |"
        )
    lines.extend(
        [
            "",
            "Interpretation: empty DNS rows mean the Kawamura profile has not been",
            "digitized yet.  Any generated plot remains a machinery overlay, not a",
            "certification result, until the C1 run is statistically meaningful and",
            "boundary-condition matched.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--c1-profile", default=str(FIGS / "c1_thermal_profile_les_channel_c1_thermal_pr5_preflight.csv"))
    ap.add_argument("--c1-summary", default=None, help="Optional C1 harness summary JSON. Defaults to sidecar path.")
    ap.add_argument("--pr", type=float, default=5.0)
    ap.add_argument("--nu", type=float, default=3.09278e-6)
    ap.add_argument("--out", default=str(FIGS / "c1_kawamura_overlay_readiness"))
    args = ap.parse_args()

    c1_profile = Path(args.c1_profile)
    outbase = Path(args.out)
    outdir = outbase.parent
    outdir.mkdir(parents=True, exist_ok=True)

    c1_rows = read_csv(c1_profile)
    summary_path = Path(args.c1_summary) if args.c1_summary else default_summary_path(c1_profile)
    harness_summary = load_harness_summary(summary_path)
    if harness_summary is not None:
        harness_summary["summary_path"] = str(summary_path)
    derived = add_derived_from_summary_or_profile(c1_rows, args.pr, args.nu, harness_summary)

    summary = {
        "c1_profile": str(c1_profile.relative_to(ROOT)) if c1_profile.is_absolute() else str(c1_profile),
        "c1_summary": str(summary_path.relative_to(ROOT)) if summary_path.is_absolute() and summary_path.exists() else str(summary_path),
        "pr": args.pr,
        "derived_normalization": derived,
        "targets": [],
    }
    for target in TARGETS:
        dns_rows = read_csv(Path(target["dns_csv"]))
        plot_status = plot_target(c1_rows, dns_rows, target, outdir)
        output = "" if plot_status.startswith("skipped_") else plot_status
        summary["targets"].append(
            {
                "id": target["id"],
                "dns_rows": len(dns_rows),
                "c1_rows": len(c1_rows),
                "normalization": target["normalization"],
                "plot_status": plot_status,
                "output": output,
            }
        )

    json_path = outbase.with_suffix(".json")
    md_path = outbase.with_suffix(".md")
    json_path.write_text(json.dumps(summary, indent=2))
    write_markdown(summary, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

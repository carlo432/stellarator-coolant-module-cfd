#!/usr/bin/env python3
"""Compare C2 curved CHT wall-temperature scale against duct correlations."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def scalar_from_thermo(text: str, key: str) -> float:
    match = re.search(rf"\b{re.escape(key)}\s+([-+0-9.eE]+)\s*;", text)
    if not match:
        raise ValueError(f"could not parse {key} from thermophysicalProperties")
    return float(match.group(1))


def patch_dt(summary: dict, patch: str) -> float:
    return summary["patch_averages"][patch]["area_average"] - summary["tin_K"]


def gnielinski_nu(reynolds: float, prandtl: float, friction_factor: float) -> float:
    return (
        (friction_factor / 8.0) * (reynolds - 1000.0) * prandtl
    ) / (
        1.0 + 12.7 * math.sqrt(friction_factor / 8.0) * (prandtl ** (2.0 / 3.0) - 1.0)
    )


def petukhov_friction(reynolds: float) -> float:
    return (0.79 * math.log(reynolds) - 1.64) ** -2


def write_markdown(result: dict, path: Path) -> None:
    inputs = result["inputs"]
    correlations = result["correlations"]
    rows = result["rows"]
    lines = [
        "# C2 Curved CHT Correlation Sanity Check",
        "",
        "This compares the C2 wall-resolved curved-slice CHT temperature scale against",
        "straight turbulent-duct heat-transfer correlations. It is an engineering",
        "sanity check, not experimental validation and not a reactor-geometry claim.",
        "",
        "## Inputs",
        "",
        "| Quantity | Value |",
        "|---|---:|",
        f"| radial width | `{inputs['radial_width_m']:.6g} m` |",
        f"| height | `{inputs['height_m']:.6g} m` |",
        f"| hydraulic diameter | `{inputs['hydraulic_diameter_m']:.6g} m` |",
        f"| nominal inlet speed | `{inputs['u_in_m_s']:.6g} m/s` |",
        f"| Reynolds number | `{inputs['reynolds']:.6g}` |",
        f"| Prandtl number | `{inputs['prandtl']:.6g}` |",
        f"| fluid conductivity from `Cp*mu/Pr` | `{inputs['fluid_k_W_mK']:.6g} W/m/K` |",
        f"| heat flux | `{inputs['heat_flux_W_m2']:.6g} W/m^2` |",
        f"| solid conduction drop `q*t/k_s` | `{inputs['solid_conduction_drop_K']:.6g} K` |",
        "",
        "## Correlation Anchors",
        "",
        "| Correlation | Nu | Film dT, K | Base dT with solid drop, K |",
        "|---|---:|---:|---:|",
    ]
    for key in ("dittus_boelter", "gnielinski"):
        item = correlations[key]
        lines.append(
            f"| {item['label']} | `{item['Nu']:.6g}` | `{item['film_dT_K']:.6g}` | "
            f"`{item['base_dT_K']:.6g}` |"
        )
    lines.extend(
        [
            "",
            "## CHT Results Versus Gnielinski",
            "",
            "| Case | Interface dT, K | Interface Nu | Interface dT vs Gnielinski | Base dT, K | Base dT vs Gnielinski+solid |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['label']} | `{row['interface_dT_K']:.6g}` | `{row['interface_Nu']:.6g}` | "
            f"`{row['interface_dT_vs_gnielinski_pct']:+.2f}%` | `{row['heatedBase_dT_K']:.6g}` | "
            f"`{row['base_dT_vs_gnielinski_plus_solid_pct']:+.2f}%` |"
        )
    lines.extend(["", "## Interpretation", "", result["interpretation"], ""])
    path.write_text("\n".join(lines))


def plot_result(result: dict, path: Path) -> None:
    labels = ["Dittus", "Gnielinski"] + [row["short_label"] for row in result["rows"]]
    interface = [
        result["correlations"]["dittus_boelter"]["film_dT_K"],
        result["correlations"]["gnielinski"]["film_dT_K"],
        *[row["interface_dT_K"] for row in result["rows"]],
    ]
    base = [
        result["correlations"]["dittus_boelter"]["base_dT_K"],
        result["correlations"]["gnielinski"]["base_dT_K"],
        *[row["heatedBase_dT_K"] for row in result["rows"]],
    ]
    x = range(len(labels))
    fig, ax = plt.subplots(figsize=(10, 5.8))
    width = 0.38
    ax.bar([i - width / 2 for i in x], interface, width=width, label="film/interface dT", color="#377eb8")
    ax.bar([i + width / 2 for i in x], base, width=width, label="base dT incl. solid", color="#ff7f00")
    ax.set_ylabel("dT over inlet [K]")
    ax.set_title("C2 curved CHT scale vs duct correlation anchors")
    ax.set_xticks(list(x), labels, rotation=18, ha="right")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="c2_curved_stellarator_slice_cht_wallres_xfine")
    ap.add_argument(
        "--summaries",
        nargs="+",
        default=[
            "figs/c2_curved_stellarator_slice_cht_wallres_longrun_summary.json",
            "figs/c2_curved_stellarator_slice_cht_wallres_fine_summary.json",
            "figs/c2_curved_stellarator_slice_cht_wallres_xfine_summary.json",
        ],
    )
    ap.add_argument("--out", default="figs/c2_curved_cht_correlation_validation")
    args = ap.parse_args()

    case = Path(args.case)
    setup = load_json(case / "geometry_c2cht.json")
    params = setup["parameters"]
    thermo = (case / "constant/bottomWater/thermophysicalProperties").read_text()
    rho = scalar_from_thermo(thermo, "rho")
    cp = scalar_from_thermo(thermo, "Cp")
    mu = scalar_from_thermo(thermo, "mu")
    pr = scalar_from_thermo(thermo, "Pr")
    fluid_k = cp * mu / pr
    width = float(params["radial_width"])
    height = float(params["height"])
    dh = 2.0 * width * height / (width + height)
    u_in = float(params["u_in"])
    heat_flux = float(params["heat_flux"])
    solid_drop = heat_flux * float(params["solid_thickness"]) / float(params["solid_kappa"])
    reynolds = rho * u_in * dh / mu
    friction = petukhov_friction(reynolds)
    nu_db = 0.023 * reynolds**0.8 * pr**0.4
    nu_gn = gnielinski_nu(reynolds, pr, friction)

    def corr(label: str, nu: float) -> dict[str, float | str]:
        film_dt = heat_flux * dh / (nu * fluid_k)
        return {
            "label": label,
            "Nu": nu,
            "film_dT_K": film_dt,
            "base_dT_K": film_dt + solid_drop,
        }

    correlations = {
        "dittus_boelter": corr("Dittus-Boelter", nu_db),
        "gnielinski": corr("Gnielinski", nu_gn),
    }
    summaries = [load_json(Path(path)) for path in args.summaries]
    rows = []
    gn = correlations["gnielinski"]
    for summary in summaries:
        interface_dt = patch_dt(summary, "heater_to_bottomWater")
        base_dt = patch_dt(summary, "heatedBase")
        interface_nu = heat_flux * dh / (interface_dt * fluid_k)
        name = Path(summary["case"]).name
        short = name.replace("c2_curved_stellarator_slice_cht_", "").replace("wallres_", "")
        rows.append(
            {
                "case": summary["case"],
                "label": name,
                "short_label": short,
                "fluid_cells": summary["mesh"]["bottomWater_cells"],
                "solid_cells": summary["mesh"]["heater_cells"],
                "interface_dT_K": interface_dt,
                "heatedBase_dT_K": base_dt,
                "interface_Nu": interface_nu,
                "interface_dT_vs_gnielinski_pct": 100.0 * (interface_dt - gn["film_dT_K"]) / gn["film_dT_K"],
                "interface_Nu_vs_gnielinski_pct": 100.0 * (interface_nu - gn["Nu"]) / gn["Nu"],
                "base_dT_vs_gnielinski_plus_solid_pct": 100.0 * (base_dt - gn["base_dT_K"]) / gn["base_dT_K"],
            }
        )

    max_abs_interface = max(abs(row["interface_dT_vs_gnielinski_pct"]) for row in rows)
    max_abs_base = max(abs(row["base_dT_vs_gnielinski_plus_solid_pct"]) for row in rows)
    interpretation = (
        "The C2 wall-resolved CHT interface dT is about 25-26% below the straight-duct "
        "Gnielinski film estimate, while the heated-base dT including the analytic "
        "solid drop is about 17% below the Gnielinski-plus-solid estimate. This is "
        "same-order agreement, not exact validation. The offset is plausible for a "
        "short curved slice with developing/curved flow, but it should be framed as "
        "a correlation sanity check only."
    )
    if max_abs_interface <= 20.0 and max_abs_base <= 20.0:
        interpretation = (
            "The C2 wall-resolved CHT scale is within about 20% of the straight-duct "
            "correlation anchors. This is useful same-order support, not validation, "
            "because the C2 geometry is short, curved, and not a developed straight duct."
        )

    result = {
        "inputs": {
            "radial_width_m": width,
            "height_m": height,
            "hydraulic_diameter_m": dh,
            "u_in_m_s": u_in,
            "rho_kg_m3": rho,
            "cp_J_kgK": cp,
            "mu_Pa_s": mu,
            "prandtl": pr,
            "fluid_k_W_mK": fluid_k,
            "heat_flux_W_m2": heat_flux,
            "solid_conduction_drop_K": solid_drop,
            "reynolds": reynolds,
            "petukhov_friction_factor": friction,
        },
        "correlations": correlations,
        "rows": rows,
        "interpretation": interpretation,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json_path = out.with_suffix(".json")
    md_path = out.with_suffix(".md")
    csv_path = out.with_suffix(".csv")
    png_path = out.with_suffix(".png")
    json_path.write_text(json.dumps(result, indent=2))
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "fluid_cells",
                "solid_cells",
                "interface_dT_K",
                "heatedBase_dT_K",
                "interface_Nu",
                "interface_dT_vs_gnielinski_pct",
                "interface_Nu_vs_gnielinski_pct",
                "base_dT_vs_gnielinski_plus_solid_pct",
            ],
        )
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    write_markdown(result, md_path)
    plot_result(result, png_path)
    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()

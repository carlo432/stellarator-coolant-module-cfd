#!/usr/bin/env python3
"""Compare C2 curved and straight-control smoke summaries."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def get_speed(data: dict) -> float:
    return data["patch_owner_speed_m_per_s"]["heated_first_wall"]["area_weighted_average"]


def get_area(data: dict) -> float:
    return data["patch_owner_speed_m_per_s"]["heated_first_wall"]["area_m2"]


def get_dp(data: dict) -> float:
    return data["pressure_drop_owner_cell_m2_per_s2"]


def get_dts(data: dict) -> tuple[float, float, float]:
    metric = data["patch_owner_temperature_K"]["heated_first_wall"]
    outlet = data["patch_owner_temperature_K"]["outlet"]
    tin = data["tin_K"]
    return (
        metric["area_weighted_average"] - tin,
        metric["max"] - tin,
        outlet["area_weighted_average"] - tin,
    )


def get_setup_value(data: dict, key: str) -> float | str | None:
    setup = data.get("thermal_smoke_setup", {})
    return setup.get(key)


def pct(new: float, ref: float) -> float:
    return 100.0 * (new - ref) / ref


def fmt_optional(value: float | None, suffix: str = "") -> str:
    if value is None:
        return "n/a"
    return f"{value:.6g}{suffix}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--curved-flow", default="figs/c2_curved_stellarator_slice_preflight_smoke_summary.json")
    ap.add_argument("--curved-thermal", default="figs/c2_curved_stellarator_slice_thermal_smoke_summary.json")
    ap.add_argument("--straight-flow", default="figs/c2_straight_control_preflight_smoke_summary.json")
    ap.add_argument("--straight-thermal", default="figs/c2_straight_control_thermal_smoke_summary.json")
    ap.add_argument("--out", default="figs/c2_curved_vs_straight_control_summary")
    ap.add_argument("--title", default="C2 Curved vs Straight Control Summary")
    args = ap.parse_args()

    curved_flow = load(Path(args.curved_flow))
    curved_thermal = load(Path(args.curved_thermal))
    straight_flow = load(Path(args.straight_flow))
    straight_thermal = load(Path(args.straight_thermal))

    curved_speed = get_speed(curved_flow)
    straight_speed = get_speed(straight_flow)
    curved_area = get_area(curved_flow)
    straight_area = get_area(straight_flow)
    curved_dp = get_dp(curved_flow)
    straight_dp = get_dp(straight_flow)
    curved_dt_avg, curved_dt_max, curved_outlet = get_dts(curved_thermal)
    straight_dt_avg, straight_dt_max, straight_outlet = get_dts(straight_thermal)
    curved_product = get_setup_value(curved_thermal, "gradient_area_product")
    straight_product = get_setup_value(straight_thermal, "gradient_area_product")
    curved_gradient = get_setup_value(curved_thermal, "applied_heated_wall_gradient")
    straight_gradient = get_setup_value(straight_thermal, "applied_heated_wall_gradient")
    product_delta_percent = None
    if isinstance(curved_product, float) and isinstance(straight_product, float) and straight_product:
        product_delta_percent = pct(curved_product, straight_product)

    equal_power_basis = product_delta_percent is not None and abs(product_delta_percent) < 0.01
    thermal_basis = (
        "matched heated-wall gradient-area product, a passive equal-total-power surrogate"
        if equal_power_basis
        else "same fixed-gradient passive-scalar smoke"
    )

    summary = {
        "thermal_basis": thermal_basis,
        "curved": {
            "heated_wall_adjacent_speed_m_per_s": curved_speed,
            "heated_wall_area_m2": curved_area,
            "applied_heated_wall_gradient": curved_gradient,
            "gradient_area_product": curved_product,
            "pressure_drop_m2_per_s2": curved_dp,
            "heated_wall_adjacent_dT_avg_K": curved_dt_avg,
            "heated_wall_adjacent_dT_max_K": curved_dt_max,
            "outlet_adjacent_dT_avg_K": curved_outlet,
        },
        "straight": {
            "heated_wall_adjacent_speed_m_per_s": straight_speed,
            "heated_wall_area_m2": straight_area,
            "applied_heated_wall_gradient": straight_gradient,
            "gradient_area_product": straight_product,
            "pressure_drop_m2_per_s2": straight_dp,
            "heated_wall_adjacent_dT_avg_K": straight_dt_avg,
            "heated_wall_adjacent_dT_max_K": straight_dt_max,
            "outlet_adjacent_dT_avg_K": straight_outlet,
        },
        "curved_minus_straight": {
            "speed_delta_m_per_s": curved_speed - straight_speed,
            "speed_delta_percent": pct(curved_speed, straight_speed),
            "heated_wall_area_delta_m2": curved_area - straight_area,
            "heated_wall_area_delta_percent": pct(curved_area, straight_area),
            "gradient_area_product_delta_percent": product_delta_percent,
            "pressure_drop_delta_m2_per_s2": curved_dp - straight_dp,
            "pressure_drop_delta_percent": pct(curved_dp, straight_dp),
            "heated_wall_adjacent_dT_avg_delta_K": curved_dt_avg - straight_dt_avg,
            "heated_wall_adjacent_dT_avg_delta_percent": pct(curved_dt_avg, straight_dt_avg),
            "heated_wall_adjacent_dT_max_delta_K": curved_dt_max - straight_dt_max,
            "outlet_adjacent_dT_avg_delta_K": curved_outlet - straight_outlet,
        },
    }

    out_stem = Path(args.out)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    json_path = out_stem.with_suffix(".json")
    md_path = out_stem.with_suffix(".md")
    json_path.write_text(json.dumps(summary, indent=2))
    setup_rows = ""
    if isinstance(curved_gradient, float) and isinstance(straight_gradient, float):
        setup_rows += (
            f"| Applied heated-wall gradient | {straight_gradient:.6g} | "
            f"{curved_gradient:.6g} | {curved_gradient - straight_gradient:.6g} |\n"
        )
    if isinstance(curved_product, float) and isinstance(straight_product, float):
        setup_rows += (
            f"| Gradient-area product | {straight_product:.6g} | {curved_product:.6g} | "
            f"{curved_product - straight_product:.6g} ({fmt_optional(product_delta_percent, '%')}) |\n"
        )
    if equal_power_basis:
        interpretation = (
            "Interpretation: this is still a short passive-scalar smoke comparison, "
            "not validated CHT or LES, but it removes the first fixed-gradient area "
            "caveat by matching the heated-wall gradient-area product. The remaining "
            "signal is therefore more useful as a curved-vs-straight diagnostic than "
            "the same-gradient comparison."
        )
    else:
        interpretation = (
            "Interpretation: this is the first C2 control comparison. It is useful "
            "for workflow and relative curved-vs-straight diagnostics, but it is "
            "still a short passive-scalar smoke comparison, not a validated CHT or "
            "LES result. Because the heated-wall areas differ, the average dT result "
            "should be read as a fixed-gradient smoke diagnostic, not a clean "
            "equal-total-power comparison."
        )
    md_path.write_text(
        f"# {args.title}\n\n"
        "This compares local smoke tests only: same-scale C2 curved/helical slice "
        "against a straight duct control. Metrics are owner-cell adjacent values, "
        "not wall-patch velocity and not CHT wall temperature.\n\n"
        f"Thermal basis: {thermal_basis}.\n\n"
        "| Metric | Straight | Curved | Curved - Straight |\n"
        "|---|---:|---:|---:|\n"
        f"| Heated-wall adjacent speed, m/s | {straight_speed:.6g} | {curved_speed:.6g} | {curved_speed - straight_speed:.6g} ({pct(curved_speed, straight_speed):+.2f}%) |\n"
        f"| Heated-wall area, m^2 | {straight_area:.6g} | {curved_area:.6g} | {curved_area - straight_area:.6g} ({pct(curved_area, straight_area):+.2f}%) |\n"
        f"{setup_rows}"
        f"| Pressure drop, m^2/s^2 | {straight_dp:.6g} | {curved_dp:.6g} | {curved_dp - straight_dp:.6g} ({pct(curved_dp, straight_dp):+.2f}%) |\n"
        f"| Heated-wall adjacent dT avg, K | {straight_dt_avg:.6g} | {curved_dt_avg:.6g} | {curved_dt_avg - straight_dt_avg:.6g} ({pct(curved_dt_avg, straight_dt_avg):+.2f}%) |\n"
        f"| Heated-wall adjacent dT max, K | {straight_dt_max:.6g} | {curved_dt_max:.6g} | {curved_dt_max - straight_dt_max:.6g} |\n"
        f"| Outlet adjacent dT avg, K | {straight_outlet:.6g} | {curved_outlet:.6g} | {curved_outlet - straight_outlet:.6g} |\n\n"
        f"{interpretation}\n"
    )
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

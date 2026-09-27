#!/usr/bin/env python3
"""Estimate diagnostic heat-transfer metrics from existing summary CSVs.

These are report-planning quantities, not validation-grade correlations. The
current module is an expansion/chamber geometry with passive scalar heat
transport, so the apparent Nusselt number should be described as a diagnostic
metric rather than a developed-pipe benchmark result.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Fluid:
    rho: float
    mu: float
    k: float
    cp: float

    @property
    def prandtl(self) -> float:
        return self.cp * self.mu / self.k


@dataclass(frozen=True)
class Channel:
    width: float
    depth: float
    length: float

    @property
    def hydraulic_diameter(self) -> float:
        return 2.0 * self.width * self.depth / (self.width + self.depth)


FLIBE_LIKE = Fluid(rho=1940.0, mu=0.006, k=1.0, cp=2400.0)
INLET = Channel(width=0.02, depth=0.04, length=0.40)
REYNOLDS = 10_000.0
INLET_SPEED = 1.1598
PRESSURE_DROP_PA = 1774.0
INPUT = Path("results/heat_flux_sensitivity/flibe_thermal_refined_heat_flux_summary.csv")
OUTPUT = Path("results/heat_transfer_metrics/apparent_heat_transfer_metrics.csv")


def dittus_boelter_nusselt(reynolds: float, prandtl: float, heating: bool = True) -> float:
    exponent = 0.4 if heating else 0.3
    return 0.023 * reynolds**0.8 * prandtl**exponent


def module_loss_coefficient(delta_p: float, fluid: Fluid, speed: float) -> float:
    return 2.0 * delta_p / (fluid.rho * speed**2)


def darcy_friction_factor_if_treated_as_straight_duct(
    delta_p: float, fluid: Fluid, speed: float, channel: Channel
) -> float:
    return 2.0 * delta_p * channel.hydraulic_diameter / (
        fluid.rho * speed**2 * channel.length
    )


def main() -> None:
    rows = []
    db_nu = dittus_boelter_nusselt(REYNOLDS, FLIBE_LIKE.prandtl)
    loss_k = module_loss_coefficient(PRESSURE_DROP_PA, FLIBE_LIKE, INLET_SPEED)
    straight_f = darcy_friction_factor_if_treated_as_straight_duct(
        PRESSURE_DROP_PA, FLIBE_LIKE, INLET_SPEED, INLET
    )

    with INPUT.open(newline="", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            q = float(row["heat_flux_kw_m2"]) * 1000.0
            wall_delta_t = float(row["wall_surface_avg_dt"])
            outlet_delta_t = float(row["outlet_avg_dt"])
            h_inlet_ref = q / wall_delta_t
            nu_inlet_ref = h_inlet_ref * INLET.hydraulic_diameter / FLIBE_LIKE.k
            h_outlet_ref = q / max(wall_delta_t - outlet_delta_t, 1e-12)
            nu_outlet_ref = h_outlet_ref * INLET.hydraulic_diameter / FLIBE_LIKE.k

            rows.append(
                {
                    "case": row["case"],
                    "heat_flux_kw_m2": row["heat_flux_kw_m2"],
                    "wall_avg_deltaT_K": f"{wall_delta_t:.6g}",
                    "outlet_avg_deltaT_K": f"{outlet_delta_t:.6g}",
                    "h_app_inlet_ref_W_m2K": f"{h_inlet_ref:.6g}",
                    "Nu_app_inlet_ref": f"{nu_inlet_ref:.6g}",
                    "h_app_outlet_ref_W_m2K": f"{h_outlet_ref:.6g}",
                    "Nu_app_outlet_ref": f"{nu_outlet_ref:.6g}",
                    "Dittus_Boelter_Nu_reference": f"{db_nu:.6g}",
                    "Pr_flibe_like": f"{FLIBE_LIKE.prandtl:.6g}",
                    "Re_reference": f"{REYNOLDS:.6g}",
                    "module_loss_coefficient_K": f"{loss_k:.6g}",
                    "straight_duct_Darcy_f_diagnostic": f"{straight_f:.6g}",
                }
            )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(OUTPUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Estimate first-pass flow and heat-transfer planning values.

This script is intentionally simple and dependency-free. It helps the team
translate target Reynolds numbers into inlet velocities and mass flow rates for
the baseline module.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Fluid:
    name: str
    rho: float  # kg/m^3
    mu: float  # Pa s
    k: float  # W/m K
    cp: float  # J/kg K


@dataclass(frozen=True)
class Channel:
    width: float  # m
    depth: float  # m

    @property
    def area(self) -> float:
        return self.width * self.depth

    @property
    def hydraulic_diameter(self) -> float:
        return 2.0 * self.width * self.depth / (self.width + self.depth)


FLUIDS = [
    Fluid("water_debug_300K", rho=997.0, mu=8.9e-4, k=0.6, cp=4180.0),
    Fluid("air_debug_300K", rho=1.18, mu=1.85e-5, k=0.026, cp=1005.0),
    Fluid("flibe_placeholder_800K", rho=1940.0, mu=6.0e-3, k=1.0, cp=2400.0),
]


def velocity_for_re(reynolds: float, fluid: Fluid, channel: Channel) -> float:
    return reynolds * fluid.mu / (fluid.rho * channel.hydraulic_diameter)


def mass_flow(fluid: Fluid, channel: Channel, velocity: float) -> float:
    return fluid.rho * channel.area * velocity


def prandtl(fluid: Fluid) -> float:
    return fluid.cp * fluid.mu / fluid.k


def main() -> None:
    inlet = Channel(width=0.02, depth=0.04)
    reynolds_targets = [1_000, 10_000, 50_000]

    print("Baseline module inlet estimates")
    print(f"inlet area: {inlet.area:.6f} m^2")
    print(f"hydraulic diameter: {inlet.hydraulic_diameter:.5f} m")
    print()

    for fluid in FLUIDS:
        print(f"fluid: {fluid.name}")
        print(f"  Pr: {prandtl(fluid):.2f}")
        for re in reynolds_targets:
            velocity = velocity_for_re(re, fluid, inlet)
            mdot = mass_flow(fluid, inlet, velocity)
            print(
                f"  Re={re:>6}: U={velocity:>8.4f} m/s, "
                f"m_dot={mdot:>8.4f} kg/s"
            )
        print()


if __name__ == "__main__":
    main()

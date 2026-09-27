# v6 Scalar Temperature Water Re10000 Heat Flux

First passive-temperature transport case.

## Purpose

Use the converged Re=10,000 water velocity field from `v5` and solve passive scalar `T` as a temperature-like field.

This is not buoyant flow or conjugate heat transfer. It is a staged thermal debug case.

## Physics

- frozen flow field from `v5`
- scalar field `T`, interpreted as temperature in K
- thermal diffusivity: `1.44e-7 m^2/s`
- inlet scalar: 300 K
- heated wall: fixed gradient equivalent to about 10 kW/m^2 with water thermal conductivity 0.6 W/m/K
- adiabatic walls: zero gradient
- outlet: inletOutlet

## Visualization

After running, open this in ParaView:

- `VTK/v6_scalar_temperature_water_re10000_heatflux_*.vtm`

# v7 Scalar Temperature Water Re10000 Heat Flux 100 kW

Passive-temperature transport case using the converged Re = 10,000 water velocity field.

## Purpose

Repeat `v6` with a stronger heated-wall signal:

- `v6`: 10 kW/m^2 equivalent
- `v7`: 100 kW/m^2 equivalent

This is still a passive scalar case, not buoyant or conjugate heat transfer.

## Physics

- frozen flow field from `v5`
- scalar field `T`, interpreted as temperature in K
- thermal diffusivity: `1.44e-7 m^2/s`
- inlet scalar: 300 K
- heated wall: fixed gradient equivalent to about 100 kW/m^2 with water thermal conductivity 0.6 W/m/K
- adiabatic walls: zero gradient
- outlet: inletOutlet

## Visualization

After running, open this in ParaView:

- `VTK/v7_scalar_temperature_water_re10000_heatflux100kw_1000.vtm`

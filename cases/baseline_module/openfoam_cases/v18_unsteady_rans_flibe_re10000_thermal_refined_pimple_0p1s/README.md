# v18 Unsteady RANS FLiBe Re10000 Thermal Refined PIMPLE 0.1 s

Short continuation of the successful `v17` transient RANS smoke test.

## Purpose

This case checks whether the `pimpleFoam` setup remains stable beyond the initial `0.02 s` smoke test before committing to a longer roughly flow-through-time run.

## Setup

- solver: `pimpleFoam`
- turbulence model: k-omega SST RANS
- mesh: copied from `v17`
- initial time directory: `0.0201038`, copied from the final saved `v17` smoke-test field
- `startFrom latestTime`
- target end time: `0.1 s`
- adaptive time step: `maxCo 0.5`, `maxDeltaT 2.5e-4 s`
- write interval: `0.01 s`

## Result

The continuation completed successfully.

- final saved time: `0.100208 s`
- final max Courant number from log: about `0.474`
- inlet areaAverage(p): `0.915882 m^2/s^2`
- outlet areaAverage(p): `0 m^2/s^2`
- pressure-drop estimate with rho = 1940 kg/m^3: about `1777 Pa`
- inlet flow rate: `-0.00092784 m^3/s`
- outlet flow rate: `0.00092784 m^3/s`

Generated outputs:

- `results/slices/v18_unsteady_0p1s_midplane_z_velocity.png`
- `results/slices/v18_unsteady_0p1s_midplane_z_velocity.csv`
- `VTK/v18_unsteady_rans_flibe_re10000_thermal_refined_pimple_0p1s_322.vtm`

## Claim Boundary

This is a continuation/stability case. It is still unsteady RANS, not LES, and should be used to decide whether a longer comparison run is worth the runtime.

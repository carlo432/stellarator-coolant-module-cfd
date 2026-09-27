# v4 Isothermal RANS Water Re1000 Clean

This is the canonical first successful OpenFOAM debug case.

## Purpose

Prove the baseline CFD workflow:

1. Gmsh geometry.
2. `gmshToFoam` conversion.
3. wall patch type correction.
4. `checkMesh`.
5. `simpleFoam`.
6. `foamToVTK`.
7. basic pressure and flow post-processing.

## Physics

- water debug fluid
- Re about 1000
- inlet velocity: 0.0335 m/s
- isothermal
- steady RANS
- k-omega SST

## Key Results

- inlet area: 0.0008 m^2
- outlet area: 0.0008 m^2
- inlet flow rate: -2.68e-05 m^3/s
- outlet flow rate: 2.68e-05 m^3/s
- kinematic pressure drop: 0.00104842 m^2/s^2
- water-equivalent pressure drop: about 1.05 Pa

## Visualization

Open this in ParaView:

- `VTK/v4_isothermal_rans_water_re1000_clean_500.vtm`

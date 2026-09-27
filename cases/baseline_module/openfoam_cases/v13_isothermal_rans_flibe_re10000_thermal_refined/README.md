# v13 Isothermal RANS FLiBe Re10000 Thermal Refined

Second refinement of the FLiBe-like Re = 10,000 isothermal RANS baseline.

## Mesh

Geometry:

- `../../geometry/baseline_module_v3_thermal_refined.geo`

Mesh quality:

- cells: 73,344 prism cells
- points: 41,675
- heated-wall faces: 3,216
- max non-orthogonality: 23.5598
- max skewness: 0.500471
- `checkMesh`: OK

## Physics

- solver: `simpleFoam`
- turbulence model: k-omega SST
- FLiBe-like placeholder properties
- inlet velocity: 1.1598 m/s
- Re: about 10,000

## Current Result

Run completed to 1000 iterations.

Final post-processing:

- inlet areaAverage(p): 0.914662 m^2/s^2
- outlet areaAverage(p): 0 m^2/s^2
- pressure drop estimate: about 1774 Pa using rho = 1940 kg/m^3
- inlet volumetric flow rate: -0.00092784 m^3/s
- outlet volumetric flow rate: 0.000927841 m^3/s

VTK:

- `VTK/v13_isothermal_rans_flibe_re10000_thermal_refined_1000.vtm`

Interpretation:

The second refinement still solves cleanly, but pressure drop is lower than `v8`/`v10`. Flow results are not wildly unstable, but the mesh-refinement sequence is not fully grid independent.

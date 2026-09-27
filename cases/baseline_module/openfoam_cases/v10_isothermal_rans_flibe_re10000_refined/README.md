# v10 Isothermal RANS FLiBe Re10000 Refined

Refined-mesh repeat of the FLiBe-like Re = 10,000 isothermal RANS baseline.

## Mesh

Geometry:

- `../../geometry/baseline_module_v2_refined.geo`

Mesh quality:

- cells: 28,976 prism cells
- points: 17,068
- max non-orthogonality: 25.2729
- max skewness: 0.453184
- `checkMesh`: OK

## Physics

- solver: `simpleFoam`
- turbulence model: k-omega SST
- FLiBe-like placeholder properties
- inlet velocity: 1.1598 m/s
- Re: about 10,000

## Current Result

Run completed to 500 iterations.

Final post-processing:

- inlet areaAverage(p): 0.957927 m^2/s^2
- outlet areaAverage(p): 0 m^2/s^2
- pressure drop estimate: about 1858 Pa using rho = 1940 kg/m^3
- inlet volumetric flow rate: -0.00092784 m^3/s
- outlet volumetric flow rate: 0.000927841 m^3/s

VTK:

- `VTK/v10_isothermal_rans_flibe_re10000_refined_500.vtm`

Interpretation:

The refined RANS pressure drop is very close to coarse `v8`, so the flow baseline is not strongly mesh-sensitive at this refinement level.

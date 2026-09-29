# C3 Curved FFHR-Like Depth-Mapped CHT

This case reuses the C2 curved CHT geometry, sets the external
heated-base gradient to zero, and injects the same total power using
mapped fluid and solid cell-zone source distributions from the
FFHR-like OpenMC region-depth profile.

- total power: `427.03204 W`
- fluid power: `225.97539 W`
- solid power: `201.05665 W`
- source stations per region: `24`
- OpenMC profile: `cases/branch_d/figs/c3_openmc_ffhr_be_stack_heating_split_region_depth.csv`
- fluid source power range: `0` to `17.1286 W`
- solid source power range: `0` to `17.4113 W`

Run:

```bash
./Allpre.c2cht
./Allrun.c2cht
```

Caveat: this maps OpenMC stack depth onto streamwise mesh stations.
It exercises profile-to-CHT plumbing and source-placement sensitivity;
it is not a physical wall-normal source map or reactor neutronics.

# Case Notes

This case successfully proved the Gmsh -> OpenFOAM -> simpleFoam pipeline, but it used an incorrect Gmsh physical-surface mapping.

Problem:

- `inlet` was assigned to the lower inlet wall instead of the actual inlet cap.
- This produced an inlet patch area of `0.004 m^2` instead of the expected `0.0008 m^2`.

Use the corrected successor case instead:

- `v3_isothermal_rans_water_re1000_corrected`

# ARC Blanket 2D OpenFOAM Case

Purpose: local OpenFOAM reproduction-style case for Leffler's ARC blanket 2D
velocity-only LES workflow.

This is not a Nek5000 reproduction and not Ferrero CAD.  It is a slide-traced
OpenFOAM approximation using the explicit deck values:

- `L_ref = 0.02 m`
- `U_ref = 2.7 m/s`
- `rho = 2005.781 kg/m3`
- `mu = 0.01047 Pa s`
- `nu = 5.2199118e-06 m2/s`
- `Re = rho*U_ref*L_ref/mu = 10345.00`
- reference velocity for nondimensionalization: `2.7 m/s`
- inlet vector mode: `ferrero-average-normal`
- outlet pressure: `0` kinematic gauge pressure in OpenFOAM's incompressible form

Inlet velocity vectors:

- `inlet_aux100`: `(-0.313 0 0) m/s`
- `inlet_ch_a`: `(1.98 0 0) m/s`
- `inlet_main40`: `(0 1.98 0) m/s`
- `inlet_ch_b`: `(0.88548292 1.7709658 0) m/s`

The full Leffler run reports fields at `t* = 7566.2`.  With
`t* = U_ref*t/L_ref`, that corresponds to about
`56.05 s`.  This local case starts as a short smoke run;
extend `system/controlDict` for production-style averaging.

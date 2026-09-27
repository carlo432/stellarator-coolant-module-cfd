# Ferrero-Dimensioned ARC 2D Geometry Manifest

Source: Ferrero thesis Fig. 2.1 dimensioned sketch, PDF page 18 / body page 10.

This is the first machine-readable scaffold. The outer tank is dimension-anchored;
the inner VV/divertor exclusion is volume-scaled from the prior trace and needs
direct Fig. 2.1/2.6 refinement before promotion.

| Patch | Intended width | Trace length [m] |
|---|---:|---:|
| outlet180 | 0.180 | 0.180000 |
| inlet_aux100 | 0.100 | 0.100000 |
| inlet_ch_a | 0.020 | 0.020000 |
| inlet_main40 | 0.040 | 0.040000 |
| inlet_ch_b | 0.020 | 0.022361 |

## First-Pass Integral Diagnostics

- outer section area: `16.0240 m2`
- inner exclusion area: `1.5404 m2`
- fluid section area: `14.4836 m2`
- outer revolved volume at R=3.5 m: `352.39 m3`
- fluid revolved volume at R=3.5 m: `318.51 m3`
- provisional inner-exclusion scale: `0.720`

Ferrero Table 2.2 checks: tank volume `353.9 m3`, FLiBe volume `319.0 m3`,
section area `16.0 m2`, mean toroidal radius `3.5 m`.

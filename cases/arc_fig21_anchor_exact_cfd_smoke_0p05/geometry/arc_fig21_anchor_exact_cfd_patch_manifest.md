# ARC Fig. 2.1 Anchor-Exact CFD Patch Manifest

Source contours: `references/digitization/ferrero_arc/arc_fig21_anchor_exact_mesh_contours.json`

| Patch | Curve IDs | Provisional basis | Initial U [m/s] |
|---|---:|---|---:|
| `outlet180` | `[75]` | 180 mm split on top tank edge, centered at R=3300 mm | `(0 0 0)` |
| `inlet_aux100` | `[72]` | 100 mm split on upper R=4605 mm tank wall | `(-0.313 0 0)` |
| `inlet_ch_a` | `[245]` | nearest VV contour to legacy upper-channel inlet target | `(1.60179 1.1639 0)` |
| `inlet_main40` | `[121]` | nearest VV contour to legacy combined lower-left inlet target | `(-1.94531 -0.369029 0)` |
| `inlet_ch_b` | `[288]` | nearest VV contour to legacy lower-channel inlet target | `(1.63333 -1.11922 0)` |

The VV inlet cuts are provisional nearest-contour cuts. Run `scripts/audit_arc2d_inlet_bc.py` after meshing.

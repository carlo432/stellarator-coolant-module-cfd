# ARC Fig. 2.1 Anchor-Exact CFD Patch Manifest

Source contours: `references/digitization/ferrero_arc/arc_fig21_anchor_exact_mesh_contours.json`
Patch target file: `cases/branch_d/arc_fig21_patch_targets_thesis_derived.json`

| Patch | Curve IDs | Basis | Initial U [m/s] |
|---|---:|---|---:|
| `outlet180` | `[76]` | 180 mm on top tank edge centered R=3300 (unchanged; thesis: outlet on upper edge) | `(0 0 0)` |
| `inlet_aux100` | `[74]` | CORRECTED: thesis line 669-670: 'tank outlet and supplementary inlet are on the UPPER EDGE of the tank, being the inlet the MOST EXTERNAL channel' -> 100 mm on the Z=+3876 top edge OUTBOARD of the outlet (center R=3800), NOT on the 4605 side wall | `(-0 -0.313 0)` |
| `inlet_ch_a` | `[230]` | circuit 4 discharge: 'outflowing at the end of the leg region in the outer blanket' -> upper lens boundary outboard of the upper-neck attach | `(0.875425 1.77596 0)` |
| `inlet_main40` | `[307]` | merged lower cluster (Leffler's two combined 20mm): circuits 1+3 of thesis sec 2.1.2 discharge at the lens bottom / lower-divertor region -> VV lens lower boundary just inboard of the lower-neck attach | `(1.62729 1.12798 0)` |
| `inlet_ch_b` | `[330]` | circuit 2 discharge: 'outflowing after the other lower divertor leg' -> outboard side of the LOWER FOOT circle | `(1.8943 -0.576223 0)` |

The VV inlet cuts are provisional nearest-contour cuts. Change the patch target JSON, regenerate, and rerun `scripts/audit_arc2d_inlet_bc.py` after meshing.

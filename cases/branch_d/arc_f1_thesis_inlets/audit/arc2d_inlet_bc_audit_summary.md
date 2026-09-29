# ARC 2D Inlet Boundary-Condition Audit

Case: `cases/branch_d/arc_f1_thesis_inlets`

This compares the fixed inlet velocity vectors in `0/U` against OpenFOAM patch face-area normals.
Negative normal speed means inflow into the domain; positive means outflow; near-zero means mostly tangential.

| Patch | Current U [m/s] | Outward normal xy | Normal speed outward [m/s] | Angle to inward normal | Recommended normal-inlet U [m/s] | Class |
|---|---:|---:|---:|---:|---:|---|
| inlet_aux100 | `(-0.000, -0.313)` | `(0.000, 1.000)` | `-0.313` | `0.0 deg` | `(-0.000, -2.700)` | inflow |
| inlet_ch_a | `(0.875, 1.776)` | `(-0.442, -0.897)` | `-1.980` | `0.0 deg` | `(1.194, 2.422)` | inflow |
| inlet_ch_b | `(1.894, -0.576)` | `(-0.957, 0.291)` | `-1.980` | `0.0 deg` | `(2.583, -0.786)` | inflow |
| inlet_main40 | `(1.627, 1.128)` | `(-0.822, -0.570)` | `-1.980` | `0.0 deg` | `(2.219, 1.538)` | inflow |

CSV: `cases/branch_d/arc_f1_thesis_inlets/audit/arc2d_inlet_bc_audit.csv`

Interpretation: use this as a boundary-condition diagnostic only. Changing inlet vectors would require a new run; it would not retroactively change the completed endpoint fields.

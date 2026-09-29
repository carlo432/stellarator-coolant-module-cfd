# C3 Curved Material-Split Volumetric CHT

This case reuses the C2 curved CHT geometry, disables the external
heated-base gradient, and injects same-total-power volumetric sources
into the FLiBe-like fluid and solid heater regions.

- target id: `c3_ms_total_blanket_55_45`
- target equivalent heat flux: `100000 W/m^2`
- heated-base reference area: `0.0042703204 m^2`
- total power: `427.03204 W`
- fluid power: `234.86762 W`
- solid power: `192.16442 W`
- target source table: `/root/revengence/results/c3_material_split_targets.csv`

Run:

```bash
./Allpre.c2cht
./Allrun.c2cht
```

This is a literature-guided source placement case, not a validated
stellarator neutronics prediction.

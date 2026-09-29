#!/usr/bin/env python3
from pathlib import Path

solid_cells = 6144
fluid_cells = 38400
total_cells = solid_cells + fluid_cells
out = Path("constant/polyMesh/cellZones")

def labels(start, stop):
    return "\n".join(str(i) for i in range(start, stop))

parts = [
    "/*--------------------------------*- C++ -*----------------------------------*/",
    "FoamFile",
    "{",
    "    version 2.0;",
    "    format ascii;",
    "    class regIOobject;",
    "    location \"constant/polyMesh\";",
    "    object cellZones;",
    "    meta { names ( bottomWater heater ); }",
    "}",
    "2",
    "(",
    "bottomWater",
    "{",
    "    type cellZone;",
    "    cellLabels List<label>",
    str(fluid_cells),
    "(",
    labels(solid_cells, total_cells),
    ");",
    "}",
    "heater",
    "{",
    "    type cellZone;",
    "    cellLabels List<label>",
    str(solid_cells),
    "(",
    labels(0, solid_cells),
    ");",
    "}",
    ")",
    "",
]
out.write_text("\n".join(parts))
print(f"wrote {out} with bottomWater={fluid_cells} heater={solid_cells}")

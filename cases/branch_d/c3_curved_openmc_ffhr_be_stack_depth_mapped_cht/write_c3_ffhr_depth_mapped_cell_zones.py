#!/usr/bin/env python3
from pathlib import Path

fluid_cells = 38400
solid_cells = 6144
fluid_zones = [{'name': 'c3FFHRFluidSource_00', 'start': 0, 'stop': 1600}, {'name': 'c3FFHRFluidSource_01', 'start': 1600, 'stop': 3200}, {'name': 'c3FFHRFluidSource_02', 'start': 3200, 'stop': 4800}, {'name': 'c3FFHRFluidSource_03', 'start': 4800, 'stop': 6400}, {'name': 'c3FFHRFluidSource_04', 'start': 6400, 'stop': 8000}, {'name': 'c3FFHRFluidSource_05', 'start': 8000, 'stop': 9600}, {'name': 'c3FFHRFluidSource_06', 'start': 9600, 'stop': 11200}, {'name': 'c3FFHRFluidSource_07', 'start': 11200, 'stop': 12800}, {'name': 'c3FFHRFluidSource_08', 'start': 12800, 'stop': 14400}, {'name': 'c3FFHRFluidSource_09', 'start': 14400, 'stop': 16000}, {'name': 'c3FFHRFluidSource_10', 'start': 16000, 'stop': 17600}, {'name': 'c3FFHRFluidSource_11', 'start': 17600, 'stop': 19200}, {'name': 'c3FFHRFluidSource_12', 'start': 19200, 'stop': 20800}, {'name': 'c3FFHRFluidSource_13', 'start': 20800, 'stop': 22400}, {'name': 'c3FFHRFluidSource_14', 'start': 22400, 'stop': 24000}, {'name': 'c3FFHRFluidSource_15', 'start': 24000, 'stop': 25600}, {'name': 'c3FFHRFluidSource_16', 'start': 25600, 'stop': 27200}, {'name': 'c3FFHRFluidSource_17', 'start': 27200, 'stop': 28800}, {'name': 'c3FFHRFluidSource_18', 'start': 28800, 'stop': 30400}, {'name': 'c3FFHRFluidSource_19', 'start': 30400, 'stop': 32000}, {'name': 'c3FFHRFluidSource_20', 'start': 32000, 'stop': 33600}, {'name': 'c3FFHRFluidSource_21', 'start': 33600, 'stop': 35200}, {'name': 'c3FFHRFluidSource_22', 'start': 35200, 'stop': 36800}, {'name': 'c3FFHRFluidSource_23', 'start': 36800, 'stop': 38400}]
solid_zones = [{'name': 'c3FFHRSolidSource_00', 'start': 0, 'stop': 256}, {'name': 'c3FFHRSolidSource_01', 'start': 256, 'stop': 512}, {'name': 'c3FFHRSolidSource_02', 'start': 512, 'stop': 768}, {'name': 'c3FFHRSolidSource_03', 'start': 768, 'stop': 1024}, {'name': 'c3FFHRSolidSource_04', 'start': 1024, 'stop': 1280}, {'name': 'c3FFHRSolidSource_05', 'start': 1280, 'stop': 1536}, {'name': 'c3FFHRSolidSource_06', 'start': 1536, 'stop': 1792}, {'name': 'c3FFHRSolidSource_07', 'start': 1792, 'stop': 2048}, {'name': 'c3FFHRSolidSource_08', 'start': 2048, 'stop': 2304}, {'name': 'c3FFHRSolidSource_09', 'start': 2304, 'stop': 2560}, {'name': 'c3FFHRSolidSource_10', 'start': 2560, 'stop': 2816}, {'name': 'c3FFHRSolidSource_11', 'start': 2816, 'stop': 3072}, {'name': 'c3FFHRSolidSource_12', 'start': 3072, 'stop': 3328}, {'name': 'c3FFHRSolidSource_13', 'start': 3328, 'stop': 3584}, {'name': 'c3FFHRSolidSource_14', 'start': 3584, 'stop': 3840}, {'name': 'c3FFHRSolidSource_15', 'start': 3840, 'stop': 4096}, {'name': 'c3FFHRSolidSource_16', 'start': 4096, 'stop': 4352}, {'name': 'c3FFHRSolidSource_17', 'start': 4352, 'stop': 4608}, {'name': 'c3FFHRSolidSource_18', 'start': 4608, 'stop': 4864}, {'name': 'c3FFHRSolidSource_19', 'start': 4864, 'stop': 5120}, {'name': 'c3FFHRSolidSource_20', 'start': 5120, 'stop': 5376}, {'name': 'c3FFHRSolidSource_21', 'start': 5376, 'stop': 5632}, {'name': 'c3FFHRSolidSource_22', 'start': 5632, 'stop': 5888}, {'name': 'c3FFHRSolidSource_23', 'start': 5888, 'stop': 6144}]

def labels(start, stop):
    return "\n".join(str(i) for i in range(start, stop))

def zone_block(name, start, stop):
    return [
        name,
        "{",
        "    type cellZone;",
        "    cellLabels List<label>",
        str(stop - start),
        "(",
        labels(start, stop),
        ");",
        "}",
    ]

def write_cell_zones(region, region_total_cells, zones):
    out = Path(f"constant/{region}/polyMesh/cellZones")
    names = [region] + [zone["name"] for zone in zones]
    parts = [
        "/*--------------------------------*- C++ -*----------------------------------*/",
        "FoamFile",
        "{",
        "    version 2.0;",
        "    format ascii;",
        "    class regIOobject;",
        f"    location \"constant/{region}/polyMesh\";",
        "    object cellZones;",
        "    meta { names ( " + " ".join(names) + " ); }",
        "}",
        str(len(names)),
        "(",
    ]
    parts.extend(zone_block(region, 0, region_total_cells))
    for zone in zones:
        parts.extend(zone_block(zone["name"], zone["start"], zone["stop"]))
    parts.extend([")", ""])
    out.write_text("\n".join(parts))
    print(f"wrote {out} with {len(zones)} mapped source zones")

write_cell_zones("bottomWater", fluid_cells, fluid_zones)
write_cell_zones("heater", solid_cells, solid_zones)

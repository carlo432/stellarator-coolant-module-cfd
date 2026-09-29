#!/usr/bin/env python3
"""Generate a small curved C2 conjugate-heat-transfer smoke case.

This is the CHT rung after the C2 passive-scalar comparisons: a curved fluid
channel sharing its heated-side interface with a thin solid slab.  It is still a
local scaffold, not a reactor blanket CAD model.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import stat
from pathlib import Path

from gen_curved_stellarator_slice_c2 import centerline, foam_header, fmt, radial_normal, vadd, vscale, write


def find_overall_ratio_for_first_cell(first_cell: float, length: float, n_cells: int) -> float:
    lo, hi = 1.000001, 100.0
    for _ in range(200):
        k = 0.5 * (lo + hi)
        c = length * (k - 1.0) / (k**n_cells - 1.0)
        if c > first_cell:
            lo = k
        else:
            hi = k
    return hi ** (n_cells - 1)


def station_vertices(
    t: float,
    radius: float,
    arc: float,
    y_amp: float,
    radial_width: float,
    height: float,
    solid_thickness: float,
) -> list[tuple[float, float, float]]:
    c = centerline(t, radius, arc, y_amp)
    n = radial_normal(t, arc)
    ey = (0.0, 1.0, 0.0)

    # A/D are the fluid-solid interface on the plasma-facing heated side.
    a = vadd(vadd(c, vscale(n, -0.5 * radial_width)), vscale(ey, -0.5 * height))
    b = vadd(vadd(c, vscale(n, 0.5 * radial_width)), vscale(ey, -0.5 * height))
    cpt = vadd(vadd(c, vscale(n, 0.5 * radial_width)), vscale(ey, 0.5 * height))
    d = vadd(vadd(c, vscale(n, -0.5 * radial_width)), vscale(ey, 0.5 * height))
    e = vadd(a, vscale(n, -solid_thickness))
    f = vadd(d, vscale(n, -solid_thickness))
    return [a, b, cpt, d, e, f]


def replace_control_dict(path: Path, end_time: int, write_interval: int) -> None:
    text = (path / "system/controlDict").read_text()
    text = text.replace("endTime         5000;", f"endTime         {end_time};")
    text = text.replace("writeInterval   5000;", f"writeInterval   {write_interval};")
    # Avoid writing large VTK trees during this preflight; render/postprocess are separate.
    if "functions" in text:
        text = text.split("functions", 1)[0].rstrip() + "\n\n// functions disabled for C2 CHT smoke preflight\n"
    (path / "system/controlDict").write_text(text)


def write_cell_zone_helper(case: Path, solid_cells: int, fluid_cells: int) -> None:
    helper = case / "write_cell_zones_c2cht.py"
    write(
        helper,
        f"""#!/usr/bin/env python3
from pathlib import Path

solid_cells = {solid_cells}
fluid_cells = {fluid_cells}
total_cells = solid_cells + fluid_cells
out = Path("constant/polyMesh/cellZones")

def labels(start, stop):
    return "\\n".join(str(i) for i in range(start, stop))

parts = [
    "/*--------------------------------*- C++ -*----------------------------------*/",
    "FoamFile",
    "{{",
    "    version 2.0;",
    "    format ascii;",
    "    class regIOobject;",
    "    location \\"constant/polyMesh\\";",
    "    object cellZones;",
    "    meta {{ names ( bottomWater heater ); }}",
    "}}",
    "2",
    "(",
    "bottomWater",
    "{{",
    "    type cellZone;",
    "    cellLabels List<label>",
    str(fluid_cells),
    "(",
    labels(solid_cells, total_cells),
    ");",
    "}}",
    "heater",
    "{{",
    "    type cellZone;",
    "    cellLabels List<label>",
    str(solid_cells),
    "(",
    labels(0, solid_cells),
    ");",
    "}}",
    ")",
    "",
]
out.write_text("\\n".join(parts))
print(f"wrote {{out}} with bottomWater={{fluid_cells}} heater={{solid_cells}}")
""",
    )
    helper.chmod(helper.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="c2_curved_stellarator_slice_cht_preflight")
    ap.add_argument("--segments", type=int, default=24)
    ap.add_argument("--stream-cells", type=int, default=4)
    ap.add_argument("--radial-cells", type=int, default=10)
    ap.add_argument("--height-cells", type=int, default=8)
    ap.add_argument("--solid-cells", type=int, default=2)
    ap.add_argument(
        "--wall-adjacent-cell",
        type=float,
        default=None,
        help=(
            "Target fluid-cell width adjacent to the fluid-solid interface. "
            "If set, radial simpleGrading is reversed so the final radial cell is this size."
        ),
    )
    ap.add_argument("--radius", type=float, default=0.18)
    ap.add_argument("--arc-deg", type=float, default=95.0)
    ap.add_argument("--y-amp", type=float, default=0.010)
    ap.add_argument("--radial-width", type=float, default=0.030)
    ap.add_argument("--height", type=float, default=0.016)
    ap.add_argument("--solid-thickness", type=float, default=0.004)
    ap.add_argument("--u-in", type=float, default=1.1598)
    ap.add_argument("--heat-flux", type=float, default=100000.0)
    ap.add_argument("--solid-kappa", type=float, default=30.0)
    ap.add_argument("--end-time", type=int, default=400)
    ap.add_argument("--write-interval", type=int, default=400)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    work = Path(__file__).resolve().parent
    template = work / "cht_wallres"
    case = work / args.case
    if case.exists():
        if not args.overwrite:
            raise FileExistsError(f"{case} exists; use --overwrite to replace it")
        shutil.rmtree(case)
    case.mkdir(parents=True)

    shutil.copytree(template / "constant/bottomWater", case / "constant/bottomWater")
    shutil.copytree(template / "constant/heater", case / "constant/heater")
    shutil.rmtree(case / "constant/bottomWater/polyMesh", ignore_errors=True)
    shutil.rmtree(case / "constant/heater/polyMesh", ignore_errors=True)
    shutil.copy2(template / "constant/g", case / "constant/g")
    shutil.copy2(template / "constant/regionProperties", case / "constant/regionProperties")
    shutil.copytree(template / "system/bottomWater", case / "system/bottomWater")
    shutil.copytree(template / "system/heater", case / "system/heater")
    shutil.copy2(template / "system/fvSchemes", case / "system/fvSchemes")
    shutil.copy2(template / "system/fvSolution", case / "system/fvSolution")
    shutil.copy2(template / "system/controlDict", case / "system/controlDict")
    replace_control_dict(case, args.end_time, args.write_interval)

    arc = math.radians(args.arc_deg)
    radial_grading = 1.0
    if args.wall_adjacent_cell is not None:
        radial_grading = 1.0 / find_overall_ratio_for_first_cell(
            args.wall_adjacent_cell, args.radial_width, args.radial_cells
        )
    stations = args.segments + 1
    verts: list[tuple[float, float, float]] = []
    for i in range(stations):
        verts.extend(
            station_vertices(
                i / args.segments,
                args.radius,
                arc,
                args.y_amp,
                args.radial_width,
                args.height,
                args.solid_thickness,
            )
        )

    def idx(station: int, corner: int) -> int:
        return 6 * station + corner

    solid_blocks: list[tuple[int, ...]] = []
    fluid_blocks: list[tuple[int, ...]] = []
    heated_base_faces: list[tuple[int, ...]] = []
    solid_side_faces: list[tuple[int, ...]] = []
    back_faces: list[tuple[int, ...]] = []
    side_faces: list[tuple[int, ...]] = []

    for i in range(args.segments):
        # Solid: interface A/D to external E/F.
        a_i, a_j = idx(i, 0), idx(i + 1, 0)
        d_i, d_j = idx(i, 3), idx(i + 1, 3)
        e_i, e_j = idx(i, 4), idx(i + 1, 4)
        f_i, f_j = idx(i, 5), idx(i + 1, 5)
        solid_blocks.append((a_i, a_j, d_j, d_i, e_i, e_j, f_j, f_i))
        heated_base_faces.append((e_i, f_i, f_j, e_j))
        solid_side_faces.append((a_i, e_i, e_j, a_j))
        solid_side_faces.append((d_i, d_j, f_j, f_i))

        # Fluid: back B/C to heated-side interface A/D.
        b_i, b_j = idx(i, 1), idx(i + 1, 1)
        c_i, c_j = idx(i, 2), idx(i + 1, 2)
        fluid_blocks.append((b_i, b_j, c_j, c_i, a_i, a_j, d_j, d_i))
        back_faces.append((b_i, b_j, c_j, c_i))
        side_faces.append((b_i, a_i, a_j, b_j))
        side_faces.append((c_i, c_j, d_j, d_i))

    solid_side_faces.append((idx(0, 0), idx(0, 3), idx(0, 5), idx(0, 4)))
    solid_side_faces.append((idx(args.segments, 0), idx(args.segments, 4), idx(args.segments, 5), idx(args.segments, 3)))
    inlet_face = (idx(0, 1), idx(0, 2), idx(0, 3), idx(0, 0))
    outlet_face = (
        idx(args.segments, 1),
        idx(args.segments, 0),
        idx(args.segments, 3),
        idx(args.segments, 2),
    )

    vertices_txt = "\n".join(f"    {fmt(v)}" for v in verts)
    solid_blocks_txt = "\n".join(
        "    hex (%s) (%d %d %d) simpleGrading (1 1 1)"
        % (" ".join(str(x) for x in block), args.stream_cells, args.height_cells, args.solid_cells)
        for block in solid_blocks
    )
    fluid_blocks_txt = "\n".join(
        "    hex (%s) (%d %d %d) simpleGrading (1 1 %.8g)"
        % (
            " ".join(str(x) for x in block),
            args.stream_cells,
            args.height_cells,
            args.radial_cells,
            radial_grading,
        )
        for block in fluid_blocks
    )

    def faces_txt(faces: list[tuple[int, ...]]) -> str:
        return "\n".join("        (%s)" % " ".join(str(x) for x in face) for face in faces)

    write(
        case / "system/blockMeshDict",
        foam_header("dictionary", "blockMeshDict")
        + f"""scale 1;
vertices
(
{vertices_txt}
);
blocks
(
{solid_blocks_txt}
{fluid_blocks_txt}
);
edges ();
boundary
(
    inlet
    {{
        type patch;
        faces (({' '.join(str(x) for x in inlet_face)}));
    }}
    outlet
    {{
        type patch;
        faces (({' '.join(str(x) for x in outlet_face)}));
    }}
    blanket_back_wall
    {{
        type wall;
        faces
        (
{faces_txt(back_faces)}
        );
    }}
    side_walls
    {{
        type wall;
        faces
        (
{faces_txt(side_faces)}
        );
    }}
    heatedBase
    {{
        type wall;
        faces
        (
{faces_txt(heated_base_faces)}
        );
    }}
    solidSides
    {{
        type wall;
        faces
        (
{faces_txt(solid_side_faces)}
        );
    }}
);
mergePatchPairs ();
""",
    )

    solid_cells = args.segments * args.stream_cells * args.height_cells * args.solid_cells
    fluid_cells = args.segments * args.stream_cells * args.height_cells * args.radial_cells
    write_cell_zone_helper(case, solid_cells, fluid_cells)

    wall_regex = "(blanket_back_wall|side_walls|bottomWater_to_heater)"
    write(
        case / "0/bottomWater/U",
        foam_header("volVectorField", "U")
        + f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform ({args.u_in:.12g} 0 0);
boundaryField
{{
    inlet {{ type fixedValue; value uniform ({args.u_in:.12g} 0 0); }}
    outlet {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    "{wall_regex}" {{ type noSlip; }}
}}
""",
    )
    write(
        case / "0/bottomWater/T",
        foam_header("volScalarField", "T")
        + """dimensions [0 0 0 1 0 0 0];
internalField uniform 800;
boundaryField
{
    inlet { type fixedValue; value uniform 800; }
    outlet { type inletOutlet; inletValue uniform 800; value uniform 800; }
    "(blanket_back_wall|side_walls)" { type zeroGradient; }
    bottomWater_to_heater { type compressible::turbulentTemperatureRadCoupledMixed; Tnbr T; kappaMethod fluidThermo; value uniform 800; }
}
""",
    )
    write(
        case / "0/bottomWater/p_rgh",
        foam_header("volScalarField", "p_rgh")
        + f"""dimensions [1 -1 -2 0 0 0 0];
internalField uniform 1e5;
boundaryField
{{
    inlet {{ type fixedFluxPressure; value uniform 1e5; }}
    outlet {{ type fixedValue; value uniform 1e5; }}
    "{wall_regex}" {{ type fixedFluxPressure; value uniform 1e5; }}
}}
""",
    )
    write(
        case / "0/bottomWater/p",
        foam_header("volScalarField", "p")
        + """dimensions [1 -1 -2 0 0 0 0];
internalField uniform 1e5;
boundaryField { ".*" { type calculated; value uniform 1e5; } }
""",
    )
    write(
        case / "0/bottomWater/alphat",
        foam_header("volScalarField", "alphat")
        + f"""dimensions [1 -1 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    "(inlet|outlet)" {{ type calculated; value uniform 0; }}
    "{wall_regex}" {{ type compressible::alphatWallFunction; Prt 0.85; value uniform 0; }}
}}
""",
    )
    write(
        case / "0/bottomWater/k",
        foam_header("volScalarField", "k")
        + f"""dimensions [0 2 -2 0 0 0 0];
internalField uniform 0.00505;
boundaryField
{{
    inlet {{ type fixedValue; value uniform 0.00505; }}
    outlet {{ type inletOutlet; inletValue uniform 0.00505; value uniform 0.00505; }}
    "{wall_regex}" {{ type kqRWallFunction; value uniform 0.00505; }}
}}
""",
    )
    write(
        case / "0/bottomWater/omega",
        foam_header("volScalarField", "omega")
        + f"""dimensions [0 0 -1 0 0 0 0];
internalField uniform 69.3;
boundaryField
{{
    inlet {{ type fixedValue; value uniform 69.3; }}
    outlet {{ type inletOutlet; inletValue uniform 69.3; value uniform 69.3; }}
    "{wall_regex}" {{ type omegaWallFunction; value uniform 69.3; }}
}}
""",
    )
    write(
        case / "0/bottomWater/nut",
        foam_header("volScalarField", "nut")
        + f"""dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    "(inlet|outlet)" {{ type calculated; value uniform 0; }}
    "{wall_regex}" {{ type nutLowReWallFunction; value uniform 0; }}
}}
""",
    )

    heater_gradient = args.heat_flux / args.solid_kappa
    write(
        case / "0/heater/T",
        foam_header("volScalarField", "T")
        + f"""dimensions [0 0 0 1 0 0 0];
internalField uniform 800;
boundaryField
{{
    heatedBase {{ type fixedGradient; gradient uniform {heater_gradient:.12g}; }}
    solidSides {{ type zeroGradient; }}
    heater_to_bottomWater {{ type compressible::turbulentTemperatureRadCoupledMixed; Tnbr T; kappaMethod solidThermo; value uniform 800; }}
}}
""",
    )
    write(
        case / "0/heater/p",
        foam_header("volScalarField", "p")
        + """dimensions [1 -1 -2 0 0 0 0];
internalField uniform 1e5;
boundaryField { ".*" { type calculated; value uniform 1e5; } }
""",
    )

    allpre = case / "Allpre.c2cht"
    write(
        allpre,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
blockMesh > log.blockMesh 2>&1
./write_cell_zones_c2cht.py > log.cellZones 2>&1
splitMeshRegions -cellZones -overwrite > log.split 2>&1
checkMesh -allRegions > log.checkMesh 2>&1
""",
    )
    allpre.chmod(allpre.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    allrun = case / "Allrun.c2cht"
    write(
        allrun,
        """#!/usr/bin/env bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -euo pipefail
cd "$(dirname "$0")"
chtMultiRegionSimpleFoam > log.cht 2>&1
""",
    )
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    write(
        case / "geometry_c2cht.json",
        json.dumps(
            {
                "vertices": verts,
                "solid_blocks": solid_blocks,
                "fluid_blocks": fluid_blocks,
                "patches": {
                    "inlet": [inlet_face],
                    "outlet": [outlet_face],
                    "blanket_back_wall": back_faces,
                    "side_walls": side_faces,
                    "heatedBase": heated_base_faces,
                    "solidSides": solid_side_faces,
                },
                "centerline": [centerline(i / args.segments, args.radius, arc, args.y_amp) for i in range(stations)],
                "solid_cells": solid_cells,
                "fluid_cells": fluid_cells,
                "total_cells": solid_cells + fluid_cells,
                "radial_grading": radial_grading,
                "parameters": vars(args),
            },
            indent=2,
        ),
    )
    write(
        case / "README.c2cht.md",
        f"""# C2 Curved CHT Preflight

Generated by `gen_c2_curved_cht.py`.

- fluid cells: `{fluid_cells}`
- solid cells: `{solid_cells}`
- total cells before region split: `{solid_cells + fluid_cells}`
- solid thickness: `{args.solid_thickness} m`
- fluid radial grading: `{radial_grading:.8g}`
- requested wall-adjacent fluid cell: `{args.wall_adjacent_cell}`
- external base heat flux: `{args.heat_flux} W/m^2`
- heater gradient: `{heater_gradient:.8g} K/m`

Run:

```bash
./Allpre.c2cht
./Allrun.c2cht
```

This is a CHT smoke/preflight scaffold on a representative curved C2 slice. It
is not a converged or validated stellarator blanket calculation.
""",
    )
    print(
        f"wrote {case}/ fluid_cells={fluid_cells} solid_cells={solid_cells} "
        f"total={solid_cells + fluid_cells}"
    )


if __name__ == "__main__":
    main()

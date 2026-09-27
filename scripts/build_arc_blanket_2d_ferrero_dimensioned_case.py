#!/usr/bin/env python3
"""Build a first Ferrero-dimensioned ARC 2D OpenFOAM case scaffold.

This is a separate geometry path from the old Leffler slide trace.  The outer
tank outline is anchored to the explicit millimeter dimensions in Ferrero
thesis Fig. 2.1.  The inner VV/divertor exclusion is still a first-pass
polyline reconstruction and must be verified against Ferrero Table 2.2 before
being promoted as source-level geometry recovery.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from build_arc_blanket_2d_case import THICKNESS
from build_arc_blanket_2d_case import dist, write, write_openfoam_case


R0_MM = 3300.0
DIM_MM = {
    "r_inner_step": 1521.0,
    "r_left_outer": 1934.0,
    "r_left_inner": 2023.0,
    "r_midplane": 3300.0,
    "r_top_width": 3708.0,
    "r_right_inner": 4605.0,
    "r_right_outer": 4731.0,
    "r_far_right": 5803.0,
    "z_left_low": 2338.0,
    "z_left_mid": 2807.0,
    "z_left_high": 3194.0,
    "z_right_step": 1859.0,
    "z_full": 3876.0,
}
INNER_EXCLUSION_SCALE = 0.72


def local_r(mm: float) -> float:
    return (mm - R0_MM) / 1000.0


def half(mm: float) -> float:
    return mm / 2000.0


def polygon_area(points: list[tuple[float, float, float]]) -> float:
    xy = [(p[0], p[1]) for p in points]
    acc = 0.0
    for (x0, y0), (x1, y1) in zip(xy, xy[1:] + xy[:1]):
        acc += x0 * y1 - x1 * y0
    return abs(acc) * 0.5


def write_dimensioned_geo(case: Path, mesh_size: float, outer_mode: str) -> None:
    geom = case / "geometry"
    geom.mkdir(parents=True, exist_ok=True)

    x1521 = local_r(DIM_MM["r_inner_step"])
    x1934 = local_r(DIM_MM["r_left_outer"])
    x2023 = local_r(DIM_MM["r_left_inner"])
    x3300 = local_r(DIM_MM["r_midplane"])
    x3708 = local_r(DIM_MM["r_top_width"])
    x4605 = local_r(DIM_MM["r_right_inner"])
    x4731 = local_r(DIM_MM["r_right_outer"])
    x5803 = local_r(DIM_MM["r_far_right"])

    z_full = half(DIM_MM["z_full"])
    z_left_high = half(DIM_MM["z_left_high"])
    z_left_mid = half(DIM_MM["z_left_mid"])
    z_left_low = half(DIM_MM["z_left_low"])
    z_right = half(DIM_MM["z_right_step"])

    # Two modes are intentionally kept.  "area-matched" reproduces the Table
    # 2.2 integral scale.  "figure-shaped" follows the visible Fig. 2.1 stepped
    # silhouette more closely, exposing the remaining digitization ambiguity.
    if outer_mode == "area-matched":
        outer = [
            (x1521, z_full, mesh_size),
            (x3300 - 0.090, z_full, mesh_size * 0.35),
            (x3300 + 0.090, z_full, mesh_size * 0.35),  # 180 mm outlet
            (x3708, z_full, mesh_size),
            (x4731, z_full, mesh_size * 0.5),
            (x4731, z_full - 0.100, mesh_size * 0.35),  # 100 mm auxiliary inlet
            (x4731, z_left_mid, mesh_size),
            (x5803, z_full, mesh_size),
            (x5803, -z_full, mesh_size),
            (x4731, -z_left_mid, mesh_size),
            (x4731, -z_full, mesh_size),
            (x3708, -z_full, mesh_size),
            (x1521, -z_full, mesh_size),
        ]
    elif outer_mode == "figure-shaped":
        outer = [
            (x2023, z_full, mesh_size),
            (x3300 - 0.090, z_full, mesh_size * 0.35),
            (x3300 + 0.090, z_full, mesh_size * 0.35),  # 180 mm outlet
            (x3708, z_full, mesh_size),
            (x4731, z_full, mesh_size * 0.5),
            (x4731, z_full - 0.100, mesh_size * 0.35),  # 100 mm auxiliary inlet
            (x4731, z_left_mid, mesh_size),
            (x5803, z_right, mesh_size),
            (x5803, -z_right, mesh_size),
            (x4731, -z_left_mid, mesh_size),
            (x4731, -z_full, mesh_size),
            (x3708, -z_full, mesh_size),
            (x2023, -z_full, mesh_size),
            (x1521, -z_left_high, mesh_size),
            (x1521, -z_right, mesh_size),
            (x1934, -z_right, mesh_size),
            (x1934, z_right, mesh_size),
            (x1521, z_right, mesh_size),
            (x1521, z_left_high, mesh_size),
        ]
    else:
        raise ValueError(f"unknown outer mode {outer_mode}")

    # First-pass VV/divertor exclusion, rescaled to the Ferrero coordinate
    # frame.  The short inlet segments preserve Leffler's explicit 20/40/20 mm
    # channel openings for controlled F1 comparisons.
    inner_raw = [
        (0.30, 1.72, mesh_size),
        (0.52, 1.60, mesh_size),
        (0.56, 1.38, mesh_size),
        (0.43, 1.18, mesh_size),
        (0.22, 1.12, mesh_size * 0.35),
        (0.22, 1.10, mesh_size * 0.35),  # 20 mm inlet_ch_a
        (0.05, 0.80, mesh_size),
        (0.40, 0.55, mesh_size),
        (0.68, 0.18, mesh_size),
        (0.78, -0.25, mesh_size),
        (0.68, -0.65, mesh_size),
        (0.43, -0.95, mesh_size),
        (0.25, -1.16, mesh_size * 0.35),
        (0.27, -1.17, mesh_size * 0.35),  # 20 mm inlet_ch_b
        (0.45, -1.20, mesh_size),
        (0.65, -1.36, mesh_size),
        (0.66, -1.60, mesh_size),
        (0.50, -1.78, mesh_size),
        (0.25, -1.82, mesh_size),
        (0.06, -1.66, mesh_size),
        (0.02, -1.42, mesh_size),
        (0.12, -1.25, mesh_size),
        (-0.39, -1.05, mesh_size * 0.5),
        (-0.35, -1.05, mesh_size * 0.5),  # 40 mm inlet_main40 group
        (-0.65, -0.65, mesh_size),
        (-0.78, -0.20, mesh_size),
        (-0.75, 0.35, mesh_size),
        (-0.55, 0.75, mesh_size),
        (-0.25, 0.98, mesh_size),
        (0.05, 1.05, mesh_size),
    ]
    inner = [(x * INNER_EXCLUSION_SCALE, y * INNER_EXCLUSION_SCALE, h) for x, y, h in inner_raw]
    inner[5] = (inner[4][0], inner[4][1] - 0.020, inner[5][2])
    inner[13] = (inner[12][0] + 0.020, inner[12][1] - 0.010, inner[13][2])
    inner[23] = (inner[22][0] + 0.040, inner[22][1], inner[23][2])

    all_points = outer + inner
    lines = []
    for i, (x, y, h) in enumerate(all_points, 1):
        lines.append(f"Point({i}) = {{{x:.8g}, {y:.8g}, 0, {h:.8g}}};")

    outer_lines = []
    for i in range(len(outer)):
        lid = i + 1
        lines.append(f"Line({lid}) = {{{i + 1}, {(i + 1) % len(outer) + 1}}};")
        outer_lines.append(lid)

    offset = len(outer)
    inner_lines = []
    for i in range(len(inner)):
        lid = len(outer) + i + 1
        a = offset + i + 1
        b = offset + ((i + 1) % len(inner)) + 1
        lines.append(f"Line({lid}) = {{{a}, {b}}};")
        inner_lines.append(lid)

    n_outer = len(outer_lines)
    outlet180 = [2]
    inlet_aux100 = [5]
    inlet_ch_a = [n_outer + 5]
    inlet_ch_b = [n_outer + 13]
    inlet_main40 = [n_outer + 23]
    inlet_or_outlet = set(outlet180 + inlet_aux100 + inlet_ch_a + inlet_main40 + inlet_ch_b)
    all_curve_indices = list(range(1, n_outer + len(inner_lines) + 1))
    wall_indices = [idx for idx in all_curve_indices if idx not in inlet_or_outlet]

    def surf_refs(curve_indices: list[int]) -> str:
        return ", ".join(f"ext[{2 + (idx - 1)}]" for idx in curve_indices)

    geo = "\n".join(
        [
            'SetFactory("Built-in");',
            f"Mesh.CharacteristicLengthMin = {mesh_size * 0.3:.8g};",
            f"Mesh.CharacteristicLengthMax = {mesh_size:.8g};",
            "Mesh.Algorithm = 6;",
            "Mesh.RecombineAll = 0;",
            *lines,
            f"Curve Loop(1) = {{{', '.join(str(i) for i in outer_lines)}}};",
            f"Curve Loop(2) = {{{', '.join(str(-i) for i in reversed(inner_lines))}}};",
            "Plane Surface(1) = {1, 2};",
            f"ext[] = Extrude {{0, 0, {THICKNESS:.8g}}} {{ Surface{{1}}; Layers{{1}}; Recombine; }};",
            'Physical Surface("frontAndBack") = {1, ext[0]};',
            f'Physical Surface("outlet180") = {{{surf_refs(outlet180)}}};',
            f'Physical Surface("inlet_aux100") = {{{surf_refs(inlet_aux100)}}};',
            f'Physical Surface("inlet_ch_a") = {{{surf_refs(inlet_ch_a)}}};',
            f'Physical Surface("inlet_main40") = {{{surf_refs(inlet_main40)}}};',
            f'Physical Surface("inlet_ch_b") = {{{surf_refs(inlet_ch_b)}}};',
            f'Physical Surface("walls") = {{{surf_refs(wall_indices)}}};',
            'Physical Volume("fluid") = {ext[1]};',
            "",
        ]
    )
    write(geom / "arc_blanket_2d.geo", geo)

    outer_area = polygon_area(outer)
    inner_area = polygon_area(inner)
    fluid_area = outer_area - inner_area
    flibe_volume_est = fluid_area * 2.0 * math.pi * 3.5
    tank_section_volume_est = outer_area * 2.0 * math.pi * 3.5
    manifest = [
        "# Ferrero-Dimensioned ARC 2D Geometry Manifest",
        "",
        "Source: Ferrero thesis Fig. 2.1 dimensioned sketch, PDF page 18 / body page 10.",
        "",
        f"Outer tank mode: `{outer_mode}`.",
        "",
        "This is the first machine-readable scaffold. The outer tank is dimension-anchored;",
        "the inner VV/divertor exclusion is volume-scaled from the prior trace and needs",
        "direct Fig. 2.1/2.6 refinement before promotion.",
        "",
        "| Patch | Intended width | Trace length [m] |",
        "|---|---:|---:|",
        f"| outlet180 | 0.180 | {dist(outer[1], outer[2]):.6f} |",
        f"| inlet_aux100 | 0.100 | {dist(outer[4], outer[5]):.6f} |",
        f"| inlet_ch_a | 0.020 | {dist(inner[4], inner[5]):.6f} |",
        f"| inlet_main40 | 0.040 | {dist(inner[22], inner[23]):.6f} |",
        f"| inlet_ch_b | 0.020 | {dist(inner[12], inner[13]):.6f} |",
        "",
        "## First-Pass Integral Diagnostics",
        "",
        f"- outer section area: `{outer_area:.4f} m2`",
        f"- inner exclusion area: `{inner_area:.4f} m2`",
        f"- fluid section area: `{fluid_area:.4f} m2`",
        f"- outer revolved volume at R=3.5 m: `{tank_section_volume_est:.2f} m3`",
        f"- fluid revolved volume at R=3.5 m: `{flibe_volume_est:.2f} m3`",
        f"- provisional inner-exclusion scale: `{INNER_EXCLUSION_SCALE:.3f}`",
        "",
        "Ferrero Table 2.2 checks: tank volume `353.9 m3`, FLiBe volume `319.0 m3`,",
        "section area `16.0 m2`, mean toroidal radius `3.5 m`.",
        "",
    ]
    write(geom / "ferrero_dimensioned_manifest.md", "\n".join(manifest))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d_ferrero_dimensioned_candidate")
    ap.add_argument("--mesh-size", type=float, default=0.022)
    ap.add_argument("--end-time", type=float, default=0.01)
    ap.add_argument("--delta-t", type=float, default=1.0e-4)
    ap.add_argument("--outer-mode", choices=("area-matched", "figure-shaped"), default="area-matched")
    ap.add_argument(
        "--inlet-mode",
        choices=("legacy-global", "normal", "ferrero-average-normal"),
        default="ferrero-average-normal",
    )
    args = ap.parse_args()

    case = Path(args.case)
    case.mkdir(parents=True, exist_ok=True)
    write_dimensioned_geo(case, args.mesh_size, args.outer_mode)
    write_openfoam_case(case, args.end_time, args.delta_t, args.inlet_mode)
    write_dimensioned_geo(case, args.mesh_size, args.outer_mode)
    write(
        case / "README.ferrero_dimensioned.md",
        f"""# Ferrero-Dimensioned ARC Blanket 2D Candidate

Purpose: first machine-readable ARC 2D scaffold after recovering Ferrero thesis
Fig. 2.1 dimensioned geometry.

This case is separate from the old slide-traced ARC branch. The outer tank
outline is anchored to the printed Fig. 2.1 mm dimensions and checked against
Ferrero Table 2.2. The inner VV/divertor exclusion is still provisional and
volume-scaled from the previous trace, so this case is ready for mesh/BC smoke,
not for final reproduction claims.

- inlet mode: `{args.inlet_mode}`
- outer tank mode: `{args.outer_mode}`
- geometry manifest: `geometry/ferrero_dimensioned_manifest.md`
- source anchor manifest: `references/digitization/ferrero_arc/fig2_1_dimension_anchor_manifest.md`
- source crop: `references/digitization/ferrero_arc/fig2_1_arc_geometric_measures_crop.png`
""",
    )
    print(f"wrote {case}")
    print(f"inlet_mode={args.inlet_mode}")
    print(f"outer_mode={args.outer_mode}")


if __name__ == "__main__":
    main()

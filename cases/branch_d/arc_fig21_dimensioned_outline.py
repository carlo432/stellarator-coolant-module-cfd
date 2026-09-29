#!/usr/bin/env python3
"""Clean, anchor-exact rebuild of Ferrero thesis Fig. 2.1 (ARC tank + VV poloidal section).

Replaces the lumpy traced-polyline scaffold: every anchored mm value from Fig. 2.1 is used
EXACTLY; the VV lens is a smooth Miller-family curve; the divertor legs/feet are analytic
(circles + tangent channels). Pixel-calibrated readings that CORRECT the anchor manifest:
  - 3708 (radial) and 3194 (vertical) dimension lines terminate at the UPPER DIVERTOR FOOT
    CENTER -> foot center = (R 3.708, Z 3.194) m, mirrored below; foot radius 0.5 m ([9]/thesis;
    Ferrero later reduced to 0.4 m in COMSOL -- we draw the sketch's 0.5 m).
  - 5803 is the APEX of a circular outboard bulge spanning steps at Z = +/-1859
    (arc center R=3.655 m on the midplane, radius 2.148 m), not a flat wall corner.
  - Top/bottom tank edges at Z = +/-3876; 4605 wall steps to 4731 at |Z| = 1938 = 3876/2.
  - Left wall R=1521 between Z = +/-2807, chamfering to the top/bottom edge at R=1934.
  - Fitting a Miller curve to the lens's anchored edges recovers R0 ~ 3.28 m ~ 3.3 m (the
    thesis's plasma major radius) -- independent consistency check.
Verification: FLiBe cross-section area / centroid / Pappus volume vs Ferrero Table 2.2.

Outputs: figs/arc_fig21_dimensioned_outline.png  (drawing, source-style)
         figs/arc_fig21_outline_coords.json      (coordinate tables for the CFD rebuild)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

MM = 1e-3
# ---- anchored values from Fig. 2.1 [mm], used exactly ----
R_LEFT_INNER = 1521
R_LEFT_CHAMF = 1934
R_LEFT_LOWER_OUT = 2023          # outer face of lower-left wall (drawn); informational
R0 = 3300                        # plasma axis (dash-dot datum)
R_RIGHT_WALL = 4605
R_RIGHT_STEP = 4731
R_BULGE_APEX = 5803
Z_LENS_SHOULDER = 2338
Z_LEFT_CHAMF = 2807
FOOT_CENTER = (3708, 3194)       # upper divertor foot center; mirrored below
Z_BULGE_STEP = 1859
Z_TOP = 3876
Z_WALL_STEP = Z_TOP // 2         # 1938: 4605->4731 step height (= Z_TOP/2 in the CAD)
WALL = 30
FOOT_R = 500
# ---- thesis analytics ----
A_MINOR = 1130                   # plasma minor radius [mm]
KAPPA = 1.8
DELTA = 0.375

# ---- tank interior outline (FLiBe-wetted, wall inner faces), upper half then mirrored ----
def tank_outline() -> np.ndarray:
    # outboard bulge arc through (4731, +/-1859) with apex (5803, 0)
    rc = ((R_BULGE_APEX**2 - R_RIGHT_STEP**2 - Z_BULGE_STEP**2)
          / (2 * (R_BULGE_APEX - R_RIGHT_STEP)))
    rad = R_BULGE_APEX - rc
    th0 = np.arctan2(Z_BULGE_STEP, R_RIGHT_STEP - rc)
    th = np.linspace(th0, -th0, 61)
    arc = np.c_[rc + rad * np.cos(th), rad * np.sin(th)]
    top_right_chamf = Z_TOP - (R_RIGHT_WALL - 4389)  # 45-deg chamfer, top edge ends at R=4389
    upper = [
        (R_LEFT_INNER, 0),
        (R_LEFT_INNER, Z_LEFT_CHAMF),
        (R_LEFT_CHAMF, Z_TOP),           # left chamfer (as drawn, steeper than 45 deg)
        (4389, Z_TOP),                   # top edge
        (R_RIGHT_WALL, top_right_chamf), # 45-deg top-right chamfer
        (R_RIGHT_WALL, Z_WALL_STEP),
        (R_RIGHT_STEP, Z_WALL_STEP),
        (R_RIGHT_STEP, Z_BULGE_STEP),
    ]
    lower = [(r, -z) for r, z in reversed(upper[1:])]
    # lower-left detail as drawn: wall at 1934 below -2807 with 45-deg chamfer to bottom edge
    lower_left = [
        (R_LEFT_CHAMF, -(Z_TOP - (2425 - R_LEFT_CHAMF))),  # chamfer start (45 deg to bottom edge)
        (2425, -Z_TOP),
    ]
    # splice: replace the naive mirrored left chamfer with the stepped lower-left
    lo = [(r, -z) for r, z in reversed(upper[3:])]         # from bottom edge back up the right side... build explicitly instead
    pts = upper + list(arc) + [
        (R_RIGHT_STEP, -Z_BULGE_STEP),
        (R_RIGHT_STEP, -Z_WALL_STEP),
        (R_RIGHT_WALL, -Z_WALL_STEP),
        (R_RIGHT_WALL, -top_right_chamf),
        (4389, -Z_TOP),
        (2425, -Z_TOP),                   # bottom edge
        (R_LEFT_CHAMF, lower_left[0][1]), # 45-deg chamfer up to the 1934 wall
        (R_LEFT_CHAMF, -Z_LEFT_CHAMF),
        (R_LEFT_INNER, -Z_LEFT_CHAMF),    # step in to the 1521 wall
        (R_LEFT_INNER, 0),
    ]
    return np.asarray(pts, float)

# ---- VV outer surface: Miller-family lens fitted to the anchors ----
# R0 fixed at 3300 (anchor); a_vv from the drawn midplane edges; kappa_vv to hit the
# anchored shoulder height 2338 exactly; delta from the thesis (0.375).
A_VV = 1210
KAPPA_VV = Z_LENS_SHOULDER / A_VV
def miller(R0_, a_, kappa_, delta_, n=241) -> np.ndarray:
    th = np.linspace(-np.pi, np.pi, n)
    return np.c_[R0_ + a_ * np.cos(th + np.arcsin(delta_) * np.sin(th)),
                 kappa_ * a_ * np.sin(th)]

def plasma() -> np.ndarray:
    return miller(R0, A_MINOR, KAPPA, DELTA)

def vv_lens() -> np.ndarray:
    return miller(R0, A_VV, KAPPA_VV, DELTA)

# ---- divertor legs: tangent channel from lens top toward the foot circle ----
def foot_circle(sign: int, n=121, r: float = FOOT_R) -> np.ndarray:
    cx, cz = FOOT_CENTER[0], sign * FOOT_CENTER[1]
    th = np.linspace(0, 2 * np.pi, n)
    return np.c_[cx + r * np.cos(th), cz + r * np.sin(th)]

NECK_HALF_WIDTH = 160  # mm, read from the drawing (~55-65 deg channel per thesis sec. 2.1.2)
def neck(sign: int) -> np.ndarray:
    # channel axis: from the lens top-inner shoulder toward the foot center
    lens = vv_lens()
    top = lens[np.argmax(sign * lens[:, 1])]
    c = np.array([FOOT_CENTER[0], sign * FOOT_CENTER[1]], float)
    ax = c - top
    ax /= np.hypot(*ax)
    nrm = np.array([-ax[1], ax[0]])
    a0, a1 = top - 300 * ax, c            # start slightly inside the lens; end at foot center
    return np.asarray([a0 + NECK_HALF_WIDTH * nrm, a1 + NECK_HALF_WIDTH * nrm,
                       a1 - NECK_HALF_WIDTH * nrm, a0 - NECK_HALF_WIDTH * nrm], float)

# ---- polygon utilities (area, centroid; union via shapely if present, else summed) ----
def poly_area_centroid(p: np.ndarray) -> tuple[float, float]:
    x, y = p[:, 0], p[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    a = cr.sum() / 2
    cx = ((x + x1) * cr).sum() / (6 * a)
    return abs(a), cx

def checks(foot_r: float = FOOT_R) -> dict:
    """Integral checks vs Ferrero Table 2.2. The sketch draws 0.5 m divertor feet, but
    Table 2.2 comes from Ferrero's COMSOL model where he reduced them to 0.4 m -- run both:
    at 0.4 m the area (+0.7%), Pappus volume (-1.9%) and VV plasma-facing lateral area
    (385.2 vs the published 381.2/388.1 inner/outer band) all close within 2%."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    tank = Polygon(tank_outline())
    vv = unary_union([Polygon(vv_lens()),
                      Polygon(neck(+1)), Polygon(foot_circle(+1, r=foot_r)),
                      Polygon(neck(-1)), Polygon(foot_circle(-1, r=foot_r))]).buffer(0)
    flibe = tank.difference(vv)
    A = flibe.area * MM * MM
    rbar = flibe.centroid.x * MM
    V = 2 * np.pi * rbar * A
    pl = Polygon(plasma())
    # lateral (revolved) areas: integral 2*pi*R ds
    def lateral(poly_pts: np.ndarray) -> float:
        p = poly_pts * MM
        ds = np.hypot(*np.diff(p, axis=0).T)
        Rm = 0.5 * (p[:-1, 0] + p[1:, 0])
        return float((2 * np.pi * Rm * ds).sum())
    return dict(
        foot_radius_mm=foot_r,
        flibe_section_area_m2=A, ferrero_section_area_m2=16.0,
        mean_toroidal_radius_m=rbar, ferrero_mean_radius_m=3.5,
        pappus_volume_m3=V, ferrero_tank_volume_m3=353.9,
        vv_plasma_facing_lateral_area_m2=lateral(np.asarray(vv.exterior.coords)),
        ferrero_plasma_vv_areas_m2=[381.2, 388.1],
        miller_boundary_lateral_area_m2=lateral(plasma()),
        plasma_section_area_m2=pl.area * MM * MM,
        vv_exclusion_area_m2=vv.area * MM * MM,
    )

def draw(out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    fig, ax = plt.subplots(figsize=(9.5, 11), facecolor="white")
    t = tank_outline()
    ax.plot(t[:, 0], t[:, 1], "-", color="#1a355e", lw=2.2, label="tank inner wall (anchor-exact)")
    ax.plot(t[:, 0] + np.where(t[:, 0] > R0, WALL, -WALL), t[:, 1], "-",
            color="#1a355e", lw=0.8, alpha=0.45)  # 30 mm wall hint
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    vv = unary_union([Polygon(vv_lens()), Polygon(neck(+1)), Polygon(foot_circle(+1)),
                      Polygon(neck(-1)), Polygon(foot_circle(-1))]).buffer(0)
    vb = np.asarray(vv.exterior.coords)
    ax.fill(vb[:, 0], vb[:, 1], color="#b3cde0", alpha=0.5, zorder=2)
    ax.plot(vb[:, 0], vb[:, 1], "-", color="#7a1f1f", lw=2.0,
            label=f"VV exclusion (Miller fit: R0=3.30, a={A_VV/1e3:.2f}, k={KAPPA_VV:.2f}, d=0.375"
                  " + legs/feet)")
    p = plasma()
    ax.plot(p[:, 0], p[:, 1], "--", color="#b05a00", lw=1.6,
            label="plasma boundary (thesis: a=1.13, k=1.8, d=0.375)")
    for s in (+1, -1):
        ax.add_patch(Circle((FOOT_CENTER[0], s * FOOT_CENTER[1]), 12, color="#7a1f1f"))
    # anchor callouts
    anno = [
        (R_LEFT_INNER, 0, "R=1521"), (R_LEFT_CHAMF, Z_TOP, "R=1934"),
        (R0, 0, "R0=3300"), (R_RIGHT_WALL, Z_WALL_STEP, "4605"),
        (R_RIGHT_STEP, Z_BULGE_STEP, "4731 @ Z=1859"), (R_BULGE_APEX, 0, "R=5803 (arc apex)"),
        (FOOT_CENTER[0], FOOT_CENTER[1], "foot (3708, 3194) r=500"),
        (R_LEFT_INNER, Z_LEFT_CHAMF, "Z=2807"), (2900, Z_LENS_SHOULDER, "Z=2338"),
        (3100, Z_TOP, "Z=3876"),
    ]
    for x, z, s in anno:
        ax.annotate(s, (x, z), textcoords="offset points", xytext=(6, 6), fontsize=8.5,
                    color="0.25", fontstyle="italic")
    ax.axhline(0, color="0.6", lw=0.7, ls="-.")
    ax.axvline(R0, color="0.6", lw=0.7, ls="-.")
    c5, c4 = checks(500), checks(400)
    txt = ("integral checks vs Ferrero Table 2.2 (feet 0.5 m as drawn | 0.4 m as in his COMSOL):\n"
           f"  FLiBe section area {c5['flibe_section_area_m2']:.2f} | {c4['flibe_section_area_m2']:.2f} m$^2$  (16.0)\n"
           f"  mean toroidal radius {c5['mean_toroidal_radius_m']:.2f} | {c4['mean_toroidal_radius_m']:.2f} m  (3.5)\n"
           f"  Pappus tank volume {c5['pappus_volume_m3']:.1f} | {c4['pappus_volume_m3']:.1f} m$^3$  (353.9)\n"
           f"  VV plasma-facing area {c5['vv_plasma_facing_lateral_area_m2']:.1f} | "
           f"{c4['vv_plasma_facing_lateral_area_m2']:.1f} m$^2$  (381.2$-$388.1)")
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=9.5, va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.4"))
    ax.set_aspect("equal"); ax.set_xlabel("R [mm]"); ax.set_ylabel("Z [mm]")
    ax.set_title("ARC tank + VV poloidal section -- anchor-exact rebuild of Ferrero Fig. 2.1\n"
                 "(smooth Miller-family VV, analytic divertor legs; every printed mm anchor exact)",
                 fontsize=11.5)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print("wrote", out)

if __name__ == "__main__":
    figs = Path("figs"); figs.mkdir(exist_ok=True)
    for fr in (500, 400):
        c = checks(fr)
        print(f"--- foot radius {fr} mm " + ("(as drawn in Fig 2.1)" if fr == 500 else "(Ferrero COMSOL, Table 2.2 basis)"))
        for k, v in c.items():
            print(f"  {k}: {v}")
    c = checks(400)
    draw(figs / "arc_fig21_dimensioned_outline.png")
    coords = dict(
        units="mm", midplane="Z=0", datum="machine axis R=0",
        anchors_exact=dict(R_left_inner=R_LEFT_INNER, R_left_chamf=R_LEFT_CHAMF, R0=R0,
                           R_right_wall=R_RIGHT_WALL, R_right_step=R_RIGHT_STEP,
                           R_bulge_apex=R_BULGE_APEX, Z_lens_shoulder=Z_LENS_SHOULDER,
                           Z_left_chamf=Z_LEFT_CHAMF, Z_bulge_step=Z_BULGE_STEP, Z_top=Z_TOP,
                           foot_center=list(FOOT_CENTER), foot_radius=FOOT_R, wall=WALL),
        interpretation_corrections=[
            "3708/3194 locate the UPPER DIVERTOR FOOT CENTER (dimension lines terminate there)",
            "5803 is the apex of a circular outboard bulge between Z=+/-1859 steps (center R=3655, radius 2148)",
            "top/bottom edges at Z=+/-3876; 4605->4731 step at |Z|=1938=3876/2",
            "VV lens is Miller-family; fitting drawn edges recovers R0~3.3 m (anchor-consistent)",
        ],
        estimated_not_anchored=["top edge right end R=4389 (45-deg chamfer)",
                                "bottom-left chamfer start (as drawn, 45-deg to bottom edge)",
                                f"neck half-width {NECK_HALF_WIDTH} mm", "neck axis = lens top -> foot center"],
        tank_inner_wall=tank_outline().tolist(),
        vv_outer=vv_lens().tolist(),
        plasma_boundary=plasma().tolist(),
        upper_neck=neck(+1).tolist(), lower_neck=neck(-1).tolist(),
        upper_foot=foot_circle(+1).tolist(), lower_foot=foot_circle(-1).tolist(),
        integral_checks={"foot_r_500_as_drawn": checks(500), "foot_r_400_ferrero_comsol": checks(400)},
    )
    (figs / "arc_fig21_outline_coords.json").write_text(json.dumps(coords, indent=2))
    print("wrote figs/arc_fig21_outline_coords.json")

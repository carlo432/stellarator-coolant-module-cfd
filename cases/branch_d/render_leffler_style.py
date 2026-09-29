#!/usr/bin/env python3
"""Leffler-style 2D field contours for the stellarator slice: flat mid-gap cross-section, white
background, rainbow colormap, clean colorbar + caption box -- matching the ARC deck aesthetic, but
with the field they listed as future work (temperature) alongside velocity.
"""
import os, numpy as np
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri

pv.OFF_SCREEN = True
open("geom_heat/case.foam", "w").close()
r = pv.OpenFOAMReader("geom_heat/case.foam"); r.set_active_time_value(0.8)
vol = r.read().combine().cell_data_to_point_data()
b = vol.bounds; ymid = (b[2] + b[3]) / 2 + 0.1 * (b[3] - b[2])   # off the degenerate mid-face plane
sl = vol.slice(normal=(0, 1, 0), origin=(0, ymid, 0)).connectivity("largest").triangulate()
print("slice points:", sl.n_points, "at y =", round(ymid, 4))

pts = sl.points
x, z = pts[:, 0], pts[:, 2]
f = sl.faces.reshape(-1, 4)[:, 1:]            # triangles
tri = mtri.Triangulation(x, z, f)
Umag = np.linalg.norm(sl["U"], axis=1)
T = sl["T"]

def panel(ax, field, cmap, clim, title, barlabel, caption):
    lv = np.linspace(clim[0], clim[1], 24)
    cf = ax.tricontourf(tri, np.clip(field, clim[0], clim[1]), levels=lv, cmap=cmap, extend="both")
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=8)
    cb = plt.colorbar(cf, ax=ax, fraction=0.05, pad=0.02, shrink=0.85)
    cb.set_label(barlabel, fontsize=11)
    cb.set_ticks(np.linspace(clim[0], clim[1], 5))
    # caption box (Leffler style)
    ax.text(0.5, -0.06, caption, transform=ax.transAxes, ha="center", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.4", lw=0.8))

fig, axes = plt.subplots(1, 2, figsize=(15, 8), facecolor="white")
panel(axes[0], Umag, "jet", (0, float(np.percentile(Umag, 99))),
      "Coolant Velocity Magnitude", "velocity_mag [m/s]",
      "Mid-Gap Velocity Field\nFLiBe coolant, curved first-wall slice, Re~9900")
panel(axes[1], T, "jet", (900, float(np.percentile(T, 98))),
      "Coolant Temperature", "T [K]",
      "Mid-Gap Temperature Field\n0.5 MW/m2 first-wall load (the ARC deck's 'future work')")
fig.suptitle("Stellarator First-Wall Coolant Channel -- 3D LES Results", fontsize=17, fontweight="bold", y=0.98)
fig.tight_layout(rect=[0, 0.02, 1, 0.95])
fig.savefig("figs/vision_leffler_style.png", dpi=130, facecolor="white", bbox_inches="tight")
print("wrote figs/vision_leffler_style.png  | Umax", round(float(Umag.max()),2), "Tmax", round(float(T.max()),1))

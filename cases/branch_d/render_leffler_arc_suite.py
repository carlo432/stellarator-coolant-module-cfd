#!/usr/bin/env python3
"""Leffler-style result images for the ARC blanket 2D replication (OpenFOAM, ferrero-split F1 case).
Mirrors his result slides: instantaneous velocity (|U|, Ux, Uy) + time-averaged velocity, as flat 2D
contours, white bg, rainbow cmap, caption boxes. Case is a one-cell-thick 2D mesh, so no slicing --
the mesh IS the plane. Honest framing: pimpleFoam FV (not his Nek5000 spectral), slide-traced geometry
(exact CAD/inlet split unrecoverable -- see docs/leffler_arc_geometry_recovery_audit.md).

Run: LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a python3 render_leffler_arc_suite.py [case_dir]
"""
import os, sys, numpy as np
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import pyvista as pv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri

CASE = sys.argv[1] if len(sys.argv) > 1 else "../../cases/arc_blanket_2d_ferrero_split_candidate"
TAG = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(os.path.normpath(CASE))

pv.OFF_SCREEN = True
foam = os.path.join(CASE, "case.foam")
open(foam, "w").close()
r = pv.OpenFOAMReader(foam)
t = max(r.time_values); r.set_active_time_value(t)
vol = r.read().combine().cell_data_to_point_data()
print(f"t={t}  n_points={vol.n_points}  fields={[k for k in vol.point_data]}")

# one-cell-thick 2D mesh: find the thin axis, take the midplane slice there, plot the other two axes
b = vol.bounds
ext = [b[1]-b[0], b[3]-b[2], b[5]-b[4]]
thin = int(np.argmin(ext))
axes2d = [i for i in range(3) if i != thin]
nrm = [0,0,0]; nrm[thin] = 1
org = [0,0,0]; org[thin] = (b[2*thin]+b[2*thin+1])/2
sl = vol.slice(normal=tuple(nrm), origin=tuple(org)).triangulate()
print(f"thin axis={'xyz'[thin]}  slice points={sl.n_points}")
a, c = axes2d
x, y = sl.points[:, a], sl.points[:, c]
tri = mtri.Triangulation(x, y, sl.faces.reshape(-1, 4)[:, 1:])
U, Um = sl["U"], sl["UMean"]
lab = ["x", "y", "z"]

def panel(ax, fld, cmap, clim, title, label):
    lv = np.linspace(clim[0], clim[1], 24)
    cf = ax.tricontourf(tri, np.clip(fld, clim[0], clim[1]), levels=lv, cmap=cmap, extend="both")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=13, fontweight="bold")
    cb = plt.colorbar(cf, ax=ax, fraction=0.045, pad=0.02, shrink=0.9); cb.set_label(label, fontsize=10)
    cb.set_ticks(np.linspace(clim[0], clim[1], 5))

def sym(arr):
    m = float(np.percentile(np.abs(arr), 99)); return (-m, m)

def make(fname, suptitle, caption, specs):
    n = len(specs)
    # tank is closer to square than the duct was -> side-by-side row reads better
    fig, axes = plt.subplots(1, n, figsize=(6.2*n, 6.5), facecolor="white")
    if n == 1: axes = [axes]
    for ax, (fld, cmap, clim, title, label) in zip(axes, specs): panel(ax, fld, cmap, clim, title, label)
    fig.suptitle(suptitle, fontsize=15, fontweight="bold", y=0.98)
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9.5,
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.4", lw=0.8))
    fig.tight_layout(rect=[0, 0.06, 1, 0.94]); fig.savefig(fname, dpi=130, facecolor="white", bbox_inches="tight")
    print("wrote", fname); plt.close(fig)

Umag = np.linalg.norm(U, axis=1); Ummag = np.linalg.norm(Um, axis=1)
vmax = float(np.percentile(Umag, 99))
cap = (("ARC liquid-immersion blanket  |  FLiBe, Ferrero-average inlet split  |  "
       "pimpleFoam LES (FV; Leffler used Nek5000 spectral)  |  t ~ 500 nd-s  |  " + TAG))

make(f"figs/leffler_arc_suite_1_instantaneous_velocity_{TAG}.png",
     "ARC Blanket -- Instantaneous Velocity Field (replication in OpenFOAM)", cap,
     [(Umag, "jet", (0, vmax), "Velocity Magnitude", "|U| [m/s]"),
      (U[:, a], "jet", sym(U[:, a]), f"{lab[a]} Velocity", f"U_{lab[a]} [m/s]"),
      (U[:, c], "jet", sym(U[:, c]), f"{lab[c]} Velocity", f"U_{lab[c]} [m/s]")])

make(f"figs/leffler_arc_suite_2_timeaveraged_velocity_{TAG}.png",
     "ARC Blanket -- Time-Averaged Velocity Field",
     cap + "  |  averaged over final window (2.704-3.704 s)",
     [(Ummag, "jet", (0, vmax), "Mean Velocity Magnitude", "|U_mean| [m/s]"),
      (Um[:, a], "jet", sym(Um[:, a]), f"Mean {lab[a]} Velocity", f"U_{lab[a]},mean [m/s]"),
      (Um[:, c], "jet", sym(Um[:, c]), f"Mean {lab[c]} Velocity", f"U_{lab[c]},mean [m/s]")])

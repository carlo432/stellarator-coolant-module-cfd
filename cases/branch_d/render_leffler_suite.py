#!/usr/bin/env python3
"""Fuller Leffler-style results SUITE for the stellarator slice -- mirrors the ARC deck's multiple
result slides (instantaneous velocity magnitude + components, time-averaged velocity, line plots),
PLUS the temperature field the deck listed as future work. Flat 2D contours, white bg, rainbow cmap.
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
b = vol.bounds; ymid = (b[2]+b[3])/2 + 0.1*(b[3]-b[2])
sl = vol.slice(normal=(0,1,0), origin=(0,ymid,0)).connectivity("largest").triangulate()
x, z = sl.points[:,0], sl.points[:,2]
tri = mtri.Triangulation(x, z, sl.faces.reshape(-1,4)[:,1:])
U, Um = sl["U"], sl["UMean"]; T, Tm = sl["T"], sl["TMean"]

def panel(ax, fld, cmap, clim, title, label):
    lv = np.linspace(clim[0], clim[1], 24)
    cf = ax.tricontourf(tri, np.clip(fld, clim[0], clim[1]), levels=lv, cmap=cmap, extend="both")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=13, fontweight="bold")
    cb = plt.colorbar(cf, ax=ax, fraction=0.045, pad=0.02, shrink=0.9); cb.set_label(label, fontsize=10)
    cb.set_ticks(np.linspace(clim[0], clim[1], 5))

def sym(a):  # symmetric clim for signed components
    m = float(np.percentile(np.abs(a), 99)); return (-m, m)

def make(fname, suptitle, caption, specs):
    n = len(specs); fig, axes = plt.subplots(n, 1, figsize=(9, 3.1*n), facecolor="white")
    if n == 1: axes = [axes]
    for ax,(fld,cmap,clim,title,label) in zip(axes, specs): panel(ax, fld, cmap, clim, title, label)
    fig.suptitle(suptitle, fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.5, 0.005, caption, ha="center", fontsize=9.5,
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.4", lw=0.8))
    fig.tight_layout(rect=[0,0.03,1,0.98]); fig.savefig(fname, dpi=130, facecolor="white", bbox_inches="tight")
    print("wrote", fname); plt.close(fig)

Umag = np.linalg.norm(U,axis=1); Ummag = np.linalg.norm(Um,axis=1)
vmax = float(np.percentile(Umag,99))
make("figs/leffler_suite_1_instantaneous_velocity.png",
     "Instantaneous Velocity Field (3D LES, t = 0.8 s)",
     "Mid-gap cross-section  |  FLiBe coolant, curved stellarator first-wall slice, Re~9900",
     [(Umag,"jet",(0,vmax),"Velocity Magnitude","|U| [m/s]"),
      (U[:,0],"jet",sym(U[:,0]),"Streamwise (x) Velocity","U_x [m/s]"),
      (U[:,2],"jet",sym(U[:,2]),"Cross-stream (z) Velocity","U_z [m/s]")])

make("figs/leffler_suite_2_timeaveraged_velocity.png",
     "Time-Averaged Velocity Field",
     "Mean over the developed window  |  recirculation/secondary flow persists in the average",
     [(Ummag,"jet",(0,vmax),"Mean Velocity Magnitude","|U_mean| [m/s]"),
      (Um[:,0],"jet",sym(Um[:,0]),"Mean Streamwise Velocity","U_x,mean [m/s]"),
      (Um[:,2],"jet",sym(Um[:,2]),"Mean Cross-stream Velocity","U_z,mean [m/s]")])

make("figs/leffler_suite_3_temperature.png",
     "Temperature Field -- the step the ARC deck left as 'future work'",
     "0.5 MW/m2 first-wall load  |  hot layer on the heated wall, cool core (locations resolved; peak = coarse-wall, pin w/ OPT-1+HPC)",
     [(T,"jet",(900,float(np.percentile(T,98))),"Instantaneous Temperature","T [K]"),
      (Tm,"jet",(900,float(np.percentile(Tm,98))),"Time-Averaged Temperature","T_mean [K]")])

#!/usr/bin/env python3
"""Headless 3D render: module geometry (flat-shaded shell) + heated wall colored by temperature.
matplotlib-only. Plot remap (x,y,z)->(x, z, y) so chamber height y is VERTICAL (natural view)."""
import os, re, glob
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm, colors

ROOT = "cases/branch_d"
REMAP = [0, 2, 1]  # physical (x,y,z) -> plot (x, z, y): y becomes vertical

def read_stl(path):
    v = np.array(re.findall(r"vertex\s+([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)", open(path).read()), float)
    return v.reshape(-1, 3, 3)

def read_wall_vtp(path):
    txt = open(path).read()
    pts = np.fromstring(re.search(r"Name=['\"]Points['\"][^>]*>(.*?)</DataArray>", txt, re.DOTALL).group(1), sep=" ").reshape(-1, 3)
    conn = np.fromstring(re.search(r"Name=['\"]connectivity['\"][^>]*>(.*?)</DataArray>", txt, re.DOTALL).group(1), sep=" ").astype(int)
    offs = np.fromstring(re.search(r"Name=['\"]offsets['\"][^>]*>(.*?)</DataArray>", txt, re.DOTALL).group(1), sep=" ").astype(int)
    T = np.fromstring(re.search(r"<CellData>.*?Name=['\"]T['\"][^>]*>(.*?)</DataArray>", txt, re.DOTALL).group(1), sep=" ")
    faces, s = [], 0
    for e in offs:
        faces.append(pts[conn[s:e]][:, REMAP]); s = e
    return faces, T

def flat_shade(tris, light=np.array([0.4, 0.6, 0.7])):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-12)
    s = 0.55 + 0.45 * np.abs(n @ (light / np.linalg.norm(light)))
    return s

def render(stl, wall_vtp, heated_y, out, title):
    tris = read_stl(stl)
    heated = np.all(np.abs(tris[:, :, 1] - heated_y) < 1e-4, axis=1)
    ptris = tris[:, :, REMAP]
    shell = ptris[~heated]
    shade = flat_shade(shell)
    fig = plt.figure(figsize=(12, 5.6), dpi=170)
    ax = fig.add_subplot(111, projection="3d")
    sc = Poly3DCollection(shell, alpha=0.32, linewidths=0)
    sc.set_facecolor(cm.Blues(0.25 + 0.5 * (shade - shade.min()) / (np.ptp(shade) + 1e-9)))
    ax.add_collection3d(sc)
    faces, T = read_wall_vtp(wall_vtp)
    dT = T - 800.0
    norm = colors.Normalize(dT.min(), dT.max())
    pc = Poly3DCollection(faces, facecolor=cm.inferno(norm(dT)), edgecolor="none")
    ax.add_collection3d(pc)
    P = ptris.reshape(-1, 3)
    ax.set_xlim(P[:,0].min(), P[:,0].max()); ax.set_ylim(P[:,1].min(), P[:,1].max()); ax.set_zlim(P[:,2].min(), P[:,2].max())
    ax.set_box_aspect((np.ptp(P[:,0]), np.ptp(P[:,1]) * 1.6, np.ptp(P[:,2]) * 1.6))
    ax.view_init(elev=20, azim=-68)
    ax.set_xlabel("x  (flow direction) [m]"); ax.set_ylabel("z (depth) [m]"); ax.set_zlabel("y (height) [m]")
    ax.set_title(title, fontsize=12)
    m = cm.ScalarMappable(norm=norm, cmap=cm.inferno); m.set_array(dT)
    cb = fig.colorbar(m, ax=ax, shrink=0.55, pad=0.01); cb.set_label("heated-wall ΔT = T−800 K  [K]")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(out); print("wrote", out, f"(peak ΔT={dT.max():.1f} K)")

v23_vtp = sorted(glob.glob(f"{ROOT}/v23_*/VTK/**/boundary/heated_wall.vtp", recursive=True))[0]
render(f"{ROOT}/v4_bent.stl", v23_vtp, -0.04, f"{ROOT}/figs/module3d_bent_hotspot.png",
       "Outlet-offset cooling module — 3D geometry with heated-wall temperature (hotspot)")

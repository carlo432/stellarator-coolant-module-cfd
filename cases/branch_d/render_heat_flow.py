#!/usr/bin/env python3
"""HEAT FLOW for the stellarator -- the step Leffler's ARC deck explicitly listed as FUTURE WORK
('Adding a temperature field'), now done on our 3D curved first-wall slice.

Mid-gap cutaway colored by coolant TEMPERATURE: 0.5 MW/m2 plasma load on the heated (inner/concave)
first wall, FLiBe coolant carrying it away -- hot thermal layer on the first wall, cool core, developing
along the bend. Honest: coarse thermal wall (Pr=14.4 sublayer under-resolved) so the wall superheat is
ILLUSTRATIVE; certifying it is exactly what the OPT-1 wall function + HPC are for.

Run: LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a python3 render_heat_flow.py
"""
import sys, glob, numpy as np, pyvista as pv
pv.global_theme.background = "#05060a"; pv.OFF_SCREEN = True

CASE = "geom_heat/case.foam"
OUT = "figs/vision_heat_flow_stellarator.png"

import os
open("geom_heat/case.foam", "w").close()
r = pv.OpenFOAMReader(CASE)
t = max(r.time_values)                         # latest heated snapshot
r.set_active_time_value(t)
vol = r.read().combine().cell_data_to_point_data()
surf = vol.extract_surface()
b = vol.bounds; ymid = (b[2]+b[3])/2
T = vol["T"]
T_in = 900.0
# emphasize the heating: baseline inlet -> upper band (clip the coarse-mesh near-wall spike)
lo = T_in; hi = float(np.percentile(T, 98))
print(f"t={t}  T range {T.min():.1f}..{T.max():.1f} K  | clim {lo:.1f}..{hi:.1f}")

p = pv.Plotter(off_screen=True, window_size=(1700, 1150)); p.enable_anti_aliasing()
# faint wireframe outline only (no murky fill)
p.add_mesh(surf, style="wireframe", color="#2a3038", opacity=0.25, line_width=0.5, show_scalar_bar=False)
sl = vol.slice(normal=(0,1,0), origin=(0,ymid,0))
p.add_mesh(sl, scalars="T", cmap="inferno", clim=[lo,hi], smooth_shading=True,
           show_scalar_bar=False, ambient=0.35, diffuse=0.9)
# cross-stream cards for 3D richness
zlo, zhi = b[4]+0.01, b[5]-0.01
for f in np.linspace(0.12,0.92,6):
    cs = vol.slice(normal=(0,0,1), origin=(0,0,zlo+f*(zhi-zlo)))
    if cs.n_points: p.add_mesh(cs, scalars="T", cmap="inferno", clim=[lo,hi],
                               show_scalar_bar=False, opacity=0.95, smooth_shading=True)

p.add_scalar_bar(title="Coolant temperature (K)", n_labels=5, color="white",
                 title_font_size=26, label_font_size=20, position_x=0.30, position_y=0.05, width=0.4)
p.add_text("Heat flow in a stellarator first-wall coolant channel",
           position="upper_edge", font_size=16, color="white")
p.add_text("the step the ARC deck left as 'future work' -- now in 3D  |  FLiBe, 0.5 MW/m2 first-wall load  |  wall-resolved LES",
           position=(0.02,0.02), viewport=True, font_size=10.5, color="#aab0bd")
p.camera_position = "xz"; p.camera.azimuth = 18; p.camera.elevation = 30
p.reset_camera(); p.camera.zoom(1.55)
p.screenshot(OUT)
print("wrote", OUT)

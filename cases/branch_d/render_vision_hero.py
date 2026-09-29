#!/usr/bin/env python3
"""VISION hero render v2: clean cutaway of the curved FLiBe stellarator first-wall coolant slice.
A mid-gap longitudinal slice colored by coolant velocity (smooth contour) inside the faint channel walls
+ a few cross-stream profile cards. Reads as professional CFD, not a bundle of streamlines.

Run: LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a python3 render_vision_hero.py
"""
import numpy as np, pyvista as pv
pv.global_theme.background = "#05060a"
pv.OFF_SCREEN = True

CASE = "geom_c11/case.foam"
OUT = "figs/vision_hero_coolant_slice.png"

r = pv.OpenFOAMReader(CASE); r.set_active_time_value(0.6)
vol = r.read().combine().cell_data_to_point_data()
vol["Umag"] = np.linalg.norm(vol["U"], axis=1)
surf = vol.extract_surface()
b = vol.bounds
ymid = (b[2] + b[3]) / 2

p = pv.Plotter(off_screen=True, window_size=(1700, 1150))
p.enable_anti_aliasing()

# faint glass channel walls for 3D context
p.add_mesh(surf, color="#5a6472", opacity=0.10, specular=0.4, smooth_shading=True, show_scalar_bar=False)

# MAIN: longitudinal mid-gap slice colored by coolant speed -- the clean cutaway contour
sl = vol.slice(normal=(0, 1, 0), origin=(0, ymid, 0))
p.add_mesh(sl, scalars="Umag", cmap="inferno", smooth_shading=True, show_scalar_bar=False,
           ambient=0.25, diffuse=0.9, clim=[0, float(vol["Umag"].max())])

# cross-stream "profile cards" along the bend (perpendicular-ish to the duct) for richness
zlo, zhi = b[4] + 0.01, b[5] - 0.01
for frac in np.linspace(0.12, 0.92, 6):
    zc = zlo + frac * (zhi - zlo)
    try:
        cs = vol.slice(normal=(0, 0, 1), origin=(0, 0, zc))
        if cs.n_points > 0:
            p.add_mesh(cs, scalars="Umag", cmap="inferno", show_scalar_bar=False,
                       clim=[0, float(vol["Umag"].max())], opacity=0.92, smooth_shading=True)
    except Exception:
        pass

p.add_scalar_bar(title="Coolant speed (m/s)", n_labels=5, color="white", italic=False,
                 title_font_size=26, label_font_size=20, position_x=0.30, position_y=0.05, width=0.4)
p.add_text("Compact FLiBe-cooled stellarator first-wall coolant channel",
           position="upper_edge", font_size=16, color="white")
p.add_text("wall-resolved LES   |   curved first-wall slice   |   cutaway velocity field   |   Re ~ 9900",
           position=(0.02, 0.02), viewport=True, font_size=11, color="#aab0bd")

p.camera_position = "xz"
p.camera.azimuth = 18; p.camera.elevation = 30
p.reset_camera(); p.camera.zoom(1.55)
p.screenshot(OUT)
print("wrote", OUT, "| Umax =", round(float(vol["Umag"].max()), 2))

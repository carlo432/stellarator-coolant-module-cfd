#!/usr/bin/env python3
"""First-wall temperature map -- the engineering money question Leffler's deck only hypothesized
('recirculating regions may lead to hotspots'): WHERE and HOW HOT does the plasma-facing wall run?
Heated first wall surface colored by T, 0.5 MW/m2 FLiBe-cooled stellarator slice.
Honest: coarse thermal wall (Pr=14.4) -> peak is illustrative; certifying it = OPT-1 WF + HPC.
"""
import numpy as np, pyvista as pv
pv.global_theme.background = "#05060a"; pv.OFF_SCREEN = True
open("geom_heat/case.foam","w").close()
r = pv.OpenFOAMReader("geom_heat/case.foam"); r.set_active_time_value(0.8)
m = r.read()
hw = m["boundary"]["heated_first_wall"].cell_data_to_point_data()
allb = m.combine().extract_surface()
T = hw["T"]
lo = 900.0; hi = float(np.percentile(T, 97))
print(f"first-wall T: min {T.min():.0f} max {T.max():.0f} mean {T.mean():.0f} | clim {lo:.0f}..{hi:.0f}")

p = pv.Plotter(off_screen=True, window_size=(1700,1150)); p.enable_anti_aliasing()
p.add_mesh(allb, style="wireframe", color="#2a3038", opacity=0.18, line_width=0.4, show_scalar_bar=False)
p.add_mesh(hw, scalars="T", cmap="inferno", clim=[lo,hi], smooth_shading=True, show_scalar_bar=False,
           ambient=0.35, diffuse=0.95, specular=0.2)
p.add_scalar_bar(title="First-wall temperature (K)", n_labels=5, color="white",
                 title_font_size=26, label_font_size=20, position_x=0.30, position_y=0.05, width=0.4)
p.add_text("How hot does the plasma-facing first wall run?",
           position="upper_edge", font_size=16, color="white")
p.add_text("FLiBe, 0.5 MW/m2 first-wall load   |   hotspots form where near-wall flow slows -- the ARC deck only hypothesized these; we locate them   |   exact peak: pin with OPT-1 WF + HPC",
           position=(0.012,0.015), viewport=True, font_size=9, color="#aab0bd")
# view the wall BROADSIDE: aim the camera down the wall's average surface normal
hwn = hw.compute_normals(cell_normals=True, point_normals=False, auto_orient_normals=True)
nrm = np.asarray(hwn.cell_data["Normals"]).mean(0); nrm /= (np.linalg.norm(nrm) + 1e-12)
ctr = np.array(hw.center)
diag = np.linalg.norm(np.array(hw.bounds[1::2]) - np.array(hw.bounds[0::2]))
cam = ctr + nrm * diag * 1.6
p.camera.position = tuple(cam); p.camera.focal_point = tuple(ctr); p.camera.up = (0, 1, 0)
p.camera.zoom(1.35)
p.screenshot("figs/vision_firstwall_temp.png")
print("wrote figs/vision_firstwall_temp.png")

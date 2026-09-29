#!/usr/bin/env python3
"""Publication-quality headless 3D render (PyVista/VTK under Xvfb + software GL):
module shell + heated-wall temperature hotspot + recirculation streamlines."""
import os, numpy as np
os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"
os.environ["GALLIUM_DRIVER"] = "llvmpipe"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import pyvista as pv
pv.OFF_SCREEN = True

VTM = ("cases/branch_d/v23_scalar_temperature_flibe_re10000_heatflux100kw_bent_recirculation_dt00025/"
       "VTK/v23_scalar_temperature_flibe_re10000_heatflux100kw_bent_recirculation_dt00025_8000.vtm")
OUT = "cases/branch_d/figs/module3d_pyvista_bent.png"

m = pv.read(VTM)
internal = m["internal"].cell_data_to_point_data()
hw = m["boundary"]["heated_wall"]
Tw = np.asarray(hw.cell_data["T"] if "T" in hw.cell_data else hw["T"]).ravel()
hw.cell_data["dT [K]"] = Tw - 800.0  # cell data preserves the true peak

# seed streamlines from a point cloud filling the inlet duct (x=0, y[-0.01,0.01], z[0,0.04])
ys = np.linspace(-0.0085, 0.0085, 11); zs = np.linspace(0.004, 0.036, 8)
Y, Z = np.meshgrid(ys, zs)
pts = np.column_stack([np.full(Y.size, 0.006), Y.ravel(), Z.ravel()])
seeds = pv.PolyData(pts)
strl = internal.streamlines_from_source(seeds, vectors="U", max_length=6.0,
                                        integration_direction="forward")
strl["|U| [m/s]"] = np.linalg.norm(np.asarray(strl["U"]), axis=1)

shell = internal.extract_surface()

p = pv.Plotter(off_screen=True, window_size=(1800, 1000))
p.set_background("white")
p.add_mesh(shell, color="#c6d4e8", opacity=0.12, smooth_shading=True, show_scalar_bar=False)
p.add_mesh(hw, scalars="dT [K]", cmap="inferno", smooth_shading=True,
           scalar_bar_args=dict(title="heated-wall dT [K]", n_labels=5, fmt="%.0f", color="black",
                                position_x=0.04, position_y=0.05, width=0.30, height=0.045))
if strl.n_cells > 0:
    p.add_mesh(strl.tube(radius=0.0007), scalars="|U| [m/s]", cmap="viridis",
               scalar_bar_args=dict(title="|U| [m/s]", n_labels=5, fmt="%.2f", color="black",
                                    position_x=0.66, position_y=0.05, width=0.30, height=0.045))
p.add_text("Outlet-offset cooling module: recirculation streamlines + heated-wall hotspot\n"
           "FLiBe-like, Re=10,000, 100 kW/m^2  -  peak wall dT = %.0f K" % (Tw.max()-800.0),
           position="upper_edge", font_size=13, color="black")
p.add_axes(color="black")
p.camera_position = "iso"
p.camera.azimuth = 20
p.camera.elevation = 15
p.reset_camera()
p.camera.zoom(1.25)
p.screenshot(OUT)
print("wrote", OUT, "| peak dT=%.1f K | streamline cells=%d" % (Tw.max()-800.0, strl.n_cells))

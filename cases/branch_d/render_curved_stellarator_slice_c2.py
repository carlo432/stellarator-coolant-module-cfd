#!/usr/bin/env python3
"""Render the C2 curved-slice geometry from `geometry_c2.json`."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"
os.environ["GALLIUM_DRIVER"] = "llvmpipe"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import pyvista as pv

pv.OFF_SCREEN = True


def make_poly(vertices: np.ndarray, faces: list[list[int]]) -> pv.PolyData:
    flat = []
    for face in faces:
        flat.extend([len(face), *face])
    return pv.PolyData(vertices, np.asarray(flat, dtype=np.int64))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="c2_curved_stellarator_slice_preflight")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    case = Path(args.case)
    data = json.loads((case / "geometry_c2.json").read_text())
    vertices = np.asarray(data["vertices"], dtype=float)
    patches = data["patches"]
    centerline = np.asarray(data["centerline"], dtype=float)
    out = Path(args.out) if args.out else case.parent / "figs" / f"{case.name}.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    heated = make_poly(vertices, patches["heated_first_wall"])
    back = make_poly(vertices, patches["blanket_back_wall"])
    sides = make_poly(vertices, patches["side_walls"])
    inlet = make_poly(vertices, patches["inlet"])
    outlet = make_poly(vertices, patches["outlet"])
    line = pv.Spline(centerline, len(centerline) * 8)

    p = pv.Plotter(off_screen=True, window_size=(1700, 1050))
    p.set_background("white")
    p.add_mesh(sides, color="#d7dde8", opacity=0.35, smooth_shading=True, show_scalar_bar=False)
    p.add_mesh(back, color="#7197b8", opacity=0.42, smooth_shading=True, show_scalar_bar=False)
    p.add_mesh(heated, color="#d45a32", opacity=0.92, smooth_shading=True, show_scalar_bar=False)
    p.add_mesh(inlet, color="#2f7f5f", opacity=0.95, show_scalar_bar=False)
    p.add_mesh(outlet, color="#6c56a3", opacity=0.95, show_scalar_bar=False)
    p.add_mesh(line.tube(radius=0.0012), color="black", show_scalar_bar=False)
    p.add_text(
        "C2 curved stellarator-relevant coolant slice\n"
        "red = heated first wall, blue = blanket back wall, green/purple = inlet/outlet",
        position="upper_edge",
        font_size=13,
        color="black",
    )
    p.add_axes(color="black")
    p.camera_position = "iso"
    p.camera.azimuth = -35
    p.camera.elevation = 18
    p.reset_camera()
    p.camera.zoom(1.25)
    p.screenshot(str(out))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

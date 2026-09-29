#!/usr/bin/env python3
"""Render a C2 passive-scalar thermal hotspot view from an OpenFOAM case."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

os.environ["LIBGL_ALWAYS_SOFTWARE"] = "1"
os.environ["GALLIUM_DRIVER"] = "llvmpipe"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import pyvista as pv

from postprocess_curved_stellarator_slice_c2 import (
    latest_time_dir,
    patch_face_ids,
    read_boundary,
    read_faces,
    read_internal_scalar,
    read_owner,
    read_points,
)

pv.OFF_SCREEN = True


def make_patch_poly(points: np.ndarray, faces: list[list[int]], face_ids: list[int]) -> pv.PolyData:
    flat: list[int] = []
    for face_id in face_ids:
        face = faces[face_id]
        flat.extend([len(face), *face])
    return pv.PolyData(points, np.asarray(flat, dtype=np.int64))


def face_center(points: np.ndarray, face: list[int]) -> np.ndarray:
    return points[np.asarray(face, dtype=int)].mean(axis=0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="c2_curved_stellarator_slice_thermal_equal_power")
    ap.add_argument("--out", default=None)
    ap.add_argument("--tin", type=float, default=800.0)
    args = ap.parse_args()

    case = Path(args.case)
    latest = latest_time_dir(case)
    points = np.asarray(read_points(case / "constant/polyMesh/points"), dtype=float)
    faces = read_faces(case / "constant/polyMesh/faces")
    owners = read_owner(case / "constant/polyMesh/owner")
    patches = read_boundary(case / "constant/polyMesh/boundary")
    temperature = read_internal_scalar(latest / "T")

    heated_ids = list(patch_face_ids(patches["heated_first_wall"]))
    heated = make_patch_poly(points, faces, heated_ids)
    dts = np.asarray([temperature[owners[face_id]] - args.tin for face_id in heated_ids], dtype=float)
    heated.cell_data["adjacent_dT_K"] = dts

    back = make_patch_poly(points, faces, list(patch_face_ids(patches["blanket_back_wall"])))
    sides = make_patch_poly(points, faces, list(patch_face_ids(patches["side_walls"])))
    inlet = make_patch_poly(points, faces, list(patch_face_ids(patches["inlet"])))
    outlet = make_patch_poly(points, faces, list(patch_face_ids(patches["outlet"])))

    max_idx = int(np.argmax(dts))
    max_face_id = heated_ids[max_idx]
    max_center = face_center(points, faces[max_face_id])
    max_dt = float(dts[max_idx])
    avg_dt = float(dts.mean())

    out = Path(args.out) if args.out else case.parent / "figs" / f"{case.name}_hotspot.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    plotter = pv.Plotter(off_screen=True, window_size=(1700, 1050))
    plotter.set_background("white")
    plotter.add_mesh(sides, color="#d7dde8", opacity=0.18, smooth_shading=True, show_scalar_bar=False)
    plotter.add_mesh(back, color="#6b8fb3", opacity=0.25, smooth_shading=True, show_scalar_bar=False)
    plotter.add_mesh(inlet, color="#2f7f5f", opacity=0.55, show_scalar_bar=False)
    plotter.add_mesh(outlet, color="#6c56a3", opacity=0.55, show_scalar_bar=False)
    plotter.add_mesh(
        heated,
        scalars="adjacent_dT_K",
        cmap="inferno",
        clim=(0.0, max(max_dt, 1e-6)),
        opacity=0.96,
        smooth_shading=True,
        scalar_bar_args={
            "title": "adjacent dT [K]",
            "title_font_size": 13,
            "label_font_size": 11,
            "color": "black",
        },
    )
    plotter.add_mesh(pv.Sphere(radius=0.003, center=max_center), color="#00b7ff", show_scalar_bar=False)
    plotter.add_text(
        f"C2 equal-power passive-scalar hotspot\n"
        f"heated-wall adjacent dT avg {avg_dt:.3f} K, max {max_dt:.3f} K",
        position="upper_edge",
        font_size=13,
        color="black",
    )
    plotter.add_axes(color="black")
    plotter.camera_position = "iso"
    plotter.camera.azimuth = -35
    plotter.camera.elevation = 18
    plotter.reset_camera()
    plotter.camera.zoom(1.25)
    plotter.screenshot(str(out))
    print(f"wrote {out}")
    print(f"heated-wall adjacent dT avg={avg_dt:.6g} K max={max_dt:.6g} K")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render the C2 wall-resolved CHT interface temperature hotspot."""
from __future__ import annotations

import argparse
import os
import re
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
    read_points,
)

pv.OFF_SCREEN = True


def patch_block(text: str, patch_name: str) -> str:
    start = re.search(rf"\b{re.escape(patch_name)}\s*\{{", text)
    if not start:
        raise ValueError(f"patch {patch_name} not found")
    idx = start.end()
    depth = 1
    while idx < len(text) and depth:
        if text[idx] == "{":
            depth += 1
        elif text[idx] == "}":
            depth -= 1
        idx += 1
    return text[start.end() : idx - 1]


def read_patch_value_list(field_path: Path, patch_name: str) -> list[float]:
    block = patch_block(field_path.read_text(), patch_name)
    matches = list(
        re.finditer(
            r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;",
            block,
            re.S,
        )
    )
    if not matches:
        raise ValueError(f"{field_path}: no nonuniform scalar value list for {patch_name}")
    match = matches[-1]
    count = int(match.group(1))
    values = [float(x) for x in re.findall(r"[-+0-9.eE]+", match.group(2))]
    if len(values) != count:
        raise ValueError(f"{field_path}: expected {count} values, parsed {len(values)}")
    return values


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
    ap.add_argument("case", nargs="?", default="c2_curved_stellarator_slice_cht_wallres_smoke")
    ap.add_argument("--region", default="heater")
    ap.add_argument("--patch", default="heater_to_bottomWater")
    ap.add_argument("--tin", type=float, default=800.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    case = Path(args.case)
    latest = latest_time_dir(case)
    region_root = case / "constant" / args.region / "polyMesh"
    time_root = latest / args.region
    points = np.asarray(read_points(region_root / "points"), dtype=float)
    faces = read_faces(region_root / "faces")
    patches = read_boundary(region_root / "boundary")
    face_ids = list(patch_face_ids(patches[args.patch]))
    values = np.asarray(read_patch_value_list(time_root / "T", args.patch), dtype=float)
    dts = values - args.tin

    interface = make_patch_poly(points, faces, face_ids)
    interface.cell_data["interface_dT_K"] = dts
    heated_base = make_patch_poly(points, faces, list(patch_face_ids(patches["heatedBase"])))
    solid_sides = make_patch_poly(points, faces, list(patch_face_ids(patches["solidSides"])))

    max_idx = int(np.argmax(dts))
    max_center = face_center(points, faces[face_ids[max_idx]])
    max_dt = float(dts[max_idx])
    avg_dt = float(dts.mean())

    out = Path(args.out) if args.out else case.parent / "figs" / f"{case.name}_cht_interface_hotspot.png"
    out.parent.mkdir(parents=True, exist_ok=True)

    plotter = pv.Plotter(off_screen=True, window_size=(1700, 1050))
    plotter.set_background("white")
    plotter.add_mesh(solid_sides, color="#c9d0d8", opacity=0.18, smooth_shading=True, show_scalar_bar=False)
    plotter.add_mesh(heated_base, color="#7c8796", opacity=0.25, smooth_shading=True, show_scalar_bar=False)
    plotter.add_mesh(
        interface,
        scalars="interface_dT_K",
        cmap="inferno",
        clim=(0.0, max(max_dt, 1e-6)),
        opacity=0.98,
        smooth_shading=True,
        scalar_bar_args={
            "title": "CHT interface dT [K]",
            "title_font_size": 13,
            "label_font_size": 11,
            "color": "black",
        },
    )
    plotter.add_mesh(pv.Sphere(radius=0.003, center=max_center), color="#00b7ff", show_scalar_bar=False)
    plotter.add_text(
        f"C2 wall-resolved CHT interface hotspot\n"
        f"interface dT avg {avg_dt:.3f} K, max {max_dt:.3f} K",
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
    print(f"interface dT avg={avg_dt:.6g} K max={max_dt:.6g} K")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit ARC 2D inlet velocity directions against OpenFOAM patch normals."""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


U_REF = 2.7
INLET_PATCHES = ("inlet_aux100", "inlet_ch_a", "inlet_ch_b", "inlet_main40")


def strip_header(text: str) -> str:
    i = text.find("\n(")
    if i < 0:
        raise ValueError("could not find OpenFOAM list body")
    return text[i + 1 :]


def read_points(path: Path) -> np.ndarray:
    txt = strip_header(path.read_text())
    vals = re.findall(r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", txt)
    return np.asarray(vals, dtype=float)


def read_faces(path: Path) -> list[list[int]]:
    txt = strip_header(path.read_text())
    faces: list[list[int]] = []
    for m in re.finditer(r"(\d+)\(([^()]*)\)", txt):
        n = int(m.group(1))
        ids = [int(x) for x in m.group(2).split()]
        if len(ids) == n:
            faces.append(ids)
    return faces


def read_boundary(path: Path) -> dict[str, dict[str, int]]:
    txt = path.read_text()
    patches: dict[str, dict[str, int]] = {}
    for m in re.finditer(r"\n\s*(\w+)\s*\{(.*?)\n\s*\}", txt, re.S):
        name = m.group(1)
        body = m.group(2)
        nf = re.search(r"nFaces\s+(\d+);", body)
        sf = re.search(r"startFace\s+(\d+);", body)
        if nf and sf:
            patches[name] = {"nFaces": int(nf.group(1)), "startFace": int(sf.group(1))}
    return patches


def read_fixed_u(path: Path) -> dict[str, np.ndarray]:
    txt = path.read_text()
    out: dict[str, np.ndarray] = {}
    for patch in INLET_PATCHES:
        m = re.search(rf"{patch}\s*\{{[^}}]*value\s+uniform\s+\(([^)]*)\)", txt, re.S)
        if m:
            out[patch] = np.asarray([float(v) for v in m.group(1).split()], dtype=float)
    return out


def face_area_vector(points: np.ndarray, face: list[int]) -> np.ndarray:
    pts = points[np.asarray(face, dtype=int)]
    acc = np.zeros(3)
    for a, b in zip(pts, np.roll(pts, -1, axis=0)):
        acc += np.cross(a, b)
    return 0.5 * acc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="cases/arc_blanket_2d")
    ap.add_argument("--outdir", default="results/arc_blanket_2d")
    args = ap.parse_args()

    case = Path(args.case)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    poly = case / "constant/polyMesh"
    points = read_points(poly / "points")
    faces = read_faces(poly / "faces")
    boundary = read_boundary(poly / "boundary")
    fixed_u = read_fixed_u(case / "0/U")

    rows = []
    for patch in INLET_PATCHES:
        meta = boundary[patch]
        face_ids = range(meta["startFace"], meta["startFace"] + meta["nFaces"])
        sf = np.sum([face_area_vector(points, faces[i]) for i in face_ids], axis=0)
        area = float(np.linalg.norm(sf))
        normal = sf / area if area else np.zeros(3)
        u = fixed_u.get(patch, np.zeros(3))
        speed = float(np.linalg.norm(u))
        outward_flux = float(np.dot(u, sf))
        normal_speed_outward = outward_flux / area if area else float("nan")
        tangential_speed = math.sqrt(max(speed**2 - normal_speed_outward**2, 0.0)) if area else float("nan")
        recommended = -U_REF * normal
        angle = float("nan")
        if speed > 0 and area > 0:
            cosang = np.clip(float(np.dot(u, -normal) / speed), -1.0, 1.0)
            angle = math.degrees(math.acos(cosang))
        rows.append(
            {
                "patch": patch,
                "nFaces": meta["nFaces"],
                "area_m2": area,
                "outward_normal_x": normal[0],
                "outward_normal_y": normal[1],
                "current_Ux": u[0],
                "current_Uy": u[1],
                "current_speed": speed,
                "normal_speed_outward": normal_speed_outward,
                "outward_flux_m3_per_s": outward_flux,
                "tangential_speed": tangential_speed,
                "angle_to_inward_normal_deg": angle,
                "recommended_Ux_for_2p7_normal_inlet": recommended[0],
                "recommended_Uy_for_2p7_normal_inlet": recommended[1],
                "classification": "inflow" if normal_speed_outward < -1e-9 else ("outflow" if normal_speed_outward > 1e-9 else "tangential"),
            }
        )

    csv_path = outdir / "arc2d_inlet_bc_audit.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md_path = outdir / "arc2d_inlet_bc_audit_summary.md"
    lines = [
        "# ARC 2D Inlet Boundary-Condition Audit",
        "",
        f"Case: `{case}`",
        "",
        "This compares the fixed inlet velocity vectors in `0/U` against OpenFOAM patch face-area normals.",
        "Negative normal speed means inflow into the domain; positive means outflow; near-zero means mostly tangential.",
        "",
        "| Patch | Current U [m/s] | Outward normal xy | Normal speed outward [m/s] | Angle to inward normal | Recommended normal-inlet U [m/s] | Class |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        lines.append(
            "| {patch} | `({current_Ux:.3f}, {current_Uy:.3f})` | `({outward_normal_x:.3f}, {outward_normal_y:.3f})` | `{normal_speed_outward:.3f}` | `{angle_to_inward_normal_deg:.1f} deg` | `({recommended_Ux_for_2p7_normal_inlet:.3f}, {recommended_Uy_for_2p7_normal_inlet:.3f})` | {classification} |".format(
                **r
            )
        )
    lines.extend(
        [
            "",
            f"CSV: `{csv_path}`",
            "",
            "Interpretation: use this as a boundary-condition diagnostic only. Changing inlet vectors would require a new run; it would not retroactively change the completed endpoint fields.",
            "",
        ]
    )
    md_path.write_text("\n".join(lines))
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

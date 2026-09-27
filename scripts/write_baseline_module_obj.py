#!/usr/bin/env python3
"""Write a simple OBJ preview of the baseline module footprint.

This is not a CFD mesh. It is a visual sanity-check model for the proposed
expanded heated channel geometry.
"""

from pathlib import Path


WIDTH_INLET = 0.02
WIDTH_CHAMBER = 0.08
WIDTH_OUTLET = 0.02
DEPTH = 0.04
L_INLET = 0.10
L_CHAMBER = 0.20
L_OUTLET = 0.10


def footprint_points() -> list[tuple[float, float]]:
    x0 = 0.0
    x1 = L_INLET
    x2 = L_INLET + L_CHAMBER
    x3 = L_INLET + L_CHAMBER + L_OUTLET

    return [
        (x0, -WIDTH_INLET / 2),
        (x1, -WIDTH_INLET / 2),
        (x1, -WIDTH_CHAMBER / 2),
        (x2, -WIDTH_CHAMBER / 2),
        (x2, -WIDTH_OUTLET / 2),
        (x3, -WIDTH_OUTLET / 2),
        (x3, WIDTH_OUTLET / 2),
        (x2, WIDTH_OUTLET / 2),
        (x2, WIDTH_CHAMBER / 2),
        (x1, WIDTH_CHAMBER / 2),
        (x1, WIDTH_INLET / 2),
        (x0, WIDTH_INLET / 2),
    ]


def main() -> None:
    out_path = Path("cases/baseline_module/geometry/baseline_module_v1.obj")
    pts = footprint_points()

    vertices: list[tuple[float, float, float]] = []
    for z in (0.0, DEPTH):
        for x, y in pts:
            vertices.append((x, y, z))

    faces: list[list[int]] = []
    n = len(pts)

    # bottom and top faces
    faces.append(list(range(1, n + 1)))
    faces.append(list(range(n + 1, 2 * n + 1)))

    # side faces
    for i in range(n):
        j = (i + 1) % n
        faces.append([i + 1, j + 1, j + 1 + n, i + 1 + n])

    with out_path.open("w", encoding="ascii") as f:
        f.write("# Baseline module v1 OBJ preview\n")
        f.write("# Not a CFD mesh\n")
        for x, y, z in vertices:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for face in faces:
            f.write("f " + " ".join(str(idx) for idx in face) + "\n")

    print(out_path)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create a passive scalar C2 thermal-smoke case from a completed C2 flow case."""
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

from postprocess_curved_stellarator_slice_c2 import (
    face_area,
    patch_face_ids,
    read_boundary,
    read_faces,
    read_points,
)


def foam_header(cls: str, obj: str) -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*/\n"
        f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
    )


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def latest_time_dir(case: Path) -> Path:
    time_dirs = []
    for path in case.iterdir():
        if path.is_dir():
            try:
                time_dirs.append((float(path.name), path))
            except ValueError:
                continue
    if not time_dirs:
        raise ValueError(f"{case}: no numeric time directories found")
    return max(time_dirs, key=lambda item: item[0])[1]


def copy_with_location_zero(src: Path, dst: Path) -> None:
    text = src.read_text()
    text = text.replace(f'location    "{src.parent.name}";', 'location    "0";')
    dst.write_text(text)


def patch_area(case: Path, patch_name: str = "heated_first_wall") -> float:
    points = read_points(case / "constant/polyMesh/points")
    faces = read_faces(case / "constant/polyMesh/faces")
    patches = read_boundary(case / "constant/polyMesh/boundary")
    return sum(face_area(faces[face_id], points) for face_id in patch_face_ids(patches[patch_name]))


def read_heated_wall_gradient(t_path: Path) -> float:
    text = t_path.read_text()
    block = re.search(r"\bheated_first_wall\s*\{(.*?)\n\s*\}", text, re.S)
    if not block:
        raise ValueError(f"{t_path}: could not find heated_first_wall block")
    gradient = re.search(r"\bgradient\s+uniform\s+([-+0-9.eE]+)\s*;", block.group(1))
    if not gradient:
        raise ValueError(f"{t_path}: heated_first_wall has no uniform gradient")
    return float(gradient.group(1))


def set_heated_wall_gradient(t_path: Path, gradient: float) -> None:
    text = t_path.read_text()
    pattern = (
        r"(\bheated_first_wall\s*\{.*?\bgradient\s+uniform\s+)"
        r"([-+0-9.eE]+)"
        r"(\s*;.*?\n\s*\})"
    )

    def repl(match: re.Match[str]) -> str:
        return f"{match.group(1)}{gradient:.12g}{match.group(3)}"

    new_text, count = re.subn(pattern, repl, text, count=1, flags=re.S)
    if count != 1:
        raise ValueError(f"{t_path}: could not replace heated_first_wall gradient")
    t_path.write_text(new_text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("flow_case")
    ap.add_argument("thermal_case")
    ap.add_argument("--end-time", type=float, default=0.2)
    ap.add_argument("--delta-t", type=float, default=0.00025)
    ap.add_argument(
        "--heated-wall-gradient",
        type=float,
        help="Override the heated_first_wall fixedGradient value after copying 0/T.",
    )
    ap.add_argument(
        "--match-total-gradient-area-to",
        help=(
            "Reference case whose heated_first_wall gradient-area product should "
            "be matched. This is a passive-scalar equal-total-power surrogate."
        ),
    )
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    if args.heated_wall_gradient is not None and args.match_total_gradient_area_to:
        raise ValueError("use either --heated-wall-gradient or --match-total-gradient-area-to, not both")

    flow_case = Path(args.flow_case)
    thermal_case = Path(args.thermal_case)
    if thermal_case.exists():
        if not args.overwrite:
            raise FileExistsError(f"{thermal_case} exists; use --overwrite to replace it")
        shutil.rmtree(thermal_case)
    thermal_case.mkdir(parents=True)

    shutil.copytree(flow_case / "constant", thermal_case / "constant")
    shutil.copytree(flow_case / "0", thermal_case / "0")
    if (flow_case / "geometry_c2.json").exists():
        shutil.copy2(flow_case / "geometry_c2.json", thermal_case / "geometry_c2.json")

    latest = latest_time_dir(flow_case)
    copy_with_location_zero(latest / "U", thermal_case / "0/U")
    copy_with_location_zero(latest / "phi", thermal_case / "0/phi")

    this_area = patch_area(flow_case)
    source_gradient = read_heated_wall_gradient(thermal_case / "0/T")
    applied_gradient = source_gradient
    setup: dict[str, float | str | None] = {
        "source_flow_case": str(flow_case),
        "source_flow_time": latest.name,
        "thermal_case": str(thermal_case),
        "end_time_s": args.end_time,
        "delta_t_s": args.delta_t,
        "heated_wall_area_m2": this_area,
        "source_heated_wall_gradient": source_gradient,
        "applied_heated_wall_gradient": applied_gradient,
        "gradient_area_product": applied_gradient * this_area,
        "gradient_mode": "copied_from_source",
        "match_reference_case": None,
    }
    if args.heated_wall_gradient is not None:
        applied_gradient = args.heated_wall_gradient
        set_heated_wall_gradient(thermal_case / "0/T", applied_gradient)
        setup.update(
            {
                "applied_heated_wall_gradient": applied_gradient,
                "gradient_area_product": applied_gradient * this_area,
                "gradient_mode": "manual_override",
            }
        )
    elif args.match_total_gradient_area_to:
        reference_case = Path(args.match_total_gradient_area_to)
        reference_area = patch_area(reference_case)
        reference_gradient = read_heated_wall_gradient(reference_case / "0/T")
        applied_gradient = reference_gradient * reference_area / this_area
        set_heated_wall_gradient(thermal_case / "0/T", applied_gradient)
        setup.update(
            {
                "reference_heated_wall_area_m2": reference_area,
                "reference_heated_wall_gradient": reference_gradient,
                "target_gradient_area_product": reference_gradient * reference_area,
                "applied_heated_wall_gradient": applied_gradient,
                "gradient_area_product": applied_gradient * this_area,
                "gradient_mode": "matched_reference_gradient_area",
                "match_reference_case": str(reference_case),
            }
        )

    write(
        thermal_case / "system/controlDict",
        foam_header("dictionary", "controlDict")
        + f"""application scalarTransportFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {args.end_time:g};
deltaT {args.delta_t:g};
writeControl timeStep;
writeInterval {max(round(args.end_time / args.delta_t), 1)};
purgeWrite 0;
writeFormat ascii;
writePrecision 6;
writeCompression off;
timeFormat general;
timePrecision 6;
runTimeModifiable true;
""",
    )
    write(
        thermal_case / "system/fvSchemes",
        foam_header("dictionary", "fvSchemes")
        + """ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { div(phi,T) Gauss limitedLinear 1; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""",
    )
    write(
        thermal_case / "system/fvSolution",
        foam_header("dictionary", "fvSolution")
        + """solvers
{
    T
    {
        solver smoothSolver;
        smoother symGaussSeidel;
        tolerance 1e-8;
        relTol 0;
    }
}
relaxationFactors
{
    fields {}
    equations {}
}
""",
    )
    write(
        thermal_case / "README.c2.md",
        f"""# C2 Passive Thermal Smoke

Generated by `make_c2_thermal_smoke_case.py`.

- source flow case: `{flow_case}`
- source flow time: `{latest.name}`
- scalar end time: `{args.end_time:g} s`
- scalar deltaT: `{args.delta_t:g} s`
- heated-wall area: `{this_area:.8g} m^2`
- applied heated-wall gradient: `{applied_gradient:.8g}`
- gradient mode: `{setup['gradient_mode']}`

Run:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
scalarTransportFoam > log.scalarTransportFoam 2>&1
```
""",
    )
    write(thermal_case / "thermal_smoke_setup.json", json.dumps(setup, indent=2) + "\n")
    print(f"wrote {thermal_case}/ from {flow_case}/ latest={latest.name}")


if __name__ == "__main__":
    main()

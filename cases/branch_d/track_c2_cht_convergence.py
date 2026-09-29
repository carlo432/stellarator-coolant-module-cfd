#!/usr/bin/env python3
"""Track C2 curved CHT temperature drift across written time directories."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


def numeric_time_dirs(case: Path) -> list[Path]:
    times: list[tuple[float, Path]] = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0.0:
            times.append((value, path))
    return [path for _, path in sorted(times, key=lambda item: item[0])]


def patch_average(case: Path, region: str, patch: str, time_name: str, field: str = "T") -> dict[str, float]:
    cmd = (
        "source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 "
        "|| source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1; "
        f"chtMultiRegionSimpleFoam -postProcess -region {region} -time {time_name} "
        f"-func 'patchAverage(name={patch},fields=({field}))'"
    )
    completed = subprocess.run(
        ["bash", "-lc", cmd],
        cwd=case,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    text = completed.stdout
    area = re.search(r"total area\s+=\s+([-+0-9.eE]+)", text)
    avg = re.search(rf"areaAverage\({re.escape(patch)}\) of {field} =\s+([-+0-9.eE]+)", text)
    if not area or not avg:
        raise ValueError(f"could not parse patchAverage output for {region}/{patch}/{field} at {time_name}")
    return {
        "area_m2": float(area.group(1)),
        "area_average": float(avg.group(1)),
    }


def parse_solver_log(case: Path) -> dict[str, dict[str, float]]:
    log = case / "log.cht"
    if not log.exists():
        return {}
    text = log.read_text()
    chunks = re.split(r"^Time =\s+", text, flags=re.M)
    by_time: dict[str, dict[str, float]] = {}
    for chunk in chunks[1:]:
        first_line, _, rest = chunk.partition("\n")
        time_name = first_line.strip()
        minmax = re.findall(r"Min/max T:([-+0-9.eE]+)\s+([-+0-9.eE]+)", rest)
        if len(minmax) >= 2:
            fluid = minmax[-2]
            solid = minmax[-1]
            by_time[time_name] = {
                "fluid_min_T_K": float(fluid[0]),
                "fluid_max_T_K": float(fluid[1]),
                "solid_min_T_K": float(solid[0]),
                "solid_max_T_K": float(solid[1]),
            }
    return by_time


def scalar_internal_minmax(path: Path) -> tuple[float, float]:
    text = path.read_text()
    uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if uniform:
        value = float(uniform.group(1))
        return value, value
    nonuniform = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\(\s*(.*?)\s*\)\s*;",
        text,
        re.S,
    )
    if not nonuniform:
        raise ValueError(f"{path}: could not parse scalar internalField")
    values = [float(x) for x in re.findall(r"[-+0-9.eE]+", nonuniform.group(1))]
    if not values:
        raise ValueError(f"{path}: no scalar internalField values parsed")
    return min(values), max(values)


def percent_change(new: float, old: float) -> float:
    return 100.0 * (new - old) / old if old else 0.0


def write_markdown(summary: dict, path: Path) -> None:
    rows = summary["time_rows"]
    lines = [
        "# C2 Curved CHT Longrun Convergence Track",
        "",
        f"- case: `{summary['case']}`",
        f"- time directories sampled: `{', '.join(row['time'] for row in rows)}`",
        f"- Tin: `{summary['tin_K']}` K",
        "",
        "| Time | Interface mean dT, K | Heated-base mean dT, K | Internal fluid max T, K | Internal solid max T, K |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['time']} | {row['interface_mean_dT_K']:.6g} | "
            f"{row['heatedBase_mean_dT_K']:.6g} | "
            f"{row.get('fluid_max_T_K', float('nan')):.6g} | "
            f"{row.get('solid_max_T_K', float('nan')):.6g} |"
        )
    lines.extend(["", "## Drift Between Written Times", ""])
    if len(rows) < 2:
        lines.append("- Not enough written times for a drift check.")
    else:
        lines.extend(
            [
                "| Interval | Interface dT change | Heated-base dT change | Solid max T change |",
                "|---|---:|---:|---:|",
            ]
        )
        for prev, curr in zip(rows, rows[1:]):
            interface_delta = curr["interface_mean_dT_K"] - prev["interface_mean_dT_K"]
            base_delta = curr["heatedBase_mean_dT_K"] - prev["heatedBase_mean_dT_K"]
            solid_delta = curr.get("solid_max_T_K", 0.0) - prev.get("solid_max_T_K", 0.0)
            lines.append(
                f"| {prev['time']} -> {curr['time']} | "
                f"{interface_delta:+.6g} K ({percent_change(curr['interface_mean_dT_K'], prev['interface_mean_dT_K']):+.3f}%) | "
                f"{base_delta:+.6g} K ({percent_change(curr['heatedBase_mean_dT_K'], prev['heatedBase_mean_dT_K']):+.3f}%) | "
                f"{solid_delta:+.6g} K |"
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            summary["interpretation"],
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default="c2_curved_stellarator_slice_cht_wallres_longrun")
    ap.add_argument("--out-dir", default="figs")
    ap.add_argument("--tin", type=float, default=800.0)
    args = ap.parse_args()

    case = Path(args.case)
    log_metrics = parse_solver_log(case)
    rows = []
    for time_dir in numeric_time_dirs(case):
        time_name = time_dir.name
        heater_interface = patch_average(case, "heater", "heater_to_bottomWater", time_name)
        heated_base = patch_average(case, "heater", "heatedBase", time_name)
        fluid_interface = patch_average(case, "bottomWater", "bottomWater_to_heater", time_name)
        row = {
            "time": time_name,
            "interface_area_m2": heater_interface["area_m2"],
            "interface_mean_T_K": heater_interface["area_average"],
            "interface_mean_dT_K": heater_interface["area_average"] - args.tin,
            "fluid_side_interface_mean_T_K": fluid_interface["area_average"],
            "fluid_side_interface_mean_dT_K": fluid_interface["area_average"] - args.tin,
            "heatedBase_area_m2": heated_base["area_m2"],
            "heatedBase_mean_T_K": heated_base["area_average"],
            "heatedBase_mean_dT_K": heated_base["area_average"] - args.tin,
        }
        row.update(log_metrics.get(time_name, {}))
        fluid_min, fluid_max = scalar_internal_minmax(time_dir / "bottomWater/T")
        solid_min, solid_max = scalar_internal_minmax(time_dir / "heater/T")
        row.update(
            {
                "fluid_min_T_K": fluid_min,
                "fluid_max_T_K": fluid_max,
                "solid_min_T_K": solid_min,
                "solid_max_T_K": solid_max,
            }
        )
        rows.append(row)

    interpretation = (
        "This is a continuation diagnostic, not a formal convergence proof. Use it "
        "to decide whether the wall-resolved CHT case has flattened enough to quote "
        "as a stable temperature scale or whether another continuation/grid point is needed."
    )
    if len(rows) >= 2:
        prev, curr = rows[-2], rows[-1]
        interface_pct = abs(percent_change(curr["interface_mean_dT_K"], prev["interface_mean_dT_K"]))
        base_pct = abs(percent_change(curr["heatedBase_mean_dT_K"], prev["heatedBase_mean_dT_K"]))
        if interface_pct <= 1.0 and base_pct <= 1.0:
            interpretation = (
                "The last written interval changes both interface and heated-base mean dT "
                "by <= 1%, so this is acceptable as a local longrun temperature-scale "
                "artifact. It is still not a grid-converged or validated CHT result."
            )
        else:
            interpretation = (
                "The last written interval still changes the wall-temperature metrics by "
                "> 1%, so keep the convergence caveat and consider another continuation "
                "before treating the C2 CHT temperature scale as stable."
            )

    summary = {
        "case": str(case),
        "tin_K": args.tin,
        "time_rows": rows,
        "interpretation": interpretation,
    }

    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = case.parent / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{case.name}_convergence_track"
    json_path = out_dir / f"{stem}.json"
    md_path = out_dir / f"{stem}.md"
    json_path.write_text(json.dumps(summary, indent=2))
    write_markdown(summary, md_path)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()

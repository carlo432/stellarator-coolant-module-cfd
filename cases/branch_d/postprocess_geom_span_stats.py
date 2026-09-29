#!/usr/bin/env python3
"""OPT-5 spanwise averaging harness for geometry-tier cases.

The curved slice has a side-to-side span direction.  This script groups cells
by the two non-span coordinates and collapses each group across span, reporting
how much sampling scatter remains after span averaging.  It is a postprocess
tool only: it does not claim statistical convergence by itself.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def read_text(path: Path) -> str:
    return path.read_text(errors="replace")


def parse_internal_block(text: str, field_type: str) -> str:
    pattern = rf"internalField\s+nonuniform\s+List<{field_type}>\s*\d+\s*\((.*?)\)\s*;"
    match = re.search(pattern, text, re.S)
    if not match:
        raise ValueError(f"no nonuniform List<{field_type}> internalField")
    return match.group(1)


def read_scalar(path: Path) -> np.ndarray:
    text = read_text(path)
    uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if uniform:
        return np.array([float(uniform.group(1))], dtype=float)
    block = parse_internal_block(text, "scalar")
    return np.array([float(x) for x in block.split()], dtype=float)


def read_vector(path: Path) -> np.ndarray:
    block = parse_internal_block(read_text(path), "vector")
    rows = re.findall(r"\(([^()]*)\)", block)
    return np.array([[float(v) for v in row.split()] for row in rows], dtype=float)


def read_field_magnitude(path: Path) -> tuple[np.ndarray, str]:
    text = read_text(path)
    if "List<vector>" in text:
        return np.linalg.norm(read_vector(path), axis=1), "vector_magnitude"
    if "List<scalar>" in text or re.search(r"internalField\s+uniform\s+[-+0-9.eE]+", text):
        return read_scalar(path), "scalar"
    raise ValueError(f"{path}: unsupported field type")


def numeric_times(case: Path) -> list[str]:
    out: list[tuple[float, str]] = []
    for entry in case.iterdir():
        if not entry.is_dir():
            continue
        try:
            out.append((float(entry.name), entry.name))
        except ValueError:
            continue
    return [name for _, name in sorted(out)]


def latest_time_with(case: Path, required: list[str]) -> str:
    for time_name in reversed(numeric_times(case)):
        tdir = case / time_name
        if all((tdir / f).exists() for f in required):
            return time_name
    raise FileNotFoundError(f"{case}: no time has required fields {required}")


def rounded_key(values: np.ndarray, decimals: int) -> np.ndarray:
    return np.round(values.astype(float), decimals=decimals)


def group_span(
    c0: np.ndarray,
    c1: np.ndarray,
    span: np.ndarray,
    field: np.ndarray,
    decimals: int,
) -> list[dict[str, float]]:
    keys = np.rec.fromarrays(
        [rounded_key(c0, decimals), rounded_key(c1, decimals)],
        names="a,b",
    )
    unique, inverse = np.unique(keys, return_inverse=True)
    rows: list[dict[str, float]] = []
    for idx in range(len(unique)):
        mask = inverse == idx
        vals = field[mask]
        spans = span[mask]
        n = int(mask.sum())
        if n < 1:
            continue
        std = float(vals.std(ddof=1)) if n > 1 else 0.0
        sem = std / math.sqrt(n) if n > 0 else float("nan")
        mean = float(vals.mean())
        rows.append(
            {
                "bin_a": float(unique[idx].a),
                "bin_b": float(unique[idx].b),
                "n_span_cells": n,
                "span_min": float(spans.min()),
                "span_max": float(spans.max()),
                "field_mean": mean,
                "field_std_across_span": std,
                "field_sem_after_span_average": sem,
                "relative_sem": abs(sem / mean) if mean else float("nan"),
            }
        )
    return rows


def summarize(rows: list[dict[str, float]], field: np.ndarray) -> dict[str, float]:
    n_groups = len(rows)
    n_span = np.array([r["n_span_cells"] for r in rows], dtype=float)
    stds = np.array([r["field_std_across_span"] for r in rows], dtype=float)
    sems = np.array([r["field_sem_after_span_average"] for r in rows], dtype=float)
    rel = np.array([r["relative_sem"] for r in rows if np.isfinite(r["relative_sem"])], dtype=float)
    point_std = float(field.std(ddof=1)) if len(field) > 1 else 0.0
    group_means = np.array([r["field_mean"] for r in rows], dtype=float)
    span_mean_std = float(group_means.std(ddof=1)) if len(group_means) > 1 else 0.0
    return {
        "n_cells": int(len(field)),
        "n_span_groups": int(n_groups),
        "median_span_cells_per_group": float(np.median(n_span)) if len(n_span) else 0.0,
        "max_span_cells_per_group": float(n_span.max()) if len(n_span) else 0.0,
        "point_field_std": point_std,
        "span_group_mean_std": span_mean_std,
        "median_within_span_std": float(np.median(stds)) if len(stds) else 0.0,
        "median_span_average_sem": float(np.median(sems)) if len(sems) else 0.0,
        "median_relative_sem": float(np.median(rel)) if len(rel) else float("nan"),
        "nominal_sample_reduction_factor": float(np.median(n_span)) if len(n_span) else 0.0,
    }


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_plot(path: Path, rows: list[dict[str, float]], title: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"plot skipped: {exc}")
        return

    x = np.array([r["bin_a"] for r in rows], dtype=float)
    y = np.array([r["bin_b"] for r in rows], dtype=float)
    sem = np.array([r["field_sem_after_span_average"] for r in rows], dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    sc = ax.scatter(x, y, c=sem, s=24, cmap="viridis")
    ax.set_xlabel("non-span coordinate A")
    ax.set_ylabel("non-span coordinate B")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("span-average standard error")
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", nargs="?", default="geom_c11")
    parser.add_argument("--time")
    parser.add_argument("--field", default="TMean")
    parser.add_argument("--span-axis", choices=("x", "y", "z"), default="y")
    parser.add_argument("--decimals", type=int, default=8)
    parser.add_argument("--outdir", default="figs")
    args = parser.parse_args()

    case = Path(args.case)
    required = ["Cx", "Cy", "Cz", args.field]
    time_name = args.time or latest_time_with(case, required)
    tdir = case / time_name
    coords = {
        "x": read_scalar(tdir / "Cx"),
        "y": read_scalar(tdir / "Cy"),
        "z": read_scalar(tdir / "Cz"),
    }
    field, field_kind = read_field_magnitude(tdir / args.field)
    if any(len(v) != len(field) for v in coords.values()):
        raise ValueError("coordinate and field lengths do not match")

    non_span = [axis for axis in ("x", "y", "z") if axis != args.span_axis]
    rows = group_span(
        coords[non_span[0]],
        coords[non_span[1]],
        coords[args.span_axis],
        field,
        args.decimals,
    )
    if not rows:
        raise ValueError("no span groups generated")
    summary = summarize(rows, field)
    summary.update(
        {
            "case": str(case),
            "time": time_name,
            "field": args.field,
            "field_kind": field_kind,
            "span_axis": args.span_axis,
            "bin_axes": non_span,
            "decimals": args.decimals,
        }
    )

    outdir = Path(args.outdir)
    if not outdir.is_absolute():
        outdir = Path(__file__).resolve().parent / outdir
    stem = f"opt5_span_stats_{case.name}_{args.field}_{time_name}".replace("/", "_")
    csv_path = outdir / f"{stem}.csv"
    json_path = outdir / f"{stem}.json"
    md_path = outdir / f"{stem}.md"
    png_path = outdir / f"{stem}.png"
    write_csv(csv_path, rows)
    json_path.write_text(json.dumps(summary, indent=2) + "\n")
    write_plot(png_path, rows, f"OPT-5 span SEM: {case.name} {args.field} t={time_name}")
    md_path.write_text(
        "# OPT-5 Spanwise Statistics Summary\n\n"
        f"- case: `{case}`\n"
        f"- time: `{time_name}`\n"
        f"- field: `{args.field}` ({field_kind})\n"
        f"- span axis: `{args.span_axis}`\n"
        f"- span groups: `{summary['n_span_groups']}`\n"
        f"- median span cells/group: `{summary['median_span_cells_per_group']:.3g}`\n"
        f"- median span-average SEM: `{summary['median_span_average_sem']:.6g}`\n"
        f"- median relative SEM: `{summary['median_relative_sem']:.6g}`\n"
        f"- nominal sample-reduction factor: `{summary['nominal_sample_reduction_factor']:.3g}`\n\n"
        "Interpretation: this is OPT-5 post-processing evidence only.  It estimates\n"
        "how much cell-to-cell scatter can be reduced by averaging over the spanwise\n"
        "direction; it does not prove time convergence or physical validation.\n"
    )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"wrote {png_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

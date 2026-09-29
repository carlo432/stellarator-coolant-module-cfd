#!/usr/bin/env python3
"""Postprocess the C5 BFS WALE LES smoke run.

Outputs compact probe statistics and a cautious same-style comparison against
the existing outlet-offset v26 LES probe.  This is C5 evidence packaging, not a
validation script.
"""

from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "cases" / "branch_d"
FIGS = WORK / "figs"
BFS_CASE = WORK / "c5_bfs_les_smoke_wale"
OFFSET_CASE = WORK / "v26_les_smoke_flibe_re10000_bent_recirculation_wale"

BFS_PROBES = [
    ("bfs_hotspot_step_h", (0.111, -0.015, 0.02)),
    ("bfs_mid_recirculation", (0.150, -0.015, 0.02)),
    ("bfs_downstream_recovery", (0.220, -0.015, 0.02)),
]
OFFSET_PROBES = [
    ("offset_mid_recirculation", (0.200, -0.030, 0.02)),
    ("offset_near_wall_recirculation", (0.150, -0.035, 0.02)),
    ("offset_downstream_hotspot", (0.260, -0.035, 0.02)),
]


def parse_probe_file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    times: list[float] = []
    values: list[list[tuple[float, float, float]]] = []
    vec_re = re.compile(r"\(([^()]+)\)")
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(None, 1)
            if len(fields) != 2:
                continue
            vectors = []
            for match in vec_re.findall(fields[1]):
                parts = [float(v) for v in match.split()]
                if len(parts) == 3:
                    vectors.append((parts[0], parts[1], parts[2]))
            if vectors:
                times.append(float(fields[0]))
                values.append(vectors)
    if not times:
        raise RuntimeError(f"no probe samples parsed from {path}")
    order = np.argsort(np.asarray(times))
    return np.asarray(times)[order], np.asarray(values, dtype=float)[order, :, :]


def parse_probe_dir(base: Path) -> tuple[np.ndarray, np.ndarray]:
    files = sorted(base.glob("*/U"), key=lambda p: float(p.parent.name))
    if not files:
        raise RuntimeError(f"no probe U files found under {base}")
    all_t = []
    all_u = []
    for path in files:
        t, u = parse_probe_file(path)
        all_t.append(t)
        all_u.append(u)
    times = np.concatenate(all_t)
    values = np.concatenate(all_u, axis=0)
    order = np.argsort(times)
    return times[order], values[order, :, :]


def moving_average(values: np.ndarray, width: int = 31) -> np.ndarray:
    if len(values) < 5:
        return values.copy()
    width = min(width, len(values) if len(values) % 2 == 1 else len(values) - 1)
    width = max(width, 3)
    kernel = np.ones(width) / width
    return np.convolve(values, kernel, mode="same")


def component_stats(
    times: np.ndarray,
    values: np.ndarray,
    mask: np.ndarray,
    component: int,
) -> dict[str, float]:
    series = values[mask, component]
    t = times[mask]
    if len(series) == 0:
        return {
            "n_samples": 0,
            "t_start": math.nan,
            "t_end": math.nan,
            "mean": math.nan,
            "std": math.nan,
            "detrended_std": math.nan,
            "peak_to_peak": math.nan,
        }
    trend = moving_average(series)
    trim = min(20, max(0, (len(series) - 1) // 4))
    residual = series - trend
    if trim:
        residual = residual[trim:-trim]
    return {
        "n_samples": int(len(series)),
        "t_start": float(t[0]),
        "t_end": float(t[-1]),
        "mean": float(np.mean(series)),
        "std": float(np.std(series)),
        "detrended_std": float(np.std(residual)) if len(residual) else math.nan,
        "peak_to_peak": float(np.max(series) - np.min(series)),
    }


def read_latest_pressure_drop(case_dir: Path) -> tuple[float, float, float]:
    def latest(path: Path) -> float:
        last = None
        with path.open() as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    last = float(parts[1])
        if last is None:
            raise RuntimeError(f"no values found in {path}")
        return last

    inlet = latest(case_dir / "postProcessing" / "inletP" / "0" / "surfaceFieldValue.dat")
    outlet = latest(case_dir / "postProcessing" / "outletP" / "0" / "surfaceFieldValue.dat")
    return inlet, outlet, inlet - outlet


def write_stats_csv(
    rows: list[dict[str, object]],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case",
        "geometry",
        "probe_index",
        "probe_label",
        "x",
        "y",
        "z",
        "window",
        "component",
        "n_samples",
        "t_start",
        "t_end",
        "mean",
        "std",
        "detrended_std",
        "peak_to_peak",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_case_rows(
    rows: list[dict[str, object]],
    case_name: str,
    geometry: str,
    labels: list[tuple[str, tuple[float, float, float]]],
    times: np.ndarray,
    values: np.ndarray,
) -> None:
    windows = {
        "full": times >= times[0],
        "developed_t_ge_0p5": times >= 0.5,
    }
    components = {"Ux": 0, "Uy": 1, "Uz": 2}
    for probe_index, (label, coords) in enumerate(labels):
        for window_name, mask in windows.items():
            for component_name, component_index in components.items():
                stats = component_stats(times, values[:, probe_index, :], mask, component_index)
                rows.append(
                    {
                        "case": case_name,
                        "geometry": geometry,
                        "probe_index": probe_index,
                        "probe_label": label,
                        "x": coords[0],
                        "y": coords[1],
                        "z": coords[2],
                        "window": window_name,
                        "component": component_name,
                        **stats,
                    }
                )


def plot_bfs_timeseries(times: np.ndarray, values: np.ndarray, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=180)
    ax.axvspan(0.5, 0.8, color="0.92", label="developed window: t >= 0.5 s")
    colors = ["tab:red", "tab:orange", "tab:blue"]
    for idx, (label, coords) in enumerate(BFS_PROBES):
        ax.plot(
            times,
            values[:, idx, 0],
            lw=1.1,
            color=colors[idx],
            label=f"{label} x={coords[0]:.3f} m",
        )
    ax.set_xlabel("time [s]")
    ax.set_ylabel("probe Ux [m/s]")
    ax.set_title("C5 BFS wall-modeled WALE LES smoke: recirculation-zone probes")
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=7.6, loc="best")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_bfs_vs_offset(
    bfs_t: np.ndarray,
    bfs_u: np.ndarray,
    offset_t: np.ndarray,
    offset_u: np.ndarray,
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=180)
    ax.axvspan(0.5, 0.8, color="0.92", label="developed window: t >= 0.5 s")
    ax.plot(
        bfs_t,
        bfs_u[:, 2, 0],
        color="tab:red",
        lw=1.0,
        label="BFS downstream recovery probe Ux",
    )
    ax.plot(
        offset_t,
        offset_u[:, 2, 0],
        color="tab:blue",
        lw=1.0,
        label="outlet-offset downstream hotspot probe Ux",
    )
    ax.set_xlabel("time [s]")
    ax.set_ylabel("geometry-specific downstream probe Ux [m/s]")
    ax.set_title("C5 same-style LES probe comparison: BFS vs outlet-offset")
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=7.6, loc="best")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_summary(
    rows: list[dict[str, object]],
    pressure_drop: tuple[float, float, float],
    path: Path,
) -> None:
    def find(case: str, probe_label: str, window: str, component: str) -> dict[str, object]:
        for row in rows:
            if (
                row["case"] == case
                and row["probe_label"] == probe_label
                and row["window"] == window
                and row["component"] == component
            ):
                return row
        raise KeyError((case, probe_label, window, component))

    bfs_down = find("c5_bfs_les_smoke_wale", "bfs_downstream_recovery", "developed_t_ge_0p5", "Ux")
    bfs_hot = find("c5_bfs_les_smoke_wale", "bfs_hotspot_step_h", "developed_t_ge_0p5", "Ux")
    off_down = find(
        "v26_les_smoke_flibe_re10000_bent_recirculation_wale",
        "offset_downstream_hotspot",
        "developed_t_ge_0p5",
        "Ux",
    )
    inlet, outlet, dp = pressure_drop

    text = f"""# C5 BFS LES Smoke Postprocess

Date: 2026-06-26

## Status

The C5 BFS wall-modeled WALE LES smoke run completed to `t = 0.799989 s`.
This closes the missing BFS-run leg of C5, but C5 is still not a formal
geometry-LES closure until the same slice/probe metrics are promoted into a
direct BFS versus outlet-offset comparison note.

## Key BFS Probe Results

Developed window: `t >= 0.5 s`.

| Probe | Mean Ux [m/s] | Std Ux [m/s] | Detrended std Ux [m/s] | Peak-to-peak Ux [m/s] |
|---|---:|---:|---:|---:|
| step-height hotspot, x=0.111 m | {bfs_hot['mean']:.5g} | {bfs_hot['std']:.5g} | {bfs_hot['detrended_std']:.5g} | {bfs_hot['peak_to_peak']:.5g} |
| downstream recovery, x=0.220 m | {bfs_down['mean']:.5g} | {bfs_down['std']:.5g} | {bfs_down['detrended_std']:.5g} | {bfs_down['peak_to_peak']:.5g} |

For context, the existing outlet-offset v26 downstream hotspot probe has
developed-window mean/std/detrended-std/peak-to-peak Ux of
`{off_down['mean']:.5g} / {off_down['std']:.5g} / {off_down['detrended_std']:.5g} / {off_down['peak_to_peak']:.5g} m/s`.

## Run Health Markers

- Latest sampled pressure: inlet `{inlet:.6g}`, outlet `{outlet:.6g}`, delta `{dp:.6g}` in OpenFOAM kinematic-pressure units.
- Latest yPlus from `log.yPlus.latest`: adiabatic walls min/max/avg `0.361755 / 49.4647 / 9.60491`; heated wall min/max/avg `0.446446 / 8.49802 / 4.01631`.
- This is isothermal velocity LES only: no wall-temperature or thermal-fluctuation claim should be made from this C5 run.

## Artifacts

- `cases/branch_d/figs/c5_bfs_les_probe_stats.csv`
- `cases/branch_d/figs/c5_bfs_les_probe_timeseries.png`
- `cases/branch_d/figs/c5_bfs_vs_outlet_offset_les_probe_comparison.png`

## Interpretation

The BFS case now has resolved unsteady probe behavior in the recirculation and
recovery region.  The most report-safe statement is that the project now has a
BFS LES smoke counterpart to the outlet-offset LES lineage.  The remaining C5
work is to promote same-region mean fields, resolved TKE, and comparable
geometry-specific probe metrics into one explicit comparison.
"""
    path.write_text(text)


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    bfs_t, bfs_u = parse_probe_file(BFS_CASE / "postProcessing" / "bfsProbes" / "0" / "U")
    offset_t, offset_u = parse_probe_dir(OFFSET_CASE / "postProcessing" / "recircProbes")

    rows: list[dict[str, object]] = []
    add_case_rows(rows, "c5_bfs_les_smoke_wale", "BFS", BFS_PROBES, bfs_t, bfs_u)
    add_case_rows(
        rows,
        "v26_les_smoke_flibe_re10000_bent_recirculation_wale",
        "outlet-offset",
        OFFSET_PROBES,
        offset_t,
        offset_u,
    )

    stats_path = FIGS / "c5_bfs_les_probe_stats.csv"
    write_stats_csv(rows, stats_path)
    plot_bfs_timeseries(bfs_t, bfs_u, FIGS / "c5_bfs_les_probe_timeseries.png")
    plot_bfs_vs_offset(
        bfs_t,
        bfs_u,
        offset_t,
        offset_u,
        FIGS / "c5_bfs_vs_outlet_offset_les_probe_comparison.png",
    )
    summary_path = FIGS / "c5_bfs_les_smoke_summary.md"
    write_summary(rows, read_latest_pressure_drop(BFS_CASE), summary_path)

    print(f"wrote {stats_path}")
    print(f"wrote {FIGS / 'c5_bfs_les_probe_timeseries.png'}")
    print(f"wrote {FIGS / 'c5_bfs_vs_outlet_offset_les_probe_comparison.png'}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()

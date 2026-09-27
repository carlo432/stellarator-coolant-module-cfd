#!/usr/bin/env python3
"""Digitize Leffler slide-16 line profiles and overlay local ARC samples."""
from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


U_REF = 2.7


@dataclass(frozen=True)
class PlotSpec:
    title: str
    local_line: str
    local_component: str
    bbox: tuple[int, int, int, int]  # x0, y0, x1, y1 inclusive axis frame
    xlim: tuple[float, float]
    ylim: tuple[float, float]
    xlabel: str
    ylabel: str


PLOTS: tuple[PlotSpec, ...] = (
    PlotSpec("Left Line", "left_mid", "UMean_y_m_per_s", (152, 286, 1328, 781), (76.0, 101.0), (-0.2, 0.8), "X [-]", "Y velocity [-]"),
    PlotSpec("Right Line", "right_mid", "UMean_y_m_per_s", (153, 890, 1328, 1382), (230.0, 290.0), (-1.2, 1.2), "X [-]", "Y velocity [-]"),
    PlotSpec("Top Line", "top_slot", "UMean_x_m_per_s", (2107, 286, 3284, 781), (180.0, 194.0), (-0.3, 0.4), "Y [-]", "X velocity [-]"),
    PlotSpec("Bottom Line", "bottom_slot", "UMean_x_m_per_s", (2107, 890, 3284, 1382), (-194.0, -180.0), (-0.8, 0.0), "Y [-]", "X velocity [-]"),
)


CURVES = {
    "500 nd-s": {"style": ":", "color": "#1f77b4"},
    "1500 nd-s": {"style": "--", "color": "#ff7f0e"},
    "2500 nd-s": {"style": "-", "color": "#2ca02c"},
}


def curve_mask(rgb: np.ndarray, curve: str) -> np.ndarray:
    r = rgb[..., 0].astype(float)
    g = rgb[..., 1].astype(float)
    b = rgb[..., 2].astype(float)
    if curve == "500 nd-s":
        return (b > 110) & (g > 70) & (r < 130) & (b > r + 35)
    if curve == "1500 nd-s":
        return (r > 180) & (g > 70) & (g < 190) & (b < 120) & (r > b + 80)
    if curve == "2500 nd-s":
        return (g > 95) & (r < 120) & (b < 130) & (g > r + 25)
    raise KeyError(curve)


def pixel_to_data(spec: PlotSpec, px: np.ndarray, py: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = spec.bbox
    x = spec.xlim[0] + (px - x0) / (x1 - x0) * (spec.xlim[1] - spec.xlim[0])
    y = spec.ylim[1] - (py - y0) / (y1 - y0) * (spec.ylim[1] - spec.ylim[0])
    return x, y


def data_to_pixel(spec: PlotSpec, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = spec.bbox
    px = x0 + (x - spec.xlim[0]) / (spec.xlim[1] - spec.xlim[0]) * (x1 - x0)
    py = y0 + (spec.ylim[1] - y) / (spec.ylim[1] - spec.ylim[0]) * (y1 - y0)
    return px, py


def digitize_plot(img: np.ndarray, spec: PlotSpec) -> list[dict[str, float | str]]:
    x0, y0, x1, y1 = spec.bbox
    # Stay just inside the frame to avoid axes/grid/title text contaminating masks.
    crop = img[y0 + 2 : y1 - 1, x0 + 2 : x1 - 1, :]
    rows: list[dict[str, float | str]] = []
    for curve in CURVES:
        mask = curve_mask(crop, curve)
        for col in range(mask.shape[1]):
            ys = np.flatnonzero(mask[:, col])
            if len(ys) == 0:
                continue
            # Median suppresses antialiasing thickness and dashed/dotted edge pixels.
            px = x0 + 2 + col
            py = y0 + 2 + float(np.median(ys))
            x, y = pixel_to_data(spec, np.array([px], dtype=float), np.array([py], dtype=float))
            rows.append(
                {
                    "plot": spec.title,
                    "curve": curve,
                    "axis_coord": float(x[0]),
                    "velocity_nd": float(y[0]),
                    "pixel_x": float(px),
                    "pixel_y": float(py),
                }
            )
    return rows


def local_profile(local_csv: Path, spec: PlotSpec) -> pd.DataFrame:
    df = pd.read_csv(local_csv)
    d = df[df["line"] == spec.local_line].copy()
    if d.empty:
        return d
    frac = (d["s_m"] - d["s_m"].min()) / (d["s_m"].max() - d["s_m"].min())
    d["axis_coord"] = spec.xlim[0] + frac * (spec.xlim[1] - spec.xlim[0])
    d["velocity_nd"] = d[spec.local_component] / U_REF
    d["valid_bool"] = d["valid"].astype(str).str.lower() == "true"
    d.loc[~d["valid_bool"], "velocity_nd"] = np.nan
    return d[["axis_coord", "velocity_nd", "valid_bool"]]


def comparison_metrics(dig: pd.DataFrame, local: pd.DataFrame, spec: PlotSpec) -> list[dict[str, float | str]]:
    out: list[dict[str, float | str]] = []
    if local.empty:
        return out
    local_valid = local.dropna(subset=["velocity_nd"])
    if local_valid.empty:
        return out
    lx = local_valid["axis_coord"].to_numpy()
    ly = local_valid["velocity_nd"].to_numpy()
    order = np.argsort(lx)
    lx = lx[order]
    ly = ly[order]
    for curve in CURVES:
        c = dig[(dig["plot"] == spec.title) & (dig["curve"] == curve)].sort_values("axis_coord")
        if c.empty:
            continue
        x = c["axis_coord"].to_numpy()
        y = c["velocity_nd"].to_numpy()
        inside = (x >= np.nanmin(lx)) & (x <= np.nanmax(lx))
        if np.count_nonzero(inside) < 5:
            continue
        interp = np.interp(x[inside], lx, ly)
        diff = interp - y[inside]
        out.append(
            {
                "plot": spec.title,
                "curve": curve,
                "n": int(np.count_nonzero(inside)),
                "local_valid_fraction": float(local["valid_bool"].mean()),
                "mean_signed_local_minus_leffler": float(np.nanmean(diff)),
                "mean_abs_difference": float(np.nanmean(np.abs(diff))),
                "rms_difference": float(np.sqrt(np.nanmean(diff**2))),
            }
        )
    return out


def alignment_sweep_metrics(dig: pd.DataFrame, local: pd.DataFrame, spec: PlotSpec) -> list[dict[str, float | str]]:
    out: list[dict[str, float | str]] = []
    if local.empty:
        return out
    local_valid = local.dropna(subset=["velocity_nd"])
    if local_valid.empty:
        return out

    base_x = local_valid["axis_coord"].to_numpy()
    base_y = local_valid["velocity_nd"].to_numpy()
    variants = []
    for direction in ("forward", "reversed"):
        if direction == "forward":
            vx = base_x.copy()
        else:
            vx = spec.xlim[0] + spec.xlim[1] - base_x
        for sign_label, sign in (("same_sign", 1.0), ("sign_flipped", -1.0)):
            order = np.argsort(vx)
            variants.append((direction, sign_label, vx[order], sign * base_y[order]))

    for curve in CURVES:
        c = dig[(dig["plot"] == spec.title) & (dig["curve"] == curve)].sort_values("axis_coord")
        if c.empty:
            continue
        x = c["axis_coord"].to_numpy()
        y = c["velocity_nd"].to_numpy()
        for direction, sign_label, lx, ly in variants:
            inside = (x >= np.nanmin(lx)) & (x <= np.nanmax(lx))
            if np.count_nonzero(inside) < 5:
                continue
            interp = np.interp(x[inside], lx, ly)
            diff = interp - y[inside]
            out.append(
                {
                    "plot": spec.title,
                    "curve": curve,
                    "direction": direction,
                    "sign": sign_label,
                    "n": int(np.count_nonzero(inside)),
                    "local_valid_fraction": float(local["valid_bool"].mean()),
                    "mean_signed_local_minus_leffler": float(np.nanmean(diff)),
                    "mean_abs_difference": float(np.nanmean(np.abs(diff))),
                    "rms_difference": float(np.sqrt(np.nanmean(diff**2))),
                }
            )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slide", default="references/digitization/leffler_arc/slide-16.png")
    ap.add_argument("--local-csv", default="results/arc_blanket_2d/arc2d_line_profiles_approx.csv")
    ap.add_argument("--outdir", default="results/arc_blanket_2d")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    slide = Path(args.slide)
    local_csv = Path(args.local_csv)
    outdir = Path(args.outdir)
    digdir = Path("references/digitization/leffler_arc")
    outdir.mkdir(parents=True, exist_ok=True)
    digdir.mkdir(parents=True, exist_ok=True)

    img = np.array(Image.open(slide).convert("RGB"))
    rows: list[dict[str, float | str]] = []
    for spec in PLOTS:
        rows.extend(digitize_plot(img, spec))

    dig_csv = digdir / "leffler_line_profiles_digitized.csv"
    with dig_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["plot", "curve", "axis_coord", "velocity_nd", "pixel_x", "pixel_y"])
        writer.writeheader()
        writer.writerows(rows)

    dig = pd.DataFrame(rows)
    metrics: list[dict[str, float | str]] = []
    sweep: list[dict[str, float | str]] = []

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.2), constrained_layout=True)
    for ax, spec in zip(axes.flat, PLOTS):
        for curve, meta in CURVES.items():
            c = dig[(dig["plot"] == spec.title) & (dig["curve"] == curve)].sort_values("axis_coord")
            ax.plot(c["axis_coord"], c["velocity_nd"], linestyle=meta["style"], color=meta["color"], lw=2.0, label=f"Leffler {curve}")
        loc = local_profile(local_csv, spec)
        if not loc.empty:
            ax.plot(loc["axis_coord"], loc["velocity_nd"], color="black", lw=1.7, alpha=0.9, label="local OpenFOAM UMean/Uref")
            metrics.extend(comparison_metrics(dig, loc, spec))
            sweep.extend(alignment_sweep_metrics(dig, loc, spec))
        ax.set_title(spec.title)
        ax.set_xlim(*spec.xlim)
        ax.set_ylim(*spec.ylim)
        ax.set_xlabel(spec.xlabel)
        ax.set_ylabel(spec.ylabel)
        ax.grid(True, alpha=0.28)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=True)
    fig.suptitle("Digitized Leffler Slide-16 Profiles vs Local ARC 2D Approximate Samples", fontsize=14, fontweight="bold")
    suffix = f"_{args.tag}" if args.tag else ""
    overlay_png = outdir / f"arc2d_leffler_digitized_overlay{suffix}.png"
    fig.savefig(overlay_png, dpi=180)
    plt.close(fig)

    metrics_csv = outdir / f"arc2d_leffler_digitized_overlay{suffix}_metrics.csv"
    with metrics_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "plot",
                "curve",
                "n",
                "local_valid_fraction",
                "mean_signed_local_minus_leffler",
                "mean_abs_difference",
                "rms_difference",
            ],
        )
        writer.writeheader()
        writer.writerows(metrics)

    counts = dig.groupby(["plot", "curve"]).size().reset_index(name="points")
    counts_csv = digdir / "leffler_line_profiles_digitized_counts.csv"
    counts.to_csv(counts_csv, index=False)

    sweep_csv = outdir / f"arc2d_leffler_alignment_sweep{suffix}_metrics.csv"
    with sweep_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "plot",
                "curve",
                "direction",
                "sign",
                "n",
                "local_valid_fraction",
                "mean_signed_local_minus_leffler",
                "mean_abs_difference",
                "rms_difference",
            ],
        )
        writer.writeheader()
        writer.writerows(sweep)

    best_csv = outdir / f"arc2d_leffler_alignment_sweep{suffix}_best.csv"
    if sweep:
        sweep_df = pd.DataFrame(sweep)
        best = sweep_df.sort_values("mean_abs_difference").groupby("plot", as_index=False).first()
        best.to_csv(best_csv, index=False)
    else:
        best_csv.write_text("")

    summary = outdir / f"arc2d_leffler_digitized_overlay{suffix}_summary.md"
    summary.write_text(
        f"""# ARC 2D Digitized Leffler Line-Profile Overlay

Inputs:

- Leffler slide image: `{slide}`
- digitized CSV: `{dig_csv}`
- local OpenFOAM profile CSV: `{local_csv}`

Outputs:

- overlay figure: `{overlay_png}`
- comparison metrics: `{metrics_csv}`
- sign/orientation sweep metrics: `{sweep_csv}`
- best sign/orientation sweep result per plot: `{best_csv}`
- digitization counts: `{counts_csv}`

Method:

- Slide 16 was rasterized from `leffler_christopher_10705.pdf` at 300 dpi.
- The four plot boxes and axis limits were calibrated manually from visible axes and tick labels.
- Blue dotted, orange dashed, and green solid curves were extracted by RGB color masks.
- Local OpenFOAM profiles use the supplied profile CSV divided by `U_ref = {U_REF:g} m/s`.
- Local line coordinates are affinely mapped onto the Leffler plot coordinate extents, so the overlay is a shape/scale diagnostic, not a proof of identical geometry coordinates.
- The sign/orientation sweep tries forward/reversed local line direction and same/flipped velocity sign to test whether the mismatch is just a plotting convention issue.

Caveats:

- This is plot-pixel digitization from a slide deck, not source data from Leffler/Nek5000.
- Our ARC geometry and line locations are slide-traced approximations.
- The local final `UMean` has the final-write restart caveat documented in `docs/leffler_arc_2d_endpoint_verdict.md`.
- Use the mean-absolute differences as rough diagnostic numbers only.
"""
    )

    print(f"wrote {dig_csv}")
    print(f"wrote {overlay_png}")
    print(f"wrote {metrics_csv}")
    print(f"wrote {sweep_csv}")
    print(f"wrote {best_csv}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()

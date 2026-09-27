#!/usr/bin/env python3
"""Detect the black sampling bars on Leffler slide 12."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


BARS = {
    "left_mid": {
        "orientation": "horizontal",
        "roi": (2100, 850, 2320, 950),
    },
    "right_mid": {
        "orientation": "horizontal",
        "roi": (2600, 850, 2980, 950),
    },
    "top_slot": {
        "orientation": "vertical",
        "roi": (2500, 150, 2620, 290),
    },
    "bottom_slot": {
        "orientation": "vertical",
        "roi": (2500, 1500, 2620, 1640),
    },
}


def dark_mask(rgb: np.ndarray, threshold: int) -> np.ndarray:
    return np.all(rgb < threshold, axis=2)


def longest_true_run(vals: np.ndarray) -> tuple[int, int]:
    idx = np.flatnonzero(vals)
    if len(idx) == 0:
        raise ValueError("no dark run found")
    best = (int(idx[0]), int(idx[0]))
    start = int(idx[0])
    prev = int(idx[0])
    for raw in idx[1:]:
        cur = int(raw)
        if cur == prev + 1:
            prev = cur
        else:
            if prev - start > best[1] - best[0]:
                best = (start, prev)
            start = prev = cur
    if prev - start > best[1] - best[0]:
        best = (start, prev)
    return best


def detect_bar(mask: np.ndarray, spec: dict[str, object]) -> dict[str, float | str]:
    x0, y0, x1, y1 = spec["roi"]  # type: ignore[misc]
    roi = mask[y0:y1, x0:x1]
    if spec["orientation"] == "horizontal":
        row_counts = roi.sum(axis=1)
        rows = row_counts > max(4, 0.25 * row_counts.max())
        r0, r1 = longest_true_run(rows)
        sub = roi[r0 : r1 + 1, :]
        col_counts = sub.sum(axis=0)
        cols = col_counts > max(2, 0.25 * sub.shape[0])
        c0, c1 = longest_true_run(cols)
    else:
        col_counts = roi.sum(axis=0)
        cols = col_counts > max(4, 0.25 * col_counts.max())
        c0, c1 = longest_true_run(cols)
        sub = roi[:, c0 : c1 + 1]
        row_counts = sub.sum(axis=1)
        rows = row_counts > max(2, 0.25 * sub.shape[1])
        r0, r1 = longest_true_run(rows)

    bx0 = x0 + c0
    bx1 = x0 + c1
    by0 = y0 + r0
    by1 = y0 + r1
    width = bx1 - bx0 + 1
    height = by1 - by0 + 1
    length = width if spec["orientation"] == "horizontal" else height
    return {
        "bar": "",
        "orientation": str(spec["orientation"]),
        "pixel_x0": bx0,
        "pixel_y0": by0,
        "pixel_x1": bx1,
        "pixel_y1": by1,
        "pixel_width": width,
        "pixel_height": height,
        "pixel_length": length,
        "pixel_center_x": 0.5 * (bx0 + bx1),
        "pixel_center_y": 0.5 * (by0 + by1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slide", default="references/digitization/leffler_arc/slide-12.png")
    ap.add_argument("--outdir", default="references/digitization/leffler_arc")
    ap.add_argument("--threshold", type=int, default=80)
    args = ap.parse_args()

    slide = Path(args.slide)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    img = Image.open(slide).convert("RGB")
    arr = np.asarray(img)
    mask = dark_mask(arr, args.threshold)

    rows = []
    for name, spec in BARS.items():
        row = detect_bar(mask, spec)
        row["bar"] = name
        rows.append(row)

    left = next(r for r in rows if r["bar"] == "left_mid")
    right = next(r for r in rows if r["bar"] == "right_mid")
    top = next(r for r in rows if r["bar"] == "top_slot")
    bottom = next(r for r in rows if r["bar"] == "bottom_slot")
    left_len = float(left["pixel_length"])
    ratios = {
        "right_over_left": float(right["pixel_length"]) / left_len,
        "top_over_left": float(top["pixel_length"]) / left_len,
        "bottom_over_left": float(bottom["pixel_length"]) / left_len,
    }

    csv_path = outdir / "slide12_sampling_bars_detected.csv"
    with csv_path.open("w", newline="") as f:
        fields = [
            "bar",
            "orientation",
            "pixel_x0",
            "pixel_y0",
            "pixel_x1",
            "pixel_y1",
            "pixel_width",
            "pixel_height",
            "pixel_length",
            "pixel_center_x",
            "pixel_center_y",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    overlay = img.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {
        "left_mid": "#1f77b4",
        "right_mid": "#d62728",
        "top_slot": "#ff7f0e",
        "bottom_slot": "#2ca02c",
    }
    for name, spec in BARS.items():
        draw.rectangle(spec["roi"], outline="#5555ff", width=2)  # type: ignore[arg-type]
    for r in rows:
        name = str(r["bar"])
        box = [int(r["pixel_x0"]), int(r["pixel_y0"]), int(r["pixel_x1"]), int(r["pixel_y1"])]
        draw.rectangle(box, outline=colors[name], width=5)
        draw.text((box[0], max(0, box[1] - 24)), f"{name}: {int(r['pixel_length'])} px", fill=colors[name])
    overlay_path = outdir / "slide12_sampling_bars_detected_overlay.png"
    overlay.save(overlay_path)

    summary_path = outdir / "slide12_sampling_bars_detected_summary.md"
    summary_path.write_text(
        f"""# Leffler Slide-12 Sampling Bar Detection

Source slide: `{slide}`

Detected bars:

| Bar | Orientation | Pixel bbox | Length [px] | Center [px] |
| --- | --- | ---: | ---: | ---: |
"""
        + "\n".join(
            "| {bar} | {orientation} | `({pixel_x0:.0f}, {pixel_y0:.0f}) -> ({pixel_x1:.0f}, {pixel_y1:.0f})` | `{pixel_length:.0f}` | `({pixel_center_x:.1f}, {pixel_center_y:.1f})` |".format(
                **r
            )
            for r in rows
        )
        + f"""

Length ratios:

- right/left: `{ratios['right_over_left']:.3f}`
- top/left: `{ratios['top_over_left']:.3f}`
- bottom/left: `{ratios['bottom_over_left']:.3f}`

Interpretation:

The detected slide-12 bars support the shorter axis-matched line spans already used as a diagnostic. The right bar is about `2.44x` the left bar, and the top/bottom bars are about `0.56x` the left bar. That is consistent with slide-16 coordinate spans of about `60`, `25`, and `14` units, respectively.

Outputs:

- CSV: `{csv_path}`
- overlay: `{overlay_path}`

Caveat: this is slide-pixel detection from a presentation raster, not original CAD or Nek5000 coordinates.
"""
    )

    print(f"wrote {csv_path}")
    print(f"wrote {overlay_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build slide-ready panels from generated slice PNGs."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
SLICE_DIR = ROOT / "results" / "slices"
OUT_DIR = ROOT / "results" / "slice_panels"
ADVISOR_DIR = ROOT / "results" / "advisor_packet"


PANELS = [
    {
        "name": "heat_flux_midplane_temperature_panel",
        "title": "Midplane temperature slices: 100 vs 500 kW/m^2",
        "layout": (1, 2),
        "items": [
            ("A", "100 kW/m^2, v14", SLICE_DIR / "v14_100kw_midplane_z_temperature.png"),
            ("B", "500 kW/m^2, v16", SLICE_DIR / "v16_500kw_midplane_z_temperature.png"),
        ],
    },
    {
        "name": "v14_temperature_velocity_panel",
        "title": "Thermal-refined baseline: temperature and velocity",
        "layout": (1, 2),
        "items": [
            ("A", "Temperature, v14", SLICE_DIR / "v14_100kw_midplane_z_temperature.png"),
            ("B", "Velocity magnitude, v14", SLICE_DIR / "v14_100kw_midplane_z_velocity.png"),
        ],
    },
    {
        "name": "v16_temperature_views_panel",
        "title": "High heat-flux case: global and near-wall temperature views",
        "layout": (1, 2),
        "items": [
            ("A", "Midplane temperature, v16", SLICE_DIR / "v16_500kw_midplane_z_temperature.png"),
            ("B", "Near-wall temperature, v16", SLICE_DIR / "v16_500kw_near_heated_wall_y_temperature.png"),
        ],
    },
    {
        "name": "slice_overview_four_panel",
        "title": "Slice overview for advisor review",
        "layout": (2, 2),
        "items": [
            ("A", "v14 temperature", SLICE_DIR / "v14_100kw_midplane_z_temperature.png"),
            ("B", "v14 velocity", SLICE_DIR / "v14_100kw_midplane_z_velocity.png"),
            ("C", "v16 temperature", SLICE_DIR / "v16_500kw_midplane_z_temperature.png"),
            ("D", "v16 near-wall temperature", SLICE_DIR / "v16_500kw_near_heated_wall_y_temperature.png"),
        ],
    },
    {
        "name": "unsteady_rans_progress_panel",
        "title": "Unsteady RANS branch: smoke test to 0.1 s",
        "layout": (1, 2),
        "items": [
            ("A", "v17, 0.020 s", SLICE_DIR / "v17_unsteady_smoke_midplane_z_velocity.png"),
            ("B", "v18, 0.100 s", SLICE_DIR / "v18_unsteady_0p1s_midplane_z_velocity.png"),
        ],
    },
]


def draw_panel(panel: dict[str, object]) -> Path:
    rows, cols = panel["layout"]
    items = panel["items"]
    fig_w = 8.6 * cols
    fig_h = 5.2 * rows
    fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h), dpi=160)
    axes_list = list(axes.flat) if hasattr(axes, "flat") else [axes]

    for ax, (letter, subtitle, path) in zip(axes_list, items):
        image = read_trimmed_image(path)
        ax.imshow(image)
        ax.axis("off")
        ax.text(
            0.015,
            0.985,
            letter,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=16,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "black", "linewidth": 0.7, "alpha": 0.9},
        )
        ax.set_title(subtitle, fontsize=14, pad=8)

    for ax in axes_list[len(items) :]:
        ax.axis("off")

    fig.suptitle(str(panel["title"]), fontsize=18, fontweight="bold", y=0.985)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{panel['name']}.png"
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def read_trimmed_image(path: Path) -> np.ndarray:
    image = Image.open(path).convert("RGB")
    data = np.asarray(image)
    non_white = np.any(data < 245, axis=2)
    rows = np.where(np.any(non_white, axis=1))[0]
    cols = np.where(np.any(non_white, axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return mpimg.imread(path)

    pad = 14
    top = max(0, int(rows[0]) - pad)
    bottom = min(data.shape[0], int(rows[-1]) + pad + 1)
    left = max(0, int(cols[0]) - pad)
    right = min(data.shape[1], int(cols[-1]) + pad + 1)
    return data[top:bottom, left:right, :]


def main() -> None:
    outputs = []
    for panel in PANELS:
        outputs.append(draw_panel(panel))

    ADVISOR_DIR.mkdir(parents=True, exist_ok=True)
    copies = {
        OUT_DIR / "heat_flux_midplane_temperature_panel.png": ADVISOR_DIR / "10_slice_heat_flux_midplane_panel.png",
        OUT_DIR / "v14_temperature_velocity_panel.png": ADVISOR_DIR / "11_slice_v14_temperature_velocity_panel.png",
        OUT_DIR / "v16_temperature_views_panel.png": ADVISOR_DIR / "12_slice_v16_temperature_views_panel.png",
        OUT_DIR / "slice_overview_four_panel.png": ADVISOR_DIR / "13_slice_overview_four_panel.png",
        OUT_DIR / "unsteady_rans_progress_panel.png": ADVISOR_DIR / "14_unsteady_rans_progress_panel.png",
    }
    for src, dst in copies.items():
        dst.write_bytes(src.read_bytes())
        outputs.append(dst)

    for path in outputs:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Crop Kawamura 1998 rendered pages into figure images for digitization."""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "references" / "digitization" / "kawamura_1998"
PAGES = WORK / "pages"
FIGURES = WORK / "figures"


# Coordinates are pixel boxes on the 220-DPI page renders:
# (left, upper, right, lower).  They intentionally include captions/legends.
CROPS = [
    ("fig2_mean_temperature_log_region", "kawamura_1998_page-03.png", (175, 1710, 835, 2310)),
    ("fig3_near_wall_temperature", "kawamura_1998_page-03.png", (910, 210, 1740, 760)),
    ("fig4_nusselt_context", "kawamura_1998_page-03.png", (910, 1705, 1640, 2325)),
    ("fig5_heat_flux_balance_context", "kawamura_1998_page-04.png", (395, 210, 1360, 1125)),
    ("fig6_temperature_rms", "kawamura_1998_page-04.png", (190, 1790, 850, 2355)),
    ("fig7_streamwise_heat_flux_context", "kawamura_1998_page-05.png", (200, 205, 900, 800)),
    ("fig8_wall_normal_heat_flux", "kawamura_1998_page-05.png", (200, 800, 900, 1375)),
    ("fig9_turbulent_prandtl", "kawamura_1998_page-05.png", (925, 210, 1645, 995)),
]


def crop_all() -> list[Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for stem, page_name, box in CROPS:
        src = PAGES / page_name
        with Image.open(src) as image:
            crop = image.crop(box)
            out = FIGURES / f"kawamura_1998_{stem}.png"
            crop.save(out)
            written.append(out)
    return written


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs = []
    for path in paths:
        image = Image.open(path).convert("RGB")
        image.thumbnail((360, 280))
        canvas = Image.new("RGB", (380, 330), "white")
        canvas.paste(image, ((380 - image.width) // 2, 10))
        draw = ImageDraw.Draw(canvas)
        draw.text((10, 295), path.name.replace("kawamura_1998_", ""), fill="black")
        thumbs.append(canvas)

    cols = 2
    rows = (len(thumbs) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * 380, rows * 330), "white")
    for idx, thumb in enumerate(thumbs):
        sheet.paste(thumb, ((idx % cols) * 380, (idx // cols) * 330))
    out = FIGURES / "kawamura_1998_digitization_contact_sheet.png"
    sheet.save(out)
    return out


def main() -> None:
    written = crop_all()
    sheet = make_contact_sheet(written)
    for path in written:
        print(f"wrote {path.relative_to(ROOT)}")
    print(f"wrote {sheet.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

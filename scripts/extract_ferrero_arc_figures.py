#!/usr/bin/env python3
"""Extract focused Ferrero thesis figure crops used for ARC geometry rebuild."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "references" / "digitization" / "ferrero_arc" / "thesis_pages"
OUT_DIR = ROOT / "references" / "digitization" / "ferrero_arc"


CROPS = {
    "fig2_1_arc_geometric_measures_crop.png": {
        "page": "page-18.png",
        "box": (260, 250, 2410, 1800),
        "caption": "Ferrero thesis Fig. 2.1: ARC geometric measures in mm",
    },
    "fig2_6_comsol_geometry_crop.png": {
        "page": "page-21.png",
        "box": (350, 300, 1800, 2350),
        "caption": "Ferrero thesis Fig. 2.6: COMSOL 2D geometry",
    },
}


def add_caption(image: Image.Image, caption: str) -> Image.Image:
    caption_h = 52
    out = Image.new("RGB", (image.width, image.height + caption_h), "white")
    out.paste(image.convert("RGB"), (0, 0))
    draw = ImageDraw.Draw(out)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 24)
    except OSError:
        font = ImageFont.load_default()
    draw.text((16, image.height + 12), caption, fill=(20, 20, 20), font=font)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, spec in CROPS.items():
        src = SRC_DIR / spec["page"]
        if not src.exists():
            raise FileNotFoundError(src)
        with Image.open(src) as page:
            crop = page.crop(spec["box"])
        out = add_caption(crop, spec["caption"])
        out.save(OUT_DIR / name)
        print(f"wrote {OUT_DIR / name}")


if __name__ == "__main__":
    main()

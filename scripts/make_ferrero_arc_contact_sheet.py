#!/usr/bin/env python3
"""Build a contact sheet for Ferrero ARC geometry/inlet recovery pages."""
from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="references/digitization/ferrero_arc")
    ap.add_argument("--out", default="references/digitization/ferrero_arc/ferrero_arc_geometry_recovery_contact_sheet.png")
    ap.add_argument("--thumb-width", type=int, default=420)
    args = ap.parse_args()

    indir = Path(args.indir)
    files = sorted(indir.glob("article_pages/*.png")) + sorted(indir.glob("thesis_pages/*.png"))
    if not files:
        raise SystemExit(f"no png files found under {indir}")

    font = ImageFont.load_default()
    padding = 18
    label_h = 30
    cols = 3
    thumbs: list[tuple[Path, Image.Image]] = []
    for path in files:
        img = Image.open(path).convert("RGB")
        scale = args.thumb_width / img.width
        thumb = img.resize((args.thumb_width, int(img.height * scale)), Image.Resampling.LANCZOS)
        thumbs.append((path, thumb))

    cell_w = args.thumb_width + 2 * padding
    cell_h = max(t.height for _, t in thumbs) + label_h + 2 * padding
    rows = (len(thumbs) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * cell_w, rows * cell_h), "white")
    draw = ImageDraw.Draw(sheet)

    for i, (path, thumb) in enumerate(thumbs):
        row, col = divmod(i, cols)
        x0 = col * cell_w + padding
        y0 = row * cell_h + padding
        label = str(path.relative_to(indir))
        draw.text((x0, y0), label, fill=(0, 0, 0), font=font)
        sheet.paste(thumb, (x0, y0 + label_h))
        draw.rectangle(
            [x0, y0 + label_h, x0 + thumb.width - 1, y0 + label_h + thumb.height - 1],
            outline=(120, 120, 120),
            width=1,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

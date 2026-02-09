#!/usr/bin/env python
"""
make_acamar_2026_card.py

Goal:
- Paste the ACMAR logo onto the 2026 (丙午马年) greeting card
- Add text following the 2024 card layout:
  * Top: HAPPY / NEW / YEAR with small horizontal lines around "NEW"
  * Left: "2026" + "丙午（马年）"
  * Right: vertical circles "新 年 快 乐"
  * Bottom banner: ACMAR logo + Chinese/English organization name

Only uses Python + Pillow (no DALL·E).

Usage (run in the folder that contains the images):
    python make_acamar_2026_card.py

Optional:
    python make_acamar_2026_card.py --bg "ACAMAR New Year 2026 6.png" --logo "ACAMAR.png" -o out.png

If Chinese characters render as □□□ on your machine, pass a CJK font file:
    python make_acamar_2026_card.py --font-cn "/System/Library/Fonts/PingFang.ttc"
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path
from typing import Optional, Sequence, Tuple

from PIL import Image, ImageDraw, ImageFont, ImageFilter


def _pick_first_existing(paths: Sequence[str]) -> Optional[str]:
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None


def _default_cn_font() -> Optional[str]:
    """
    Try common CJK font locations across Linux/Mac/Windows.
    You can override with --font-cn.
    """
    # macOS: KaiTi is often shipped via MobileAsset and not in /System/Library/Fonts.
    mac_kaiti = glob.glob(
        "/System/Library/AssetsV2/com_apple_MobileAsset_Font*/**/AssetData/Kaiti.ttc",
        recursive=True,
    )

    candidates = [
        *mac_kaiti,
        # Linux (Noto)
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/opentype/noto/NotoSerifCJK-Bold.ttc",
        "/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc",
        # macOS
        "/System/Library/Fonts/PingFang.ttc",
        "/System/Library/Fonts/STHeiti Light.ttc",
        "/System/Library/Fonts/STHeiti Medium.ttc",
        # Windows
        r"C:\Windows\Fonts\msyh.ttc",   # Microsoft YaHei
        r"C:\Windows\Fonts\simhei.ttf", # SimHei
        r"C:\Windows\Fonts\simsun.ttc", # SimSun
    ]
    return _pick_first_existing(candidates)


def _default_en_font() -> Optional[str]:
    """
    Try common Latin font locations.
    You can override with --font-en.
    """
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        # macOS (system fonts)
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
        "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        # macOS (installed fonts)
        "/Library/Fonts/Arial Bold.ttf",
        "/Library/Fonts/Arial.ttf",
        "/Library/Fonts/Arial Unicode.ttf",
        os.path.expanduser("~/Library/Fonts/Arial Bold.ttf"),
        os.path.expanduser("~/Library/Fonts/Arial.ttf"),
        r"C:\Windows\Fonts\arialbd.ttf",
        r"C:\Windows\Fonts\arial.ttf",
    ]
    return _pick_first_existing(candidates)


def _default_title_en_font() -> Optional[str]:
    """
    Try a less-curvy serif font for the top title.
    """
    candidates = [
        # macOS
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
        "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        # Windows
        r"C:\Windows\Fonts\timesbd.ttf",
        r"C:\Windows\Fonts\times.ttf",
    ]
    return _pick_first_existing(candidates)


def _ft(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


def _draw_text(
    draw: ImageDraw.ImageDraw,
    xy: Tuple[int, int],
    text: str,
    font: ImageFont.FreeTypeFont,
    fill=(255, 244, 220, 255),
    shadow=(0, 0, 0, 140),
    shadow_offset=(3, 3),
    stroke_width: int = 2,
    stroke_fill=(60, 15, 0, 220),
    anchor: Optional[str] = None,
) -> None:
    x, y = xy
    draw.text(
        (x + shadow_offset[0], y + shadow_offset[1]),
        text,
        font=font,
        fill=shadow,
        stroke_width=stroke_width,
        stroke_fill=stroke_fill,
        anchor=anchor,
    )
    draw.text(
        (x, y),
        text,
        font=font,
        fill=fill,
        stroke_width=stroke_width,
        stroke_fill=stroke_fill,
        anchor=anchor,
    )


def _text_wh(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont) -> Tuple[int, int]:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def make_card(
    bg_path: Path,
    logo_path: Path,
    out_path: Path,
    year: str = "2026",
    ganzhi: str = "丙午",
    animal: str = "马",
    cn_org: str = "中澳天文联合研究中心",
    en_org: str = "Australia-China Consortium for Astrophysical Research",
    font_cn: Optional[str] = None,
    font_en: Optional[str] = None,
    banner_frac: float = 0.17,
    logo_box_frac: float = 0.20,
    margin_x_frac: float = 0.055,
) -> Path:
    bg = Image.open(bg_path).convert("RGBA")
    w, h = bg.size

    overlay = Image.new("RGBA", bg.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # ---------------------------
    # Bottom banner
    # ---------------------------
    # Intentionally NOT drawing a shaded banner: keep the background untouched.
    banner_h = int(banner_frac * h)
    banner_y0 = h - banner_h

    # ---------------------------
    # Logo box (white rounded square with soft shadow)
    # ---------------------------
    margin_x = int(margin_x_frac * w)
    margin_y = int(0.04 * banner_h)
    box_size = int(logo_box_frac * w)
    box = [margin_x, banner_y0 + margin_y, margin_x + box_size, banner_y0 + margin_y + box_size]
    radius = max(8, int(0.07 * box_size))

    shadow = Image.new("RGBA", bg.size, (0, 0, 0, 0))
    sdraw = ImageDraw.Draw(shadow)
    sdraw.rounded_rectangle([box[0] + 6, box[1] + 6, box[2] + 6, box[3] + 6], radius=radius, fill=(0, 0, 0, 120))
    shadow = shadow.filter(ImageFilter.GaussianBlur(8))
    overlay = Image.alpha_composite(shadow, overlay)
    draw = ImageDraw.Draw(overlay)
    draw.rounded_rectangle(box, radius=radius, fill=(255, 255, 255, 255), outline=(255, 255, 255, 255), width=2)

    # Paste logo (scaled to fit)
    logo = Image.open(logo_path).convert("RGBA")
    pad = int(0.08 * box_size)
    avail = box_size - 2 * pad
    lw, lh = logo.size
    scale = min(avail / lw, avail / lh)
    new_sz = (max(1, int(lw * scale)), max(1, int(lh * scale)))
    logo_r = logo.resize(new_sz, Image.Resampling.LANCZOS)
    lx = box[0] + (box_size - new_sz[0]) // 2
    ly = box[1] + (box_size - new_sz[1]) // 2
    overlay.alpha_composite(logo_r, (lx, ly))

    # ---------------------------
    # Fonts
    # ---------------------------
    if font_cn is None:
        font_cn = _default_cn_font()
    if font_en is None:
        font_en = _default_en_font()

    font_top = _default_title_en_font() or font_en

    if font_en is None:
        raise RuntimeError("Cannot find a default English font. Please pass --font-en /path/to/font.ttf")
    if font_cn is None:
        # Chinese can still run, but may render as tofu squares if the chosen font lacks CJK glyphs.
        # We'll fall back to the English font so the script still works.
        font_cn = font_en

    # Size scales (tuned for 1365x2048; scale with image height)
    big_size = int(0.055 * h)     # Top HAPPY/NEW/YEAR
    cn_size = int(0.032 * h)      # Bottom Chinese org name
    en_size = int(0.032 * h)      # Bottom English org name
    side_size = int(0.023 * h) + 2    # Left year label

    f_big = _ft(font_top, big_size)
    f_cn  = _ft(font_cn, cn_size)
    f_en  = _ft(font_en, en_size)
    f_side = _ft(font_cn, side_size)

    # Helper for measuring
    draw = ImageDraw.Draw(overlay)

    # ---------------------------
    # Top: HAPPY / NEW / YEAR (centered)
    # ---------------------------
    lines = ["HAPPY", "NEW", "YEAR"]
    # shrink if too wide
    max_line_w = max(_text_wh(draw, t, f_big)[0] for t in lines)
    if max_line_w > 0.62 * w:
        scale = (0.62 * w) / max_line_w
        f_big = _ft(font_en, max(10, int(big_size * scale)))

    line_h = _text_wh(draw, "Hg", f_big)[1]
    top_y = int(0.025 * h)
    gap = int(0.010 * h)

    for i, t in enumerate(lines):
        tw, th = _text_wh(draw, t, f_big)
        x = (w - tw) // 2
        y = top_y + i * (line_h + gap)
        _draw_text(draw, (x, y), t, f_big, stroke_width=3, shadow_offset=(4, 4))
        if t == "NEW":
            yline = y + th // 2
            padx = int(0.03 * w)
            line_len = int(0.11 * w)
            draw.line([(x - padx - line_len, yline), (x - padx, yline)], fill=(255, 244, 220, 200), width=4)
            draw.line([(x + tw + padx, yline), (x + tw + padx + line_len, yline)], fill=(255, 244, 220, 200), width=4)

    # ---------------------------
    # Left label: year + ganzhi
    # ---------------------------
    left_text = f"{year}\n{ganzhi}（{animal}年）"
    xL = int(0.06 * w)
    yL = int(0.46 * h)
    _draw_text(draw, (xL, yL), left_text, f_side, stroke_width=1, shadow_offset=(2, 2))

    # ---------------------------
    # Right vertical circles: 新年快乐
    # ---------------------------
    chars = list("新年快乐")
    r = int(0.030 * w) + 2
    cx = int(0.92 * w)
    cy0 = int(0.46 * h)
    step = int(1.95 * r)
    circle_fill = (80, 10, 0, 90)
    circle_outline = (255, 244, 220, 160)
    f_char = _ft(font_cn, int(r * 0.95) + 2)

    for i, ch in enumerate(chars):
        cy = cy0 + i * step
        draw.ellipse([cx - r, cy - r, cx + r, cy + r], fill=circle_fill, outline=circle_outline, width=3)
        _draw_text(draw, (cx, cy), ch, f_char, stroke_width=1, shadow_offset=(2, 2), anchor="mm")

    # ---------------------------
    # Bottom org text (centered in remaining space right of logo)
    # ---------------------------
    text_x0 = box[2] + int(0.04 * w)
    text_x1 = w - margin_x
    area_w = text_x1 - text_x0

    # Reduce English font if needed (single line)
    en_w, _ = _text_wh(draw, en_org, f_en)
    if en_w > area_w:
        scale = area_w / en_w
        f_en = _ft(font_en, max(10, int(en_size * scale)))

    cn_w, cn_h = _text_wh(draw, cn_org, f_cn)
    en_w, en_h = _text_wh(draw, en_org, f_en)

    total_h = cn_h + en_h + int(0.012 * h)
    start_y = banner_y0 + (banner_h - total_h) // 2 + int(0.12 * banner_h)

    cn_x = text_x0 + (area_w - cn_w) // 2
    en_x = text_x0 + (area_w - en_w) // 2

    _draw_text(draw, (cn_x, start_y), cn_org, f_cn, stroke_width=2, shadow_offset=(3, 3))
    _draw_text(draw, (en_x, start_y + cn_h + int(0.012 * h)), en_org, f_en, stroke_width=1, shadow_offset=(2, 2))

    # Composite and save
    out = Image.alpha_composite(bg, overlay).convert("RGB")
    out.save(out_path, quality=95)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bg", default="ACAMAR New Year 2026 6.png", help="Background image (2026 card).")
    ap.add_argument("--logo", default="ACAMAR.png", help="ACMAR logo image.")
    ap.add_argument("-o", "--out", default="ACAMAR_New_Year_2026_ACAMAR.png", help="Output image path.")
    ap.add_argument("--year", default="2026")
    ap.add_argument("--ganzhi", default="丙午")
    ap.add_argument("--animal", default="马")
    ap.add_argument("--cn-org", default="中澳天文联合研究中心")
    ap.add_argument("--en-org", default="Australia-China ConsortiuM for Astrophysical Research")
    ap.add_argument("--font-cn", default=None, help="Path to a CJK font (ttf/ttc) to render Chinese.")
    ap.add_argument("--font-en", default=None, help="Path to an English font (ttf/ttc).")
    args = ap.parse_args()

    bg = Path(args.bg)
    logo = Path(args.logo)
    out = Path(args.out)

    if not bg.exists():
        raise FileNotFoundError(f"Cannot find background image: {bg.resolve()}")
    if not logo.exists():
        raise FileNotFoundError(f"Cannot find logo image: {logo.resolve()}")

    out_path = make_card(
        bg_path=bg,
        logo_path=logo,
        out_path=out,
        year=args.year,
        ganzhi=args.ganzhi,
        animal=args.animal,
        cn_org=args.cn_org,
        en_org=args.en_org,
        font_cn=args.font_cn,
        font_en=args.font_en,
    )

    print(f"Saved: {out_path.resolve()}")


if __name__ == "__main__":
    main()

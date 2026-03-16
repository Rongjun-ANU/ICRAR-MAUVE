#!/usr/bin/env python
"""
Build ACAMAR11 combined QR posters by date.

Input filename pattern in working directory:
	ACAMAR11_wechat_YYYYMMDD.png
	ACAMAR11_slack_YYYYMMDD.png

For each date that has both files, this script produces:
	ACAMAR11_QRcodes_YYYYMMDD.png

Layout:
- A4 landscape canvas
- Title on top: "Welcome to ACAMAR 11!"
- WeChat QR on the left, Slack QR on the right
- Slack guidance text under the right QR
- Optional ACAMAR logo (ACAMAR.png) near top-right
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Optional, Tuple

from PIL import Image, ImageDraw, ImageFont


QR_PATTERN = re.compile(r"^ACAMAR11_(wechat|slack)_(\d{8})\.png$", re.IGNORECASE)


def _pick_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
	candidates = [
		"/System/Library/Fonts/Supplemental/Arial.ttf",
		"/System/Library/Fonts/Supplemental/Arial Bold.ttf",
		"/Library/Fonts/Arial.ttf",
		"/Library/Fonts/Arial Bold.ttf",
		"/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
	]
	for font_path in candidates:
		try:
			return ImageFont.truetype(font_path, size=size)
		except OSError:
			continue
	return ImageFont.load_default()


def _fit_image(img: Image.Image, box_w: int, box_h: int) -> Image.Image:
	w, h = img.size
	scale = min(box_w / w, box_h / h)
	new_size = (max(1, int(w * scale)), max(1, int(h * scale)))
	return img.resize(new_size, Image.Resampling.LANCZOS)


def _discover_pairs(input_dir: Path) -> Dict[str, Dict[str, Path]]:
	grouped: Dict[str, Dict[str, Path]] = {}
	for p in input_dir.iterdir():
		if not p.is_file():
			continue
		m = QR_PATTERN.match(p.name)
		if not m:
			continue
		kind = m.group(1).lower()
		date = m.group(2)
		grouped.setdefault(date, {})[kind] = p
	return grouped


def _draw_centered_text(
	draw: ImageDraw.ImageDraw,
	text: str,
	center_x: int,
	y: int,
	font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
	fill: Tuple[int, int, int] = (25, 25, 25),
) -> None:
	bbox = draw.textbbox((0, 0), text, font=font)
	text_w = bbox[2] - bbox[0]
	draw.text((center_x - text_w // 2, y), text, font=font, fill=fill)


def make_poster_for_date(
	date: str,
	wechat_path: Path,
	slack_path: Path,
	output_path: Path,
	logo_path: Optional[Path] = None,
	a4_size_px: Tuple[int, int] = (3508, 2480),
) -> None:
	canvas = Image.new("RGB", a4_size_px, (250, 252, 255))
	draw = ImageDraw.Draw(canvas)
	w, h = canvas.size

	# Light header/footer bands to structure content on A4 landscape.
	draw.rectangle([(0, 0), (w, int(0.17 * h))], fill=(232, 241, 252))
	draw.rectangle([(0, int(0.84 * h)), (w, h)], fill=(241, 247, 255))

	title_font = _pick_font(size=int(0.045 * h))
	label_font = _pick_font(size=int(0.022 * h))
	body_font = _pick_font(size=int(0.019 * h))
	date_font = _pick_font(size=int(0.018 * h))

	_draw_centered_text(
		draw,
		text="Welcome to ACAMAR 11!",
		center_x=w // 2,
		y=int(0.05 * h),
		font=title_font,
		fill=(16, 50, 94),
	)

	_draw_centered_text(
		draw,
		text=f"QR codes for {date}",
		center_x=w // 2,
		y=int(0.115 * h),
		font=date_font,
		fill=(50, 70, 95),
	)

	if logo_path is not None and logo_path.exists():
		logo = Image.open(logo_path).convert("RGBA")
		logo_fit = _fit_image(logo, box_w=int(0.13 * w), box_h=int(0.11 * h))
		lx = w - logo_fit.size[0] - int(0.03 * w)
		ly = int(0.02 * h)
		canvas.paste(logo_fit, (lx, ly), logo_fit)

	wechat_raw = Image.open(wechat_path).convert("RGB")
	slack_raw = Image.open(slack_path).convert("RGB")

	qr_box_w = int(0.31 * w)
	qr_box_h = int(0.50 * h)
	qy = int(0.23 * h)
	wx = int(0.12 * w)
	sx = w - qr_box_w - int(0.12 * w)

	wechat_qr = _fit_image(wechat_raw, qr_box_w, qr_box_h)
	slack_qr = _fit_image(slack_raw, qr_box_w, qr_box_h)

	# White frames around each QR for cleaner print results.
	frame_pad = int(0.012 * h)
	for x in (wx, sx):
		draw.rectangle(
			[
				(x - frame_pad, qy - frame_pad),
				(x + qr_box_w + frame_pad, qy + qr_box_h + frame_pad),
			],
			fill=(255, 255, 255),
			outline=(180, 190, 205),
			width=3,
		)

	wx_img = wx + (qr_box_w - wechat_qr.size[0]) // 2
	wy_img = qy + (qr_box_h - wechat_qr.size[1]) // 2
	sx_img = sx + (qr_box_w - slack_qr.size[0]) // 2
	sy_img = qy + (qr_box_h - slack_qr.size[1]) // 2

	canvas.paste(wechat_qr, (wx_img, wy_img))
	canvas.paste(slack_qr, (sx_img, sy_img))

	_draw_centered_text(
		draw,
		text="WeChat",
		center_x=wx + qr_box_w // 2,
		y=qy + qr_box_h + int(0.03 * h),
		font=label_font,
		fill=(20, 90, 35),
	)
	_draw_centered_text(
		draw,
		text="Slack",
		center_x=sx + qr_box_w // 2,
		y=qy + qr_box_h + int(0.03 * h),
		font=label_font,
		fill=(90, 40, 20),
	)

	# Wrapped guidance text below the Slack QR.
	line1 = "After joining the Slack group, please search"
	line2 = '"acamar-11" channel and join it.'
	_draw_centered_text(
		draw,
		text=line1,
		center_x=sx + qr_box_w // 2,
		y=qy + qr_box_h + int(0.08 * h),
		font=body_font,
		fill=(40, 40, 40),
	)
	_draw_centered_text(
		draw,
		text=line2,
		center_x=sx + qr_box_w // 2,
		y=qy + qr_box_h + int(0.11 * h),
		font=body_font,
		fill=(40, 40, 40),
	)

	output_path.parent.mkdir(parents=True, exist_ok=True)
	canvas.save(output_path, format="PNG")


def main() -> None:
	parser = argparse.ArgumentParser(description="Combine ACAMAR11 WeChat/Slack QR codes by date.")
	parser.add_argument(
		"--input-dir",
		default=".",
		help="Directory containing ACAMAR11_wechat/slack_YYYYMMDD.png files.",
	)
	parser.add_argument(
		"--output-dir",
		default=".",
		help="Directory to save ACAMAR11_QRcodes_YYYYMMDD.png files.",
	)
	parser.add_argument(
		"--logo",
		default="ACAMAR.png",
		help="Path to ACAMAR logo image (optional).",
	)
	args = parser.parse_args()

	input_dir = Path(args.input_dir).resolve()
	output_dir = Path(args.output_dir).resolve()
	logo_path = Path(args.logo)
	if not logo_path.is_absolute():
		logo_path = (input_dir / logo_path).resolve()

	if not input_dir.exists():
		raise FileNotFoundError(f"Input directory not found: {input_dir}")

	grouped = _discover_pairs(input_dir)
	if not grouped:
		print("No matching QR files found.")
		return

	done = 0
	skipped = 0
	for date in sorted(grouped.keys()):
		entry = grouped[date]
		if "wechat" not in entry or "slack" not in entry:
			print(f"Skip {date}: need both wechat and slack QR files.")
			skipped += 1
			continue

		out_name = f"ACAMAR11_QRcodes_{date}.png"
		out_path = output_dir / out_name
		make_poster_for_date(
			date=date,
			wechat_path=entry["wechat"],
			slack_path=entry["slack"],
			output_path=out_path,
			logo_path=logo_path if logo_path.exists() else None,
		)
		done += 1
		print(f"Saved: {out_path}")

	print(f"Finished. Generated {done} file(s), skipped {skipped} date(s).")


if __name__ == "__main__":
	main()

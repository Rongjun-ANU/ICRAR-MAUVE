#!/usr/bin/env python

from __future__ import annotations

import glob
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


try:
	from PIL import Image
except ModuleNotFoundError as e:  # pragma: no cover
	raise SystemExit(
		"Missing dependency: Pillow (PIL).\n"
		f"Interpreter: {sys.executable}\n"
		"Install it into this interpreter (e.g. `python -m pip install pillow`) or run via the "
		"conda env's `python` where Pillow is installed."
	) from e


@dataclass(frozen=True)
class Rect:
	x: int
	y: int
	w: int
	h: int

	@property
	def right(self) -> int:
		return self.x + self.w

	@property
	def bottom(self) -> int:
		return self.y + self.h


def _intersects(a: Rect, b: Rect) -> bool:
	return not (a.right <= b.x or b.right <= a.x or a.bottom <= b.y or b.bottom <= a.y)


def _contains(outer: Rect, inner: Rect) -> bool:
	return (
		inner.x >= outer.x
		and inner.y >= outer.y
		and inner.right <= outer.right
		and inner.bottom <= outer.bottom
	)


def _split_free_rect(free: Rect, used: Rect) -> list[Rect]:
	if not _intersects(free, used):
		return [free]

	pieces: list[Rect] = []

	if used.x > free.x:
		pieces.append(Rect(free.x, free.y, used.x - free.x, free.h))

	if used.right < free.right:
		pieces.append(Rect(used.right, free.y, free.right - used.right, free.h))

	if used.y > free.y:
		pieces.append(Rect(free.x, free.y, free.w, used.y - free.y))

	if used.bottom < free.bottom:
		pieces.append(Rect(free.x, used.bottom, free.w, free.bottom - used.bottom))

	# Filter degenerate pieces
	return [r for r in pieces if r.w > 0 and r.h > 0]


def _prune_free_rects(free_rects: list[Rect]) -> list[Rect]:
	pruned: list[Rect] = []
	for r in free_rects:
		contained = False
		for other in free_rects:
			if r is other:
				continue
			if _contains(other, r):
				contained = True
				break
		if not contained:
			pruned.append(r)
	return pruned


def _best_short_side_fit(free_rects: list[Rect], w: int, h: int) -> Rect | None:
	best: Rect | None = None
	best_short = math.inf
	best_long = math.inf

	for fr in free_rects:
		if w > fr.w or h > fr.h:
			continue
		leftover_h = fr.h - h
		leftover_w = fr.w - w
		short_side = min(leftover_w, leftover_h)
		long_side = max(leftover_w, leftover_h)
		if (short_side < best_short) or (short_side == best_short and long_side < best_long):
			best_short = short_side
			best_long = long_side
			best = Rect(fr.x, fr.y, w, h)
	return best


def pack_maxrects_bin(
	width: int, height: int, sizes: list[tuple[int, int]]
) -> list[tuple[int, int]] | None:
	"""Pack rectangles into a width x height bin.

	Returns (x,y) placements in the same order as `sizes`, or None if impossible.
	"""

	# Sort by area/side to get more stable packing, but return placements in original order.
	indexed = list(enumerate(sizes))
	indexed.sort(key=lambda it: (max(it[1]), it[1][0] * it[1][1]), reverse=True)

	free_rects: list[Rect] = [Rect(0, 0, width, height)]
	placements: dict[int, tuple[int, int]] = {}

	for idx, (w, h) in indexed:
		node = _best_short_side_fit(free_rects, w, h)
		if node is None:
			return None

		placements[idx] = (node.x, node.y)

		new_free: list[Rect] = []
		for fr in free_rects:
			new_free.extend(_split_free_rect(fr, node))
		free_rects = _prune_free_rects(new_free)

	return [placements[i] for i in range(len(sizes))]


def pack_maxrects_square(side: int, sizes: list[tuple[int, int]]) -> list[tuple[int, int]] | None:
	return pack_maxrects_bin(side, side, sizes)


def _common_suffix(strings: list[str]) -> str:
	if not strings:
		return ""
	rev = [s[::-1] for s in strings]
	prefix = os.path.commonprefix(rev)
	return prefix[::-1]


def _expand_args_to_files(args: Iterable[str]) -> list[Path]:
	files: list[Path] = []
	for a in args:
		# Support both shell-expanded arguments and literal globs.
		matches = glob.glob(a)
		if matches:
			files.extend(Path(m) for m in matches)
		else:
			files.append(Path(a))

	# Normalize, de-dup, and keep stable order
	seen: set[Path] = set()
	out: list[Path] = []
	for p in files:
		p = p.expanduser()
		if p in seen:
			continue
		seen.add(p)
		out.append(p)
	return out


def _minimal_square_side(sizes: list[tuple[int, int]]) -> int:
	if not sizes:
		return 0
	total_area = sum(w * h for w, h in sizes)
	max_w = max(w for w, _ in sizes)
	max_h = max(h for _, h in sizes)
	lower = max(max_w, max_h, int(math.ceil(math.sqrt(total_area))))
	return lower


def _parse_ratio(args: list[str]) -> tuple[int, int, bool, list[str]]:
	"""Parse optional trailing ratio args.

	Accepts:
	  - no ratio provided -> 1:1
	  - trailing two integers X Y -> X:Y
	  - special '-1 -1' -> 1:1
	Returns (x, y, ratio_explicitly_given, remaining_args).
	"""

	if len(args) >= 3:
		a, b = args[-2], args[-1]
		try:
			x = int(a)
			y = int(b)
		except ValueError:
			return 1, 1, False, args

		rest = args[:-2]
		if x == -1 and y == -1:
			return 1, 1, False, rest

		if x <= 0 or y <= 0:
			raise SystemExit(f"Invalid ratio {x}:{y}. Use positive integers (e.g. 16 9) or -1 -1.")
		return x, y, True, rest

	return 1, 1, False, args


def _ceil_div(a: int, b: int) -> int:
	return (a + b - 1) // b


def _minimal_scale_for_ratio(sizes: list[tuple[int, int]], x: int, y: int) -> int:
	if not sizes:
		return 0
	total_area = sum(w * h for w, h in sizes)
	max_w = max(w for w, _ in sizes)
	max_h = max(h for _, h in sizes)

	k_dim = max(_ceil_div(max_w, x), _ceil_div(max_h, y))
	k_area = int(math.ceil(math.sqrt(total_area / (x * y))))
	return max(1, k_dim, k_area)


def main(argv: list[str]) -> int:
	if len(argv) < 2:
		print(
			"Usage: ./auto_arange_and_combine.py <images...> [X Y]\n"
			"  - <images...> can be shell-expanded or globs\n"
			"  - [X Y] (optional) sets output aspect ratio X:Y (e.g. 16 9)\n"
			"  - omit [X Y] for default 1:1, or pass -1 -1 for 1:1\n"
			"Example: ./auto_arange_and_combine.py *combined.png\n"
			"Example: ./auto_arange_and_combine.py *combined.png 16 9",
			file=sys.stderr,
		)
		return 2

	ratio_x, ratio_y, ratio_given, image_args = _parse_ratio(argv[1:])
	paths = _expand_args_to_files(image_args)
	paths = [p for p in paths if p.is_file() and not (p.name.startswith("ALL_") or p.name.startswith("All_"))]
	if not paths:
		print("No input files found (or all were skipped because they start with 'ALL_' or 'All_').", file=sys.stderr)
		return 2

	basenames = [p.name for p in paths]
	suffix = _common_suffix(basenames)
	if not suffix:
		# Fallback: use extension of the first input
		suffix = Path(basenames[0]).suffix
	if suffix.startswith("_") or suffix.startswith("."):
		out_name = f"All{suffix}"
	else:
		out_name = f"All_{suffix}"
	
	out_path = Path(out_name)
	if ratio_given:
		out_path = out_path.with_name(f"{out_path.stem}_{ratio_x}_{ratio_y}{out_path.suffix}")

	images: list[Image.Image] = []
	sizes: list[tuple[int, int]] = []
	for p in paths:
		img = Image.open(p)
		img.load()
		img = img.convert("RGBA")
		images.append(img)
		sizes.append(img.size)

	k0 = _minimal_scale_for_ratio(sizes, ratio_x, ratio_y)
	k = k0
	width = ratio_x * k
	height = ratio_y * k

	placements = pack_maxrects_bin(width, height, sizes)
	while placements is None:
		k *= 2
		width = ratio_x * k
		height = ratio_y * k
		placements = pack_maxrects_bin(width, height, sizes)

	lo, hi = k0, k
	while lo < hi:
		mid = (lo + hi) // 2
		trial = pack_maxrects_bin(ratio_x * mid, ratio_y * mid, sizes)
		if trial is None:
			lo = mid + 1
		else:
			hi = mid
			placements = trial

	k = lo
	width = ratio_x * k
	height = ratio_y * k

	canvas = Image.new("RGBA", (width, height), (0, 0, 0, 255))
	for img, (x, y) in zip(images, placements, strict=True):
		canvas.paste(img, (x, y), img)

	ext = out_path.suffix.lower()
	save_kwargs: dict[str, object] = {}

	# Important: we never rescale any input image; we only paste them at native pixel size.
	# Here we pick save settings that avoid quality loss where possible.
	if ext in {".png"}:
		# PNG is lossless; `compress_level` only affects file size and CPU.
		save_kwargs.update({"compress_level": 0, "optimize": False})
	elif ext in {".jpg", ".jpeg"}:
		# JPEG is inherently lossy. Use settings that minimize additional loss.
		canvas = canvas.convert("RGB")
		save_kwargs.update({"quality": 100, "subsampling": 0, "optimize": False})
	elif ext in {".webp"}:
		# Prefer lossless output for WebP if requested by extension.
		save_kwargs.update({"lossless": True, "quality": 100})

	canvas.save(out_path, **save_kwargs)
	print(
		f"Wrote {out_path} ({width}x{height}, ratio {ratio_x}:{ratio_y}) from {len(images)} images"
	)
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))


#!/usr/bin/env python
"""Pack observed VRI footprints and horizontal yellow names on a black canvas.

Examples:
    python auto_arrange_and_combine_named.py
    python auto_arrange_and_combine_named.py extended
    python auto_arrange_and_combine_named.py '*_observed_VRI.png' 16 9
    python auto_arrange_and_combine_named.py --rotation-step 90 --no-pa-alignment

Only exactly black pixels are free space. Rotated image rectangles can overlap,
but visible pixels and padded label rectangles cannot. The search uses a
conservative occupancy grid and validates the final layout at full resolution.
This is a heuristic layout, not an OR-Tools optimality proof. Input images are
never resized; non-quarter-turn rotations use bicubic interpolation.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont
    from scipy.ndimage import maximum_filter
    from scipy.signal import fftconvolve
except ModuleNotFoundError as exc:
    raise SystemExit(f"Missing dependency: {exc.name}. Install numpy, pillow and scipy in {sys.executable}.") from exc


EXTENDED_GALAXY_IDS = (
    "IC3392",
    "NGC4064",
    "NGC4192",
    "NGC4293",
    "NGC4298",
    "NGC4330",
    "NGC4383",
    "NGC4396",
    "NGC4419",
    "NGC4457",
    "NGC4501",
    "NGC4522",
    "NGC4694",
    "NGC4698",
)


@dataclass
class Source:
    path: Path
    label: str
    image: Image.Image
    sha256: str


@dataclass
class Variant:
    angle: float
    rotated_crop_box: tuple[int, int, int, int]
    image_offset: tuple[int, int]
    label_box: tuple[int, int, int, int]
    size: tuple[int, int]
    mask: np.ndarray


@dataclass
class Placement:
    variant: Variant
    x: int
    y: int


def rotate_image(image: Image.Image, angle: float) -> Image.Image:
    """Use exact pixel permutations for quarter turns, interpolation otherwise."""
    angle %= 360
    if angle == 0:
        return image.copy()
    if angle in (90, 180, 270):
        method = {90: Image.Transpose.ROTATE_90, 180: Image.Transpose.ROTATE_180,
                  270: Image.Transpose.ROTATE_270}[angle]
        return image.transpose(method)
    return image.rotate(angle, resample=Image.Resampling.BICUBIC, expand=True, fillcolor="black")


def general_rotation_angles(step: int) -> list[int]:
    """Return rotations around the original orientation, limited to +/-90 deg."""
    angles = {0, -90, 90}
    for angle in range(step, 90, step):
        angles.update((-angle, angle))
    return sorted(angles)


def _least_equivalent_axis_rotation(angle: float) -> float:
    """Map a major-axis rotation to its least-magnitude 180-deg equivalent."""
    return ((angle + 90) % 180) - 90


def coarsen_mask(mask: np.ndarray, grid: int) -> np.ndarray:
    """Every cell containing even one occupied pixel stays occupied."""
    height, width = mask.shape
    h, w = math.ceil(height / grid), math.ceil(width / grid)
    padded = np.zeros((h * grid, w * grid), dtype=bool)
    padded[:height, :width] = mask
    return padded.reshape(h, grid, w, grid).any(axis=(1, 3))


def load_font(size: int, path: Path | None) -> ImageFont.FreeTypeFont:
    if path is not None:
        return ImageFont.truetype(str(path), size)
    for candidate in ("DejaVuSans.ttf", "/System/Library/Fonts/Supplemental/Arial.ttf",
                      "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    # Modern Pillow ships a scalable default font.
    try:
        return ImageFont.load_default(size=size)
    except TypeError as exc:
        raise ValueError("Supply --font PATH or use Pillow >= 10.1 for scalable labels.") from exc


def load_source(path: Path) -> Source:
    with Image.open(path) as original:
        rgba = original.convert("RGBA")
        background = Image.new("RGBA", rgba.size, (0, 0, 0, 255))
        image = Image.alpha_composite(background, rgba).convert("RGB")
    if not np.asarray(image).any():
        raise ValueError(f"Input is entirely black: {path}")
    return Source(path.resolve(), path.name[:-len("_observed_VRI.png")], image,
                  hashlib.sha256(path.read_bytes()).hexdigest())


def read_pa_table(path: Path) -> dict[str, dict[str, float]]:
    """Read the local tab-separated Brown table: its fifth column is i, not a/b."""
    catalogue = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        columns = line.split("\t")
        if len(columns) != 6:
            continue
        match = re.match(r"^(IC|NGC|VCC)\s*(\d+)(?:\s|$)", columns[0])
        if match:
            inclination, pa = float(columns[4]), float(columns[5])
            if not (0 <= inclination <= 90 and math.isfinite(pa)):
                raise ValueError(f"Invalid inclination/PA in {path}: {columns[0]}")
            catalogue["".join(match.groups())] = {"pa_deg": pa % 360, "inclination_deg": inclination}
    if not catalogue:
        raise ValueError(f"No galaxy PA rows found in {path}")
    return catalogue


def pa_alignment(source: Source, catalogue: dict[str, dict[str, float]]) -> tuple[list[float], dict]:
    entry = catalogue.get(source.label)
    if entry is None:
        return [], {"status": "no_unique_catalogue_entry"}
    metadata = dict(entry)
    try:
        from astropy import units as u
        from astropy.io import fits
        from astropy.wcs import WCS
    except ModuleNotFoundError:
        return [], dict(metadata, status="astropy_unavailable")
    candidates = []
    for stem in (f"{source.label}_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits",
                 f"{source.label}_PHANGS_DATACUBE_native_VRI.fits"):
        for suffix in ("", ".gz"):
            path = source.path.with_name(stem + suffix)
            if path.is_file():
                candidates.append(path)
    if not candidates:
        return [], dict(metadata, status="matching_vri_fits_missing")
    # Header-only access: the PNG is the image input, FITS supplies orientation.
    with fits.open(candidates[0], memmap=True) as hdus:
        if "V_FLUX" not in hdus:
            return [], dict(metadata, status="V_FLUX_header_missing")
        header = hdus["V_FLUX"].header
        if (header.get("NAXIS1"), header.get("NAXIS2")) != source.image.size:
            return [], dict(metadata, status="fits_png_dimensions_differ")
        wcs = WCS(header).celestial
    if not wcs.has_celestial:
        return [], dict(metadata, status="celestial_wcs_missing")
    x, y = (source.image.width - 1) / 2, (source.image.height - 1) / 2
    centre = wcs.pixel_to_world(x, y)
    tip = centre.directional_offset_by(entry["pa_deg"] * u.deg, 1 * u.arcsec)
    px, py = wcs.world_to_pixel(tip)
    # The observed renderer uses np.flipud before PNG export. Hence FITS +y
    # maps to screen-up, and atan2(dy, dx) is the visual CCW angle from right.
    theta = math.degrees(math.atan2(float(py) - y, float(px) - x))
    # A major axis is unchanged by a 180-degree turn. Keep only the least
    # rotation needed for horizontal and vertical alignment.
    angles = sorted({round(_least_equivalent_axis_rotation(-theta + quarter * 90), 6)
                     for quarter in range(2)})
    metadata.update(status="available", wcs_source=str(candidates[0]),
                    major_axis_angle_from_right_deg=theta, alignment_angles=angles)
    return angles, metadata


def make_variant(source: Source, angle: float, font: ImageFont.FreeTypeFont,
                 grid: int, gap: int) -> Variant:
    rotated = rotate_image(source.image, angle)
    support = np.asarray(rotated).max(axis=2) > 0
    ys, xs = np.nonzero(support)
    crop = (int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1)
    support = support[crop[1]:crop[3], crop[0]:crop[2]]
    image_h, image_w = support.shape
    left, top, right, bottom = font.getbbox(source.label)
    text_pad = max(3, font.size // 8)
    label_w, label_h = right - left + 2 * text_pad, bottom - top + 2 * text_pad
    width = max(image_w, label_w) + 2 * gap
    height = image_h + label_h + 3 * gap
    image_x, image_y = (width - image_w) // 2, gap
    mask = np.zeros((height, width), dtype=bool)
    mask[image_y:image_y + image_h, image_x:image_x + image_w] = support

    # Find a fully black rectangle near the galaxy, including black corners.
    # A footer is always available when the observed footprint fills its frame.
    guarded = maximum_filter(mask, size=2 * gap + 1, mode="constant") if gap else mask
    integral = np.pad(guarded.astype(np.int32), ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    candidates_x = np.unique(np.append(np.arange(gap, width - label_w - gap + 1, max(1, grid)),
                                       width - label_w - gap))
    candidates_y = np.unique(np.append(np.arange(gap, height - label_h - gap + 1, max(1, grid)),
                                       height - label_h - gap))
    yy, xx = np.meshgrid(candidates_y, candidates_x, indexing="ij")
    sums = (integral[yy + label_h, xx + label_w] - integral[yy, xx + label_w]
            - integral[yy + label_h, xx] + integral[yy, xx])
    # Prefer a nearby black corner, with a small preference for below the centre.
    distance = ((xx + label_w / 2 - (image_x + image_w / 2)) ** 2
                + (yy + label_h / 2 - (image_y + image_h * 0.65)) ** 2)
    distance[sums != 0] = np.inf
    if not np.isfinite(distance).any():
        raise ValueError(f"Cannot reserve a label for {source.label}")
    best = np.unravel_index(np.argmin(distance), distance.shape)
    label_x, label_y = int(xx[best]), int(yy[best])
    mask[label_y:label_y + label_h, label_x:label_x + label_w] = True
    if gap:
        mask = maximum_filter(mask, size=2 * gap + 1, mode="constant")
    ys, xs = np.nonzero(mask)
    x0, y0, x1, y1 = int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1
    mask = mask[y0:y1, x0:x1]
    coarse = coarsen_mask(mask, grid)
    return Variant(angle, crop, (image_x - x0, image_y - y0),
                   (label_x - x0, label_y - y0, label_w, label_h),
                   (coarse.shape[1] * grid, coarse.shape[0] * grid), coarse)


def pack_at_size(variants: list[list[Variant]], order: list[int], width: int,
                 height: int, grid: int) -> list[Placement] | None:
    occupied = np.zeros((height // grid, width // grid), dtype=bool)
    placements: list[Placement | None] = [None] * len(variants)
    for index in order:
        best = None
        for variant in variants[index]:
            h, w = variant.mask.shape
            if h > occupied.shape[0] or w > occupied.shape[1]:
                continue
            if occupied.any():
                # Integer-valued correlations: < 0.5 absorbs floating-point FFT error.
                collisions = fftconvolve(occupied.astype(np.float32),
                                         variant.mask[::-1, ::-1].astype(np.float32), mode="valid")
                valid = collisions < 0.5
                rows = np.flatnonzero(valid.any(axis=1))
                if not len(rows):
                    continue
                y = int(rows[0])
                x = int(np.flatnonzero(valid[y])[0])
            else:
                x = y = 0
            # At a tested canvas size, keep each galaxy as close as possible to
            # its original PNG orientation. Rotate further only when the lower-
            # angle variants cannot fit without a protected-pixel collision.
            score = (abs(variant.angle), y + h, x, y, w, variant.angle)
            if best is None or score < best[0]:
                best = (score, variant, x, y)
        if best is None:
            return None
        _, variant, x, y = best
        h, w = variant.mask.shape
        region = occupied[y:y + h, x:x + w]
        if np.any(region & variant.mask):
            raise ValueError("Packing collision: floating-point occupancy search failed validation")
        region |= variant.mask
        placements[index] = Placement(variant, x * grid, y * grid)
    return [placement for placement in placements if placement is not None]


def find_layout(variants: list[list[Variant]], ratio: tuple[int, int], grid: int,
                attempts: int, seed: int) -> tuple[int, int, list[Placement]]:
    rx, ry = ratio
    areas = [min(int(v.mask.sum()) for v in choices) * grid ** 2 for choices in variants]
    min_k = max(1, math.ceil(math.sqrt(sum(areas) / (rx * ry))))
    min_k = max(min_k, max(min(max(math.ceil(v.size[0] / rx), math.ceil(v.size[1] / ry))
                              for v in choices) for choices in variants))
    order = sorted(range(len(variants)), key=lambda i: (-areas[i], i))
    rng = random.Random(seed)
    best = None
    best_k = math.ceil(min_k * 1.25)
    # Establish a feasible layout before trying smaller canvases. Growth must
    # eventually fit even a shelf of all panels; no solver or stored proof used.
    while best is None:
        best = pack_at_size(variants, order, rx * best_k, ry * best_k, grid)
        print(f"Initial search: {rx * best_k}x{ry * best_k}: {'fits' if best else 'growing'}", flush=True)
        if best is None:
            best_k = max(best_k + grid, math.ceil(best_k * 1.18))
    for attempt in range(attempts):
        if attempt:
            # Largest pieces first, with reproducible order changes to escape
            # one greedy packing. Each order gets its own search interval.
            priority = [area * rng.uniform(0.7, 1.3) for area in areas]
            order = sorted(range(len(variants)), key=lambda i: (-priority[i], i))
        low, high = min_k - 1, best_k
        # Heuristic failures are not infeasibility proofs. Binary search merely
        # samples useful canvas sizes; another order may fit a smaller canvas.
        for _ in range(8):
            if high - low <= 1:
                break
            candidate = (low + high) // 2
            result = pack_at_size(variants, order, rx * candidate, ry * candidate, grid)
            if result is None:
                low = candidate
            else:
                high = candidate
                if candidate < best_k:
                    best_k, best = candidate, result
        print(f"Search {attempt + 1}/{attempts}: best canvas {rx * best_k}x{ry * best_k}", flush=True)
    width, height = rx * best_k, ry * best_k
    used_right = max(p.x + p.variant.size[0] for p in best)
    used_bottom = max(p.y + p.variant.size[1] for p in best)
    dx, dy = (width - used_right) // 2, (height - used_bottom) // 2
    return width, height, [Placement(p.variant, p.x + dx, p.y + dy) for p in best]


def render_variant(source: Source, variant: Variant, font: ImageFont.FreeTypeFont) -> tuple[Image.Image, np.ndarray]:
    image = rotate_image(source.image, variant.angle).crop(variant.rotated_crop_box)
    sprite = Image.new("RGB", variant.size, "black")
    sprite.paste(image, variant.image_offset)
    reserved = np.asarray(sprite).max(axis=2) > 0
    x, y, w, h = variant.label_box
    if reserved[y:y + h, x:x + w].any():
        raise ValueError(f"Label overlaps visible pixels for {source.label}")
    reserved[y:y + h, x:x + w] = True
    left, top, right, bottom = font.getbbox(source.label)
    ImageDraw.Draw(sprite).text((x + (w - (right - left)) // 2 - left,
                                y + (h - (bottom - top)) // 2 - top),
                               source.label, font=font, fill=(255, 255, 0))
    return sprite, reserved


def paste_checked(canvas: Image.Image, occupied: np.ndarray, sprite: Image.Image,
                  reserved: np.ndarray, position: tuple[int, int]) -> None:
    x, y = position
    w, h = sprite.size
    if x < 0 or y < 0 or x + w > canvas.width or y + h > canvas.height:
        raise ValueError("Placement extends outside the canvas")
    region = occupied[y:y + h, x:x + w]
    if np.any(region & reserved):
        raise ValueError("Visible pixels or label rectangles overlap")
    canvas.paste(sprite, (x, y), Image.fromarray(reserved.astype(np.uint8) * 255))
    region |= reserved


def save_layout(sources: list[Source], placements: list[Placement], width: int, height: int,
                font: ImageFont.FreeTypeFont, output: Path, report_path: Path,
                settings: dict, alignments: list[dict]) -> dict:
    canvas = Image.new("RGB", (width, height), "black")
    occupied = np.zeros((height, width), dtype=bool)
    rows = []
    for source, placement, alignment in zip(sources, placements, alignments, strict=True):
        variant = placement.variant
        sprite, reserved = render_variant(source, variant, font)
        paste_checked(canvas, occupied, sprite, reserved, (placement.x, placement.y))
        ix, iy = variant.image_offset
        lx, ly, lw, lh = variant.label_box
        rows.append({"source": str(source.path), "source_sha256": source.sha256,
                     "label": source.label, "original_size": list(source.image.size),
                     "angle_degrees": variant.angle, "rotated_crop_box": list(variant.rotated_crop_box),
                     "sprite_position": [placement.x, placement.y], "sprite_size": list(variant.size),
                     "image_position": [placement.x + ix, placement.y + iy],
                     "label_box": [placement.x + lx, placement.y + ly, lw, lh],
                     "catalogue_alignment": alignment})
    overlaps = 0
    for i, a in enumerate(rows):
        ax, ay = a["image_position"]
        ac = a["rotated_crop_box"]
        for b in rows[:i]:
            bx, by = b["image_position"]
            bc = b["rotated_crop_box"]
            overlaps += int(ax < bx + bc[2] - bc[0] and bx < ax + ac[2] - ac[0]
                            and ay < by + bc[3] - bc[1] and by < ay + ac[3] - ac[1])
    report = {"status": "HEURISTIC_ONLY", "canvas": [width, height], "settings": settings,
              "note": "No optimality proof. Only exact black background is free; no brightness threshold or resizing. "
                      "Angles not divisible by 90 use bicubic interpolation. Positions and boxes are in output pixels; "
                      "rotated_crop_box is left, top, right, bottom in the expanded rotated source.",
              "validation": {"overlapping_reserved_pixels": 0, "image_rectangle_overlap_pairs": overlaps,
                             "reserved_fraction": float(occupied.mean()),
                             "all_placements_inside_canvas": True},
              "placements": rows}
    # Complete full-resolution validation before writing either deliverable.
    canvas.save(output, compress_level=6)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("items", nargs="*", help="Use no items for all galaxies, 'extended' for the fixed 14, or provide observed PNG files/globs; optionally followed by ratio X Y [label]")
    parser.add_argument("--output", type=Path, help="Output PNG (default All_observed_VRI_named[_X_Y].png)")
    parser.add_argument("--report-file", type=Path, help="Layout JSON (default output stem + .layout.json)")
    parser.add_argument("--rotation-step", type=int, default=15, help="General angle step within -90..+90 (default 15); use 90 with --no-pa-alignment for exact pixels")
    parser.add_argument("--no-rotate", action="store_true", help="Keep the original orientations")
    parser.add_argument("--pa-table", type=Path, help="Brown PA table (default Brown2021Table1.txt beside this script, if present)")
    parser.add_argument("--no-pa-alignment", action="store_true", help="Do not add PA-based horizontal/vertical candidates")
    parser.add_argument("--grid-size", type=int, default=12, help="Conservative packing grid in pixels (default 12; smaller is tighter/slower)")
    parser.add_argument("--gap", type=int, default=6, help="Clearance radius in pixels around footprints and labels (default 6)")
    parser.add_argument("--font-size", type=int, default=36, help="Horizontal yellow label size in pixels (default 36)")
    parser.add_argument("--font", type=Path, help="Optional TrueType/OpenType font path")
    parser.add_argument("--attempts", type=int, default=4, help="Number of deterministic packing orders (default 4)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible packing orders")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.rotation_step <= 0 or args.rotation_step > 90:
        parser.error("--rotation-step must be between 1 and 90 degrees")
    if min(args.grid_size, args.font_size, args.attempts) < 1 or args.gap < 0:
        parser.error("grid size, font size and attempts must be positive; gap must be nonnegative")
    items = args.items.copy()
    if items and items[-1] == "label":
        items.pop()  # Accepted for familiarity; labels are always enabled here.
    ratio, explicit_ratio = (1, 1), False
    if len(items) >= 2:
        try:
            x, y = int(items[-2]), int(items[-1])
        except ValueError:
            pass
        else:
            del items[-2:]
            if (x, y) != (-1, -1):
                if x <= 0 or y <= 0:
                    parser.error("The aspect ratio must use two positive integers")
                ratio, explicit_ratio = (x, y), True
    paths = []
    selection_mode = "custom"
    if not items:
        selection_mode = "all"
        if not explicit_ratio:
            ratio = (16, 9)
        paths = sorted(Path.cwd().glob("*_observed_VRI.png"))
    elif items == ["extended"]:
        selection_mode = "extended"
        if not explicit_ratio:
            ratio = (16, 9)
        paths = [Path.cwd() / f"{galaxy_id}_observed_VRI.png"
                 for galaxy_id in EXTENDED_GALAXY_IDS]
        missing = [path.name for path in paths if not path.is_file()]
        if missing:
            parser.error(
                f"extended selection is missing {len(missing)} required observed_VRI file(s): "
                + ", ".join(missing)
            )
    else:
        for item in items:
            matches = sorted(glob.glob(str(Path(item).expanduser())))
            if not matches:
                parser.error(f"No files match: {item}")
            paths.extend(Path(match).resolve() for match in matches)
    paths = list(dict.fromkeys(path for path in paths if path.is_file()
                              and path.name.endswith("_observed_VRI.png")
                              and not path.name.lower().startswith("all_")))
    if not paths:
        parser.error("No individual *_observed_VRI.png inputs found; combined_VRI inputs are excluded")
    suffix = f"_{ratio[0]}_{ratio[1]}" if explicit_ratio else ""
    if selection_mode == "all":
        default_output = Path(f"All_observed_VRI_named__{ratio[0]}_{ratio[1]}.png")
    elif selection_mode == "extended":
        default_output = Path(f"14_observed_VRI_named_{ratio[0]}_{ratio[1]}.png")
    else:
        default_output = Path(f"All_observed_VRI_named{suffix}.png")
    output = (args.output or default_output).resolve()
    report_path = (args.report_file or output.with_suffix(".layout.json")).resolve()
    if output.suffix.lower() != ".png":
        parser.error("--output must end in .png for lossless output")
    if report_path.suffix.lower() != ".json":
        parser.error("--report-file must end in .json")
    if output.name.endswith("_combined_VRI.png") and not output.name.lower().startswith("all_"):
        parser.error("Output cannot replace an individual combined_VRI input product")
    if output == report_path or output in paths or report_path in paths:
        parser.error("Output/report paths must be distinct and cannot overwrite input images")
    if output.name.endswith("_observed_VRI.png") and not output.name.lower().startswith("all_"):
        parser.error("Output name would be rediscovered as a galaxy; use an All_ prefix or a different suffix")
    if not output.parent.is_dir() or not report_path.parent.is_dir():
        parser.error("Output/report parent directories must already exist")
    try:
        font = load_font(args.font_size, args.font)
        sources = [load_source(path) for path in paths]
        angles = [0] if args.no_rotate else general_rotation_angles(args.rotation_step)
        table_path = args.pa_table or Path(__file__).with_name("Brown2021Table1.txt")
        catalogue = {}
        if not args.no_pa_alignment and not args.no_rotate:
            if args.pa_table or table_path.is_file():
                catalogue = read_pa_table(table_path)
        print(f"Preparing {len(sources)} observed images: {len(angles)} general rotation(s), plus available PA alignments", flush=True)
        variants = []
        alignments = []
        for source in sources:
            extra, alignment = pa_alignment(source, catalogue) if catalogue else ([], {"status": "disabled_or_table_missing"})
            source_angles = sorted(set(angles + extra))
            variants.append([make_variant(source, angle, font, args.grid_size, args.gap) for angle in source_angles])
            alignments.append(alignment)
            print(f"Prepared {source.label}: {len(source_angles)} angles; PA {alignment['status']}", flush=True)
        divisor = math.gcd(*ratio)
        width, height, placements = find_layout(variants, (ratio[0] // divisor, ratio[1] // divisor),
                                                args.grid_size, args.attempts, args.seed)
        settings = {"ratio": list(ratio), "angles": angles,
                    "selection_mode": selection_mode,
                    "selected_galaxy_ids": [source.label for source in sources],
                    "rotation_limit_degrees": [-90, 90],
                    "rotation_preference": "least absolute angle that fits the tested canvas",
                    "grid_size": args.grid_size,
                    "gap": args.gap, "font_size": args.font_size, "font": str(getattr(font, "path", "Pillow default")),
                    "label_color": [255, 255, 0], "attempts": args.attempts, "seed": args.seed,
                    "pa_table": str(table_path.resolve()) if catalogue else None,
                    "pa_table_sha256": hashlib.sha256(table_path.read_bytes()).hexdigest() if catalogue else None}
        report = save_layout(sources, placements, width, height, font, output, report_path, settings, alignments)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"Saved {output} ({width}x{height})\nSaved {report_path}\n"
          f"Validated: no visible-pixel/label overlap; {report['validation']['image_rectangle_overlap_pairs']} "
          "image-rectangle pairs share black space. HEURISTIC_ONLY.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

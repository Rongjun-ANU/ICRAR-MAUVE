from __future__ import annotations

from pathlib import Path


def fits_stem(path: Path) -> str:
	name = path.name
	lower_name = name.lower()
	if lower_name.endswith(".fits.gz"):
		return name[: -len(".fits.gz")]
	if lower_name.endswith(".fits"):
		return name[: -len(".fits")]
	return path.stem


def fits_pattern_variants(pattern: str) -> list[str]:
	lower_pattern = pattern.lower()
	if lower_pattern.endswith(".fits.gz"):
		return [pattern, pattern[:-3]]
	if lower_pattern.endswith(".fits"):
		return [pattern, f"{pattern}.gz"]
	return [pattern]


def expand_fits_glob(input_dir: Path, pattern: str) -> list[Path]:
	matches: list[Path] = []
	for expanded_pattern in fits_pattern_variants(pattern):
		matches.extend(sorted(input_dir.glob(expanded_pattern)))

	seen: set[Path] = set()
	unique: list[Path] = []
	for path in matches:
		resolved = path.resolve()
		if resolved in seen:
			continue
		seen.add(resolved)
		unique.append(path)
	return unique


def strip_optional_gzip_suffix(path: Path) -> str:
	name = path.name
	if name.lower().endswith(".fits.gz"):
		return name[:-3]
	return name

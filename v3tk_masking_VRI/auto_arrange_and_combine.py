#!/usr/bin/env python

from __future__ import annotations

import argparse
import glob
import math
import os
import random
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


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


@dataclass(frozen=True)
class HeuristicResult:
	k: int
	width: int
	height: int
	placements: list[tuple[int, int]]
	method: str


@dataclass(frozen=True)
class OptimizeResult:
	k: int
	width: int
	height: int
	placements: list[tuple[int, int]]
	status: str
	objective_bound: float | None
	wall_time: float | None
	best_objective: float | None


@dataclass(frozen=True)
class StoredLayout:
	proof_path: Path
	output_path: Path | None
	status: str
	ratio_x: int
	ratio_y: int
	k: int
	width: int
	height: int
	paths: list[Path]
	sizes: list[tuple[int, int]]
	placements: list[tuple[int, int]]


def _elapsed(start_time: float) -> str:
	return f"{time.monotonic() - start_time:8.1f}s"


def log_progress(start_time: float, message: str) -> None:
	print(f"[{_elapsed(start_time)}] {message}", flush=True)


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


def _placement_score(mode: str, free_rect: Rect, w: int, h: int, width: int, height: int) -> tuple[int, ...] | None:
	if w > free_rect.w or h > free_rect.h:
		return None

	leftover_w = free_rect.w - w
	leftover_h = free_rect.h - h
	short_side = min(leftover_w, leftover_h)
	long_side = max(leftover_w, leftover_h)
	area_fit = free_rect.w * free_rect.h - w * h

	if mode == "bssf":
		return (short_side, long_side, free_rect.y, free_rect.x)
	if mode == "baf":
		return (area_fit, short_side, long_side, free_rect.y, free_rect.x)
	if mode == "bl":
		return (free_rect.y, free_rect.x, short_side, long_side)
	if mode == "tight":
		return (long_side, short_side, area_fit, free_rect.y, free_rect.x)
	if mode == "contact":
		contact = 0
		if free_rect.x == 0:
			contact += h
		if free_rect.y == 0:
			contact += w
		if free_rect.x + w == width:
			contact += h
		if free_rect.y + h == height:
			contact += w
		return (-contact, short_side, long_side, free_rect.y, free_rect.x)

	raise ValueError(f"Unknown packing score mode: {mode}")


def _choose_free_rect(free_rects: list[Rect], w: int, h: int, width: int, height: int, mode: str) -> Rect | None:
	best: Rect | None = None
	best_score: tuple[int, ...] | None = None

	for fr in free_rects:
		score = _placement_score(mode, fr, w, h, width, height)
		if score is None:
			continue
		if best_score is None or score < best_score:
			best_score = score
			best = Rect(fr.x, fr.y, w, h)
	return best


def pack_maxrects_bin(
	width: int,
	height: int,
	sizes: list[tuple[int, int]],
	order: list[int] | None = None,
	score_mode: str = "bssf",
) -> list[tuple[int, int]] | None:
	"""Pack rectangles into a width x height bin.

	Returns (x,y) placements in the same order as `sizes`, or None if impossible.
	"""

	if order is None:
		# Preserve the original script's default ordering for callers that do not
		# request the broader heuristic search.
		indexed = list(enumerate(sizes))
		indexed.sort(key=lambda it: (max(it[1]), it[1][0] * it[1][1]), reverse=True)
	else:
		indexed = [(i, sizes[i]) for i in order]

	free_rects: list[Rect] = [Rect(0, 0, width, height)]
	placements: dict[int, tuple[int, int]] = {}

	for idx, (w, h) in indexed:
		node = _choose_free_rect(free_rects, w, h, width, height, score_mode)
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


def _heuristic_orders(sizes: list[tuple[int, int]]) -> list[tuple[str, list[int]]]:
	orders: list[tuple[str, list[int]]] = []
	seen: set[tuple[int, ...]] = set()

	def add(name: str, order: list[int]) -> None:
		key = tuple(order)
		if key in seen:
			return
		seen.add(key)
		orders.append((name, order))

	indices = list(range(len(sizes)))
	add("input", indices[:])

	def sorted_order(name: str, key_fn: Callable[[int], float | int | tuple[float | int, ...]], reverse: bool) -> None:
		add(name, sorted(indices, key=key_fn, reverse=reverse))

	sorted_order("maxside_desc_area_desc", lambda i: (max(sizes[i]), sizes[i][0] * sizes[i][1]), True)
	sorted_order("area_desc", lambda i: sizes[i][0] * sizes[i][1], True)
	sorted_order("area_asc", lambda i: sizes[i][0] * sizes[i][1], False)
	sorted_order("width_desc", lambda i: sizes[i][0], True)
	sorted_order("height_desc", lambda i: sizes[i][1], True)
	sorted_order("perimeter_desc", lambda i: sizes[i][0] + sizes[i][1], True)
	sorted_order("aspect_wide_first", lambda i: sizes[i][0] / sizes[i][1], True)
	sorted_order("aspect_tall_first", lambda i: sizes[i][1] / sizes[i][0], True)
	sorted_order("squares_first", lambda i: -abs(sizes[i][0] - sizes[i][1]), True)

	area_order = sorted(indices, key=lambda i: sizes[i][0] * sizes[i][1], reverse=True)
	for seed in range(32):
		rng = random.Random(seed)
		order = area_order[:]
		for start in range(0, len(order), 5):
			chunk = order[start : start + 5]
			rng.shuffle(chunk)
			order[start : start + 5] = chunk
		add(f"area_chunk_shuffle_{seed}", order)

	for seed in range(16):
		rng = random.Random(1000 + seed)
		order = indices[:]
		rng.shuffle(order)
		add(f"shuffle_{seed}", order)

	return orders


def _try_heuristic_at_k(
	sizes: list[tuple[int, int]], ratio_x: int, ratio_y: int, k: int
) -> tuple[list[tuple[int, int]], str] | None:
	width = ratio_x * k
	height = ratio_y * k
	score_modes = ["bssf", "baf", "bl", "tight", "contact"]
	for order_name, order in _heuristic_orders(sizes):
		for score_mode in score_modes:
			placements = pack_maxrects_bin(width, height, sizes, order=order, score_mode=score_mode)
			if placements is not None:
				return placements, f"{order_name}+{score_mode}"
	return None


def find_heuristic_solution(sizes: list[tuple[int, int]], ratio_x: int, ratio_y: int) -> HeuristicResult:
	lower = _minimal_scale_for_ratio(sizes, ratio_x, ratio_y)
	k = max(1, lower)

	found = _try_heuristic_at_k(sizes, ratio_x, ratio_y, k)
	while found is None:
		k *= 2
		found = _try_heuristic_at_k(sizes, ratio_x, ratio_y, k)

	best_k = k
	best_placements, best_method = found
	for trial_k in range(lower, k):
		trial = _try_heuristic_at_k(sizes, ratio_x, ratio_y, trial_k)
		if trial is None:
			continue
		best_k = trial_k
		best_placements, best_method = trial
		break

	return HeuristicResult(
		k=best_k,
		width=ratio_x * best_k,
		height=ratio_y * best_k,
		placements=best_placements,
		method=best_method,
	)


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
		matches = sorted(glob.glob(a))
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


def prepare_ortools_runtime(argv: list[str]) -> None:
	"""On macOS conda envs, prefer OR-Tools' bundled dylibs before importing CP-SAT."""
	if sys.platform != "darwin" or os.environ.get("AUTO_ARANGE_ORTOOLS_DYLD") == "1":
		return
	try:
		import ortools
	except ModuleNotFoundError:
		return

	libs_dir = Path(ortools.__file__).resolve().parent / ".libs"
	if not libs_dir.is_dir():
		return

	current = [p for p in os.environ.get("DYLD_LIBRARY_PATH", "").split(os.pathsep) if p]
	if str(libs_dir) in current:
		return

	env = os.environ.copy()
	env["DYLD_LIBRARY_PATH"] = os.pathsep.join([str(libs_dir), *current])
	env["AUTO_ARANGE_ORTOOLS_DYLD"] = "1"
	os.execve(sys.executable, [sys.executable, *argv], env)


def _status_name(cp_model: object, solver: object, status: int) -> str:
	if hasattr(solver, "StatusName"):
		return str(solver.StatusName(status))
	status_names = {
		getattr(cp_model, "OPTIMAL"): "OPTIMAL",
		getattr(cp_model, "FEASIBLE"): "FEASIBLE",
		getattr(cp_model, "INFEASIBLE"): "INFEASIBLE",
		getattr(cp_model, "MODEL_INVALID"): "MODEL_INVALID",
		getattr(cp_model, "UNKNOWN"): "UNKNOWN",
	}
	return status_names.get(status, f"STATUS_{status}")


def _call_model_method(model: object, snake_name: str, camel_name: str, *args: object) -> object:
	if hasattr(model, snake_name):
		return getattr(model, snake_name)(*args)
	return getattr(model, camel_name)(*args)


def _configure_solver(
	solver: object,
	cp_model: object,
	time_limit: float,
	progress: Callable[[str], None] | None,
	verbose_solver: bool,
) -> None:
	solver.parameters.max_time_in_seconds = float(time_limit)
	solver.parameters.num_search_workers = max(1, min(os.cpu_count() or 1, 8))
	solver.parameters.random_seed = 0
	if progress is not None and verbose_solver:
		solver.parameters.log_search_progress = True
		if hasattr(solver.parameters, "log_to_stdout"):
			solver.parameters.log_to_stdout = False
		if hasattr(solver, "log_callback"):
			def _solver_log_callback(msg: str) -> None:
				text = msg.strip()
				if text:
					progress(f"CP-SAT: {text}")

			solver.log_callback = _solver_log_callback


def _add_redundant_cumulative_constraints(
	model: object,
	x_intervals: list[object],
	y_intervals: list[object],
	sizes: list[tuple[int, int]],
	width: int,
	height: int,
) -> None:
	if not (hasattr(model, "add_cumulative") or hasattr(model, "AddCumulative")):
		return
	heights = [h for _, h in sizes]
	widths = [w for w, _ in sizes]
	_call_model_method(model, "add_cumulative", "AddCumulative", x_intervals, heights, height)
	_call_model_method(model, "add_cumulative", "AddCumulative", y_intervals, widths, width)


def _add_decision_strategy(model: object, cp_model: object, vars_to_branch: list[object]) -> None:
	if not vars_to_branch:
		return
	if not (hasattr(model, "add_decision_strategy") or hasattr(model, "AddDecisionStrategy")):
		return
	if not hasattr(cp_model, "CHOOSE_MIN_DOMAIN_SIZE") or not hasattr(cp_model, "SELECT_MIN_VALUE"):
		return
	_call_model_method(
		model,
		"add_decision_strategy",
		"AddDecisionStrategy",
		vars_to_branch,
		cp_model.CHOOSE_MIN_DOMAIN_SIZE,
		cp_model.SELECT_MIN_VALUE,
	)


def solve_with_ortools(
	sizes: list[tuple[int, int]],
	ratio_x: int,
	ratio_y: int,
	upper: HeuristicResult,
	time_limit: float,
	progress: Callable[[str], None] | None = None,
	verbose_solver: bool = False,
) -> OptimizeResult:
	try:
		from ortools.sat.python import cp_model
	except ModuleNotFoundError as e:  # pragma: no cover
		raise SystemExit(
			"Missing dependency: ortools.\n"
			f"Interpreter: {sys.executable}\n"
			"Certified mode requires OR-Tools CP-SAT. Install it into this interpreter with:\n"
			f"  {sys.executable} -m pip install ortools\n"
			"Use --fast for heuristic-only output without a mathematical optimality proof."
		) from e
	except ImportError as e:  # pragma: no cover
		raise SystemExit(
			"Could not import OR-Tools CP-SAT.\n"
			f"Interpreter: {sys.executable}\n"
			f"Import error: {e}\n"
			"Certified mode requires a working OR-Tools CP-SAT install. Use --fast for heuristic-only "
			"output, or reinstall OR-Tools in this environment."
		) from e

	lower = _minimal_scale_for_ratio(sizes, ratio_x, ratio_y)
	model = cp_model.CpModel()

	new_int_var = lambda lb, ub, name: _call_model_method(model, "new_int_var", "NewIntVar", lb, ub, name)
	new_fixed_interval = lambda start, size, name: _call_model_method(
		model, "new_fixed_size_interval_var", "NewFixedSizeIntervalVar", start, size, name
	)
	add_constraint = lambda expr: _call_model_method(model, "add", "Add", expr)

	k_var = new_int_var(lower, upper.k, "k")
	x_vars = []
	y_vars = []
	x_intervals = []
	y_intervals = []

	for i, (w, h) in enumerate(sizes):
		x_var = new_int_var(0, upper.width - w, f"x_{i}")
		y_var = new_int_var(0, upper.height - h, f"y_{i}")
		add_constraint(x_var + w <= ratio_x * k_var)
		add_constraint(y_var + h <= ratio_y * k_var)
		x_vars.append(x_var)
		y_vars.append(y_var)
		x_intervals.append(new_fixed_interval(x_var, w, f"x_interval_{i}"))
		y_intervals.append(new_fixed_interval(y_var, h, f"y_interval_{i}"))

	_call_model_method(model, "add_no_overlap_2d", "AddNoOverlap2D", x_intervals, y_intervals)
	_call_model_method(model, "minimize", "Minimize", k_var)
	_add_decision_strategy(model, cp_model, [k_var, *x_vars, *y_vars])

	if hasattr(model, "add_hint"):
		model.add_hint(k_var, upper.k)
		for var, (x, _) in zip(x_vars, upper.placements, strict=True):
			model.add_hint(var, x)
		for var, (_, y) in zip(y_vars, upper.placements, strict=True):
			model.add_hint(var, y)
	elif hasattr(model, "AddHint"):
		model.AddHint(k_var, upper.k)
		for var, (x, _) in zip(x_vars, upper.placements, strict=True):
			model.AddHint(var, x)
		for var, (_, y) in zip(y_vars, upper.placements, strict=True):
			model.AddHint(var, y)

	solver = cp_model.CpSolver()
	_configure_solver(solver, cp_model, time_limit, progress, verbose_solver)

	status = solver.Solve(model)
	status_name = _status_name(cp_model, solver, status)
	objective_bound = float(solver.BestObjectiveBound()) if hasattr(solver, "BestObjectiveBound") else None
	wall_time = float(solver.WallTime()) if hasattr(solver, "WallTime") else None
	best_objective = float(solver.ObjectiveValue()) if status in {cp_model.OPTIMAL, cp_model.FEASIBLE} else None

	if status_name == "FEASIBLE":
		k = int(solver.Value(k_var))
		placements = [(int(solver.Value(x_var)), int(solver.Value(y_var))) for x_var, y_var in zip(x_vars, y_vars, strict=True)]
		return OptimizeResult(
			k=k,
			width=ratio_x * k,
			height=ratio_y * k,
			placements=placements,
			status=status_name,
			objective_bound=objective_bound,
			wall_time=wall_time,
			best_objective=best_objective,
		)

	if status_name != "OPTIMAL":
		return OptimizeResult(
			k=upper.k,
			width=upper.width,
			height=upper.height,
			placements=upper.placements,
			status=status_name,
			objective_bound=objective_bound,
			wall_time=wall_time,
			best_objective=best_objective,
		)

	k = int(solver.Value(k_var))
	placements = [(int(solver.Value(x_var)), int(solver.Value(y_var))) for x_var, y_var in zip(x_vars, y_vars, strict=True)]
	return OptimizeResult(
		k=k,
		width=ratio_x * k,
		height=ratio_y * k,
		placements=placements,
		status=status_name,
		objective_bound=objective_bound,
		wall_time=wall_time,
		best_objective=best_objective,
	)


def _solve_fixed_k_candidate(
	sizes: list[tuple[int, int]],
	ratio_x: int,
	ratio_y: int,
	k: int,
	hint_placements: list[tuple[int, int]],
	time_limit: float,
	progress: Callable[[str], None] | None = None,
	verbose_solver: bool = False,
) -> tuple[str, list[tuple[int, int]] | None, float | None]:
	try:
		from ortools.sat.python import cp_model
	except ModuleNotFoundError as e:  # pragma: no cover
		raise SystemExit(
			"Missing dependency: ortools.\n"
			f"Interpreter: {sys.executable}\n"
			"Certified mode requires OR-Tools CP-SAT. Install it into this interpreter with:\n"
			f"  {sys.executable} -m pip install ortools\n"
			"Use --fast for heuristic-only output without a mathematical optimality proof."
		) from e
	except ImportError as e:  # pragma: no cover
		raise SystemExit(
			"Could not import OR-Tools CP-SAT.\n"
			f"Interpreter: {sys.executable}\n"
			f"Import error: {e}\n"
			"Certified mode requires a working OR-Tools CP-SAT install. Use --fast for heuristic-only "
			"output, or reinstall OR-Tools in this environment."
		) from e

	width = ratio_x * k
	height = ratio_y * k
	model = cp_model.CpModel()

	new_int_var = lambda lb, ub, name: _call_model_method(model, "new_int_var", "NewIntVar", lb, ub, name)
	new_fixed_interval = lambda start, size, name: _call_model_method(
		model, "new_fixed_size_interval_var", "NewFixedSizeIntervalVar", start, size, name
	)

	x_vars = []
	y_vars = []
	x_intervals = []
	y_intervals = []
	for i, (w, h) in enumerate(sizes):
		x_var = new_int_var(0, width - w, f"x_{i}")
		y_var = new_int_var(0, height - h, f"y_{i}")
		x_vars.append(x_var)
		y_vars.append(y_var)
		x_intervals.append(new_fixed_interval(x_var, w, f"x_interval_{i}"))
		y_intervals.append(new_fixed_interval(y_var, h, f"y_interval_{i}"))

	_call_model_method(model, "add_no_overlap_2d", "AddNoOverlap2D", x_intervals, y_intervals)
	_add_redundant_cumulative_constraints(model, x_intervals, y_intervals, sizes, width, height)
	_add_decision_strategy(model, cp_model, [*x_vars, *y_vars])

	if hasattr(model, "add_hint"):
		for var, (x_hint, _), (w, _) in zip(x_vars, hint_placements, sizes, strict=True):
			model.add_hint(var, max(0, min(int(x_hint), width - w)))
		for var, (_, y_hint), (_, h) in zip(y_vars, hint_placements, sizes, strict=True):
			model.add_hint(var, max(0, min(int(y_hint), height - h)))
	elif hasattr(model, "AddHint"):
		for var, (x_hint, _), (w, _) in zip(x_vars, hint_placements, sizes, strict=True):
			model.AddHint(var, max(0, min(int(x_hint), width - w)))
		for var, (_, y_hint), (_, h) in zip(y_vars, hint_placements, sizes, strict=True):
			model.AddHint(var, max(0, min(int(y_hint), height - h)))

	solver = cp_model.CpSolver()
	_configure_solver(solver, cp_model, time_limit, progress, verbose_solver)
	status = solver.Solve(model)
	status_name = _status_name(cp_model, solver, status)
	wall_time = float(solver.WallTime()) if hasattr(solver, "WallTime") else None

	if status_name in {"OPTIMAL", "FEASIBLE"}:
		placements = [(int(solver.Value(x_var)), int(solver.Value(y_var))) for x_var, y_var in zip(x_vars, y_vars, strict=True)]
		validate_placements(width, height, sizes, placements)
		return status_name, placements, wall_time
	return status_name, None, wall_time


def solve_with_fixed_k_search(
	sizes: list[tuple[int, int]],
	ratio_x: int,
	ratio_y: int,
	upper: HeuristicResult,
	time_limit: float,
	progress: Callable[[str], None] | None = None,
	verbose_solver: bool = False,
) -> OptimizeResult:
	lower = _minimal_scale_for_ratio(sizes, ratio_x, ratio_y)
	total_wall_time = 0.0
	best = OptimizeResult(
		k=upper.k,
		width=upper.width,
		height=upper.height,
		placements=upper.placements,
		status="FEASIBLE",
		objective_bound=float(lower),
		wall_time=total_wall_time,
		best_objective=float(upper.k),
	)

	if upper.k <= lower:
		if progress is not None:
			progress(f"Current feasible k={upper.k} equals the lower bound; optimality is immediate")
		return OptimizeResult(
			k=upper.k,
			width=upper.width,
			height=upper.height,
			placements=upper.placements,
			status="OPTIMAL",
			objective_bound=float(upper.k),
			wall_time=total_wall_time,
			best_objective=float(upper.k),
		)

	deadline = time.monotonic() + float(time_limit)
	candidate = upper.k - 1
	while candidate >= lower:
		remaining = deadline - time.monotonic()
		if remaining <= 0:
			break

		width = ratio_x * candidate
		height = ratio_y * candidate
		if progress is not None:
			progress(
				f"Fixed-k search: testing k={candidate} ({width}x{height}); "
				f"current best k={best.k}; remaining {remaining:.1f}s"
			)

		try:
			validate_placements(width, height, sizes, best.placements)
		except ValueError:
			pass
		else:
			if progress is not None:
				progress(f"Existing placement already fits k={candidate}; lowering current best without solving")
			best = OptimizeResult(
				k=candidate,
				width=width,
				height=height,
				placements=best.placements,
				status="FEASIBLE",
				objective_bound=float(lower),
				wall_time=total_wall_time,
				best_objective=float(candidate),
			)
			candidate = best.k - 1
			continue

		status_name, placements, wall_time = _solve_fixed_k_candidate(
			sizes,
			ratio_x,
			ratio_y,
			candidate,
			best.placements,
			remaining,
			progress=progress,
			verbose_solver=verbose_solver,
		)
		if wall_time is not None:
			total_wall_time += wall_time

		if status_name in {"OPTIMAL", "FEASIBLE"} and placements is not None:
			best = OptimizeResult(
				k=candidate,
				width=width,
				height=height,
				placements=placements,
				status="FEASIBLE",
				objective_bound=float(lower),
				wall_time=total_wall_time,
				best_objective=float(candidate),
			)
			if progress is not None:
				progress(f"Found feasible layout at k={candidate}; continuing search for smaller k")
			candidate = best.k - 1
			continue

		if status_name == "INFEASIBLE":
			if progress is not None:
				progress(f"Proved k={candidate} infeasible; current best k={best.k} is optimal")
			return OptimizeResult(
				k=best.k,
				width=best.width,
				height=best.height,
				placements=best.placements,
				status="OPTIMAL",
				objective_bound=float(best.k),
				wall_time=total_wall_time,
				best_objective=float(best.k),
			)

		if progress is not None:
			progress(f"Fixed-k search was inconclusive at k={candidate} with status {status_name}")
		break

	return OptimizeResult(
		k=best.k,
		width=best.width,
		height=best.height,
		placements=best.placements,
		status="FEASIBLE",
		objective_bound=float(lower),
		wall_time=total_wall_time,
		best_objective=float(best.k),
	)


def validate_placements(width: int, height: int, sizes: list[tuple[int, int]], placements: list[tuple[int, int]]) -> None:
	if len(placements) != len(sizes):
		raise ValueError(f"Expected {len(sizes)} placements, got {len(placements)}")

	rects: list[Rect] = []
	for idx, ((x, y), (w, h)) in enumerate(zip(placements, sizes, strict=True)):
		rect = Rect(int(x), int(y), w, h)
		if rect.x < 0 or rect.y < 0 or rect.right > width or rect.bottom > height:
			raise ValueError(
				f"Placement {idx} is outside canvas: ({rect.x},{rect.y},{rect.w},{rect.h}) "
				f"not inside {width}x{height}"
			)
		for prev_idx, prev in enumerate(rects):
			if _intersects(prev, rect):
				raise ValueError(f"Placements {prev_idx} and {idx} overlap")
		rects.append(rect)


def center_placements(
	width: int, height: int, sizes: list[tuple[int, int]], placements: list[tuple[int, int]]
) -> list[tuple[int, int]]:
	if not placements:
		return placements
	min_x = min(x for x, _ in placements)
	min_y = min(y for _, y in placements)
	max_x = max(x + w for (x, _), (w, _) in zip(placements, sizes, strict=True))
	max_y = max(y + h for (_, y), (_, h) in zip(placements, sizes, strict=True))
	offset_x = (width - (max_x - min_x)) // 2 - min_x
	offset_y = (height - (max_y - min_y)) // 2 - min_y
	centered = [(x + offset_x, y + offset_y) for x, y in placements]
	validate_placements(width, height, sizes, centered)
	return centered


def _default_proof_path(out_path: Path) -> Path:
	return out_path.with_suffix(".proof.txt")


def _report_mode(fast_mode: bool, status: str) -> str:
	if fast_mode:
		return "heuristic-only (--fast), not certified"
	if status == "OPTIMAL":
		return "OR-Tools CP-SAT certified"
	return "OR-Tools CP-SAT fallback, not certified"


def _paths_equal(a: Path, b: Path) -> bool:
	try:
		return a.resolve() == b.resolve()
	except OSError:
		return a.absolute() == b.absolute()


def _output_path_for_proof(proof_path: Path) -> Path | None:
	if not proof_path.name.endswith(".proof.txt"):
		return None
	base = proof_path.name[: -len(".proof.txt")]
	for suffix in (".png", ".jpg", ".jpeg", ".webp", ".tif", ".tiff"):
		candidate = proof_path.with_name(f"{base}{suffix}")
		if candidate.is_file():
			return candidate
	return None


def _parse_int(value: str | None) -> int | None:
	if value is None:
		return None
	try:
		return int(value)
	except ValueError:
		return None


def _read_layout_from_proof(proof_path: Path) -> StoredLayout | None:
	try:
		lines = proof_path.read_text(encoding="utf-8").splitlines()
	except OSError:
		return None

	values: dict[str, str] = {}
	for line in lines:
		if "\t" in line or ":" not in line:
			continue
		key, value = line.split(":", 1)
		values[key.strip()] = value.strip()

	ratio_text = values.get("ratio")
	if not ratio_text or ":" not in ratio_text:
		return None
	try:
		ratio_x, ratio_y = (int(part.strip()) for part in ratio_text.split(":", 1))
	except ValueError:
		return None

	k = _parse_int(values.get("final_k"))
	canvas_text = values.get("final_canvas")
	if k is None or not canvas_text or "x" not in canvas_text:
		return None
	try:
		width, height = (int(part.strip()) for part in canvas_text.split("x", 1))
	except ValueError:
		return None

	paths: list[Path] = []
	sizes: list[tuple[int, int]] = []
	placements: list[tuple[int, int]] = []
	in_table = False
	for line in lines:
		if line == "index\tfile\tw\th\tx\ty":
			in_table = True
			continue
		if not in_table:
			continue
		if not line.strip():
			break
		parts = line.split("\t")
		if len(parts) != 6:
			return None
		_, file_name, w_text, h_text, x_text, y_text = parts
		try:
			w, h = int(w_text), int(h_text)
			x, y = int(x_text), int(y_text)
		except ValueError:
			return None
		path = Path(file_name)
		if not path.is_absolute():
			path = proof_path.parent / path
		paths.append(path)
		sizes.append((w, h))
		placements.append((x, y))

	if not placements:
		return None

	if width != ratio_x * k or height != ratio_y * k:
		return None

	try:
		validate_placements(width, height, sizes, placements)
	except ValueError:
		return None

	return StoredLayout(
		proof_path=proof_path,
		output_path=_output_path_for_proof(proof_path),
		status=values.get("status", "UNKNOWN"),
		ratio_x=ratio_x,
		ratio_y=ratio_y,
		k=k,
		width=width,
		height=height,
		paths=paths,
		sizes=sizes,
		placements=placements,
	)


def _layout_sort_key(layout: StoredLayout) -> tuple[int, int, str]:
	return (layout.k, 0 if layout.status == "OPTIMAL" else 1, layout.proof_path.name)


def find_compatible_existing_layouts(
	search_dir: Path,
	sizes: list[tuple[int, int]],
	ratio_x: int,
	ratio_y: int,
) -> list[StoredLayout]:
	layouts: list[StoredLayout] = []
	for proof_path in sorted(search_dir.glob("*.proof.txt")):
		layout = _read_layout_from_proof(proof_path)
		if layout is None:
			continue
		if layout.ratio_x != ratio_x or layout.ratio_y != ratio_y:
			continue
		if layout.sizes != sizes:
			continue
		layouts.append(layout)
	return sorted(layouts, key=_layout_sort_key)


def _layout_as_heuristic(layout: StoredLayout) -> HeuristicResult:
	return HeuristicResult(
		k=layout.k,
		width=layout.width,
		height=layout.height,
		placements=layout.placements,
		method=f"existing:{layout.proof_path.name}",
	)


def write_proof_report(
	proof_path: Path,
	paths: list[Path],
	sizes: list[tuple[int, int]],
	ratio_x: int,
	ratio_y: int,
	lower: int,
	heuristic: HeuristicResult,
	result: OptimizeResult | HeuristicResult,
	status: str,
	placements: list[tuple[int, int]],
	time_limit: float,
	fast_mode: bool,
	notes: list[str] | None = None,
) -> None:
	width = result.width
	height = result.height
	k = result.k
	total_area = sum(w * h for w, h in sizes)
	canvas_area = width * height
	density = total_area / canvas_area if canvas_area else 0.0

	lines = [
		"Best dense image arrangement report",
		"===================================",
		"",
		f"mode: {_report_mode(fast_mode, status)}",
		f"status: {status}",
		f"input_count: {len(sizes)}",
		f"ratio: {ratio_x}:{ratio_y}",
		f"lower_bound_k: {lower}",
		f"heuristic_upper_k: {heuristic.k}",
		f"heuristic_method: {heuristic.method}",
		f"final_k: {k}",
		f"final_canvas: {width}x{height}",
		f"total_image_area: {total_area}",
		f"canvas_area: {canvas_area}",
		f"density: {density:.6f}",
		f"time_limit_seconds: {time_limit:g}",
	]

	if isinstance(result, OptimizeResult):
		lines.extend(
			[
				f"solver_objective_bound: {result.objective_bound}",
				f"solver_best_objective: {result.best_objective}",
				f"solver_wall_time_seconds: {result.wall_time}",
			]
		)
	if notes:
		for note in notes:
			lines.append(f"note: {note}")

	lines.extend(
		[
			"",
			"Lower-bound formula",
			"-------------------",
			"k >= max(max_i ceil(w_i / ratio_x), max_i ceil(h_i / ratio_y),",
			"         ceil(sqrt(sum_i(w_i*h_i) / (ratio_x*ratio_y))))",
			"",
			"Mathematical proof statement",
			"----------------------------",
		]
	)
	if fast_mode:
		lines.append(
			"This run used heuristic-only mode, so it proves only that the listed placement is feasible; "
			"it does not prove global optimality."
		)
	elif status == "OPTIMAL":
		lines.append(
			"Any valid integer-pixel layout at this ratio must satisfy the lower-bound formula and the "
			"containment/non-overlap constraints encoded in the CP-SAT model. OR-Tools returned OPTIMAL "
			"for final_k, proving no feasible assignment exists with a smaller k in the searched integer domain."
		)
	else:
		lines.append(
			"OR-Tools did not return OPTIMAL, so no mathematical optimality proof is claimed and the output "
			"image is a feasible fallback layout rather than a certified best layout."
		)

	lines.extend(["", "Placements", "----------", "index\tfile\tw\th\tx\ty"])
	for i, (path, (w, h), (x, y)) in enumerate(zip(paths, sizes, placements, strict=True)):
		lines.append(f"{i}\t{path.name}\t{w}\t{h}\t{x}\t{y}")

	proof_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_canvas(
	out_path: Path,
	images: list[Image.Image],
	width: int,
	height: int,
	placements: list[tuple[int, int]],
) -> None:
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


def sync_compatible_outputs(
	layouts: list[StoredLayout],
	current_out_path: Path,
	current_proof_path: Path,
	result: OptimizeResult | HeuristicResult,
	placements: list[tuple[int, int]],
	ratio_x: int,
	ratio_y: int,
	lower: int,
	upper: HeuristicResult,
	status: str,
	time_limit: float,
	fast_mode: bool,
	start_time: float,
) -> int:
	synced = 0
	seen_outputs: set[Path] = set()
	for layout in layouts:
		if layout.output_path is None:
			continue
		if _paths_equal(layout.output_path, current_out_path) or _paths_equal(layout.proof_path, current_proof_path):
			continue
		if any(_paths_equal(layout.output_path, seen) for seen in seen_outputs):
			continue

		missing = [path for path in layout.paths if not path.is_file()]
		if missing:
			log_progress(
				start_time,
				f"Skipping compatible sync target {layout.output_path.name}; missing input {missing[0].name}",
			)
			continue

		images: list[Image.Image] = []
		sizes: list[tuple[int, int]] = []
		for path in layout.paths:
			img = Image.open(path)
			img.load()
			img = img.convert("RGBA")
			images.append(img)
			sizes.append(img.size)

		if sizes != layout.sizes:
			log_progress(
				start_time,
				f"Skipping compatible sync target {layout.output_path.name}; input image sizes changed",
			)
			continue

		validate_placements(result.width, result.height, sizes, placements)
		save_canvas(layout.output_path, images, result.width, result.height, placements)
		write_proof_report(
			layout.proof_path,
			layout.paths,
			sizes,
			ratio_x,
			ratio_y,
			lower,
			upper,
			result,
			status,
			placements,
			time_limit,
			fast_mode=fast_mode,
			notes=[f"synchronized_layout_from: {current_out_path.name}"],
		)
		seen_outputs.add(layout.output_path)
		synced += 1
		log_progress(
			start_time,
			f"Synchronized compatible mosaic {layout.output_path.name} to k={result.k} "
			f"({result.width}x{result.height})",
		)
	return synced


def _build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description="Arrange images into the densest fixed-ratio non-overlapping mosaic.",
		epilog=(
			"Ratio is still accepted as optional trailing integers X Y, e.g. "
			"./auto_arrange_and_combine.py '*combined.png' 16 9"
		),
	)
	parser.add_argument("items", nargs="*", help="Input images/globs, optionally followed by ratio X Y")
	parser.add_argument("--time-limit", type=float, default=300.0, help="OR-Tools solve time limit in seconds")
	parser.add_argument("--fast", action="store_true", help="Use heuristic packing only; no optimality proof")
	parser.add_argument("--proof-file", type=Path, help="Path for the proof/report text file")
	parser.add_argument("--verbose-solver", action="store_true", help="Print detailed OR-Tools CP-SAT search logs")
	parser.add_argument(
		"--search-mode",
		choices=("fixed-k", "objective"),
		default="fixed-k",
		help=(
			"OR-Tools search strategy: fixed-k tests k=current_best-1 first; "
			"objective uses one minimize-k model"
		),
	)
	parser.add_argument(
		"--no-reuse-existing",
		action="store_true",
		help="Ignore compatible proof reports in the current directory when choosing the warm start",
	)
	parser.add_argument(
		"--no-sync-compatible",
		action="store_true",
		help="Do not rerender other compatible mosaics in the current directory with the final layout",
	)
	return parser


def main(argv: list[str]) -> int:
	start_time = time.monotonic()
	parser = _build_parser()
	args = parser.parse_args(argv[1:])

	if not args.items:
		parser.print_help(sys.stderr)
		return 2

	if args.time_limit <= 0:
		print("--time-limit must be positive.", file=sys.stderr)
		return 2

	if not args.fast:
		prepare_ortools_runtime(argv)

	ratio_x, ratio_y, ratio_given, image_args = _parse_ratio(args.items)
	log_progress(start_time, f"Expanding input files for ratio {ratio_x}:{ratio_y}")
	paths = _expand_args_to_files(image_args)
	paths = [p for p in paths if p.is_file() and not (p.name.startswith("ALL_") or p.name.startswith("All_"))]
	if not paths:
		print("No input files found (or all were skipped because they start with 'ALL_' or 'All_').", file=sys.stderr)
		return 2
	log_progress(start_time, f"Found {len(paths)} input image(s)")

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
	proof_path = args.proof_file or _default_proof_path(out_path)

	images: list[Image.Image] = []
	sizes: list[tuple[int, int]] = []
	log_progress(start_time, "Loading images at native pixel size")
	for p in paths:
		img = Image.open(p)
		img.load()
		img = img.convert("RGBA")
		images.append(img)
		sizes.append(img.size)
	log_progress(start_time, "Finished loading images")

	lower = _minimal_scale_for_ratio(sizes, ratio_x, ratio_y)
	log_progress(start_time, f"Mathematical lower bound: k >= {lower} ({ratio_x * lower}x{ratio_y * lower})")
	log_progress(start_time, "Searching deterministic heuristic upper bound")
	heuristic = find_heuristic_solution(sizes, ratio_x, ratio_y)
	validate_placements(heuristic.width, heuristic.height, sizes, heuristic.placements)
	log_progress(
		start_time,
		f"Heuristic feasible layout: k={heuristic.k}, canvas={heuristic.width}x{heuristic.height}, "
		f"method={heuristic.method}",
	)

	existing_layouts: list[StoredLayout] = []
	best_existing: StoredLayout | None = None
	effective_upper = heuristic
	if not args.no_reuse_existing:
		existing_layouts = find_compatible_existing_layouts(Path.cwd(), sizes, ratio_x, ratio_y)
		if existing_layouts:
			best_existing = existing_layouts[0]
			log_progress(
				start_time,
				f"Found {len(existing_layouts)} compatible existing layout(s) in {Path.cwd()}; "
				f"best is {best_existing.proof_path.name} with k={best_existing.k} "
				f"({best_existing.width}x{best_existing.height}, status={best_existing.status})",
			)
			if best_existing.k <= heuristic.k:
				effective_upper = _layout_as_heuristic(best_existing)
				validate_placements(effective_upper.width, effective_upper.height, sizes, effective_upper.placements)
				log_progress(start_time, f"Using existing layout as warm start: {effective_upper.method}")
			else:
				log_progress(start_time, "Existing layouts are larger than the heuristic layout; keeping heuristic warm start")

	if args.fast:
		log_progress(start_time, "Fast mode selected; saving best available heuristic/existing layout")
		placements = center_placements(effective_upper.width, effective_upper.height, sizes, effective_upper.placements)
		save_canvas(out_path, images, effective_upper.width, effective_upper.height, placements)
		write_proof_report(
			proof_path,
			paths,
			sizes,
			ratio_x,
			ratio_y,
			lower,
			effective_upper,
			effective_upper,
			"HEURISTIC_ONLY",
			placements,
			args.time_limit,
			fast_mode=True,
			notes=[f"warm_start_source: {effective_upper.method}"],
		)
		if not args.no_sync_compatible:
			sync_compatible_outputs(
				existing_layouts,
				out_path,
				proof_path,
				effective_upper,
				placements,
				ratio_x,
				ratio_y,
				lower,
				effective_upper,
				"HEURISTIC_ONLY",
				args.time_limit,
				fast_mode=True,
				start_time=start_time,
			)
		print(
			f"Wrote {out_path} ({effective_upper.width}x{effective_upper.height}, ratio {ratio_x}:{ratio_y}) "
			f"from {len(images)} images [heuristic-only; proof report: {proof_path}; runtime {_elapsed(start_time)}]"
		)
		return 0

	if best_existing is not None and best_existing.status == "OPTIMAL":
		log_progress(start_time, "Compatible existing layout is already OPTIMAL; reusing it without a new solve")
		placements = center_placements(effective_upper.width, effective_upper.height, sizes, effective_upper.placements)
		validate_placements(effective_upper.width, effective_upper.height, sizes, placements)
		reused_result = OptimizeResult(
			k=effective_upper.k,
			width=effective_upper.width,
			height=effective_upper.height,
			placements=placements,
			status="OPTIMAL",
			objective_bound=float(effective_upper.k),
			wall_time=0.0,
			best_objective=float(effective_upper.k),
		)
		save_canvas(out_path, images, effective_upper.width, effective_upper.height, placements)
		write_proof_report(
			proof_path,
			paths,
			sizes,
			ratio_x,
			ratio_y,
			lower,
			effective_upper,
			reused_result,
			"OPTIMAL",
			placements,
			args.time_limit,
			fast_mode=False,
			notes=[f"reused_existing_optimal_layout: {best_existing.proof_path.name}"],
		)
		if not args.no_sync_compatible:
			sync_compatible_outputs(
				existing_layouts,
				out_path,
				proof_path,
				reused_result,
				placements,
				ratio_x,
				ratio_y,
				lower,
				effective_upper,
				"OPTIMAL",
				args.time_limit,
				fast_mode=False,
				start_time=start_time,
			)
		print(
			f"Wrote {out_path} ({effective_upper.width}x{effective_upper.height}, ratio {ratio_x}:{ratio_y}) "
			f"from {len(images)} images [reused existing OPTIMAL layout; proof report: {proof_path}; "
			f"runtime {_elapsed(start_time)}]"
		)
		return 0

	if args.search_mode == "fixed-k":
		log_progress(
			start_time,
			f"Starting OR-Tools fixed-k search with {args.time_limit:g}s time limit "
			f"(first target k={effective_upper.k - 1})",
		)
	else:
		log_progress(start_time, f"Starting OR-Tools minimize-k search with {args.time_limit:g}s time limit")
	stop_heartbeat = threading.Event()

	def _heartbeat() -> None:
		while not stop_heartbeat.wait(30.0):
			log_progress(
				start_time,
				f"Still optimizing; current guaranteed fallback is k={effective_upper.k} "
				f"({effective_upper.width}x{effective_upper.height})",
			)

	heartbeat_thread = threading.Thread(target=_heartbeat, daemon=True)
	heartbeat_thread.start()
	try:
		if args.search_mode == "fixed-k":
			result = solve_with_fixed_k_search(
				sizes,
				ratio_x,
				ratio_y,
				effective_upper,
				args.time_limit,
				progress=lambda message: log_progress(start_time, message),
				verbose_solver=args.verbose_solver,
			)
		else:
			result = solve_with_ortools(
				sizes,
				ratio_x,
				ratio_y,
				effective_upper,
				args.time_limit,
				progress=lambda message: log_progress(start_time, message),
				verbose_solver=args.verbose_solver,
			)
	finally:
		stop_heartbeat.set()
		heartbeat_thread.join(timeout=1.0)
	log_progress(
		start_time,
		f"OR-Tools finished with status {result.status}; best k={result.k}, "
		f"bound={result.objective_bound}, solver_wall_time={result.wall_time}",
	)
	if result.status != "OPTIMAL":
		log_progress(start_time, "No optimality proof was returned; saving fallback feasible layout")
		placements = center_placements(result.width, result.height, sizes, result.placements)
		validate_placements(result.width, result.height, sizes, placements)
		save_canvas(out_path, images, result.width, result.height, placements)
		write_proof_report(
			proof_path,
			paths,
			sizes,
			ratio_x,
			ratio_y,
			lower,
			effective_upper,
			result,
			result.status,
			placements,
			args.time_limit,
			fast_mode=False,
			notes=[f"warm_start_source: {effective_upper.method}", f"solver_search_mode={args.search_mode}"],
		)
		if not args.no_sync_compatible:
			sync_compatible_outputs(
				existing_layouts,
				out_path,
				proof_path,
				result,
				placements,
				ratio_x,
				ratio_y,
				lower,
				effective_upper,
				result.status,
				args.time_limit,
				fast_mode=False,
				start_time=start_time,
			)
		print(
			f"Wrote {out_path} ({result.width}x{result.height}, ratio {ratio_x}:{ratio_y}) "
			f"from {len(images)} images [fallback after OR-Tools {result.status}; "
			f"not certified optimal; report: {proof_path}; runtime {_elapsed(start_time)}]"
		)
		return 0

	log_progress(start_time, "Optimality proved; saving certified layout")
	placements = center_placements(result.width, result.height, sizes, result.placements)
	validate_placements(result.width, result.height, sizes, placements)
	save_canvas(out_path, images, result.width, result.height, placements)
	write_proof_report(
		proof_path,
		paths,
		sizes,
		ratio_x,
		ratio_y,
		lower,
		effective_upper,
		result,
		result.status,
		placements,
		args.time_limit,
		fast_mode=False,
		notes=[f"warm_start_source: {effective_upper.method}", f"solver_search_mode={args.search_mode}"],
	)
	if not args.no_sync_compatible:
		sync_compatible_outputs(
			existing_layouts,
			out_path,
			proof_path,
			result,
			placements,
			ratio_x,
			ratio_y,
			lower,
			effective_upper,
			result.status,
			args.time_limit,
			fast_mode=False,
			start_time=start_time,
		)
	print(
		f"Wrote {out_path} ({result.width}x{result.height}, ratio {ratio_x}:{ratio_y}) "
		f"from {len(images)} images [OR-Tools {result.status}; proof report: {proof_path}; runtime {_elapsed(start_time)}]"
	)
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))

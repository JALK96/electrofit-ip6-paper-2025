"""Gaussian adapter (skeleton).

Purpose: Encapsulate creation & execution of Gaussian input (.gcrt) and related
caching / log handling behind a narrow API so domain logic does not shell out
directly.

Planned functions (initial minimal surface):
- build_input(mol2_file: str, name: str, net_charge: int, cwd: str) -> str
- run(input_file: str, name: str, cwd: str) -> None
- maybe_use_cache(name: str, scratch_dir: str, cache_env: str | None) -> bool

Current state: Domain still uses cli.run_commands.* wrappers. Migration will
wrap those here and adjust domain imports.
"""
from __future__ import annotations
import logging
import os
import shutil
import time
from .results import GaussianResult

__all__ = ["build_input", "run", "maybe_use_cache", "run_with_result"]


def build_input(mol2_file: str, name: str, net_charge: int, cwd: str) -> str:
	"""Create Gaussian input file (delegates for now to existing wrapper).

	Returns path to created .gcrt relative to cwd.
	"""
	from electrofit.cli.run_commands import create_gaussian_input
	from electrofit.io.files import modify_gaussian_input
	# If an input file already exists in the target cwd, assume the
	# preparation step (which may have run earlier) already created it and
	# skip re-running antechamber. This prevents duplicate antechamber
	# invocations when both `prepare_inputs` and
	# `build_input` call the assembler.
	gcrt = f"{name}.gcrt"
	gcrt_path = os.path.join(cwd, gcrt)
	if os.path.isfile(gcrt_path):
		logging.info("build_input: found existing %s, skipping antechamber", gcrt_path)
	else:
		create_gaussian_input(mol2_file, name, net_charge, cwd)

	# Ensure we modify the actual file on disk (pass full path).
	modify_gaussian_input(gcrt_path)
	# return relative filename (legacy contract)
	return gcrt


def maybe_use_cache(name: str, scratch_dir: str, cache_env: str | None = None) -> bool:
	"""Attempt to hydrate Gaussian artifacts from a debug cache.

	Returns True if a full skip (valid gesp present) is possible.
	"""
	cache = cache_env or os.environ.get("ELECTROFIT_DEBUG_GAUSSIAN_CACHE")
	if not cache:
		return False
	cache_base = os.path.join(cache, name)
	copied = []
	for ext in ("gcrt", "gcrt.log", "gesp"):
		src = f"{cache_base}.{ext}"
		if os.path.isfile(src):
			dst = os.path.join(scratch_dir, f"{name}.{ext}")
			shutil.copy2(src, dst)
			copied.append(dst)
	if any(p.endswith(".gesp") for p in copied):
		logging.info("[gaussian-cache] hydrated %s", copied)
		return True
	if copied:
		logging.warning("[gaussian-cache] partial cache (missing .gesp) -> still running Gaussian")
	return False


def run(input_file: str, name: str, cwd: str) -> None:
	"""Execute Gaussian calculation (delegation wrapper)."""
	from electrofit.cli.run_commands import run_gaussian_calculation
	run_gaussian_calculation(input_file, name, cwd)


def run_with_result(mol2_file: str, name: str, net_charge: int, cwd: str, try_cache: bool = True) -> GaussianResult:
	"""High-level convenience: optional cache hydration, build (if needed), run & collect timing.

	Returns GaussianResult; keeps side-effects identical to legacy path.
	"""
	t0 = time.perf_counter()
	used_cache = False
	created: list[str] = []
	warnings: list[str] = []
	gesp_path: str | None = None
	log_path: str | None = None
	if try_cache:
		used_cache = maybe_use_cache(name, cwd)
		if used_cache:
			g = os.path.join(cwd, f"{name}.gesp")
			if os.path.isfile(g):
				gesp_path = g
	if not used_cache:
		gcrt = build_input(mol2_file, name, net_charge, cwd)
		created.append(os.path.join(cwd, gcrt))
		run(input_file=gcrt, name=name, cwd=cwd)
		log_candidate = os.path.join(cwd, f"{gcrt}.log")
		if os.path.isfile(log_candidate):
			log_path = log_candidate
		g = os.path.join(cwd, f"{name}.gesp")
		if os.path.isfile(g):
			gesp_path = g
		else:
			warnings.append("Gaussian run completed but .gesp missing")
	t1 = time.perf_counter()
	return GaussianResult(
		name=name,
		input_file=os.path.join(cwd, f"{name}.gcrt"),
		gesp_file=gesp_path,
		log_file=log_path,
		used_cache=used_cache,
		wall_time_s=t1 - t0,
		created_files=created,
		warnings=warnings,
	)

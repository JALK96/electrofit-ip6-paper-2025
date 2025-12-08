from __future__ import annotations

from typing import Tuple, Iterable, Set, Optional
from pathlib import Path


def normalize_micro_name(name: str) -> str:
    """Normalize a microstate/molecule name to include the 'IP_' prefix.

    Accepts either 'IP_101101_50' or '101101_50' and returns 'IP_101101_50'.
    """
    n = (name or "").strip()
    if not n:
        return n
    return n if n.startswith("IP_") else f"IP_{n}"


def parse_molecule_args(values: Optional[Iterable[str]]) -> Optional[Set[str]]:
    """Parse --molecule flags which may be repeated or comma-separated.

    Returns a set of normalized names with 'IP_' prefix, or None if not provided.
    """
    if not values:
        return None
    result: Set[str] = set()
    for v in values:
        if not v:
            continue
        parts = [p.strip() for p in v.split(',') if p.strip()]
        for p in parts:
            result.add(normalize_micro_name(p))
    return result or None


def resolve_stage(stage: str) -> Tuple[str, str]:
    """Map a human-friendly stage to (run_dir_name, analyze_base).

    - 'final'  -> ('run_final_gmx_simulation', 'analyze_final_sim')
    - 'sample' -> ('run_gmx_simulation',       'analyze_sample_sim')
    - 'remd'   -> ('run_remd_gmx_simulation',  'analyze_remd_sim')
    """
    s = (stage or "final").strip().lower()
    if s == "sample":
        return "run_gmx_simulation", "analyze_sample_sim"
    if s == "remd":
        return "run_remd_gmx_simulation", "analyze_remd_sim"
    # default/fallback
    return "run_final_gmx_simulation", "analyze_final_sim"


def resolve_run_and_analyze_dirs(
    folder_path: Path,
    stage: str,
    run_dir_name: str,
    analyze_base: str,
    rep: Optional[int],
) -> Tuple[Path, Path]:
    """Return (run_dir, analyze_base_dir) for a given stage and optional replica.

    For 'final' and 'sample' stages, behaves as before:
      run_dir = folder_path / run_dir_name
      analyze = folder_path / analyze_base

    For 'remd', a replica index must be provided and paths are:
      run_dir = folder_path / run_dir_name / f"rep_{rep}"
      analyze = folder_path / analyze_base / f"rep_{rep}"
    """
    stage_norm = (stage or "final").strip().lower()
    run_dir = folder_path / run_dir_name
    analyze_dir = folder_path / analyze_base
    if stage_norm == "remd":
        if rep is None:
            raise ValueError("Replica index 'rep' is required when stage='remd'.")
        run_dir = run_dir / f"rep_{rep}"
        analyze_dir = analyze_dir / f"rep_{rep}"
    return run_dir, analyze_dir

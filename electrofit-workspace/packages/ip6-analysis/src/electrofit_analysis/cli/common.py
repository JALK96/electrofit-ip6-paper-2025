from __future__ import annotations

from typing import Tuple, Iterable, Set, Optional


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
    """
    s = (stage or "final").strip().lower()
    if s == "sample":
        return "run_gmx_simulation", "analyze_sample_sim"
    # default/fallback
    return "run_final_gmx_simulation", "analyze_final_sim"

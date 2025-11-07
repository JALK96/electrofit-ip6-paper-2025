from __future__ import annotations

from typing import Tuple


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


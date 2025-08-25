"""Result dataclasses for structured tool execution outcomes (experimental).

Contract-first introduction; domain orchestration does not yet consume them.
Subsequent refactors will start returning / aggregating these from higher-level
functions.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

__all__ = [
    "GaussianResult",
    "RespStageResult",
    "ConformerChargeResult",
]


@dataclass(slots=True)
class GaussianResult:
    """Outcome of a Gaussian stage (build + run or cache hydration)."""
    name: str
    input_file: str
    gesp_file: str | None
    log_file: str | None
    used_cache: bool
    wall_time_s: float
    created_files: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass(slots=True)
class RespStageResult:
    """Outcome of two-stage RESP fitting for a single conformer."""
    name: str
    resp1_in: str
    resp2_in: str
    resp1_out: str
    resp2_out: str
    final_mol2: str
    charges_file: str
    symmetry_applied: bool
    wall_time_s: float
    warnings: List[str] = field(default_factory=list)


@dataclass(slots=True)
class ConformerChargeResult:
    """High-level aggregation (Gaussian + RESP) for future orchestration use.

    Contract: symmetry_file (if not None) WAS written before two-stage RESP started.
    """
    gaussian: GaussianResult
    resp: RespStageResult
    protocol: str
    residue_name: str
    scratch_dir: str
    symmetry_file: str | None = None
    # Future metadata fields (hashes, timing breakdowns, etc.) could be added here

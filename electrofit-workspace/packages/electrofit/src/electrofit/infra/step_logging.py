"""Utilities for unified step logging and selective configuration summarisation.

English only (project guideline). Provides small helpers to:
    * Emit a consistent step prefix.
    * Log only the *relevant* subset of the loaded configuration for a step.

Historical note: A former helper ``log_consequences`` (protocol/symmetry rationale)
has been superseded by the structured DecisionModel in ``infra.decisions`` and
was removed to avoid divergent logic. Use the DecisionModel builders for any
semantic / consequence logging going forward.

These helpers are intentionally lightweight to avoid importing large config
objects or causing circular dependencies.
"""
from __future__ import annotations
from typing import Any, Iterable, Callable
import os
import logging

StepLogFn = Callable[[str], None]


def _extract(obj: Any, path: str) -> Any:
    cur = obj
    for part in path.split('.'):
        if cur is None:
            return None
        if hasattr(cur, part):
            cur = getattr(cur, part)
        elif isinstance(cur, dict):
            cur = cur.get(part)
        else:
            return None
    return cur


def log_step_header(step: str, msg: str, log_fn: StepLogFn | None = None) -> None:
    (log_fn or logging.info)(f"[{step}]{msg}")


def log_relevant_config(step: str, cfg: Any, fields: Iterable[str], log_fn: StepLogFn | None = None) -> dict[str, Any]:
    """Log only selected dotted attribute paths from cfg.

    Returns the mapping for potential further reasoning.
    """
    summary: dict[str, Any] = {}
    for f in fields:
        summary[f] = _extract(cfg, f)
    # Human readability: add space after commas already provided by join, and stable key order
    formatted = ", ".join(f"{k}={summary[k]!r}" for k in summary)
    mode = os.environ.get('ELECTROFIT_LOG_TABLE_MODE', 'table').lower()
    # modes: 'line' -> only single line, 'table' -> only table, 'both' (default) -> both
    emit_line = mode in {'both', 'line', ''}
    emit_table = mode in {'both', 'table'} and os.environ.get('ELECTROFIT_LOG_TABLE', '1') not in {'0', 'false', 'False'}

    if emit_line:
        (log_fn or logging.info)(f"[{step}][cfg] {formatted}")

    if emit_table:
        try:
            rows = list(summary.items())
            if rows:
                k_width = max(len(k) for k, _ in rows)
                k_width = min(k_width, 40)
                header = f"[{step}][cfg] ── configuration summary ──"
                (log_fn or logging.info)(header)
                for k, v in rows:
                    key = (k[:37] + '...') if len(k) > 40 else k
                    (log_fn or logging.info)(f"[{step}][cfg] {key.ljust(k_width)} : {v!r}")
        except Exception:
            logging.debug(f"[{step}][cfg] pretty table logging failed", exc_info=True)
    return summary


__all__ = [
    'log_step_header',
    'log_relevant_config',
]

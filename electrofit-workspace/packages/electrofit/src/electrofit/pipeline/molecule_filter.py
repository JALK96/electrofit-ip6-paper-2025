"""Utilities for per‑molecule step filtering.

Provides a single small helper that, given a collection of paths somewhere
under a project ``process/`` tree, returns only those belonging to a target
molecule name. This centralises the logic so every step (0‑8) can expose a
``--molecule <name>`` CLI flag without duplicating path heuristics.

Rationale:
    The various step orchestrators iterate different *levels* of the process
    hierarchy (``process/<mol>``, ``process/<mol>/run_gmx_simulation``,
    ``process/<mol>/extracted_conforms/<conf>``, etc.). We consistently define
    the *molecule* component as the immediate path element that follows the
    ``process`` segment. This helper extracts that element in a robust way.

Design notes:
    * Pure functions (no side effects / logging) to keep reuse trivial.
    * Gracefully ignore paths that do not contain a ``process`` segment.
    * Accept generic ``Iterable[Path]`` to cover lists, generators, etc.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

__all__ = [
    "molecule_component",
    "filter_paths_for_molecule",
]


def molecule_component(p: Path) -> str | None:
    """Return molecule directory name for a path inside ``process/``.

    The molecule name is defined as the path segment immediately following the
    first occurrence of ``process``. Returns ``None`` if the pattern does not
    match (e.g. path outside a process tree).
    """
    parts = list(p.parts)
    try:
        idx = parts.index("process")
    except ValueError:
        return None
    if idx + 1 >= len(parts):
        return None
    return parts[idx + 1]


def filter_paths_for_molecule(paths: Iterable[Path], molecule: str | None) -> List[Path]:
    """Return subset of ``paths`` whose molecule component equals ``molecule``.

    If ``molecule`` is ``None`` or empty, the original order is preserved and
    all paths are returned as a list. Paths without a resolvable molecule
    component are excluded when a target molecule is provided.
    """
    if not molecule:
        return list(paths)
    target = molecule.strip()
    return [p for p in paths if molecule_component(p) == target]

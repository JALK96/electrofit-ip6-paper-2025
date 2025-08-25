"""Sampling domain logic (frame selection, conformer extraction helpers).

Public API:
	* select_frame_indices(traj, n, method, seed) -> list[int]
	* prepare_conformer_directory(...)

Implementation details (e.g. RMSD farthest‑point selection) live here to keep
pipeline step modules thin and side‑effect free on import.
"""
from .frames import select_frame_indices  # noqa: F401
from .conformer_io import prepare_conformer_directory  # noqa: F401

# Backwards compatibility alias for legacy tests expecting internal name.
_select_indices = select_frame_indices  # noqa: F401

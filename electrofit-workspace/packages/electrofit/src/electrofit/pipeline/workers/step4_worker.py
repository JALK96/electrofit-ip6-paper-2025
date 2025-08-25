"""Isolated worker function for step4 to ensure stable module path for multiprocessing pickling.

This module is intentionally minimal: it imports the heavy logic lazily via the public
API of `electrofit.pipeline.steps.step4` to avoid circular import at interpreter
startup while still giving the ProcessPoolExecutor a resolvable reference.
"""
from __future__ import annotations
from pathlib import Path
from typing import Tuple

# Import the extraction helper from the dedicated extract module (no orchestrator side effects).
try:
    from electrofit.pipeline.workers.step4_extract import _extract_for_molecule  # type: ignore
except Exception as e:  # pragma: no cover
    _extract_for_molecule = None  # type: ignore
    _IMPORT_ERROR = e
else:
    _IMPORT_ERROR = None


def worker(args_tuple) -> Tuple[str, bool, str, str | None]:  # pragma: no cover
    mol_dir_str, project_root_str, sample, method, seed, override_cfg, multi_mol, verbose = args_tuple
    if _extract_for_molecule is None:
        return (mol_dir_str, False, "import_error", f"step4 extraction helper import failed: {_IMPORT_ERROR}")
    try:
        ok, msg = _extract_for_molecule(Path(mol_dir_str), Path(project_root_str), sample, method, seed, override_cfg, multi_mol, verbose)
        return (mol_dir_str, ok, msg, None)
    except Exception as e:  # pragma: no cover
        return (mol_dir_str, False, "exception", str(e))

__all__ = ["worker"]

"""Utilities to load and programmatically validate histogram manifest files.

The histogram manifest (``hist_manifest.json``) is an auxiliary diagnostic file
emitted by step6 when histogram plotting is requested. It records, per logical
histogram stage (``initial``, ``after_outlier``, ``adjusted``, ``group_combined``),
if the artefact was *expected* under the active flags and whether it was
actually *created*. A textual ``reason`` explains skipped or failed generation.

Consumer Guidelines:
* Treat unknown keys as forward-compatible additions.
* Only assert that artefacts for which ``expected`` is true eventually have
  ``created`` true; otherwise surface ``reason``.

This module centralises schema validation so individual tests remain concise.
"""
from __future__ import annotations

from pathlib import Path
import json
from typing import Dict, Any

__all__ = ["load_hist_manifest", "validate_hist_manifest"]

REQUIRED_BOOLEAN_FIELDS = {"expected", "created"}
OPTIONAL_FIELDS = {"path", "reason"}

class HistManifestError(RuntimeError):
    """Raised when the histogram manifest does not meet minimal schema guarantees."""


def load_hist_manifest(path: Path) -> Dict[str, Any]:
    """Load the histogram manifest JSON.

    Parameters
    ----------
    path: Path
        Path to ``hist_manifest.json``.
    Raises
    ------
    HistManifestError if file missing or JSON invalid.
    """
    if not path.is_file():  # pragma: no cover - defensive
        raise HistManifestError(f"Manifest file missing: {path}")
    try:
        return json.loads(path.read_text())
    except Exception as e:  # pragma: no cover - defensive
        raise HistManifestError(f"Invalid JSON in {path}: {e}") from e


def validate_hist_manifest(manifest: Dict[str, Any]) -> None:
    """Validate manifest schema.

    Checks each entry is a dict containing the required boolean fields.
    Unknown histogram keys are ignored (forward compatibility).
    Raises HistManifestError on structural violations.
    """
    for key, entry in manifest.items():
        if not isinstance(entry, dict):
            raise HistManifestError(f"Entry '{key}' must be an object, got {type(entry).__name__}")
        missing = REQUIRED_BOOLEAN_FIELDS - entry.keys()
        if missing:
            raise HistManifestError(f"Entry '{key}' missing fields: {sorted(missing)}")
        for bf in REQUIRED_BOOLEAN_FIELDS:
            if not isinstance(entry[bf], bool):
                raise HistManifestError(f"Entry '{key}.{bf}' must be bool, got {type(entry[bf]).__name__}")
    # extraneous keys allowed silently; they may appear in future evolution
        # Validate types of optional fields if present
        if "path" in entry and entry["path"] is not None and not isinstance(entry["path"], str):
            raise HistManifestError(f"Entry '{key}.path' must be string or null")
        if "reason" in entry and entry["reason"] is not None and not isinstance(entry["reason"], str):
            raise HistManifestError(f"Entry '{key}.reason' must be string or null")


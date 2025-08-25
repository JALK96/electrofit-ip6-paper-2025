"""Central logging infrastructure for electrofit.

This module supersedes the earlier top-level :mod:`electrofit.logging` module,
which is now deprecated and slated for removal. All workflow steps and
adapters should import from ``electrofit.infra.logging``.

Design goals:
    * Single active file handler per currently processed run directory.
    * Optional console (stderr) emission without duplication.
    * Deterministic suppression of stale handlers to avoid cross-run leakage.
    * Header lines and configuration tables emitted only once (no raw
        unformatted file appends, which previously caused duplicate header lines).

Environment variables:
    ELECTROFIT_LOG_LEVEL   Override root log level (default: INFO).

Public API:
    setup_logging(path, also_console=True, suppress_initial_message=False)
    log_run_header(step_name)
    reset_logging()

All public functions have concise, side‑effect focused semantics; errors in
logging setup are deliberately non-fatal to avoid interfering with scientific
pipeline execution.
"""
from __future__ import annotations

import logging
import os
import subprocess
from threading import RLock
from logging.handlers import WatchedFileHandler
from pathlib import Path

try:
    from electrofit import __version__ as _electrofit_version
except Exception:  # pragma: no cover
    _electrofit_version = "unknown"

# --- Header de-duplication state -------------------------------------------------
# Optional suppression of identical header lines within a single logical logging
# session. Enabled via environment variable ELECTROFIT_DEDUP_HEADERS=1 / true or
# programmatically via enable_header_dedup(). Always defined irrespective of
# version import success above.
_header_cache: set[str] = set()
_dedup_headers_enabled: bool = str(os.getenv("ELECTROFIT_DEDUP_HEADERS", "0")).lower() in {"1", "true", "yes", "on"}
_header_lock = RLock()

def enable_header_dedup(enable: bool = True) -> None:
    """Enable/disable in-process header de-duplication.

    When enabled, repeated emission of the exact same header line via
    log_run_header() will be suppressed until the next reset_logging() call
    (which clears the cache). This keeps global ``step.log`` concise when
    orchestration code re-initialises logging for summaries.
    """
    global _dedup_headers_enabled
    with _header_lock:
        _dedup_headers_enabled = bool(enable)
        if not enable:
            _header_cache.clear()


class ResilientWatchedFileHandler(WatchedFileHandler):
    """WatchedFileHandler that auto-recovers if the log file (or its directory)
    was deleted between tests / runs.

    pytest's TemporaryDirectory cleanup deletes the directory *after* the test
    while the handler stays registered on the root logger. Subsequent log
    records would raise FileNotFoundError inside ``emit``. We intercept that,
    recreate the parent directory, and retry exactly once. Remaining failures
    are swallowed (logging must never break the pipeline).
    """

    def emit(self, record):  # type: ignore[override]
        try:  # first attempt (normal path)
            super().emit(record)
            return
        except FileNotFoundError:
            try:
                parent = Path(getattr(self, "baseFilename", ".")).parent
                parent.mkdir(parents=True, exist_ok=True)
            except Exception:
                pass
            try:  # retry once after recreation
                super().emit(record)
            except Exception:
                # Final suppression – avoid cascading errors in tests
                pass


def setup_logging(log_path, also_console: bool = True, suppress_initial_message: bool = False) -> None:
    """Configure root logger for the current execution context.

    Parameters
    ----------
    log_path : str | Path
        Destination log file. A ``WatchedFileHandler`` is (re)used so external
        rotation/truncation tools are respected.
    also_console : bool, default True
        When True ensure exactly one stream handler to stderr. When False any
        existing non-file stream handlers are removed to silence console
        output (useful for multi-process workers).
    suppress_initial_message : bool, default False
        Suppress the standardized "Logging initialized" line if the caller
        already produced context or wants a silent re-bind of the file handler.

    Behaviour
    ---------
    * Removes all previous ``FileHandler`` instances whose target path differs
      from ``log_path``.
    * Leaves a pre-existing handler for the same file intact (idempotent).
    * Does not downgrade an already more verbose root level (e.g. DEBUG).
    * Emits one initialization INFO line unless suppressed.
    """
    path = Path(log_path).resolve()
    root = logging.getLogger()
    env_level = os.getenv("ELECTROFIT_LOG_LEVEL", "INFO").upper()
    desired_level = getattr(logging, env_level, logging.INFO)
    if root.level > desired_level:
        root.setLevel(desired_level)
    effective_level = logging.getLevelName(root.level)
    to_remove = []
    existing_same = False
    for h in root.handlers:
        if isinstance(h, logging.FileHandler):
            try:
                existing_path = Path(getattr(h, "baseFilename", ""))
                # If the directory vanished, treat as stale regardless of same path
                if not existing_path.parent.exists():
                    to_remove.append(h)
                    continue
                if existing_path.resolve() == path:
                    existing_same = True
                else:
                    to_remove.append(h)
            except Exception:
                to_remove.append(h)
    for h in to_remove:
        root.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass
    if also_console:
        if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler) for h in root.handlers):
            fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            ch = logging.StreamHandler()
            ch.setLevel(root.level)
            ch.setFormatter(fmt)
            root.addHandler(ch)
    else:
        for h in [h for h in root.handlers if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)]:
            root.removeHandler(h)
            try:
                h.close()
            except Exception:
                pass
    if not existing_same:
        # Ensure parent directory exists (temporary dirs in tests may be removed/recreated)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:  # pragma: no cover - defensive
            logging.warning("Could not ensure parent directory exists for %s: %s", path, e)
        fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        fh = ResilientWatchedFileHandler(path, mode="a", encoding="utf-8", delay=False)
        fh.setLevel(root.level)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        if not suppress_initial_message:
            root.info(f"Logging initialized. Log file: {path} (level={effective_level})")


def _git_commit_short() -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL, timeout=1)
        return out.decode().strip()
    except Exception:  # pragma: no cover
        return ""


def log_run_header(step_name: str):
    """Emit a standardized single header line for a workflow step.

    Format: ``electrofit <version> | step=<step_name> | git=<short-hash>``
    (the git part is omitted if the repository metadata cannot be resolved).

    The line is intentionally only passed through the logging formatter to
    avoid any out-of-band file writes that could create duplicates.
    """
    commit = _git_commit_short()
    parts = [f"electrofit { _electrofit_version }", f"step={step_name}"]
    if commit:
        parts.append(f"git={commit}")
    header = " | ".join(parts)
    logger = logging.getLogger()
    # Late-binding env switch: allow enabling via env var *after* module import
    # (tests may set ELECTROFIT_DEDUP_HEADERS after another import already loaded us).
    global _dedup_headers_enabled
    if not _dedup_headers_enabled and str(os.getenv("ELECTROFIT_DEDUP_HEADERS", "0")).lower() in {"1", "true", "yes", "on"}:
        _dedup_headers_enabled = True
    if _dedup_headers_enabled:
        with _header_lock:
            if header in _header_cache:
                return
            _header_cache.add(header)
    logger.info(header)


def reset_logging():
    """Fully remove all handlers from the root logger (and known children).

    Use this before switching to a new per-run log file to prevent future log
    records from leaking into earlier run files. This performs a clean
    shutdown, closes handlers, and clears filters. It is deliberately robust
    and silent upon errors.
    """
    logging.shutdown()
    root_logger = logging.getLogger()
    handlers = root_logger.handlers[:]
    for handler in handlers:
        root_logger.removeHandler(handler)
        handler.close()
    logger_dict = logging.Logger.manager.loggerDict
    for logger_name in logger_dict:
        logger = logging.getLogger(logger_name)
        if hasattr(logger, "handlers"):
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
                handler.close()
            logger.filters = []
    # Clearing header cache upon full reset allows headers to be emitted again
    # for new log files (prevents accidental permanent suppression across steps).
    with _header_lock:
        _header_cache.clear()

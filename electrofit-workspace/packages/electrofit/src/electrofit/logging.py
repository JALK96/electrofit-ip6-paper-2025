import logging
import os
import subprocess
from logging.handlers import WatchedFileHandler
from pathlib import Path

try:
    from electrofit import __version__ as _electrofit_version
except Exception:  # pragma: no cover
    _electrofit_version = "unknown"


def setup_logging(log_path, also_console: bool = True, suppress_initial_message: bool = False) -> None:
    """Configure root logging to write ONLY to the specified file (plus optional console).

    Previous behaviour accumulated FileHandlers for each molecule/run, causing later
    log messages to be duplicated into earlier process.log files. We now enforce a
    single active FileHandler (the newest target). Console handler is preserved / added once.
    """
    path = Path(log_path).resolve()
    root = logging.getLogger()
    # Determine desired level (env overrides, default INFO)
    env_level = os.getenv("ELECTROFIT_LOG_LEVEL", "INFO").upper()
    desired_level = getattr(logging, env_level, logging.INFO)
    # Do not downgrade if a more verbose (lower numeric) level already set (e.g., DEBUG=10 < INFO=20)
    if root.level > desired_level:
        root.setLevel(desired_level)
    effective_level = logging.getLevelName(root.level)

    # Remove any existing FileHandlers pointing to different files
    to_remove = []
    existing_same = False
    for h in root.handlers:
        if isinstance(h, logging.FileHandler):
            try:
                existing_path = Path(getattr(h, "baseFilename", "")).resolve()
                if existing_path == path:
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

    # Manage console handler(s) BEFORE emitting any initialization message
    if also_console:
        # Ensure exactly one console handler (pure StreamHandler) if requested
        if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler) for h in root.handlers):
            fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            ch = logging.StreamHandler()
            ch.setLevel(root.level)
            ch.setFormatter(fmt)
            root.addHandler(ch)
    else:
        # Remove existing non-file StreamHandlers to silence console output for this context
        to_remove_console = [
            h for h in root.handlers
            if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        ]
        for h in to_remove_console:
            root.removeHandler(h)
            try:
                h.close()
            except Exception:
                pass

    # If already have a handler for this exact path, keep it (respect suppress flag logic below)
    if not existing_same:
        fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        fh = WatchedFileHandler(path, mode="a", encoding="utf-8", delay=False)
        fh.setLevel(root.level)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        if not suppress_initial_message:
            root.info(f"Logging initialized. Log file: {path} (level={effective_level})")


def _git_commit_short() -> str:
    """Return short git commit hash if available (else empty)."""
    try:
        out = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL, timeout=1)
        return out.decode().strip()
    except Exception:  # pragma: no cover - fallback on any failure
        return ""


def log_run_header(step_name: str):
    """Emit a standardized header line (version, git commit, step)."""
    commit = _git_commit_short()
    parts = [f"electrofit { _electrofit_version }", f"step={step_name}"]
    if commit:
        parts.append(f"git={commit}")
    header = " | ".join(parts)
    logger = logging.getLogger()
    logger.info(header)


def reset_logging():
    """
    Reset the logging configuration by removing all handlers.
    This allows reinitialization of logging as if starting fresh.
    """
    # Shutdown logging and flush any pending messages
    logging.shutdown()

    # Remove handlers from the root logger
    root_logger = logging.getLogger()
    handlers = root_logger.handlers[:]
    for handler in handlers:
        root_logger.removeHandler(handler)
        handler.close()

    # Optionally, reset configurations of all existing loggers
    # This is useful if you have child loggers that also need resetting
    logger_dict = logging.Logger.manager.loggerDict
    for logger_name in logger_dict:
        logger = logging.getLogger(logger_name)
        if hasattr(logger, "handlers"):
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
                handler.close()
            # Remove filters if any
            logger.filters = []

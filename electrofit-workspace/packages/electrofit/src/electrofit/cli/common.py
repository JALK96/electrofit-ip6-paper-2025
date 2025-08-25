"""Common CLI helpers for electrofit steps."""
from __future__ import annotations
import argparse
import os

def add_standard_flags(ap: argparse.ArgumentParser, project: bool = True, config: bool = True, log_console: bool = True):
    if project:
        ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()), help="Project root path (contains process/)")
    if config:
        ap.add_argument("--config", help="Optional electrofit.toml overrides merged into snapshot")
    if log_console:
        ap.add_argument("--log-console", action="store_true", help="Echo logs to console in addition to file")
    return ap

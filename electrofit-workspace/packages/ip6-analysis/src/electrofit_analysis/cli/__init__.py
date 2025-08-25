"""CLI subpackage for ip6-analysis (formerly electrofit-analysis).

This package contains command-line interfaces for various analysis tasks.
The console scripts should point to :mod:`electrofit_analysis.cli.app:main`.
This module provides a small compatibility shim so ``electrofit_analysis.cli:main``
continues to work for older entry points.
"""

from __future__ import annotations

from .app import main  # re-export for backward compatibility

__all__ = ["main"]

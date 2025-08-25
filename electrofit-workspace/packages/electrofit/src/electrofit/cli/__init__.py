"""Electrofit CLI package.

Expose the ``main`` entry point for ``python -m electrofit`` execution.
"""
from .app import main as main  # explicit re-export for linters

__all__ = ["main"]

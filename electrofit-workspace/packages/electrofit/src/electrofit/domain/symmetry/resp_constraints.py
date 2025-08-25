"""Shared helpers for applying RESP symmetry constraints.

This module centralises logic that was duplicated between initial structure
preparation and conformer charge derivation. It intentionally keeps a small
surface so future adapter refactors (e.g. alternative RESP engines) only need
to adjust one place.

All file operations are delegated through existing domain/cli wrappers via
`run_python` (callers supply it) to keep this module pure and side‑effect free
except for logging and explicit function returns.
"""
from __future__ import annotations

from dataclasses import dataclass
import os
import logging
from typing import Protocol, Callable


class PythonRunner(Protocol):  # narrow protocol for run_python
    def __call__(self, func: Callable, *args, cwd: str, **kwargs): ...  # noqa: D401


@dataclass(slots=True)
class RespSymmetryConfig:
    scratch_dir: str
    adjust_sym: bool
    ignore_sym: bool
    # filenames (relative to scratch_dir)
    resp1_in: str = "ANTECHAMBER_RESP1.IN"
    resp1_mod: str = "ANTECHAMBER_RESP1_MOD.IN"
    symmetry_txt_base: str = "symmetry_resp.txt"


def apply_and_optionally_modify(run_python: PythonRunner, find_file_with_extension, edit_resp_input, write_symmetry, cfg: RespSymmetryConfig) -> str:
    """Apply symmetry editing if requested; return RESP1 file to use.

    Returns path (relative or absolute) of the RESP1 input that downstream RESP
    stage should consume (modified or original). Raises FileNotFoundError if a
    required JSON file is missing when adjustment requested.
    """
    if not cfg.adjust_sym:
        return os.path.join(cfg.scratch_dir, cfg.resp1_in)

    json_filename = run_python(find_file_with_extension, "json", cwd=cfg.scratch_dir)
    if not json_filename:
        raise FileNotFoundError("No *.json symmetry file found in scratch.")
    json_symmetry_file = json_filename if os.path.isabs(json_filename) else os.path.join(cfg.scratch_dir, json_filename)

    # First variant: possibly honour constraints or ignore (two passes handled externally if needed)
    run_python(
        edit_resp_input,
        input_file=cfg.resp1_in,
        equiv_groups_file=json_symmetry_file,
        output_file=cfg.resp1_mod,
        ignore_sym=cfg.ignore_sym,
        cwd=cfg.scratch_dir,
    )
    run_python(
        write_symmetry,
        os.path.join(cfg.scratch_dir, cfg.resp1_mod),
        ("symmetry_resp_MOD.txt" if cfg.adjust_sym else cfg.symmetry_txt_base),
        cwd=cfg.scratch_dir,
    )
    if cfg.adjust_sym and not cfg.ignore_sym:
        logging.debug("Applied RESP symmetry constraints (ignore_sym=%s)", cfg.ignore_sym)
    return os.path.join(cfg.scratch_dir, cfg.resp1_mod)

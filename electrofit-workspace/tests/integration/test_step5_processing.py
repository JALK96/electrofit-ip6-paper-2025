"""Integration tests for step5 native conformer processing.

We simulate the presence of extracted conformer directories each with a PDB file.
For speed we monkeypatch the heavy processing function `process_conform` to only
drop a marker file so we can assert batching / execution order without running
Gaussian/RESP.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
import types
import os

from tests.helpers.project import make_project_tree


def _make_conformer(dir_: Path, idx: int):
    dir_.mkdir(parents=True, exist_ok=True)
    pdb = dir_ / f"LIGc{idx}.pdb"
    pdb.write_text("MODEL\nENDMDL\nEND\n")
    return pdb


def test_step5_dry_run_and_execute(tmp_path: Path):
    proj = tmp_path / "proj"
    make_project_tree(proj)  # sets up project skeleton with process dir

    # Create 5 fake conformer dirs with pdb under process/<mol>/extracted_conforms/*
    base = proj / "process" / "IP_011101" / "extracted_conforms"
    conf_dirs = []
    for i in range(5):
        _make_conformer(base / f"IP_011101c{i}", i)
        conf_dirs.append(base / f"IP_011101c{i}")

    # Using --mock to avoid heavy Gaussian/RESP pipeline

    # Dry run should list but not execute (no executed.txt files)
    out = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "electrofit",
            "step5",
            "--project",
            str(proj),
            "--batch-size",
            "2",
            "--interval",
            "0",
            "--dry-run",
            "--mock",
        ],
        text=True,
    )
    assert "[step5] Found 5 conformer" in out
    assert "Executing batch 1" in out
    # no executed.txt created
    assert not any((p / "executed.txt").exists() for p in conf_dirs)

    # Real run (interval 0 for speed). Expect marker files.
    out2 = subprocess.check_output(
        [
            sys.executable,
            "-m",
            "electrofit",
            "step5",
            "--project",
            str(proj),
            "--batch-size",
            "3",
            "--interval","0","--mock",
        ],
        text=True,
    )
    assert "Executing batch 1" in out2
    assert "Executing batch 2" in out2  # 5 -> 3 + 2
    markers = [(p / "executed.txt").read_text().strip() for p in conf_dirs]
    # Each marker contains run<name>; ensure count matches
    assert len(markers) == 5

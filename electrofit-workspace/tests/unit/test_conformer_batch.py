from __future__ import annotations

import os
from pathlib import Path
import pytest

from electrofit.domain.charges.conformer_batch import discover_conformer_dirs, process_conformer_dir
from tests.helpers.project import make_project_tree


def _make_conformer_dir(base: Path, name: str, with_pdb: bool = True):
    d = base / name
    d.mkdir(parents=True, exist_ok=True)
    if with_pdb:
        (d / f"{name}.pdb").write_text("MODEL\nENDMDL\nEND\n")
    return d


def test_discover_conformer_dirs_basic(tmp_path: Path):
    proj = make_project_tree(tmp_path / "proj")
    ec_root = proj["root"] / "process" / "IP_011101" / "extracted_conforms"
    # Create three conformers, one missing pdb (should be excluded)
    c0 = _make_conformer_dir(ec_root, "IP_011101c0")
    c1 = _make_conformer_dir(ec_root, "IP_011101c1")
    c2 = _make_conformer_dir(ec_root, "IP_011101c2", with_pdb=False)
    found = discover_conformer_dirs(proj["root"])
    assert c0 in found and c1 in found and c2 not in found
    # Deterministic ordering
    assert found == sorted(found)


def test_process_conformer_dir_no_pdb(tmp_path: Path):
    proj = make_project_tree(tmp_path / "proj")
    ec_root = proj["root"] / "process" / "IP_011101" / "extracted_conforms"
    bad = _make_conformer_dir(ec_root, "IP_011101cX", with_pdb=False)
    rel, ok, msg = process_conformer_dir(bad, proj["root"], None, False, mock=True, verbose=False)
    assert ok is False
    assert "PDB" in msg or "pdb" in msg


def test_process_conformer_dir_snapshot_warning(tmp_path: Path, monkeypatch):
    proj = make_project_tree(tmp_path / "proj")
    ec_root = proj["root"] / "process" / "IP_011101" / "extracted_conforms"
    good = _make_conformer_dir(ec_root, "IP_011101c1")
    # Remove upstream sources so compose_snapshot returns None, forcing warning path
    for p in [
        proj["root"] / "process" / "IP_011101" / "electrofit.toml",
        proj["root"] / "data" / "input" / "IP_011101" / "electrofit.toml",
        proj["root"] / "electrofit.toml",
    ]:
        if p.exists():
            p.unlink()
    rel, ok, msg = process_conformer_dir(good, proj["root"], None, False, mock=True, verbose=False)
    assert ok is True  # still succeeds in mock mode
    # process.log should exist
    assert (good / "process.log").is_file()


def test_disable_advanced_diag_env(tmp_path: Path, monkeypatch):
    proj = make_project_tree(tmp_path / "proj")
    ec_root = proj["root"] / "process" / "IP_011101" / "extracted_conforms"
    good = _make_conformer_dir(ec_root, "IP_011101c5")
    monkeypatch.setenv("ELECTROFIT_DISABLE_ADV_DIAG", "1")
    rel, ok, msg = process_conformer_dir(good, proj["root"], None, False, mock=True, verbose=False)
    assert ok is True
    assert rel.endswith("IP_011101c5")

"""Tests for electrofit.domain.sampling utilities.

Focus:
  * select_frame_indices edge cases (n=0, n>len, maxmin determinism with seed)
  * prepare_conformer_directory snapshot refresh & symmetry JSON copy
"""
from __future__ import annotations
import time
from pathlib import Path
import json

import numpy as np
import pytest

import mdtraj as md

from electrofit.domain.sampling import select_frame_indices, prepare_conformer_directory


def _make_traj(n_frames: int, n_atoms: int = 3) -> md.Trajectory:
    """Build a tiny trajectory with linear coordinate changes for testing."""
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue('LIG', chain)
    from mdtraj.core.element import carbon
    atoms = [top.add_atom(f'C{i}', carbon, res) for i in range(n_atoms)]  # noqa: F841
    # Coordinates: frame index encoded in x coordinate to ensure non-zero RMSD
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=float)
    for f in range(n_frames):
        xyz[f, :, 0] = f  # shift all atoms in x for distinct frames
    return md.Trajectory(xyz, top)


@pytest.mark.parametrize("n", [0, 1])
def test_select_frame_indices_n_zero_or_one(n):
    traj = _make_traj(5)
    indices = select_frame_indices(traj, n, method="linear", seed=123)
    assert indices == ([] if n == 0 else [0])


def test_select_frame_indices_n_exceeds_length():
    traj = _make_traj(4)
    indices = select_frame_indices(traj, 10, method="linear", seed=None)
    assert indices == list(range(4))


def test_select_frame_indices_maxmin_deterministic_seed():
    traj = _make_traj(6)
    a = select_frame_indices(traj, 4, method="maxmin", seed=42)
    b = select_frame_indices(traj, 4, method="maxmin", seed=42)
    c = select_frame_indices(traj, 4, method="maxmin", seed=7)
    assert a == b  # same seed -> deterministic
    assert a != c  # different seed -> likely different start frame
    assert len(a) == 4
    assert len(set(a)) == 4  # unique


def test_prepare_conformer_directory_snapshot_and_symmetry(tmp_path: Path):
    # Directory scaffold
    extracted = tmp_path / "extracted_conforms"
    extracted.mkdir()
    pis_dir = tmp_path / "run_gau_create_gmx_in"
    pis_dir.mkdir()
    # Parent snapshot
    parent_snapshot = tmp_path / "parent.toml"
    parent_snapshot.write_text("[project]\nname='x'\n")
    # Symmetry JSON (will be discovered because equiv_groups_file=None)
    (pis_dir / "equiv_groups.json").write_text(json.dumps({"groups": [[1,2]]}))

    # Fake save function writes simple marker content
    def save_fn(p: str):
        Path(p).write_text("PDB\n")

    pdb_path = prepare_conformer_directory(
        conform_index=0,
        molecule_name="MOL",
        parent_cfg_target=parent_snapshot,
        override_cfg=None,
        protocol="bcc",
        respin1_file=None,
        respin2_file=None,
        equiv_groups_file=None,
        pis_dir=pis_dir,
        extracted_conforms_dir=extracted,
        input_mol2_file=tmp_path / "missing.mol2",  # not present -> skip mol2 conversion
        traj_frame_save_fn=save_fn,
        verbose=False,
    )
    conform_dir = pdb_path.parent
    # Snapshot copied
    snap = conform_dir / "electrofit.toml"
    assert snap.is_file()
    # Symmetry JSON copied
    assert (conform_dir / "equiv_groups.json").is_file()

    # Modify parent snapshot and force refresh via override_cfg flag path
    time.sleep(0.01)  # ensure mtime advances
    parent_snapshot.write_text("[project]\nname='y'\n")
    override_cfg = tmp_path / "override.toml"
    override_cfg.write_text("[sampling]\ncount=10\n")
    prepare_conformer_directory(
        conform_index=0,
        molecule_name="MOL",
        parent_cfg_target=parent_snapshot,
        override_cfg=override_cfg,
        protocol="bcc",
        respin1_file=None,
        respin2_file=None,
        equiv_groups_file=None,
        pis_dir=pis_dir,
        extracted_conforms_dir=extracted,
        input_mol2_file=tmp_path / "missing.mol2",
        traj_frame_save_fn=save_fn,
        verbose=False,
    )
    # After refresh snapshot content should reflect new parent text
    assert "name='y'" in snap.read_text()

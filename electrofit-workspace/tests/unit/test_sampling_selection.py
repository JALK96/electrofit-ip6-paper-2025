import pathlib
import sys
import numpy as np
import mdtraj as md

# Ensure package import path
sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))

from electrofit.domain.sampling import _select_indices  # type: ignore  # noqa: E402


def _make_linear_traj(n_frames=10, n_atoms=3):
    # simple coordinate progression so RMSD has structure
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    for i in range(n_frames):
        xyz[i, :, 0] = i  # shift along x per frame
    # Build a minimal topology (chain->residue->atoms)
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue('RES', chain)
    for a in range(n_atoms):
        top.add_atom(f'C{a}', element=md.element.carbon, residue=res)
    return md.Trajectory(xyz, top)


def test_linear_selection_exact():
    traj = _make_linear_traj(9)
    idx = _select_indices(traj, 3, 'linear', None)
    assert idx == [0, 4, 8]


def test_random_seed_determinism():
    traj = _make_linear_traj(20)
    idx1 = _select_indices(traj, 5, 'random', 123)
    idx2 = _select_indices(traj, 5, 'random', 123)
    assert idx1 == idx2
    assert len(set(idx1)) == 5


def test_maxmin_diversity():
    traj = _make_linear_traj(12)
    idx = _select_indices(traj, 4, 'maxmin', 2)
    # uniqueness
    assert len(idx) == len(set(idx)) == 4
    # indices within range
    assert all(0 <= i < 12 for i in idx)


def test_request_more_than_available():
    traj = _make_linear_traj(5)
    idx = _select_indices(traj, 10, 'random', 42)
    assert idx == [0,1,2,3,4]

import pathlib, sys
import numpy as np
import mdtraj as md

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))
from electrofit.domain.sampling import _select_indices  # noqa: E402


def _traj(n=3, atoms=2):
    xyz = np.zeros((n, atoms, 3), dtype=np.float32)
    top = md.Topology(); chain = top.add_chain(); res = top.add_residue('R', chain)
    for i in range(atoms):
        top.add_atom(f'C{i}', element=md.element.carbon, residue=res)
    return md.Trajectory(xyz, top)


def test_zero_requested_frames():
    t = _traj(5)
    idx = _select_indices(t, 0, 'linear', None)
    assert idx == []


def test_maxmin_small_traj_single_frame():
    t = _traj(1)
    idx = _select_indices(t, 1, 'maxmin', 1)
    assert idx == [0]


def test_maxmin_small_traj_two_frames():
    t = _traj(2)
    idx = _select_indices(t, 2, 'maxmin', 1)
    assert set(idx) == {0,1}

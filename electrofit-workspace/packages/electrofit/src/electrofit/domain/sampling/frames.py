"""Frame selection algorithms for conformer sampling.

Currently supported methods:
  * linear  : Evenly spaced indices using linspace
  * random  : Uniform random without replacement (seeded)
  * maxmin  : Farthest point sampling in RMSD space

The API is intentionally stateless and deterministic for a given seed to aid
reproducibility and testability.
"""
from __future__ import annotations
from typing import List
import random
import numpy as np
import mdtraj as md

__all__ = ["select_frame_indices"]


def select_frame_indices(traj: md.Trajectory, n: int, method: str, seed: int | None) -> List[int]:
    """Return a list of frame indices according to the chosen sampling method.

    Parameters
    ----------
    traj : md.Trajectory
        Loaded trajectory.
    n : int
        Desired number of frames.
    method : str
        Sampling strategy name.
    seed : int | None
        Random seed (only used for stochastic strategies / starting point selection).
    """
    total = len(traj)
    if total == 0 or n <= 0:
        return []
    if n >= total:
        return list(range(total))
    method_norm = (method or "").lower()
    if method_norm in {"linear", "even", "linspace"}:
        return list(np.linspace(0, total - 1, num=n, dtype=int))
    if method_norm in {"random", "rand"}:
        rng = random.Random(seed)
        return rng.sample(range(total), n)
    if method_norm in {"maxmin", "minmax", "farthest"}:
        if seed is not None:
            rng = random.Random(seed)
            first = rng.randrange(total)
        else:
            first = 0
        selected = [first]
        min_d = md.rmsd(traj, traj, first)
        min_d[first] = -np.inf
        for _ in range(1, n):
            next_idx = int(np.argmax(min_d))
            selected.append(next_idx)
            if len(selected) == n:
                break
            d_new = md.rmsd(traj, traj, next_idx)
            min_d = np.minimum(min_d, d_new)
            min_d[next_idx] = -np.inf
        return selected
    # Fallback -> linear
    return list(np.linspace(0, total - 1, num=n, dtype=int))

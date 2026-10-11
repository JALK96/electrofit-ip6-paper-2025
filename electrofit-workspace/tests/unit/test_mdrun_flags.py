import pytest

from electrofit.adapters.gromacs import mdrun_flags


def test_cpu_flags_use_thread_mpi():
    assert mdrun_flags(threads=16, pin=False, gpu=False) == "-nt 16 -pin off -nobackup"
    assert mdrun_flags(threads=1, pin=None, gpu=False) == "-nt 1 -ntmpi 1 -nobackup"


def test_gpu_dynamics_use_one_rank_and_leave_placement_to_gromacs():
    flags = mdrun_flags(threads=8, pin=True, gpu=True)
    assert flags == "-ntmpi 1 -ntomp 8 -pin on -nb gpu -nobackup"


def test_gpu_explicit_placement_applies_to_dynamics_only():
    offload = dict(gpu_id=0, pme="gpu", bonded="cpu", update="cpu")
    dynamics = mdrun_flags(threads=8, pin=None, gpu=True, **offload)
    assert dynamics == "-ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded cpu -update cpu -gpu_id 0 -nobackup"
    minimization = mdrun_flags(threads=8, pin=None, gpu=True, minimization=True, **offload)
    assert minimization == "-ntmpi 1 -ntomp 8 -nb gpu -gpu_id 0 -nobackup"


def test_gpu_placement_ignored_without_gpu():
    assert mdrun_flags(threads=4, pin=None, gpu=False, pme="gpu", gpu_id="1") == "-nt 4 -nobackup"


def test_invalid_placement_is_rejected():
    with pytest.raises(ValueError):
        mdrun_flags(threads=4, pin=None, gpu=True, update="yes")


def test_checkpoint_interval_applies_to_dynamics_only():
    assert mdrun_flags(threads=8, pin=None, gpu=True, checkpoint_minutes=5) == "-ntmpi 1 -ntomp 8 -nb gpu -cpt 5 -nobackup"
    assert "-cpt" not in mdrun_flags(threads=8, pin=None, gpu=True, checkpoint_minutes=5, minimization=True)

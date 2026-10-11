"""Continuation of interrupted production runs with a shimmed gmx."""
import os
import types

import pytest

from electrofit.adapters.gromacs import resume_production, set_up_production
from electrofit.infra.scratch_manager import scratch_dir_for

GRO = "mol_GMX.gro"
SHIM = """#!/usr/bin/env bash
echo "gmx $*" >> "$GMX_CALLS"
if [[ $1 == dump ]]; then
  echo "   step = $(cat "${3%/*}/fake_step" 2>/dev/null || echo 100)"
  exit 0
fi
if [[ $1 == mdrun ]]; then
  [[ " $* " == *" -cpi md.cpt "* ]] || exit 3
  echo continued >> md.xtc
  echo continued >> md.log
  echo final > md.cpt
fi
exit 0
"""


@pytest.fixture
def env(tmp_path, monkeypatch):
    shim = tmp_path / "bin"
    shim.mkdir()
    (shim / "gmx").write_text(SHIM)
    (shim / "gmx").chmod(0o755)
    calls = tmp_path / "calls.txt"
    monkeypatch.setenv("PATH", f"{shim}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("GMX_CALLS", str(calls))
    run_dir = tmp_path / "proj" / "process" / "mol" / "run_final_gmx_simulation"
    (run_dir / "MDP").mkdir(parents=True)
    for name in (GRO, "mol_GMX.itp", "mol_GMX.top", "posre_mol.itp", "MDP/Production.mdp"):
        (run_dir / name).write_text(f"; {name}\n")
    (run_dir / "ff_x.ff").mkdir()
    (run_dir / "ff_x.ff" / "forcefield.itp").write_text("; ff\n")
    base = tmp_path / "scratch"
    monkeypatch.chdir(run_dir)
    return types.SimpleNamespace(run_dir=run_dir, base=base, calls=calls,
                                 scratch=scratch_dir_for(str(run_dir), str(base)))


def _resume(env):
    resume_production(GRO, "MDP", str(env.base), ff="ff_x.ff", threads=4, pin=None, gpu=True,
                      postprocess_on_scratch=False)


def _write_state(directory):
    os.makedirs(directory, exist_ok=True)
    for name, text in (("md.tpr", "tpr"), ("md.cpt", "cpt"), ("md.xtc", "frames\n"), ("md.log", "   nsteps    = 1000\nlog\n")):
        with open(os.path.join(directory, name), "w") as handle:
            handle.write(text)


def test_resume_from_scratch_after_hard_stop(env):
    _write_state(env.scratch)
    _resume(env)
    assert "gmx mdrun -deffnm md -cpi md.cpt -ntmpi 1 -ntomp 4 -nb gpu -nobackup" in env.calls.read_text()
    assert (env.run_dir / "md.xtc").read_text() == "frames\ncontinued\n"
    assert (env.run_dir / "md.cpt").read_text() == "final\n"
    assert (env.run_dir / "md.log").read_text().endswith("log\ncontinued\n")
    assert not os.path.exists(env.scratch)


def test_resume_from_run_directory_keeps_one_copy(env):
    _write_state(env.run_dir)
    _resume(env)
    assert (env.run_dir / "md.xtc").read_text() == "frames\ncontinued\n"
    leftovers = [p.name for p in env.run_dir.iterdir() if "_copy" in p.name or ".input_file" in p.name]
    assert leftovers == []
    assert not os.path.exists(env.scratch)


def test_resume_refuses_two_checkpoints(env):
    _write_state(env.scratch)
    _write_state(env.run_dir)
    with pytest.raises(RuntimeError):
        _resume(env)


def test_resume_without_checkpoint_fails(env):
    with pytest.raises(FileNotFoundError):
        _resume(env)


def test_new_run_refuses_to_overwrite_a_resumable_scratch(env):
    _write_state(env.scratch)
    with pytest.raises(RuntimeError):
        set_up_production(GRO, "MDP", str(env.base), "mol", simulation=None,
                          derived=types.SimpleNamespace(mode="target"), ff="ff_x.ff")


def test_resume_skips_a_completed_run(env):
    _write_state(env.run_dir)
    (env.run_dir / "fake_step").write_text("1000")
    assert resume_production(GRO, "MDP", str(env.base), ff="ff_x.ff", threads=4, pin=None, gpu=False,
                             postprocess_on_scratch=False) is False
    assert (env.run_dir / "md.xtc").read_text() == "frames\n"
    assert not os.path.exists(env.scratch)

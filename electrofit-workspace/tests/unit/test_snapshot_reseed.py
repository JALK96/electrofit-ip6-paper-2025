"""Regression tests for clean snapshot reseed semantics.

These tests ensure that changes in upstream or other source layers propagate
on reruns and that overrides persist while fills are re-applied.

Project language guideline: Keep docstring in English.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import tomllib

from electrofit.infra.config_snapshot import compose_snapshot


def write(path: Path, content: str) -> None:
    path.write_text(textwrap.dedent(content).lstrip())


def read_toml(path: Path) -> dict:
    with path.open("rb") as f:
        return tomllib.load(f)


def test_reseed_upstream_propagates(tmp_path: Path) -> None:
    project_root = tmp_path
    # project defaults
    write(project_root / "electrofit.toml", "[defaults]\np=1\n")

    # molecule input layer
    mol_input = project_root / "data" / "input" / "MOLX"
    mol_input.mkdir(parents=True)
    write(mol_input / "electrofit.toml", "[defaults]\nq=2\n")

    # process cfg (will hold overrides)
    process_cfg_dir = project_root / "process" / "MOLX"
    process_cfg_dir.mkdir(parents=True)
    write(process_cfg_dir / "electrofit.toml", "[overrides]\nz=9\n")

    # upstream result (initial)
    results_dir = process_cfg_dir / "results"
    results_dir.mkdir()
    write(results_dir / "electrofit.toml", "[defaults]\na=1\n[overrides]\nz=3\n")

    run_dir = process_cfg_dir / "run"
    run_dir.mkdir()

    snap1 = compose_snapshot(
        run_dir,
        project_root,
        molecule="MOLX",
        multi_molecule=False,
        log_fn=lambda s: None,
        upstream=results_dir / "electrofit.toml",
        process_cfg=process_cfg_dir / "electrofit.toml",
        molecule_input=mol_input / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
    )
    assert snap1 is not None
    data1 = read_toml(snap1)
    # upstream a=1 should seed; overrides z from process_cfg supersedes upstream z=3
    assert data1["defaults"]["a"] == 1
    assert data1["overrides"]["z"] == 9

    # mutate upstream and project defaults; add new key to defaults
    write(results_dir / "electrofit.toml", "[defaults]\na=42\n[overrides]\nz=5\n")
    write(project_root / "electrofit.toml", "[defaults]\np=7\n[new]\nadded=1\n")

    snap2 = compose_snapshot(
        run_dir,
        project_root,
        molecule="MOLX",
        multi_molecule=False,
        log_fn=lambda s: None,
        upstream=results_dir / "electrofit.toml",
        process_cfg=process_cfg_dir / "electrofit.toml",
        molecule_input=mol_input / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
    )
    assert snap2 == snap1  # same path reused
    data2 = read_toml(snap2)
    # a updated from upstream; override z persists (9); new default key p updated; fill adds [new]
    assert data2["defaults"]["a"] == 42
    assert data2["overrides"]["z"] == 9
    assert data2["defaults"]["p"] == 7
    assert data2["new"]["added"] == 1


def test_reseed_missing_upstream_fallback(tmp_path: Path) -> None:
    project_root = tmp_path
    write(project_root / "electrofit.toml", "[defaults]\nroot=1\n")
    mol_input = project_root / "data" / "input" / "MOLY"
    mol_input.mkdir(parents=True)
    write(mol_input / "electrofit.toml", "[defaults]\nlocal=2\n")
    proc_dir = project_root / "process" / "MOLY"
    proc_dir.mkdir(parents=True)
    write(proc_dir / "electrofit.toml", "[overrides]\nalpha=3\n")
    run_dir = proc_dir / "run"
    run_dir.mkdir()

    # No upstream provided -> should seed from molecule_input (next layer) then apply overrides
    snap = compose_snapshot(
        run_dir,
        project_root,
        molecule="MOLY",
        multi_molecule=False,
        log_fn=lambda s: None,
        upstream=None,
        process_cfg=proc_dir / "electrofit.toml",
        molecule_input=mol_input / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
    )
    assert snap is not None
    data = read_toml(snap)
    assert data["defaults"]["local"] == 2  # came from molecule input
    assert data["defaults"]["root"] == 1  # filled from project defaults
    assert data["overrides"]["alpha"] == 3

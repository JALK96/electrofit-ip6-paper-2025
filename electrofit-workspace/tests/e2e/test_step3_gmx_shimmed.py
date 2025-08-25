# tests/e2e/test_step3_gmx_shimmed.py
import os, sys, subprocess
from pathlib import Path
from tests.helpers.project import make_project_tree, install_minimal_gmx_fixture

def test_step3_gmx_shimmed(tmp_path, shim_bin, monkeypatch):
    monkeypatch.setenv("PATH", str(shim_bin)+os.pathsep+os.environ["PATH"])
    proj = tmp_path
    make_project_tree(proj)

    # step1
    r1 = subprocess.run([sys.executable,"-m","electrofit","step1","--project",str(proj)],
                        cwd=proj, text=True, capture_output=True)
    assert r1.returncode == 0, r1.stderr + r1.stdout

    # step2
    r2 = subprocess.run([sys.executable,"-m","electrofit","step2","--project",str(proj)],
                        cwd=proj, text=True, capture_output=True)
    assert r2.returncode == 0, r2.stderr + r2.stdout

    # Inject minimal GMX core files (fixture) to ensure presence for step3
    run_dir = proj/"process"/"IP_011101"/"run_gmx_simulation"
    install_minimal_gmx_fixture(run_dir)

    # step3 (gmx shim will fabricate further files)
    r3 = subprocess.run([sys.executable,"-m","electrofit","step3","--project",str(proj)],
                        cwd=proj, text=True, capture_output=True)
    assert r3.returncode == 0, r3.stderr + r3.stdout

    # check a few expected outputs from the workflow
    assert (run_dir/"IP_011101_GMX.top").exists()
    assert any(p.name.endswith(".gro") for p in run_dir.iterdir())
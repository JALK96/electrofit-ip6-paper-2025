import os, sys, subprocess
from pathlib import Path
from tests.helpers.project import make_project_tree, install_minimal_gmx_fixture

def _prepare_project(root: Path):
    # create two molecule input folders by duplicating the scaffold
    make_project_tree(root, name="IP_011101")
    make_project_tree(root, name="IP_111111")

    # run step1 + step2 for both (step1/2 iterate over all under data/input/process)
    for step in ("step1","step2"):
        r = subprocess.run([sys.executable, "-m", "electrofit", step, "--project", str(root)],
                            cwd=root, text=True, capture_output=True)
        assert r.returncode == 0, r.stderr + r.stdout

    # inject minimal GMX fixture into both run_gmx_simulation dirs
    for mol in ("IP_011101","IP_111111"):
        run_dir = root/"process"/mol/"run_gmx_simulation"
        install_minimal_gmx_fixture(run_dir)


def test_step3_writes_isolated_logs(tmp_path, shim_bin, monkeypatch):
    # ensure shim path first
    monkeypatch.setenv("PATH", str(shim_bin)+os.pathsep+os.environ['PATH'])
    proj = tmp_path
    _prepare_project(proj)

    r3 = subprocess.run([sys.executable, "-m", "electrofit", "step3", "--project", str(proj)],
                        cwd=proj, text=True, capture_output=True)
    assert r3.returncode == 0, r3.stderr + r3.stdout

    log1 = (proj/"process"/"IP_011101"/"run_gmx_simulation"/"process.log").read_text().splitlines()
    log2 = (proj/"process"/"IP_111111"/"run_gmx_simulation"/"process.log").read_text().splitlines()

    # Heuristik: Letzte Zeile aus log2 darf nicht in log1 vorkommen, sonst Leakage
    if log2:
        assert log2[-1] not in log1, "Log leakage: molecule2 final line present in molecule1 log"

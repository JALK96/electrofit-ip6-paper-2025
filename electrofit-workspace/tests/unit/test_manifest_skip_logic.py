import os, sys, json, pathlib, subprocess
from tests.helpers.project import make_project_tree


def test_step3_skips_incomplete_manifest(tmp_path):
    proj = tmp_path
    make_project_tree(proj)

    # Manually craft run_gmx_simulation with incomplete manifest (missing gro value)
    run_dir = proj/'process'/'IP_011101'/'run_gmx_simulation'
    run_dir.mkdir(parents=True, exist_ok=True)
    # Create minimal MDP dir referenced by mdp_dir
    mdp = run_dir/'MDP'; mdp.mkdir(exist_ok=True)
    for f in ("em_steep.mdp","NVT.mdp","NPT.mdp","Production.mdp"):
        (mdp/f).write_text(';mdp')
    manifest = {
        "molecule": "IP_011101",
        "gro": "",  # intentionally blank
        "top": "IP_011101_GMX.top",
        "itp": "IP_011101_GMX.itp",
        "posres": "posre_IP_011101.itp",
        "mdp_dir": "MDP",
    }
    (run_dir/'run.json').write_text(json.dumps(manifest, indent=2))

    r = subprocess.run([sys.executable,'-m','electrofit','step3','--project',str(proj)], cwd=proj, text=True, capture_output=True)
    assert r.returncode == 0
    assert 'Skip: incomplete manifest' in r.stdout or 'Skip: GRO file missing' in r.stdout

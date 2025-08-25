# tests/helpers/project.py
import os, json
from pathlib import Path           # <-- add this line
from .mol2 import write_minimal_mol2
from shutil import copy2

def make_project_tree(root: Path, name: str = "IP_011101") -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    (root / "process").mkdir(exist_ok=True, parents=True)
    (root / "data" / "MDP").mkdir(parents=True, exist_ok=True)

    # dummy MDPs
    for f in ("em_steep.mdp", "NVT.mdp", "NPT.mdp", "Production.mdp"):
        (root / "data" / "MDP" / f).write_text("; dummy\n")

    # molecule package under data/input/<name>/
    mol_dir = root / "data" / "input" / name
    mol2 = write_minimal_mol2(mol_dir / f"{name}.mol2", title=name)
    (mol_dir / "equiv_groups.json").write_text(json.dumps({}, indent=2))
    (mol_dir / "electrofit.toml").write_text(
        f"""[project]
molecule_name = "{name}"
residue_name  = "LIG"
charge        = 0
protocol      = "bcc"
adjust_symmetry = true
ignore_symmetry = false
atom_type = "gaff2"

[paths]
mdp_dir = "data/MDP"
base_scratch_dir = "{root}/_scratch"

[simulation]
box_type = "dodecahedron"
box_edge_distance = 1.2
ion_concentration = 0.15
cation = "NA"
anion  = "CL"
forcefield = "amber14sb.ff"
"""
    )

    # project-wide defaults (kept minimal)
    (root / "electrofit.toml").write_text(
        f"""[paths]
mdp_dir = "data/MDP"
base_scratch_dir = "{root}/_scratch"
"""
    )

    # --- NEW: also create process/<name>/run_gau_create_gmx_in with mol2 and related files ---
    run_dir = root / "process" / name / "run_gau_create_gmx_in"
    run_dir.mkdir(parents=True, exist_ok=True)
    for fname in (f"{name}.mol2", "equiv_groups.json", "electrofit.toml"):
        src = mol_dir / fname
        if src.exists():
            (run_dir / fname).write_text(src.read_text())

    return {
        "root": root,
        "mol_dir": mol_dir,
        "mol2": mol2,
        "run_dir": run_dir,
    }

FIXTURE_MIN_GMX = Path(__file__).resolve().parent.parent / "fixtures" / "minimal_gmx"

def install_minimal_gmx_fixture(run_sim_dir: Path):
    """Copy minimal GMX core files into an existing run_gmx_simulation directory.

    Used in tests to avoid patching production code. Assumes directory exists.
    """
    needed = [
        "IP_011101_GMX.gro",
        "IP_011101_GMX.itp",
        "IP_011101_GMX.top",
        "posre_IP_011101.itp",
    ]
    for name in needed:
        src = FIXTURE_MIN_GMX / name
        if src.is_file():
            copy2(src, run_sim_dir / name)


def install_fake_gmxdata(root: Path, ff_name: str = "amber14sb.ff") -> Path:
    """Create a minimal fake GMXDATA layout so forcefield validation passes in tests.

    Returns the GMXDATA root path that should be exported via the GMXDATA env var.
    """
    gmxdata = root / "fake_gmxdata"
    top = gmxdata / "top" / ff_name
    top.mkdir(parents=True, exist_ok=True)
    # Optional placeholder files (not required for validation, but useful for realism)
    for fname in ("forcefield.itp", "tip3p.itp", "ions.itp"):
        fpath = top / fname
        if not fpath.exists():
            fpath.write_text("; placeholder for tests\n")
    return gmxdata
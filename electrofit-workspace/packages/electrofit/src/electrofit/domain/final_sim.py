"""Domain functions for final simulation setup & launch (Steps7 & 8)."""
from __future__ import annotations

from pathlib import Path
import fnmatch
import logging
import os
import shutil
from typing import Tuple, Iterable

from electrofit.config.loader import load_config, dump_config
from electrofit.infra.config_snapshot import compose_snapshot
from electrofit.adapters.gromacs import set_up_production

__all__ = [
    "prepare_final_sim_directory",
    "iter_final_sim_dirs",
    "launch_final_sim_run",
]


def prepare_final_sim_directory(
    mol_dir: Path,
    project_root: Path,
    override_cfg: Path | None,
    multi_mol: bool,
) -> Tuple[bool, str]:
    """Create /run_final_gmx_simulation for a molecule and copy topology + MDP.

    Returns (ok, message) for backwards compatible console/status logging.
    """
    results_dir = mol_dir / "results"
    if not results_dir.is_dir():
        return False, "no results dir"
    dest_dir = mol_dir / "run_final_gmx_simulation"
    dest_dir.mkdir(exist_ok=True)
    acpype_dir = next((p for p in results_dir.iterdir() if p.is_dir() and p.name.endswith(".acpype")), None)
    if not acpype_dir:
        return False, "no .acpype dir in results"
    compose_snapshot(
        dest_dir,
        project_root,
        mol_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        upstream=results_dir / "electrofit.toml",
        process_cfg=mol_dir / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol_dir.name / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )
    cfg = load_config(project_root, context_dir=dest_dir, molecule_name=mol_dir.name)
    try:
        dump_config(cfg, log_fn=logging.info)
    except Exception:  # pragma: no cover
        logging.debug("[step7] config dump failed", exc_info=True)
    file_patterns = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]
    for fn in os.listdir(acpype_dir):
        for pat in file_patterns:
            if fnmatch.fnmatch(fn, pat):
                shutil.copy(acpype_dir / fn, dest_dir)
                logging.info(f"[step7] Copied {fn} -> {dest_dir}")
                break
    mdp_source_dir = project_root / "data" / "MDP"
    md_dest = dest_dir / "MDP"
    if mdp_source_dir.is_dir():
        if md_dest.exists():
            shutil.rmtree(md_dest)
        shutil.copytree(mdp_source_dir, md_dest)
        logging.info(f"[step7] Copied MDP -> {md_dest}")
    else:
        logging.warning(f"[step7][warn] no MDP source {mdp_source_dir}")
    bash_script_source = project_root / "scripts" / "gmx.sh"
    if bash_script_source.is_file():
        shutil.copy(bash_script_source, dest_dir / "gmx.sh")
    return True, "ok"


def iter_final_sim_dirs(project_root: Path) -> Iterable[Path]:
    process = project_root / "process"
    if not process.is_dir():
        return []
    for mol in process.iterdir():
        run_dir = mol / "run_final_gmx_simulation"
        if run_dir.is_dir():
            yield run_dir


def launch_final_sim_run(
    run_dir: Path,
    project_root: Path,
    multi_mol: bool,
    override_cfg: Path | None,
) -> Tuple[bool, str]:
    """Launch final production simulation in an existing prepared directory."""
    mol = run_dir.parent.name
    compose_snapshot(
        run_dir,
        project_root,
        mol,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        upstream=run_dir.parent / "results" / "electrofit.toml",
        process_cfg=run_dir.parent / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )
    cfg = load_config(project_root, context_dir=run_dir, molecule_name=mol)
    dump_config(cfg, log_fn=logging.info)
    gro = next((f.name for f in run_dir.iterdir() if f.suffix == ".gro"), None)
    if not gro:
        logging.info(f"[step8][skip] {run_dir}: no .gro file")
        return False, "no gro"
    sim = cfg.simulation
    box = sim.box
    ions = sim.ions
    runtime = cfg.gmx.runtime
    ff = getattr(sim, 'forcefield', None) or "amber14sb.ff"
    import os as _os
    prev = _os.getcwd()
    try:
        _os.chdir(run_dir)
        set_up_production(
            m_gro=gro,
            MDP_dir="MDP",
            base_scratch_dir=cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch",
            molecule_name=cfg.project.molecule_name or mol,
            box_type=box.type,
            cation=ions.cation,
            anion=ions.anion,
            d=str(box.edge_nm),
            conc=str(ions.concentration),
            ff=ff,
            threads=runtime.threads,
            pin=runtime.pin,
        )
    finally:
        _os.chdir(prev)
    return True, "ok"

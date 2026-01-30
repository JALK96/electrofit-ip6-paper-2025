"""Domain functions for REMD simulation setup & launch (Steps 7remd & 8remd)."""
from __future__ import annotations

from pathlib import Path
import fnmatch
import logging
import os
import shutil
from typing import Iterable, Tuple

from electrofit.config.loader import load_config, dump_config
from electrofit.infra.config_snapshot import compose_snapshot
from electrofit.io.ff import validate_forcefield, ensure_forcefield_installed

__all__ = [
    "prepare_remd_sim_directory",
    "iter_remd_sim_dirs",
    "launch_remd_run",
]


def prepare_remd_sim_directory(
    mol_dir: Path,
    project_root: Path,
    override_cfg: Path | None,
    multi_mol: bool,
) -> Tuple[bool, str]:
    """Create /run_remd_gmx_simulation for a molecule and copy topology + MDP.

    Does not run GROMACS. Builds a per-run config snapshot next to the directory.
    """
    results_dir = mol_dir / "results"
    if not results_dir.is_dir():
        return False, "no results dir"
    dest_dir = mol_dir / "run_remd_gmx_simulation"
    dest_dir.mkdir(exist_ok=True)

    # Compose run-time config snapshot into dest_dir
    compose_snapshot(
        dest_dir,
        project_root,
        mol_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        step="step7remd",
        upstream=results_dir / "electrofit.toml",
        process_cfg=mol_dir / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol_dir.name / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )

    # Copy ACPYPE ligand files into the run dir
    acpype_dir = next((p for p in results_dir.iterdir() if p.is_dir() and p.name.endswith(".acpype")), None)
    if not acpype_dir:
        return False, "no .acpype dir in results"
    file_patterns = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]
    for fn in os.listdir(acpype_dir):
        for pat in file_patterns:
            if fnmatch.fnmatch(fn, pat):
                shutil.copy(acpype_dir / fn, dest_dir)
                logging.info(f"[step7remd] Copied {fn} -> {dest_dir}")
                break

    # Copy MDP directory
    mdp_source_dir = project_root / "data" / "MDP"
    md_dest = dest_dir / "MDP"
    if mdp_source_dir.is_dir():
        if md_dest.exists():
            shutil.rmtree(md_dest)
        shutil.copytree(mdp_source_dir, md_dest)
        logging.info(f"[step7remd] Copied MDP -> {md_dest}")
    else:
        logging.warning(f"[step7remd][warn] no MDP source {mdp_source_dir}")

    return True, "ok"


def iter_remd_sim_dirs(project_root: Path) -> Iterable[Path]:
    process = project_root / "process"
    if not process.is_dir():
        return []
    for mol in process.iterdir():
        run_dir = mol / "run_remd_gmx_simulation"
        if run_dir.is_dir():
            yield run_dir


def launch_remd_run(
    run_dir: Path,
    project_root: Path,
    multi_mol: bool,
    override_cfg: Path | None,
    *,
    prepare_only: bool = False,
) -> Tuple[bool, str]:
    """Equilibrate and generate REMD inputs; optionally launch REMD.

    Returns (ok, message).
    """
    mol = run_dir.parent.name
    compose_snapshot(
        run_dir,
        project_root,
        mol,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        step="step8remd",
        upstream=run_dir.parent / "results" / "electrofit.toml",
        process_cfg=run_dir.parent / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )
    cfg = load_config(project_root, context_dir=run_dir, molecule_name=mol)
    dump_config(cfg, log_fn=logging.debug)

    gro = next((f.name for f in run_dir.iterdir() if f.suffix == ".gro"), None)
    if not gro:
        logging.info(f"[step8remd][skip] {run_dir}: no .gro file")
        return False, "no gro"

    # Resolve forcefield specification, supporting both short names
    # (e.g. 'amber14sb.ff') and explicit paths (installed if needed),
    # mirroring the behaviour used in step3.
    sim = getattr(cfg, "simulation", None)
    ff_attr = getattr(sim, "forcefield", None) if sim else None
    ff_raw = ff_attr or "amber14sb.ff"

    if isinstance(ff_raw, str) and (
        os.path.isabs(ff_raw) or os.sep in ff_raw or ff_raw.startswith(".")
    ):
        # Path-like: ensure it is installed under the active GMX top and
        # then use the short directory name for includes.
        ff_abs = ff_raw if os.path.isabs(ff_raw) else os.path.join(str(project_root), ff_raw)
        try:
            ff = ensure_forcefield_installed(ff_abs)
        except Exception as e:
            msg = f"[step8remd][abort] {e}"
            logging.error(msg)
            return False, msg
        explicit_ff = True
    else:
        ff = ff_raw
        explicit_ff = ff_attr is not None

    # Validate that the resolved forcefield directory exists so we fail
    # early with a clear diagnostic instead of a grompp include error.
    try:
        validate_forcefield(ff, explicit=explicit_ff)
    except FileNotFoundError as e:
        msg = f"[step8remd][abort] {e}"
        logging.error(msg)
        return False, msg

    from electrofit.adapters.gromacs_remd import set_up_remd

    import os as _os
    # Capture previous CWD defensively: in some HPC setups getcwd() can fail
    # (e.g., if the starting directory was removed). In that case we simply
    # avoid restoring it and stay in run_dir.
    try:
        prev = _os.getcwd()
        restore_prev = True
    except FileNotFoundError:  # pragma: no cover - environment dependent
        prev = None
        restore_prev = False

    try:
        _os.chdir(run_dir)
        set_up_remd(
            m_gro=gro,
            MDP_dir="MDP",
            base_scratch_dir=cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch",
            molecule_name=cfg.project.molecule_name or mol,
            simulation=cfg.simulation,
            remd_cfg=cfg.remd,
            ff=ff,
            threads=cfg.gmx.runtime.threads,
            pin=cfg.gmx.runtime.pin,
            gpu=cfg.gmx.runtime.gpu,
        )
    finally:
        if restore_prev and prev is not None:
            try:
                _os.chdir(prev)
            except FileNotFoundError:  # pragma: no cover
                pass

    # Optionally launch the REMD run now via the generated helper script
    if not prepare_only and cfg.remd.enabled:
        from electrofit.cli.run_commands import run_command
        run_command(["bash", "run_remd.sh"], cwd=run_dir)

    return True, "ok"

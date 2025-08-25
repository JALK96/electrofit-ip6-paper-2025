"""Pipeline Step 7: Prepare final GROMACS production simulation directories.

Thin orchestration wrapper that iterates molecule directories and invokes
`prepare_final_sim_directory`. All heavy lifting lives in the domain layer.
"""
from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

from electrofit.domain.final_sim import prepare_final_sim_directory
from electrofit.infra.logging import setup_logging, reset_logging, log_run_header, enable_header_dedup
from electrofit.infra.config_snapshot import CONFIG_ARG_HELP
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

__all__ = ["main", "run_step7"]


def run_step7(project: Path, override_cfg: Path | None, only_molecule: str | None = None) -> int:
    process_dir = project / "process"
    if not process_dir.is_dir():
        print("[step7] No process directory.")
        return 0
    mol_dirs_all = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    mol_dirs = filter_paths_for_molecule(mol_dirs_all, only_molecule)
    if only_molecule and not mol_dirs:
        print(f"[step7][warn] molecule '{only_molecule}' not found; nothing to do")
        return 0
    multi_mol = len(mol_dirs) > 1
    done = 0
    for mol_dir in mol_dirs:
        dest_dir = mol_dir / "run_final_gmx_simulation"
        log_path = dest_dir / "process.log"
        dest_dir.mkdir(exist_ok=True)
        try:
            reset_logging()
            setup_logging(str(log_path), also_console=False)
            log_run_header("step8")
        except Exception:  # pragma: no cover
            pass
        ok, msg = prepare_final_sim_directory(mol_dir, project, override_cfg, multi_mol)
        status = "[step7]" if ok else "[step7][skip]"
        print(f"{status} {mol_dir.name}: {msg}")
        if ok:
            done += 1

    summary = f"[step7] Prepared final sim dirs for {done}/{len(mol_dirs)} molecules."
    print(summary)
    try:
        reset_logging()
        setup_logging(str(project / "step.log"), also_console=True, suppress_initial_message=True)
        logging.info(summary)
    except Exception:
        pass
    return 0


def main(argv: list[str] | None = None):  # pragma: no cover
    ap = argparse.ArgumentParser(description="Step7: Prepare final simulation directories (copy topology, MDP, scripts).")
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--molecule", help="Prepare only this molecule name")
    args = ap.parse_args(argv)
    project_root = Path(args.project).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    # Global Header vor Ausführung
    try:
        enable_header_dedup(True)
        setup_logging(str(project_root / "step.log"), also_console=True, suppress_initial_message=True)
        log_run_header("step7")
    except Exception:
        pass
    rc = run_step7(project_root, override_cfg, getattr(args, 'molecule', None))
    raise SystemExit(rc)


if __name__ == "__main__":  # pragma: no cover
    main()

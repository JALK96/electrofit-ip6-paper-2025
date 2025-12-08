"""Pipeline Step 8 (REMD): Launch REMD preparation and optionally run REMD.

Enumerates prepared REMD directories and dispatches domain.remd_sim.launch_remd_run.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import logging

from electrofit.domain.remd_sim import iter_remd_sim_dirs, launch_remd_run
from electrofit.infra.logging import (
    setup_logging,
    reset_logging,
    log_run_header,
    enable_header_dedup,
)
from electrofit.infra.config_snapshot import CONFIG_ARG_HELP
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

__all__ = ["main", "run_step8remd"]


def run_step8remd(project: Path, override_cfg: Path | None, log_console: bool, only_molecule: str | None = None, prepare_only: bool = False) -> int:
    enable_header_dedup(True)
    setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
    log_run_header("step8remd")
    process_root = project / "process"
    multi_mol = False
    if process_root.is_dir():
        pmols = [p for p in process_root.iterdir() if p.is_dir()]
        multi_mol = len(pmols) > 1
    ran = 0
    all_dirs = list(iter_remd_sim_dirs(project))
    run_dirs = filter_paths_for_molecule(all_dirs, only_molecule)
    if only_molecule and not run_dirs:
        print(f"[step8remd][warn] molecule '{only_molecule}' not found; nothing to run")
        return 0
    for run_dir in run_dirs:
        reset_logging()
        setup_logging(str(run_dir / "process.log"), also_console=False)
        log_run_header("step8remd")
        ok, msg = launch_remd_run(run_dir, project, multi_mol, override_cfg, prepare_only=prepare_only)
        if ok:
            ran += 1
        else:
            # Mirror step3 behaviour: print per-run abort/skip messages to console
            # so configuration issues (e.g. missing forcefield) are immediately visible.
            print(msg)
    summary = f"[step8remd] Completed preparation{' and launch' if not prepare_only else ''} for {ran} REMD run(s)."
    print(summary)
    reset_logging()
    setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
    logging.info(summary)
    return 0


def main(argv: list[str] | None = None):  # pragma: no cover
    ap = argparse.ArgumentParser(description="Step8 (REMD): Prepare and optionally launch REMD runs.")
    # Top-level CLI sets ELECTROFIT_PROJECT_PATH and passes --project explicitly;
    # avoid os.getcwd() in the default to prevent FileNotFoundError when CWD is invalid.
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH"))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    ap.add_argument("--molecule", help="Only run the specified molecule name", required=False)
    ap.add_argument("--prepare-only", action="store_true", help="Prepare inputs without launching mdrun")
    args = ap.parse_args(argv)
    project_value = args.project or os.environ.get("ELECTROFIT_PROJECT_PATH")
    if not project_value:
        raise SystemExit("--project must be provided (or ELECTROFIT_PROJECT_PATH set)")
    project_root = Path(project_value).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    rc = run_step8remd(project_root, override_cfg, args.log_console, getattr(args, "molecule", None), prepare_only=args.prepare_only)
    raise SystemExit(rc)


if __name__ == "__main__":  # pragma: no cover
    main()

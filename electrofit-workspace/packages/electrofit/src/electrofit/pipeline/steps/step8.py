"""Pipeline Step 8: Launch final production GROMACS simulations.

Thin orchestrator: enumerates prepared final simulation directories and
dispatches `launch_final_sim_run` for each.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import logging

from electrofit.domain.final_sim import iter_final_sim_dirs, launch_final_sim_run
from electrofit.infra.logging import (
    setup_logging,
    reset_logging,
    log_run_header,
    enable_header_dedup,
)
from electrofit.infra.config_snapshot import CONFIG_ARG_HELP
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

__all__ = ["main", "run_step8"]


def run_step8(project: Path, override_cfg: Path | None, log_console: bool, only_molecule: str | None = None) -> int:
    enable_header_dedup(True)
    setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
    log_run_header("step8")
    process_root = project / "process"
    multi_mol = False
    if process_root.is_dir():
        pmols = [p for p in process_root.iterdir() if p.is_dir()]
        multi_mol = len(pmols) > 1
    ran = 0
    all_dirs = list(iter_final_sim_dirs(project))
    run_dirs = filter_paths_for_molecule(all_dirs, only_molecule)
    if only_molecule and not run_dirs:
        print(f"[step8][warn] molecule '{only_molecule}' not found; nothing to run")
        return 0
    for run_dir in run_dirs:
        reset_logging()
        setup_logging(str(run_dir / "process.log"), also_console=False)
        log_run_header("step8")
        ok, msg = launch_final_sim_run(run_dir, project, multi_mol, override_cfg)
        if ok:
            ran += 1
    summary = f"[step8] Completed {ran} final simulation run(s)."
    print(summary)
    reset_logging()
    setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
    logging.info(summary)
    return 0


def main(argv: list[str] | None = None):  # pragma: no cover
    ap = argparse.ArgumentParser(description="Step8: Launch final GROMACS production runs.")
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    ap.add_argument("--molecule", help="Only run the specified molecule name", required=False)
    args = ap.parse_args(argv)
    project_root = Path(args.project).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    rc = run_step8(project_root, override_cfg, args.log_console, getattr(args, "molecule", None))
    raise SystemExit(rc)


if __name__ == "__main__":  # pragma: no cover
    main()

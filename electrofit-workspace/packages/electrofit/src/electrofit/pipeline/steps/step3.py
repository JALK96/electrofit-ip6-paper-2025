"""Pipeline Step 3: Set up production GROMACS runs from manifests.

Reads each ``process/<mol>/run_gmx_simulation`` directory, re-seeds its snapshot,
loads per-run config, and invokes the GROMACS adapter to generate production inputs.
"""

from __future__ import annotations
import argparse
import json
import os
import logging
from pathlib import Path
from electrofit.config.loader import load_config, dump_config
from electrofit.adapters.gromacs import set_up_production
from electrofit.io.ff import validate_forcefield
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, reset_logging, log_run_header, enable_header_dedup
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

__all__ = ["main"]


def _iter_run_dirs(project_root: Path):
    process = project_root / "process"
    if not process.is_dir():
        return
    for mol_dir in process.iterdir():
        run_dir = mol_dir / "run_gmx_simulation"
        if run_dir.is_dir():
            yield run_dir


def main():  # pragma: no cover
    ap = argparse.ArgumentParser(
        description="Step3: Setup GROMACS production directories from run.json + config"
    )
    ap.add_argument(
        "--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
    )
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--molecule", help="Only run for this molecule name")
    args = ap.parse_args()
    project_root = Path(args.project).resolve()
    process_root = project_root / "process"
    multi_mol = (
        process_root.is_dir()
        and len([p for p in process_root.iterdir() if p.is_dir()]) > 1
    )
    ran = 0
    # Enable header de-duplication and emit initial global header (single header policy)
    try:
        enable_header_dedup(True)
        setup_logging(str(project_root / "step.log"), also_console=True, suppress_initial_message=True)
        log_run_header("step3")
    except Exception:
        pass
    run_dirs_all = list(_iter_run_dirs(project_root))
    run_dirs = filter_paths_for_molecule(run_dirs_all, args.molecule)
    if args.molecule and not run_dirs:
        print(f"[step3][warn] molecule '{args.molecule}' not found; nothing to do")
        return
    for run_dir in run_dirs:
        reset_logging()
        setup_logging(str(run_dir / "process.log"), also_console=False)
        try:
            log_run_header("step3")
        except Exception:
            logging.info("electrofit unknown | step=step3")
        mol = run_dir.parent.name
        upstream = run_dir.parent / "run_gau_create_gmx_in" / "electrofit.toml"
        molecule_input = project_root / "data" / "input" / mol / "electrofit.toml"
        process_cfg = run_dir.parent / "electrofit.toml"
        project_defaults = project_root / "electrofit.toml"
        compose_snapshot(
            run_dir,
            project_root,
            mol,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=upstream,
            process_cfg=process_cfg,
            molecule_input=molecule_input,
            project_defaults=project_defaults,
            extra_override=Path(args.config) if args.config else None,
        )
        cfg_run = load_config(project_root, context_dir=run_dir, molecule_name=mol)
        try:
            dump_config(cfg_run, log_fn=logging.info)
        except Exception:
            logging.debug("[step3] dump per-run config failed", exc_info=True)

        # Pull simulation knobs from run config (fallback to defaults if absent)
        sim = getattr(cfg_run, "simulation", None)
        ff_attr = getattr(sim, "forcefield", None)
        ff       = ff_attr or "amber14sb.ff"
        box_type = getattr(sim, "box_type", None)           or "dodecahedron"
        d        = getattr(sim, "box_edge_distance", None)  or 1.2
        cation   = getattr(sim, "cation", None)             or "NA"
        anion    = getattr(sim, "anion", None)              or "CL"
        conc     = getattr(sim, "ion_concentration", None)  or 0.15
        runtime_local = getattr(getattr(cfg_run, "gmx", None), "runtime", None)
        threads     = getattr(runtime_local, "threads", None) if runtime_local else None
        pin         = getattr(runtime_local, "pin", None) if runtime_local else None
        base_scratch_run = (
            getattr(getattr(cfg_run, "paths", None), "base_scratch_dir", None)
            or "/tmp/electrofit_scratch"
        )
        manifest = run_dir / "run.json"
        if not manifest.is_file():
            logging.info(f"[step3][skip] no run.json in {run_dir}")
            print(f"[step3][skip] no run.json in {run_dir}")
            continue
        meta = json.loads(manifest.read_text())
        molecule = meta.get("molecule")
        gro_file = meta.get("gro")
        mdp_dir = meta.get("mdp_dir") or (getattr(cfg_run.paths, "mdp_dir", None) or "MDP")
        if not molecule or not gro_file:
            msg = f"[step3][skip] Skip: incomplete manifest ({manifest})"
            logging.info(msg)
            print(msg)
            continue
        # Ensure files exist in the run directory
        gro_path = run_dir / gro_file
        mdp_path = run_dir / mdp_dir
        top_path = gro_path.with_suffix(".top")
        if not gro_path.is_file():
            msg = f"[step3][skip] Skip: GRO file missing {gro_path}"
            logging.info(msg)
            print(msg)
            continue
        if not top_path.is_file():
            msg = f"[step3][skip] Skip: TOP file missing {top_path}"
            logging.info(msg)
            print(msg)
            continue
        if not mdp_path.is_dir():
            msg = f"[step3][skip] MDP dir missing: {mdp_path}"
            logging.info(msg)
            print(msg)
            continue
        if ff_attr is None:
            logging.warning(
                "[step3][ff] No forcefield specified in config; falling back to default '%s'", ff
            )
        # Validate (raises if user explicitly set a non-existent ff; warns otherwise)
        try:
            validate_forcefield(ff, explicit=ff_attr is not None)
        except FileNotFoundError as e:
            logging.error("[step3][abort] %s", e)
            print(f"[step3][abort] {e}")
            continue
        logging.info(f"[step3] Starting GROMACS production for {molecule} in {run_dir}")
        prev = os.getcwd()
        try:
            os.chdir(run_dir)
            set_up_production(
                m_gro=gro_file,
                MDP_dir=mdp_dir,
                base_scratch_dir=base_scratch_run,
                molecule_name=molecule,
                box_type=box_type,
                cation=cation,
                anion=anion,
                d=str(d),
                conc=str(conc),
                ff=ff,
                threads=threads,
                pin=pin,
            )
            ran += 1

        finally:
            os.chdir(prev)
    summary = (
        f"[step3] Completed {ran} run(s)."
        if ran
        else "[step3] No runs executed (no manifests found)."
    )
    print(summary)
    reset_logging()
    setup_logging(str(project_root / "step.log"), also_console=True, suppress_initial_message=True)
    logging.info(summary)


if __name__ == "__main__":  # pragma: no cover
    main()

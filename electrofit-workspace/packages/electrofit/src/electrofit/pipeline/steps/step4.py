"""Pipeline Step 4: Extract conformers from production GROMACS trajectories.

Orchestration only:
    * Parse CLI arguments (sampling overrides, workers, clean mode).
    * Derive sampling defaults from layered project / override TOMLs.
    * Iterate molecule process directories and invoke per‑molecule extraction.
    * Manage optional parallel execution and progress reporting.
    * Emit per‑molecule logging into that molecule's ``extracted_conforms/process.log``.
    * Summarise results to stdout and project-level ``step.log``.

Heavy domain logic (frame selection, snapshot layering, file copies) will
move into a dedicated service module in a future refactor to further reduce
side‑effects at import time.
"""
from __future__ import annotations
import argparse
import os
import shutil
import multiprocessing
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import mdtraj as md  # noqa: F401  (imported for side-effects: topology/trajectory readers registration)
from tqdm import tqdm

try:  # Python 3.11+
    import tomllib as _toml  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml  # type: ignore

from electrofit.infra.config_snapshot import CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, reset_logging, log_run_header
from electrofit.pipeline.workers.step4_extract import _extract_for_molecule  # type: ignore
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule
from electrofit.pipeline.workers.step4_worker import worker as _worker  # type: ignore

"""Multiprocessing note: executed via runpy -> cannot rely on __main__ pickling.

We avoid fragile __module__ rewrites by moving the actual worker function into
`electrofit.pipeline.workers.step4_worker`. That module is imported by its
fully qualified name both in the parent and child processes, guaranteeing a
stable import path for pickle. This keeps this orchestrator file focused and
minimises side effects at import time.
"""

 


__all__ = ["main"]


### Extraction helper moved to `pipeline.workers.step4_extract` to avoid circular import.


# _worker kommt jetzt aus eigenem Modul; keine lokale Definition nötig.


def main():  # pragma: no cover
    ap = argparse.ArgumentParser(
        description=(
            "Step4: Extract conformers from GROMACS trajectories (md_center.xtc + md.gro).\n"
            "Sampling methods: linear, random, maxmin. Defaults read from layered electrofit.toml unless overridden."
        )
    )
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--sample", type=int, default=None, help="Number of conformers to sample (override)")
    ap.add_argument("--sampling-method", default=None, help="Sampling method: linear|random|maxmin (override)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for sampling (override)")
    ap.add_argument("--workers", type=int, default=0, help="Parallel workers (0=auto, 1=sequential)")
    ap.add_argument("--no-progress", action="store_true", help="Disable progress bar output")
    ap.add_argument("--clean", action="store_true", help="Remove existing extracted_conforms dirs before extraction")
    ap.add_argument("--molecule", help="Process only this molecule (directory under process/)")
    ap.add_argument("--verbose", action="store_true", help="Verbose bond insertion logs")
    args = ap.parse_args()

    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    project_root = Path(args.project).resolve()
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        print("No process directory found; nothing to do.")
        return
    mol_dirs_all = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    mol_dirs = filter_paths_for_molecule(mol_dirs_all, args.molecule)
    if args.molecule and not mol_dirs:
        print(f"[step4][warn] molecule '{args.molecule}' not found; nothing to do")
        return
    multi_mol = len(mol_dirs) > 1
    if not mol_dirs:
        print("[step4] No molecule directories found.")
        return

    sampling_cfg = {}
    candidate_cfgs: list[str] = []
    if getattr(args, "config", None):
        candidate_cfgs.append(args.config)
    candidate_cfgs.append(str(project_root / "electrofit.toml"))
    for c in candidate_cfgs:
        p = Path(c)
        if p.is_file():
            try:
                with p.open("rb") as fh:
                    data = _toml.load(fh)
                sampling_cfg = data.get("sampling", {}) or {}
            except Exception:
                sampling_cfg = {}
            break
    method = args.sampling_method or sampling_cfg.get("method") or "linear"
    sample_count = args.sample or sampling_cfg.get("count") or 20
    seed = args.seed if args.seed is not None else sampling_cfg.get("seed")

    if args.clean:
        for mdir in mol_dirs:
            ec_dir = mdir / "extracted_conforms"
            if ec_dir.is_dir():
                try:
                    shutil.rmtree(ec_dir)
                    print(f"[step4][clean] removed {ec_dir}")
                except Exception as e:
                    print(f"[step4][clean][warn] failed to remove {ec_dir}: {e}")

    if args.workers <= 0:
        cpu_count = max(1, multiprocessing.cpu_count() - 1)
        workers = min(cpu_count, len(mol_dirs))
    else:
        workers = min(args.workers, len(mol_dirs))

    extracted = 0
    results = []
    if workers == 1:
        iterator = mol_dirs
        if not args.no_progress:
            iterator = tqdm(iterator, desc="step4", unit="mol")
        for sub in iterator:
            ok, msg = _extract_for_molecule(sub, project_root, sample_count, method, seed, override_cfg, multi_mol, args.verbose)
            results.append((sub.name, ok, msg, None))
            if ok:
                extracted += 1

    else:
        tasks = [
            (str(p), str(project_root), sample_count, method, seed, override_cfg, multi_mol, args.verbose)
            for p in mol_dirs
        ]
        pbar = None
        if not args.no_progress:
            pbar = tqdm(total=len(mol_dirs), desc="step4", unit="mol")
        with ProcessPoolExecutor(max_workers=workers) as ex:
            future_map = {ex.submit(_worker, t): t[0] for t in tasks}
            for fut in as_completed(future_map):
                name, ok, msg, err = fut.result()
                results.append((Path(name).name, ok, msg, err))
                if ok:
                    extracted += 1
                if pbar:
                    pbar.update(1)
        if pbar:
            pbar.close()

    for name, ok, msg, err in sorted(results, key=lambda x: x[0]):
        prefix = "[step4]" if ok else "[step4][skip]"
        detail = msg if err is None else f"{msg}: {err}"
        print(f"{prefix} {name}: {detail}")
    summary = f"[step4] Conformer extraction complete for {extracted}/{len(mol_dirs)} molecules." if extracted else "[step4] No conformers extracted."
    print(summary)

    reset_logging()
    setup_logging(str(project_root / "step.log"), also_console=True, suppress_initial_message=True)
    try:
        log_run_header("step4")
    except Exception:
        logging.info("electrofit unknown | step=step4")
    logging.info(summary)


if __name__ == "__main__":  # pragma: no cover
    main()

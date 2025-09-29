"""
Run Electrofit pipeline steps via direct imports (not the CLI).

This script wires the step orchestrations programmatically so you can:
  - import and call run_pipeline(project, steps=...)
  - or execute it directly with a lightweight CLI (still using direct imports under the hood)

Covered steps: 0–8

Notes:
  - We intentionally avoid calling the step modules' CLI entry-points.
  - For steps that already expose helpers (run_step5/6/7/8), we call them directly.
  - For steps 0–4 we re-use internal helpers or inline their orchestration code.
  - We set ELECTROFIT_PROJECT_PATH to the given project so downstream modules behave consistently.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
import logging
from typing import Iterable, List, Tuple


def _ensure_import_path() -> None:
    """If the electrofit package isn't importable, add the local src path."""
    try:
        import electrofit  # noqa: F401
        return
    except Exception:
        pass
    # Try local editable layout: packages/electrofit/src
    repo_root = Path(__file__).resolve().parents[1]
    cand = repo_root / "packages" / "electrofit" / "src"
    if cand.is_dir():
        sys.path.insert(0, str(cand))


_ensure_import_path()


# Common env setup
def _set_project_env(project: Path) -> None:
    os.environ["ELECTROFIT_PROJECT_PATH"] = str(project)


# ---- Step 0 (inline from step0 logic) ---------------------------------------------------------
def run_step0(project: Path, molecule: str | None = None, log_console: bool = False) -> None:
    """Create process/ by copying data/input; optionally only one molecule.

    This mirrors step0 without invoking its CLI.
    """
    from electrofit.io.files import copy_and_rename_folders
    from electrofit.infra.logging import setup_logging, log_run_header

    project = project.resolve()
    _set_project_env(project)
    src = project / "data" / "input"
    dst = project / "process"
    dst.mkdir(exist_ok=True)
    setup_logging(str(project / "step.log"), also_console=log_console)
    log_run_header("step0")
    logging.info("[step0] Source: %s", src)
    logging.info("[step0] Destination: %s", dst)
    if not src.exists():
        raise FileNotFoundError(f"Source directory '{src}' does not exist.")
    if molecule:
        # Copy only a single molecule (mirrors step0's single-molecule code, simplified)
        mol = molecule
        sub = src / mol
        if not sub.is_dir():
            raise SystemExit(f"[step0][abort] requested molecule '{mol}' not found under {src}")
        nested = dst / mol / "run_gau_create_gmx_in"
        nested.mkdir(parents=True, exist_ok=True)
        from shutil import copy2, copytree
        for item in sub.iterdir():
            dest_path = nested / item.name
            if item.is_file():
                copy2(item, dest_path)
            elif item.is_dir():
                copytree(item, dest_path, dirs_exist_ok=True)
        logging.info("[step0] Copied molecule '%s' -> %s", mol, nested)
    else:
        copy_and_rename_folders(source=str(src), destination=str(dst))
    logging.info("[step0] Copy complete.")


# ---- Step 1 (reuse internal helpers) ----------------------------------------------------------
def run_step1(project: Path, config: Path | None = None, molecule: str | None = None) -> None:
    from electrofit.pipeline.steps import step1 as s1
    from electrofit.infra.logging import setup_logging, log_run_header, reset_logging, enable_header_dedup

    project = project.resolve()
    _set_project_env(project)
    process_root = project / "process"
    multi_mol = process_root.is_dir() and len([p for p in process_root.iterdir() if p.is_dir()]) > 1

    try:
        enable_header_dedup(True)
        setup_logging(str(project / "step.log"), also_console=False, suppress_initial_message=True)
        log_run_header("step1")
    except Exception:
        pass
    discovered = list(s1._iter_run_dirs(str(project), only_molecule=molecule))  # type: ignore[attr-defined]
    if molecule and not discovered:
        print(f"[step1][warn] requested molecule '{molecule}' not found / no run dirs")
        return
    if not discovered:
        # Fallback: assume current dir is a run dir
        discovered = [os.getcwd()]
    ran = 0
    for run_dir in discovered:
        s1._run_one_dir(run_dir, str(project), str(config) if config else None, multi_mol)  # type: ignore[attr-defined]
        ran += 1
        reset_logging()
    try:
        setup_logging(str(project / "step.log"), also_console=False, suppress_initial_message=True)
        logging.info("[step1] Completed %d run directory(ies).", ran)
    except Exception:
        pass


# ---- Step 2 (inline from step2) --------------------------------------------------------------
def run_step2(project: Path, config: Path | None = None, molecule: str | None = None, log_console: bool = False) -> None:
    import fnmatch
    import json
    import shutil
    from electrofit.infra.config_snapshot import compose_snapshot
    from electrofit.infra.logging import setup_logging, log_run_header, reset_logging, enable_header_dedup
    from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

    project = project.resolve()
    _set_project_env(project)
    process_dir = project / "process"
    mdp_source_dir = project / "data" / "MDP"
    if not process_dir.is_dir():
        print("[step2] No process directory.")
        return
    mol_dirs_all = [p for p in process_dir.iterdir() if p.is_dir()]
    mol_dirs = filter_paths_for_molecule(mol_dirs_all, molecule)
    if molecule and not mol_dirs:
        print(f"[step2][warn] molecule '{molecule}' not found; nothing to do")
        return
    multi_mol = len(mol_dirs) > 1
    enable_header_dedup(True)
    setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
    log_run_header("step2")
    _FILE_PATTERNS = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]
    def _write_manifest(dest_dir: Path, files: dict[str,str], mdp_subdir: str = "MDP") -> None:
        gro = files.get("gro")
        molecule_name = None
        if gro:
            base = os.path.splitext(os.path.basename(gro))[0]
            molecule_name = base[:-4] if base.endswith("_GMX") else base
        def _bn(v):
            return os.path.basename(v) if isinstance(v, str) and v else ""
        manifest = {
            "molecule": molecule_name or "unknown",
            "gro": _bn(files.get("gro")),
            "top": _bn(files.get("top")),
            "itp": _bn(files.get("itp")),
            "posres": _bn(files.get("posres")),
            "mdp_dir": mdp_subdir,
        }
        (dest_dir / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
        logging.info("[step2] Wrote manifest: %s", dest_dir / "run.json")

    prepared = 0
    for mol_dir in mol_dirs:
        run_gau_dir = mol_dir / "run_gau_create_gmx_in"
        if not run_gau_dir.is_dir():
            logging.info("[step2][skip] %s: no run_gau_create_gmx_in dir", mol_dir.name)
            continue
        dest_dir = mol_dir / "run_gmx_simulation"
        dest_dir.mkdir(exist_ok=True)
        logging.info("[step2] Preparing %s -> %s", mol_dir.name, dest_dir)
        reset_logging()
        setup_logging(str(dest_dir / "process.log"), also_console=False)
        log_run_header("step2")
        upstream_snap = run_gau_dir / "electrofit.toml"
        molecule_input = project / "data" / "input" / mol_dir.name / "electrofit.toml"
        process_cfg = mol_dir / "electrofit.toml"
        project_defaults = project / "electrofit.toml"
        compose_snapshot(dest_dir, project, mol_dir.name, multi_molecule=multi_mol, log_fn=logging.info,
                         upstream=upstream_snap, process_cfg=process_cfg, molecule_input=molecule_input,
                         project_defaults=project_defaults, extra_override=config)
        acpype_dir = None
        for sub in run_gau_dir.iterdir():
            if sub.is_dir() and sub.name.endswith(".acpype"):
                acpype_dir = sub
                break
        if not acpype_dir:
            logging.warning("[step2][skip] %s: no .acpype dir", mol_dir.name)
            continue
        for pattern in _FILE_PATTERNS:
            for fn in os.listdir(acpype_dir):
                if fnmatch.fnmatch(fn, pattern):
                    shutil.copy(acpype_dir / fn, dest_dir)
                    logging.info("[step2] Copied %s -> %s", fn, dest_dir)
        md_dest = dest_dir / "MDP"
        if mdp_source_dir.is_dir():
            if md_dest.exists():
                shutil.rmtree(md_dest)
            shutil.copytree(mdp_source_dir, md_dest)
            logging.info("[step2] Copied MDP -> %s", md_dest)
        else:
            logging.warning("[step2][warn] no MDP source %s", mdp_source_dir)
        selected = {"gro": None, "itp": None, "top": None, "posres": None}  # type: ignore[var-annotated]
        for name in os.listdir(dest_dir):
            if name.endswith("GMX.gro"):
                selected["gro"] = str(dest_dir / name)
            elif name.endswith("GMX.itp") and not name.startswith("posre_"):
                selected["itp"] = str(dest_dir / name)
            elif name.endswith(".top"):
                selected["top"] = str(dest_dir / name)
            elif name.startswith("posre_") and name.endswith(".itp"):
                selected["posres"] = str(dest_dir / name)
        if not selected["gro"]:
            for name in os.listdir(dest_dir):
                if name.endswith(".gro"):
                    selected["gro"] = str(dest_dir / name)
                    break
        if not selected["itp"]:
            for name in os.listdir(dest_dir):
                if name.endswith(".itp") and not name.startswith("posre_"):
                    selected["itp"] = str(dest_dir / name)
                    break
        if not selected["top"]:
            for name in os.listdir(dest_dir):
                if name.endswith(".top"):
                    selected["top"] = str(dest_dir / name)
                    break
        if not selected["posres"]:
            for name in os.listdir(dest_dir):
                if name.startswith("posre_") and name.endswith(".itp"):
                    selected["posres"] = str(dest_dir / name)
                    break
        _write_manifest(dest_dir, selected)  # type: ignore[arg-type]
        prepared += 1
        reset_logging()
    summary = f"[step2] Prepared simulation dirs for {prepared}/{len(mol_dirs)} molecules."
    print(summary)
    try:
        reset_logging()
        setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
        logging.info(summary)
    except Exception:
        pass


# ---- Step 3 (inline from step3) --------------------------------------------------------------
def run_step3(project: Path, config: Path | None = None, molecule: str | None = None, log_console: bool = False) -> None:
    import json as _json
    from electrofit.config.loader import load_config, dump_config
    from electrofit.adapters.gromacs import set_up_production
    from electrofit.io.ff import validate_forcefield
    from electrofit.infra.config_snapshot import compose_snapshot
    from electrofit.infra.logging import setup_logging, reset_logging, log_run_header, enable_header_dedup
    from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

    project = project.resolve()
    _set_project_env(project)
    process_root = project / "process"
    multi_mol = process_root.is_dir() and len([p for p in process_root.iterdir() if p.is_dir()]) > 1
    try:
        enable_header_dedup(True)
        setup_logging(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
        log_run_header("step3")
    except Exception:
        pass
    # Discover run dirs
    def _iter_run_dirs(project_root: Path):
        process = project_root / "process"
        if not process.is_dir():
            return []
        return [p / "run_gmx_simulation" for p in process.iterdir() if (p / "run_gmx_simulation").is_dir()]
    run_dirs_all = list(_iter_run_dirs(project))
    run_dirs = filter_paths_for_molecule(run_dirs_all, molecule)
    if molecule and not run_dirs:
        print(f"[step3][warn] molecule '{molecule}' not found; nothing to do")
        return
    ran = 0
    for run_dir in run_dirs:
        reset_logging()
        setup_logging(str(run_dir / "process.log"), also_console=False)
        try:
            log_run_header("step3")
        except Exception:
            logging.info("electrofit unknown | step=step3")
        mol = run_dir.parent.name
        upstream = run_dir.parent / "run_gau_create_gmx_in" / "electrofit.toml"
        molecule_input = project / "data" / "input" / mol / "electrofit.toml"
        process_cfg = run_dir.parent / "electrofit.toml"
        project_defaults = project / "electrofit.toml"
        compose_snapshot(
            run_dir,
            project,
            mol,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=upstream,
            process_cfg=process_cfg,
            molecule_input=molecule_input,
            project_defaults=project_defaults,
            extra_override=config,
        )
        cfg_run = load_config(project, context_dir=run_dir, molecule_name=mol)
        try:
            dump_config(cfg_run, log_fn=logging.debug)
        except Exception:
            logging.debug("[step3] dump per-run config failed", exc_info=True)

        sim = getattr(cfg_run, "simulation", None)
        ff_attr = getattr(sim, "forcefield", None) if sim else None
        ff = ff_attr or "amber14sb.ff"
        runtime_local = getattr(getattr(cfg_run, "gmx", None), "runtime", None)
        threads = getattr(runtime_local, "threads", None) if runtime_local else None
        pin = getattr(runtime_local, "pin", None) if runtime_local else None
        base_scratch_run = (
            getattr(getattr(cfg_run, "paths", None), "base_scratch_dir", None)
            or "/tmp/electrofit_scratch"
        )
        manifest = run_dir / "run.json"
        if not manifest.is_file():
            logging.info(f"[step3][skip] no run.json in {run_dir}")
            print(f"[step3][skip] no run.json in {run_dir}")
            continue
        meta = _json.loads(manifest.read_text())
        molecule_name = meta.get("molecule")
        gro_file = meta.get("gro")
        mdp_dir = meta.get("mdp_dir") or (getattr(cfg_run.paths, "mdp_dir", None) or "MDP")
        if not molecule_name or not gro_file:
            msg = f"[step3][skip] Skip: incomplete manifest ({manifest})"
            logging.info(msg)
            print(msg)
            continue
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
        try:
            validate_forcefield(ff, explicit=ff_attr is not None)
        except FileNotFoundError as e:
            logging.error("[step3][abort] %s", e)
            print(f"[step3][abort] {e}")
            continue
        logging.info(f"[step3] Starting GROMACS production for {molecule_name} in {run_dir}")
        prev = os.getcwd()
        try:
            os.chdir(run_dir)
            set_up_production(
                m_gro=gro_file,
                MDP_dir=mdp_dir,
                base_scratch_dir=base_scratch_run,
                molecule_name=molecule_name,
                simulation=cfg_run.simulation,
                derived=getattr(cfg_run.simulation, 'derived', None),
                ff=ff,
                threads=threads,
                pin=pin,
            )
            ran += 1
        finally:
            os.chdir(prev)
    summary = (
        f"[step3] Completed {ran} run(s)." if ran else "[step3] No runs executed (no manifests found)."
    )
    print(summary)
    from electrofit.infra.logging import setup_logging as _setup
    _setup(str(project / "step.log"), also_console=log_console, suppress_initial_message=True)
    logging.info(summary)


# ---- Step 4 (reuse worker helper; serial by default) -----------------------------------------
def run_step4(
    project: Path,
    sample: int | None = None,
    method: str | None = None,
    seed: int | None = None,
    workers: int = 1,
    clean: bool = False,
    molecule: str | None = None,
    verbose: bool = False,
    config: Path | None = None,
) -> None:
    import shutil
    from tqdm import tqdm
    from electrofit.infra.config_snapshot import CONFIG_ARG_HELP  # noqa: F401  # for compatibility docstrings
    from electrofit.infra.logging import setup_logging, reset_logging, log_run_header
    from electrofit.pipeline.molecule_filter import filter_paths_for_molecule
    from electrofit.pipeline.workers.step4_extract import _extract_for_molecule

    project = project.resolve()
    _set_project_env(project)
    process_dir = project / "process"
    if not process_dir.is_dir():
        print("[step4] No process directory found; nothing to do.")
        return
    mol_dirs_all = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    mol_dirs = filter_paths_for_molecule(mol_dirs_all, molecule)
    if molecule and not mol_dirs:
        print(f"[step4][warn] molecule '{molecule}' not found; nothing to do")
        return
    multi_mol = len(mol_dirs) > 1

    # Minimal sampling defaults from project config are handled inside the worker via compose_snapshot
    method = method or "linear"
    sample = sample or 20

    if clean:
        for mdir in mol_dirs:
            ec_dir = mdir / "extracted_conforms"
            if ec_dir.is_dir():
                try:
                    shutil.rmtree(ec_dir)
                    print(f"[step4][clean] removed {ec_dir}")
                except Exception as e:
                    print(f"[step4][clean][warn] failed to remove {ec_dir}: {e}")

    extracted = 0
    iterator: Iterable[Path] = mol_dirs
    iterator = tqdm(iterator, desc="step4", unit="mol") if workers == 1 else iterator
    for sub in iterator:
        ok, msg = _extract_for_molecule(sub, project, sample, method, seed, config, multi_mol, verbose)
        status = "[step4]" if ok else "[step4][skip]"
        print(f"{status} {sub.name}: {msg}")
        if ok:
            extracted += 1
    summary = f"[step4] Conformer extraction complete for {extracted}/{len(mol_dirs)} molecules." if extracted else "[step4] No conformers extracted."
    print(summary)
    # Record into per-step logs
    for sub in mol_dirs:
        ec_dir = sub / "extracted_conforms"
        if ec_dir.is_dir():
            reset_logging()
            setup_logging(str(ec_dir / "process.log"), also_console=False)
            try:
                log_run_header("step4")
            except Exception:
                pass


# ---- Step 5/6/7/8 (use provided helpers) -----------------------------------------------------
def run_step5(
    project: Path,
    batch_size: int = 6,
    interval: float = 0,
    max_batches: int | None = None,
    dry_run: bool = False,
    keep_going: bool = False,
    mock: bool = False,
    verbose: bool = False,
    no_parallel: bool = False,
    isolate_conformer: bool = False,
    config: Path | None = None,
    molecule: str | None = None,
) -> int:
    from electrofit.pipeline.steps.step5 import run_step5 as _run
    project = project.resolve()
    _set_project_env(project)
    return _run(
        project=project,
        batch_size=batch_size,
        interval=interval,
        max_batches=max_batches,
        dry_run=dry_run,
        keep_going=keep_going,
        mock=mock,
        verbose=verbose,
        no_parallel=no_parallel,
        isolate_conformer=isolate_conformer,
        config_override=config,
        only_molecule=molecule,
)  # type: ignore[return-value]


def run_step6(
    project: Path,
    config: Path | None = None,
    remove_outlier: bool = False,
    plot_histograms: bool = False,
    hist_combine_groups: bool = False,
    hist_bins: int = 20,
    outlier_iqr_factor: float = 1.5,
    molecule: str | None = None,
) -> int:
    from electrofit.pipeline.steps.step6 import run_step6 as _run
    project = project.resolve()
    _set_project_env(project)
    return _run(
        project=project,
        override_cfg=config,
        remove_outlier=remove_outlier,
        plot_histograms=plot_histograms,
        hist_combine_groups=hist_combine_groups,
        hist_bins=hist_bins,
        outlier_iqr_factor=outlier_iqr_factor,
        only_molecule=molecule,
    )


def run_step7(project: Path, config: Path | None = None, molecule: str | None = None) -> int:
    from electrofit.pipeline.steps.step7 import run_step7 as _run
    project = project.resolve()
    _set_project_env(project)
    return _run(project, config, molecule)


def run_step8(project: Path, config: Path | None = None, log_console: bool = False, molecule: str | None = None) -> int:
    from electrofit.pipeline.steps.step8 import run_step8 as _run
    project = project.resolve()
    _set_project_env(project)
    return _run(project, config, log_console, molecule)


# ---- High-level driver -----------------------------------------------------------------------
def run_pipeline(
    project: Path,
    steps: Iterable[int] = (0, 1, 2, 3),
    molecule: str | None = None,
    config: Path | None = None,
) -> None:
    project = project.resolve()
    _set_project_env(project)
    for s in steps:
        if s == 0:
            run_step0(project, molecule=molecule)
        elif s == 1:
            run_step1(project, config=config, molecule=molecule)
        elif s == 2:
            run_step2(project, config=config, molecule=molecule)
        elif s == 3:
            run_step3(project, config=config, molecule=molecule)
        elif s == 4:
            run_step4(project, molecule=molecule, config=config)
        elif s == 5:
            run_step5(project, config=config, molecule=molecule)
        elif s == 6:
            run_step6(project, config=config, molecule=molecule)
        elif s == 7:
            run_step7(project, config=config, molecule=molecule)
        elif s == 8:
            run_step8(project, config=config, molecule=molecule)
        else:
            raise ValueError(f"Unknown step: {s}")


def _parse_steps(expr: str) -> List[int]:
    """Parse step expressions like "0-3,5,7-8" into a sorted unique list."""
    out: set[int] = set()
    for chunk in expr.split(','):
        chunk = chunk.strip()
        if not chunk:
            continue
        if '-' in chunk:
            a, b = chunk.split('-', 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(chunk))
    return sorted(out)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Run Electrofit steps via direct imports (no CLI entrypoints).")
    ap.add_argument("--project", "-p", required=True, help="Path to project root")
    ap.add_argument("--steps", default="0-3", help="Steps to run, e.g. '0-3,5,7-8' (default: 0-3)")
    ap.add_argument("--molecule", help="Restrict to this molecule under process/", default=None)
    ap.add_argument("--config", help="Path to override TOML for snapshot layering", default=None)
    args = ap.parse_args(argv)
    project = Path(args.project)
    steps = _parse_steps(args.steps)
    cfg = Path(args.config).resolve() if args.config else None
    run_pipeline(project, steps=steps, molecule=args.molecule, config=cfg)


if __name__ == "__main__":
    main()

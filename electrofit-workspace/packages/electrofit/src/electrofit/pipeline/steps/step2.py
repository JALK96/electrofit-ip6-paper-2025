"""Pipeline Step 2: Prepare per-molecule GROMACS simulation input directories.

Creates ``process/<mol>/run_gmx_simulation`` directories, layers config snapshot,
copies topology + MDP templates, and writes a deterministic ``run.json`` manifest
for Step 3. Pure orchestration: filesystem scanning, snapshot composition,
logging setup per molecule.
"""
from __future__ import annotations
import argparse
import fnmatch
import os
import shutil
import json
import logging
from pathlib import Path
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, log_run_header, reset_logging, enable_header_dedup
from electrofit.infra.step_logging import log_relevant_config
from electrofit.config.loader import load_config, dump_config
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

__all__ = ["main"]

_FILE_PATTERNS = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]

def _write_manifest(dest_dir: Path, files: dict[str,str], mdp_subdir: str = "MDP") -> None:
    gro = files.get("gro")
    molecule = None
    if gro:
        base = os.path.splitext(os.path.basename(gro))[0]
        molecule = base[:-4] if base.endswith("_GMX") else base
    def _bn(v):
        return os.path.basename(v) if isinstance(v, str) and v else ""
    manifest = {
        "molecule": molecule or "unknown",
        "gro": _bn(files.get("gro")),
        "top": _bn(files.get("top")),
        "itp": _bn(files.get("itp")),
        "posres": _bn(files.get("posres")),
        "mdp_dir": mdp_subdir,
    }
    (dest_dir / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    logging.info("[step2] Wrote manifest: %s", dest_dir / "run.json")

def main():  # pragma: no cover
    ap = argparse.ArgumentParser(description="Step2: Prepare run_gmx_simulation directories")
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    ap.add_argument("--molecule", help="Limit to single molecule (directory under process/)")
    args = ap.parse_args()
    project_path = Path(args.project).resolve()
    process_dir = project_path / "process"
    if not process_dir.is_dir():
        print("[step2] No process directory.")
        return
    mol_dirs_all = [p for p in process_dir.iterdir() if p.is_dir()]
    mol_dirs = filter_paths_for_molecule(mol_dirs_all, args.molecule)
    if args.molecule and not mol_dirs:
        print(f"[step2][warn] molecule '{args.molecule}' not found; nothing to do")
        return
    multi_mol = len(mol_dirs) > 1
    enable_header_dedup(True)
    setup_logging(str(project_path / "step.log"), also_console=args.log_console, suppress_initial_message=True)
    log_run_header("step2")  # single logical header; later calls suppressed
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
        molecule_input = project_path / "data" / "input" / mol_dir.name / "electrofit.toml"
        process_cfg = mol_dir / "electrofit.toml"
        project_defaults = project_path / "electrofit.toml"
        compose_snapshot(dest_dir, project_path, mol_dir.name, multi_molecule=multi_mol, log_fn=logging.info,
                         step="step2",
                         upstream=upstream_snap, process_cfg=process_cfg, molecule_input=molecule_input,
                         project_defaults=project_defaults, extra_override=Path(args.config) if args.config else None)
        cfg = load_config(project_path, context_dir=dest_dir, molecule_name=mol_dir.name)
        try:
            dump_config(cfg, log_fn=logging.debug)
        except Exception:  # pragma: no cover
            logging.debug("[step2] config dump failed", exc_info=True)
        try:
            log_relevant_config(
                "step2",
                cfg,
                [
                    "project.molecule_name",
                    "project.protocol",
                    "paths.mdp_dir",
                    "paths.base_scratch_dir",
                ],
            )
        except Exception:  # pragma: no cover
            logging.debug("[step2][cfg] selective config logging failed", exc_info=True)
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
        mdp_dir_setting = getattr(getattr(cfg, "paths", None), "mdp_dir", None) or "data/MDP"
        mdp_source_dir = Path(mdp_dir_setting)
        if not mdp_source_dir.is_absolute():
            mdp_source_dir = project_path / mdp_source_dir
        if mdp_source_dir.is_dir():
            if md_dest.exists():
                shutil.rmtree(md_dest)
            shutil.copytree(mdp_source_dir, md_dest)
            logging.info("[step2] Copied MDP -> %s", md_dest)
        else:
            logging.warning("[step2][warn] no MDP source %s", mdp_source_dir)
        selected = {"gro": None, "itp": None, "top": None, "posres": None}
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
        _write_manifest(dest_dir, selected)
        prepared += 1
        reset_logging()
    summary = f"[step2] Prepared simulation dirs for {prepared}/{len(mol_dirs)} molecules."
    print(summary)
    try:
        reset_logging()
        setup_logging(str(project_path / "step.log"), also_console=args.log_console, suppress_initial_message=True)
        logging.info(summary)
    except Exception:
        pass

if __name__ == "__main__":  # pragma: no cover
    main()

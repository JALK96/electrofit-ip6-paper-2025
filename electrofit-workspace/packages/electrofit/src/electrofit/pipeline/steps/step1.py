"""Pipeline Step 1: Initial molecule processing orchestration.

Discovers run directories (``run_gau_create_gmx_in`` / future alias) under
``process/`` and executes the initial structure preparation logic directly via
the domain layer (``electrofit.domain.prep.process_initial``) while layering
config snapshots. Legacy workflow fallback removed.

Responsibility: orchestrate discovery, per-run logging setup, snapshot
composition, and delegate the actual processing.
"""
from __future__ import annotations

import argparse
import os
import logging
from pathlib import Path
from contextlib import contextmanager

from electrofit.domain.prep.process_initial import InitialPrepConfig, process_initial
from electrofit.config.loader import load_config, dump_config, resolve_symmetry_flags
from electrofit.io.files import find_file_with_extension
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, log_run_header, reset_logging, enable_header_dedup
from electrofit.pipeline.molecule_filter import molecule_component

__all__ = ["main"]

_TARGET_RUN_DIRS = {
    "run_gau_create_gmx_in",
    "init_system_create_gmx_in",
}

@contextmanager
def _pushd(path: str | None):
    prev = os.getcwd()
    if path:
        os.chdir(path)
    try:
        yield
    finally:
        if path:
            os.chdir(prev)

def _iter_run_dirs(project_root: str, only_molecule: str | None = None):
    process_root = os.path.join(project_root, "process")
    if not os.path.isdir(process_root):
        return
    for root, dirs, _files in os.walk(process_root):
        for d in dirs:
            if d in _TARGET_RUN_DIRS:
                full = os.path.join(root, d)
                if only_molecule is None or molecule_component(Path(full)) == only_molecule:
                    yield full

def _run_one_dir(run_dir: str, project_root: str, override_cfg: str | None, multi_mol: bool):
    os.environ["ELECTROFIT_PROJECT_PATH"] = project_root
    mol = Path(run_dir).parent.name
    reset_logging()
    setup_logging(str(Path(run_dir) / "process.log"), also_console=False)
    log_run_header("step1")
    molecule_input = Path(project_root) / "data" / "input" / mol / "electrofit.toml"
    project_defaults = Path(project_root) / "electrofit.toml"
    process_defaults = Path(project_root) / "process" / mol / "electrofit.toml"
    compose_snapshot(
        Path(run_dir),
        Path(project_root),
        mol,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        process_cfg=process_defaults,
        molecule_input=molecule_input,
        project_defaults=project_defaults,
        extra_override=Path(override_cfg) if override_cfg else None,
    )
    with _pushd(run_dir):
        # Direct domain integration path
        run_dir_abs = os.getcwd()
        run_cfg_path = os.path.join(run_dir_abs, "electrofit.toml")
        cfg = load_config(project_root, config_path=run_cfg_path if os.path.isfile(run_cfg_path) else None)
        dump_config(cfg, log_fn=logging.debug)
        base_scratch_dir = (
            getattr(cfg.paths, "base_scratch_dir", None)
            or os.environ.get("ELECTROFIT_SCRATCH_DIR")
            or "/tmp/electrofit_scratch"
        )
        mol2_from_name = None
        if cfg.project.molecule_name:
            cand = os.path.join(run_dir_abs, f"{cfg.project.molecule_name}.mol2")
            if os.path.isfile(cand):
                mol2_from_name = cand
        mol2_file = mol2_from_name or find_file_with_extension("mol2")
        if not mol2_file:
            logging.info("[step1][skip] no mol2 file detected in run directory; nothing to process")
            return
        molecule_name = os.path.splitext(os.path.basename(mol2_file))[0]
        adjust_initial, ignore_initial = resolve_symmetry_flags(cfg, "initial")
        cfg.project.adjust_symmetry = adjust_initial
        cfg.project.ignore_symmetry = ignore_initial
        additional_input: list[str] = []
        if adjust_initial:
            try:
                json_file = find_file_with_extension("json")
                if json_file:
                    additional_input.append(json_file)
            except Exception:  # pragma: no cover
                logging.debug("[step1] json detection failed", exc_info=True)
        residue_name = cfg.project.residue_name or molecule_name[:3].upper()
        net_charge = cfg.project.charge if cfg.project.charge is not None else 0
        atom_type = cfg.project.atom_type or "gaff2"
        protocol = cfg.project.protocol or "bcc"
    # Setup scratch (scratch manager)
        from electrofit.infra.scratch_manager import setup_scratch_directory
        input_files = [os.path.basename(mol2_file)] + list(additional_input)
        scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
    # Ensure mol2 is present inside scratch root (scratch manager copies)
        cfg_obj = InitialPrepConfig(
            molecule_name=molecule_name,
            mol2_file=os.path.basename(mol2_file),
            net_charge=net_charge,
            residue_name=residue_name,
            atom_type=atom_type,
            adjust_sym=adjust_initial,
            ignore_sym=ignore_initial,
            protocol=protocol,
        )
        # Use same ensure_finalized context to keep identical side effects
        from electrofit.cli.safe_run import ensure_finalized
        try:
            with ensure_finalized(original_dir=original_dir, scratch_dir=scratch_dir, input_files=input_files):
                result = process_initial(cfg_obj, scratch_dir, original_dir, input_files)
        except Exception:
            logging.error("[step1] error during initial processing", exc_info=True)
            raise
        # Emit manifest json (relative paths only; avoid embedding transient scratch paths)
        try:
            import json
            manifest_path = os.path.join(run_dir_abs, "initial_manifest.json")
            def _rel(p: str | None):
                return os.path.basename(p) if p else None
            # After finalize, artifacts were copied back (potentially with *_copyN). Re-scan to pick actual names.
            gro_name = None
            itp_name = None
            top_name = None
            posre_names: list[str] = []
            for fname in os.listdir(run_dir_abs):
                if fname.endswith('_GMX.gro'):
                    gro_name = fname
                elif fname.endswith('_GMX.itp'):
                    itp_name = fname
                elif fname.endswith('_GMX.top'):
                    top_name = fname
                elif fname.startswith('posre_') and fname.endswith('.itp'):
                    posre_names.append(fname)
            data = {
                'protocol': result.protocol,
                'residue_name': result.residue_name,
                'charges_source': result.charges_source,
                'mol2_final': cfg_obj.mol2_file,  # modified input now in place under this name
                'acpype_dir': _rel(result.acpype_dir),
                'gmx_gro': gro_name,
                'gmx_itp': itp_name,
                'gmx_top': top_name,
                'posre_itps': sorted(posre_names),
            }
            with open(manifest_path, 'w') as fh:
                json.dump(data, fh, indent=2)
            logging.info("[step1] wrote initial_manifest.json")
        except Exception:  # pragma: no cover
            logging.debug("[step1] manifest emission failed", exc_info=True)


def main():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description=(
            "Step1: Initial molecule processing (prepare Gaussian/Antechamber inputs). "
            "Discovers run_gau_create_gmx_in directories under process/ and merges electrofit.toml snapshots."
        )
    )
    parser.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    parser.add_argument("--config", help=CONFIG_ARG_HELP)
    parser.add_argument("--log-console", action="store_true", help="Also echo logging to console")
    parser.add_argument("--molecule", help="Process only this molecule name (directory under process/)")
    args = parser.parse_args()

    project_root = os.path.abspath(args.project)
    process_root = Path(project_root) / "process"
    multi_mol = False
    if process_root.is_dir():
        mol_dirs = [p for p in process_root.iterdir() if p.is_dir()]
        multi_mol = len(mol_dirs) > 1

    # Global header (Projektweite step.log) – auch wenn keine Runs gefunden werden.
    try:
        enable_header_dedup(True)
        setup_logging(str(Path(project_root) / "step.log"), also_console=args.log_console, suppress_initial_message=True)
        log_run_header("step1")  # single header per step (dedup suppresses later repeats)
    except Exception:
        pass

    discovered = list(_iter_run_dirs(project_root, only_molecule=args.molecule))
    if args.molecule and not discovered:
        print(f"[step1][warn] requested molecule '{args.molecule}' not found / no run dirs")
        return
    if discovered:
        ran = 0
        for run_dir in discovered:
            _run_one_dir(run_dir, project_root, args.config, multi_mol)
            ran += 1
        # Zusammenfassung im globalen Log festhalten
        try:
            reset_logging()
            setup_logging(str(Path(project_root) / "step.log"), also_console=args.log_console, suppress_initial_message=True)
            logging.info("[step1] Completed %d run directory(ies).", ran)
        except Exception:
            pass
        return
    run_dir = os.getcwd()
    _run_one_dir(run_dir, project_root, args.config, multi_mol)
    try:
        reset_logging()
        setup_logging(str(Path(project_root) / "step.log"), also_console=args.log_console, suppress_initial_message=True)
        logging.info("[step1] Completed 1 run directory (fallback current working directory).")
    except Exception:
        pass

if __name__ == "__main__":  # pragma: no cover
    main()

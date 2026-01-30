"""Domain helpers for batch processing of conformer directories (Step5).

Encapsulates discovery + per‑conformer processing so the step orchestrator
remains minimal.

Responsibilities:
    * Discover conformer directories (heuristic: directories containing a PDB).
        * Execute per‑conformer processing (config snapshot, logging, mock mode,
            Gaussian/RESP pipeline via domain ``process_conformer``).

Notes:
        * Heavy Gaussian/RESP work now uses domain ``process_conformer``; scratch
            lifecycle is managed inline here until a dedicated service is extracted.
    * Advanced diagnostic monkey‑patching is retained but initialised only once.
"""
from __future__ import annotations

from pathlib import Path
import os
import logging
import traceback
from electrofit.config.loader import load_config, dump_config, resolve_symmetry_flags
from electrofit.infra.config_snapshot import compose_snapshot
from electrofit.infra.logging import log_run_header, reset_logging, setup_logging
from electrofit.domain.charges.process_conformer import (
    ConformerConfig as _ConformerConfig,
    process_conformer as _process_conformer_impl,
)
from electrofit.io.files import find_file_with_extension, strip_extension

__all__ = [
    "discover_conformer_dirs",
    "process_conformer_dir",
]


def discover_conformer_dirs(project_root: Path) -> list[Path]:
    """Return conformer directories containing at least one PDB file.

    Mirrors historical logic from ``step5_process_conforms._discover_conformer_dirs``.
    Empty list if project/process layout not present.
    """
    process_dir = project_root / "process"
    if not process_dir.is_dir():
        return []
    dirs: list[Path] = []
    for mol_dir in sorted(process_dir.iterdir()):
        ec = mol_dir / "extracted_conforms"
        if not ec.is_dir():
            continue
        for conf in sorted(ec.iterdir()):
            if not conf.is_dir():
                continue
            if any(conf.glob("*.pdb")):
                dirs.append(conf)
    return dirs


_ADV_DIAG_INITIALISED = False


def _init_advanced_diagnostics():  # pragma: no cover - diagnostic path
    global _ADV_DIAG_INITIALISED
    if _ADV_DIAG_INITIALISED:
        return
    if os.environ.get("ELECTROFIT_DISABLE_ADV_DIAG") == "1":  # opt-out switch
        _ADV_DIAG_INITIALISED = True
        return
    try:
        import sys as _sys_d
        import faulthandler as _fh_d
        import signal as _sig_d
        import atexit as _ax_d
        import gc as _gc_d
        try:
            _fh_d.enable(file=_sys_d.stderr)
        except Exception:
            pass

        def _ax():
            try:
                print('[debug-atexit] worker normal shutdown', flush=True)
            except Exception:
                pass

        try:
            _ax_d.register(_ax)
        except Exception:
            pass

        def _sig_handler(sig, frame):  # noqa: ARG001
            try:
                print(f'[debug-signal] sig={sig}', flush=True)
                _fh_d.dump_traceback(file=_sys_d.stderr)
            except Exception:
                pass

        for _s in [getattr(_sig_d, n, None) for n in ('SIGTERM','SIGINT','SIGSEGV')]:
            if _s is not None:
                try:
                    _sig_d.signal(_s, _sig_handler)
                except Exception:
                    pass
        if hasattr(_gc_d, 'callbacks'):
            try:
                _gc_d.callbacks.append(lambda phase, info: phase=='stop' and print('[debug-gc] cycle', flush=True))  # type: ignore
            except Exception:
                pass
    except Exception:
        pass
    _ADV_DIAG_INITIALISED = True


def process_conformer_dir(
    conf_dir: Path,
    project_root: Path,
    override_cfg_path: Path | None,
    multi_mol: bool,
    mock: bool,
    verbose: bool,
) -> tuple[str, bool, str]:
    """Process a single conformer directory.

    Returns (relative_path, ok_flag, message).
    This is a near‑lift of the previous worker implementation with minor
    cleanups and a direct call to the domain conformer pipeline.
    """
    prev = os.getcwd()
    try:
        os.chdir(conf_dir)
        _init_advanced_diagnostics()
        # Monkeypatch sys/os exit for diagnostics (kept verbatim)
        import sys as _sys
        import os as _os  # local
        _orig_sys_exit = _sys.exit
        _orig_os_exit = _os._exit

        def _dbg_sys_exit(code=0):  # noqa: D401
            try:
                print(f"[debug-sys-exit] code={code}", flush=True)
            except Exception:
                pass
            return _orig_sys_exit(code)

        def _dbg_os_exit(code=0):  # noqa: D401
            try:
                print(f"[debug-os-exit] code={code}", flush=True)
            except Exception:
                pass
            return _orig_os_exit(code)

        _sys.exit = _dbg_sys_exit  # type: ignore
        _os._exit = _dbg_os_exit   # type: ignore

        # per‑conformer logging
        reset_logging()
        log_path = Path(os.getcwd()) / "process.log"
        setup_logging(str(log_path), also_console=verbose)
        log_run_header("step5")

        pdb_file = find_file_with_extension("pdb")
        if not pdb_file:
            raise FileNotFoundError("No PDB file found in conformer directory")
        mol_name = conf_dir.parent.parent.name if len(conf_dir.parts) >= 2 else conf_dir.parent.name
        snap_path = compose_snapshot(
            conf_dir,
            project_root,
            mol_name,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=conf_dir.parent / "electrofit.toml",  # parent is extracted_conforms root
            process_cfg=conf_dir.parent.parent / "electrofit.toml" if len(conf_dir.parts) >= 3 else None,
            molecule_input=project_root / "data" / "input" / mol_name / "electrofit.toml",
            project_defaults=project_root / "electrofit.toml",
            extra_override=override_cfg_path,
        )
        if not snap_path:
            logging.warning(f"[step5][warn] snapshot could not be created in {conf_dir}")
        cfg = load_config(project_root, context_dir=conf_dir)
        dump_config(cfg, header=True, log_fn=logging.debug)
        proj = cfg.project
        molecule_name = proj.molecule_name or strip_extension(pdb_file)
        if mock:
            with open("executed.txt", "w") as f:
                f.write(f"run{molecule_name}")
        else:
            # ── Protocol-aware, resilient scratch setup ────────────────────────────────────
            from electrofit.infra.scratch_manager import setup_scratch_directory, finalize_scratch_directory  # noqa: F401

            base_scratch = getattr(cfg.paths, "base_scratch_dir", None) or "/tmp/electrofit_scratch"

            pdb_path = Path(pdb_file)
            input_files: list[str] = []

            # Always include the PDB (and warn if missing)
            if pdb_path.exists():
                input_files.append(pdb_path.name)
            else:
                logging.warning("PDB not found at %s; scratch setup will proceed without it.", pdb_path)

            protocol   = (proj.protocol or "bcc").lower()
            adjust_sym, ignore_sym = resolve_symmetry_flags(cfg, "ensemble")
            proj.adjust_symmetry = adjust_sym  # type: ignore[attr-defined]
            proj.ignore_symmetry = ignore_sym  # type: ignore[attr-defined]

            if protocol == "opt":
                # RESP inputs (only include if present next to the PDB)
                resp1 = "ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN"
                resp2 = "ANTECHAMBER_RESP2.IN"
                for fname in (resp1, resp2):
                    fpath = pdb_path.parent / fname
                    if fpath.exists():
                        input_files.append(fname)
                    else:
                        logging.info("RESP file not found (skipping): %s", fpath)

            elif protocol == "bcc" and adjust_sym:
                # Prefer symmetry JSONs from conf_dir if available, else from PDB dir
                json_candidates = []
                try:
                    if 'conf_dir' in locals() and conf_dir:
                        json_candidates += list(Path(conf_dir).glob("*.json"))
                except Exception:
                    pass
                json_candidates += list(pdb_path.parent.glob("*.json"))

                if json_candidates:
                    # Prefer a JSON whose name contains the molecule name; then shortest name
                    json_candidates.sort(key=lambda p: (molecule_name not in p.name, len(p.name)))
                    input_files.append(json_candidates[0].name)
                else:
                    logging.info("No symmetry JSON found; proceeding without one (adjust_sym=True).")

            # Ensure PDB is first in the list for clarity
            if pdb_path.name in input_files:
                input_files = [pdb_path.name] + [x for x in input_files if x != pdb_path.name]

            # Create scratch directory with the selected inputs
            scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch)

            # Build the conformer config (keeps your robust defaults)
            cfg_obj = _ConformerConfig(
                molecule_name=molecule_name,
                pdb_file=pdb_file,
                net_charge=proj.charge or 0,
                residue_name=proj.residue_name or "LIG",
                adjust_sym=adjust_sym,
                ignore_sym=ignore_sym,
                symmetry_ensemble_mode=getattr(getattr(cfg, "symmetry", None), "ensemble", None),
                protocol=protocol,
            )
            try:
                _process_conformer_impl(cfg_obj, scratch_dir, original_dir, input_files, defer_finalize=False)
            finally:
                # Always finalize to replicate previous behaviour (non-deferred)
                try:
                    from electrofit.infra.scratch_manager import finalize_scratch_directory
                    finalize_scratch_directory(
                        original_dir,
                        scratch_dir,
                        input_files,
                        output_files=None,
                        overwrite=True,
                        remove_parent_if_empty=False,
                        reason="conformer-batch",
                    )
                except Exception:  # pragma: no cover
                    logging.debug("[step5] finalize scratch failed", exc_info=True)
        rel = conf_dir.relative_to(project_root)
        print(f"[worker-return] {rel} ok", flush=True)
        return (str(rel), True, "ok")
    except Exception as e:  # pragma: no cover - defensive
        try:
            tb = traceback.format_exc()
            logging.error(f"Worker failure {conf_dir}: {e}\n{tb}")
        except Exception:
            pass
        try:
            rel = conf_dir.relative_to(project_root)
            rel_str = str(rel)
        except Exception:
            rel_str = str(conf_dir)
        print(f"[worker-exc] {rel_str}: {e}", flush=True)
        return (rel_str, False, str(e))
    finally:
        try:  # restore monkeypatches
            import sys as _sys
            import os as _os
            if '_orig_sys_exit' in locals():
                _sys.exit = _orig_sys_exit  # type: ignore
            if '_orig_os_exit' in locals():
                _os._exit = _orig_os_exit   # type: ignore
        except Exception:
            pass
        try:
            os.chdir(prev)
        except Exception:
            pass

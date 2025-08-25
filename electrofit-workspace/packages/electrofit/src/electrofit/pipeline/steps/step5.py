"""Pipeline Step 5: Gaussian/RESP processing of extracted conformers.

Thin orchestration layer for conformer batch processing.

Responsibilities:
  * Discover conformer directories.
  * Batch + (optionally) parallelise processing.
  * Support isolation mode via legacy isolate runner (subprocess per conformer).
  * Maintain backwards compatible console output for existing tests & users.

Heavy logic of per‑conformer processing lives in
``domain.charges.conformer_batch.process_conformer_dir``.
"""
from __future__ import annotations

import argparse
import os
import time
import multiprocessing
from pathlib import Path
from typing import Iterable, List
from concurrent.futures import ProcessPoolExecutor, as_completed, Future
import logging

from electrofit.domain.charges.conformer_batch import (
    discover_conformer_dirs,
    process_conformer_dir,
)
from electrofit.infra.config_snapshot import CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, reset_logging, log_run_header, enable_header_dedup
from electrofit.pipeline.molecule_filter import molecule_component

__all__ = ["main", "run_step5"]


def _batched(seq: List[Path], batch_size: int) -> Iterable[list[Path]]:
    if batch_size <= 0:
        batch_size = len(seq) or 1
    for i in range(0, len(seq), batch_size):
        yield seq[i : i + batch_size]


def run_step5(
    project: Path,
    batch_size: int,
    interval: float,
    max_batches: int | None,
    dry_run: bool,
    keep_going: bool,
    mock: bool = False,
    verbose: bool = False,
    no_parallel: bool = False,
    isolate_conformer: bool = False,
    config_override: Path | None = None,
    only_molecule: str | None = None,
) -> int:
    logging.getLogger().setLevel(logging.DEBUG)
    logging.debug("[step5] Root log level set to DEBUG")
    try:
        conf_dirs = discover_conformer_dirs(project)
        if only_molecule:
            conf_dirs = [c for c in conf_dirs if molecule_component(c) == only_molecule]
            if not conf_dirs:
                print(f"[step5][warn] molecule '{only_molecule}' produced no conformer dirs.")
                return 0
        if not conf_dirs:
            print("[step5] No conformer directories found.")
            return 0
        print(f"[step5] Found {len(conf_dirs)} conformer dir(s). Batch size={batch_size}.")
        processed = 0
        batches_run = 0
        override_cfg = config_override if config_override else None
        multi_mol = len({c.parts[-3] for c in conf_dirs if len(c.parts) >= 3}) > 1
        for batch in _batched(conf_dirs, batch_size):
            batches_run += 1
            print(f"[step5] Executing batch {batches_run} containing {len(batch)} conformer(s) (parallel).")
            if dry_run:
                for conf_dir in batch:
                    rel = conf_dir.relative_to(project)
                    print(f"[step5] (dry-run) Would process {rel}")
            else:
                if no_parallel or len(batch) == 1 or isolate_conformer:
                    mode_reason = "--no-parallel" if no_parallel else "single-item batch"
                    if isolate_conformer:
                        mode_reason = "isolation" if len(batch)==1 else mode_reason+"+isolation"
                    print(f"[step5] Using direct serial execution for this batch ({mode_reason}).")
                    for conf_dir in batch:
                        rel = conf_dir.relative_to(project)
                        print(f"[step5][debug-call] starting process_one({rel})", flush=True)
                        if isolate_conformer:
                            import subprocess
                            import json  # noqa: F401
                            cmd = [
                                os.environ.get("PYTHON", "python"),
                                "-m","electrofit.pipeline.steps._step5_isolate_runner",  # pipeline isolate runner
                                "--conf", str(conf_dir),
                                "--project", str(project),
                            ]
                            if override_cfg:
                                cmd += ["--override", str(override_cfg)]
                            if multi_mol:
                                cmd.append("--multi-mol")
                            if mock:
                                cmd.append("--mock")
                            if verbose:
                                cmd.append("--verbose")
                            try:
                                cp = subprocess.run(cmd, capture_output=True, text=True)
                                stdout = cp.stdout.splitlines()
                                parsed = None
                                for line in stdout[::-1]:
                                    if line.startswith("RESULT:"):
                                        try:
                                            import json as _json
                                            parsed = _json.loads(line[len("RESULT:"):])
                                        except Exception:
                                            pass
                                        break
                                if parsed:
                                    result = (parsed.get("rel"), bool(parsed.get("ok")), parsed.get("msg"))
                                else:
                                    result = (str(rel), False, f"no RESULT line (rc={cp.returncode})")
                            except Exception as e:  # launch failure
                                result = (str(rel), False, f"isolation failure: {e}")
                        else:
                            result = process_conformer_dir(
                                conf_dir,
                                project,
                                override_cfg,
                                multi_mol,
                                mock,
                                verbose,
                            )
                        print(f"[step5][debug-result] {result!r}", flush=True)
                        if not (isinstance(result, tuple) and len(result) == 3):
                            print(f"[step5][error] {rel}: unexpected result shape {result!r}")
                            if not keep_going:
                                return 1
                            continue
                        rel_s, ok, msg = result
                        if ok:
                            processed += 1
                            if verbose:
                                print(f"[step5][done] {rel_s}")
                        else:
                            print(f"[step5][error] {rel_s}: {msg}")
                            if not keep_going:
                                return 1
                else:
                    futures: list[Future] = []
                    future_map: dict[Future, str] = {}
                    try:
                        mp_ctx = multiprocessing.get_context("spawn")
                    except ValueError:  # pragma: no cover
                        mp_ctx = None
                    if mp_ctx:
                        print("[step5] Using spawn context for process pool.")
                    else:
                        print("[step5][warn] Spawn context unavailable; falling back to default multiprocessing context.")
                    with ProcessPoolExecutor(max_workers=len(batch), mp_context=mp_ctx) as ex:
                        for conf_dir in batch:
                            rel = conf_dir.relative_to(project)
                            print(f"[step5] Queued {rel}")
                            fut = ex.submit(
                                process_conformer_dir,
                                conf_dir,
                                project,
                                override_cfg,
                                multi_mol,
                                mock,
                                verbose,
                            )
                            futures.append(fut)
                            future_map[fut] = str(rel)
                        for fut in as_completed(futures):
                            rel_lookup = future_map.get(fut, "<unknown>")
                            try:
                                result = fut.result()
                            except Exception as e:  # pragma: no cover
                                import traceback as _tb
                                print(f"[step5][future-exc] {rel_lookup}: {type(e).__name__}: {e}")
                                tb_txt = _tb.format_exc()
                                for line in tb_txt.rstrip().splitlines():
                                    print(f"[step5][trace] {line}")
                                if not keep_going:
                                    print(f"[step5] Stopping after future exception (processed {processed} / {len(conf_dirs)}).")
                                    for other in futures:
                                        if not other.done():
                                            other.cancel()
                                    return 1
                                else:
                                    continue
                            if not (isinstance(result, tuple) and len(result) == 3):
                                print(f"[step5][error] {rel_lookup}: unexpected result shape {result!r}")
                                if not keep_going:
                                    print(f"[step5] Stopping after malformed result (processed {processed} / {len(conf_dirs)}).")
                                    for other in futures:
                                        if not other.done():
                                            other.cancel()
                                    return 1
                                continue
                            relr, ok, msg = result
                            if ok:
                                processed += 1
                                if verbose:
                                    print(f"[step5][done] {relr}")
                            else:
                                print(f"[step5][error] {relr}: {msg}")
                                if not keep_going:
                                    print(
                                        f"[step5] Stopping after failure (processed {processed} / {len(conf_dirs)})."
                                    )
                                    for other in futures:
                                        if not other.done():
                                            other.cancel()
                                    return 1
            remaining = len(conf_dirs) - processed if not dry_run else len(conf_dirs) - (batches_run * batch_size)
            if remaining <= 0:
                break
            if max_batches and batches_run >= max_batches:
                print(f"[step5] Reached max-batches={max_batches}; {remaining} conformer(s) left unprocessed.")
                break
            if interval > 0:
                print(f"[step5] Sleeping {interval} seconds before next batch (remaining ~{remaining}).")
                time.sleep(interval)
        print(f"[step5] Completed processing of {processed if not dry_run else processed or (batches_run*batch_size)} conformer(s).")
        return 0
    except Exception as fatal:  # pragma: no cover
        import traceback as _tb
        print(f"[step5][fatal] Unhandled exception: {type(fatal).__name__}: {fatal}")
        for line in _tb.format_exc().rstrip().splitlines():
            print(f"[step5][trace] {line}")
        return 1


def main(argv: list[str] | None = None):  # pragma: no cover - CLI exercised in tests
    class _FlagOnlyParser(argparse.ArgumentParser):
        def error(self, message):  # override to inject isolate-conformer hint
            if '--isolate-conformer' in message and 'unrecognized arguments' in message:
                message += "\nHint: --isolate-conformer is a flag-only option. Example: electrofit step5 --project <path> --isolate-conformer"
            super().error(message)

    ap = _FlagOnlyParser(
        description=(
            "Step5: Gaussian/RESP processing of extracted conformers.\n"
            "Batching & parallelisation: --batch-size controls conformers per batch; --no-parallel forces serial execution.\n"
            "Diagnostics: --isolate-conformer runs each conformer in its own subprocess (crash isolation).\n"
            "Limit progress: --max-batches <N>. Dry run: --dry-run.\n"
            "Environment: ELECTROFIT_DEBUG_GAUSSIAN_CACHE (Gaussian cache), ELECTROFIT_DEBUG_DEFER_FINALIZE (defer scratch finalize), ELECTROFIT_LOG_LEVEL."
        )
    )
    ap.add_argument(
        "--project",
        default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()),
        help="Path to project root (contains 'process/')",
    )
    ap.add_argument(
        "--batch-size",
        type=int,
        default=6,
        help="Conformer parallelism per batch (default 6; 0=all).",
    )
    ap.add_argument(
        "--interval",
        type=float,
        default=0,
        help="Seconds to sleep between batches (default 0; 0 disables).",
    )
    ap.add_argument(
        "--max-batches",
        type=int,
        help="Optional limit of batches to run (useful for partial processing).",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="List batches/scripts without executing them.")
    ap.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue executing remaining scripts even if one fails.")
    ap.add_argument(
        "--mock",
        action="store_true",
        help="Fast mock mode (skip heavy Gaussian/RESP; create marker executed.txt).",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Also emit per-conformer log lines to console (default: only write process.log).",
    )
    ap.add_argument(
        "--no-parallel",
        action="store_true",
        help="Disable process pool; run batches serially (pool still auto-disabled for single-item batches).",
    )
    ap.add_argument(
        "--isolate-conformer",
        action="store_true",
        help="Run each conformer in an isolated subprocess (flag-only; do NOT supply a value).",
    )
    ap.add_argument(
        "--config",
        help=CONFIG_ARG_HELP,
    )
    ap.add_argument(
        "--molecule",
        help="Restrict to conformers belonging to this molecule (name under process/)",
    )
    args = ap.parse_args(argv)
    project = Path(args.project).resolve()
    config_override = Path(args.config).resolve() if getattr(args, "config", None) else None
    # Global Header vor Ausführung
    try:
        enable_header_dedup(True)
        setup_logging(str(project / "step.log"), also_console=True, suppress_initial_message=True)
        log_run_header("step5")
    except Exception:
        pass

    rc = run_step5(
        project=project,
        batch_size=args.batch_size,
        interval=args.interval,
        max_batches=args.max_batches,
        dry_run=args.dry_run,
        keep_going=args.keep_going,
        mock=args.mock,
        verbose=args.verbose,
        no_parallel=args.no_parallel,
        isolate_conformer=args.isolate_conformer,
        config_override=config_override,
        only_molecule=getattr(args, "molecule", None),
    )
    print(f"[step5][debug-end] rc={rc}", flush=True)
    # Summary Logging
    try:
        reset_logging()
        setup_logging(str(project / "step.log"), also_console=True, suppress_initial_message=True)
        logging.info(f"[step5] Exit code {rc}; 0 = success")
    except Exception:
        pass
    raise SystemExit(rc)


if __name__ == "__main__":  # pragma: no cover
    main()

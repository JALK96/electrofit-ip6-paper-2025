"""Pipeline Step 6: Aggregate RESP ensemble charges & prepare averaged MOL2.

Thin orchestration wrapper around domain aggregation logic.
"""
from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from electrofit.domain.aggregation.average_charges import process_molecule_average_charges
from electrofit.infra.logging import setup_logging, reset_logging, log_run_header, enable_header_dedup
from electrofit.infra.config_snapshot import CONFIG_ARG_HELP
from electrofit.pipeline.molecule_filter import filter_paths_for_molecule

__all__ = ["main", "run_step6"]


def run_step6(
    project: Path,
    override_cfg: Path | None,
    remove_outlier: bool,
    plot_histograms: bool,
    hist_combine_groups: bool,
    hist_bins: int,
    outlier_iqr_factor: float,
    only_molecule: str | None = None,
) -> int:
    process_dir = project / "process"
    if not process_dir.is_dir():
        print("[step6] No process directory.")
        return 0
    mol_dirs_all = [p for p in sorted(process_dir.iterdir()) if p.is_dir()]
    mol_dirs = filter_paths_for_molecule(mol_dirs_all, only_molecule)
    if only_molecule and not mol_dirs:
        print(f"[step6][warn] molecule '{only_molecule}' not found; nothing to do")
        return 0
    multi_mol = len(mol_dirs) > 1
    done = 0
    for mol_dir in mol_dirs:
        results_dir = mol_dir / "results"
        results_dir.mkdir(exist_ok=True)
        log_path = results_dir / "process.log"
        try:
            reset_logging()
            setup_logging(str(log_path), also_console=False)
        except Exception:  # pragma: no cover
            pass
        ok, msg = process_molecule_average_charges(
            mol_dir=mol_dir,
            project_root=project,
            override_cfg=override_cfg,
            multi_mol=multi_mol,
            remove_outlier=remove_outlier,
            plot_histograms=plot_histograms,
            hist_combine_groups=hist_combine_groups,
            hist_bins=hist_bins,
            outlier_iqr_factor=outlier_iqr_factor,
        )
        status = "[step6]" if ok else "[step6][skip]"
        print(f"{status} {mol_dir.name}: {msg}")
        if ok:
            done += 1
    summary = f"[step6] Completed {done}/{len(mol_dirs)} molecules."
    print(summary)
    try:
        reset_logging()
        setup_logging(str(project / "step.log"), also_console=True, suppress_initial_message=True)
        logging.info(summary)
    except Exception:
        pass
    return 0


def main(argv: list[str] | None = None):  # pragma: no cover
    ap = argparse.ArgumentParser(description="Step6: Aggregate RESP ensemble charges and update MOL2 (user charges mode).")
    ap.add_argument("--project", required=False, default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--remove-outlier", action="store_true", help="Enable outlier removal (IQR based) before group averaging")
    ap.add_argument("--plot-histograms", action="store_true", help="Generate per-atom histogram PDFs (hist*.pdf) for diagnostic inspection")
    ap.add_argument("--hist-combine-groups", action="store_true", help="When plotting adjusted histograms, overlay group-averaged distributions")
    ap.add_argument("--hist-bins", type=int, default=20, help="Number of bins for histograms (default: 20)")
    ap.add_argument("--outlier-iqr-factor", type=float, default=1.5, help="IQR factor threshold for defining outliers (default: 1.5)")
    ap.add_argument("--molecule", help="Only aggregate for this molecule name")
    args = ap.parse_args(argv)
    project_root = Path(args.project).resolve()
    override_cfg = Path(args.config).resolve() if getattr(args, "config", None) else None
    # Global Header vor Ausführung (nur einmal)
    try:
        enable_header_dedup(True)
        setup_logging(str(project_root / "step.log"), also_console=True, suppress_initial_message=True)
        log_run_header("step6")
    except Exception:
        pass
    rc = run_step6(
        project=project_root,
        override_cfg=override_cfg,
        remove_outlier=args.remove_outlier,
        plot_histograms=args.plot_histograms,
        hist_combine_groups=args.hist_combine_groups,
        hist_bins=args.hist_bins,
        outlier_iqr_factor=args.outlier_iqr_factor,
        only_molecule=getattr(args, 'molecule', None),
    )
    raise SystemExit(rc)


if __name__ == "__main__":  # pragma: no cover
    main()

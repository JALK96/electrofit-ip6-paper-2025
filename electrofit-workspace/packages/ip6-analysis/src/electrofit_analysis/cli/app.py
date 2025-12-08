"""
Unified CLI entry point for electrofit-analysis.

Provides subcommands that forward to existing CLI modules in this package.

Subcommands
-----------
- distance:     Compute Na–P distances for IP6 microstates.
- count:        Compute Na–P distances, Na+ counts, buried volume, and excess.
- summarize-nap:    Summarize distances/counts/coordination across patterns (requires out dir).
- plot-2d:      Render 2D molecule depictions for microstates.
- coordination: Analyze Na–IP6 coordination (counts, RDFs, optional projections).

Usage examples
--------------
python -m electrofit_analysis.cli.app distance -p /path/to/project
python -m electrofit_analysis.cli.app dist-count -p /path/to/project
python -m electrofit_analysis.cli.app summarize-nap -p /path/to/project -o /path/to/out
python -m electrofit_analysis.cli.app plot-2d -p /path/to/project [--subdir process]
python -m electrofit_analysis.cli.app coordination -p /path/to/project [--subdir process] [--determine-global-y] [--rdf-y-max 1800] [--plot-projection]
"""

from __future__ import annotations

import argparse
import os
from .common import parse_molecule_args


def _cmd_distance(args: argparse.Namespace) -> None:
    from .coordination.na_p_distance_ip6 import main as distance_main

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    distance_main(
        project,
        stage=getattr(args, 'stage', 'final'),
        only=only,
        rep=getattr(args, 'rep', None),
    )


def _cmd_count(args: argparse.Namespace) -> None:
    from .coordination.na_p_dist_count_ip6 import main as count_main

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    count_main(
        project,
        stage=getattr(args, 'stage', 'final'),
        only=only,
        rep=getattr(args, 'rep', None),
    )


def _cmd_summarize_nap(args: argparse.Namespace) -> None:
    # Import the updated summarizer (keep the module path you use in your package)
    from .coordination.summerize_nap_dist_count_ip6 import main as summarize_main

    project = os.path.abspath(args.project)
    root_dir = os.path.join(project, "process")
    out_dir = os.path.abspath(args.output_dir)

    # helpers to parse comma-separated lists
    def _parse_list(opt_str, allowed, default):
        if opt_str is None:
            return set(default)
        items = {s.strip().lower() for s in opt_str.split(",") if s.strip()}
        bad = items - set(allowed)
        if bad:
            raise SystemExit(
                f"Unknown option(s) {sorted(bad)}; allowed: {sorted(allowed)}"
            )
        return items

    # distance/count plot types
    plots = _parse_list(
        getattr(args, "plots", None),
        allowed={"means", "box", "overlaid", "stats"},
        default={"means", "box", "overlaid", "stats"},
    )
    # coordination plot types
    coord_plots = _parse_list(
        getattr(args, "coord_plots", None),
        allowed={"box", "stats"},
        default={"box"},
    )

    # Forward all required params to the updated summarizer.
    # (If none of --distance/--count/--coordination are set, the summarizer's
    #  own default is to run all three, so we just pass the booleans as-is.)
    summarize_main(
        root_dir=root_dir,
        output_dir=out_dir,
        do_distance=getattr(args, "distance", False),
        do_count=getattr(args, "count", False),
        do_coord=getattr(args, "coordination", False),
        plots=plots,
        coord_plots=coord_plots,
        stage=getattr(args, 'stage', 'final'),
    )


def _cmd_plot_2d(args: argparse.Namespace) -> None:
    from .plot_molecule_2d import main as plot2d_main

    project = os.path.abspath(args.project)
    plot2d_main(project, args.subdir)


def _cmd_coordination(args: argparse.Namespace) -> None:
    from .coordination.Na_IP6_coordination import main as coord_main

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    coord_main(
        project_dir=project,
        subdir=args.subdir,
        stage=getattr(args, 'stage', 'final'),
        only=only,
        determine_global_y=args.determine_global_y,
        rdf_y_max=args.rdf_y_max,
        plot_projection=args.plot_projection,
        rdf_data_path=args.rdf_data,
        rep=getattr(args, 'rep', None),
        boxplot=getattr(args, 'boxplot', True),
    )


def _cmd_charges(args: argparse.Namespace) -> None:
    from .charges.plot_charges_ip6 import main as charges_main

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    charges_main(project, remove_outlier=args.remove_outlier, only=only)


def _cmd_dihedrals(args: argparse.Namespace) -> None:
    from .dihedrals.dihedrals_ip6 import main as dihedrals_main

    project = os.path.abspath(args.project) if args.project else None
    only = parse_molecule_args(getattr(args, 'molecule', None))
    dihedrals_main(
        project,
        stage=getattr(args, 'stage', 'final'),
        only=only,
        rep=getattr(args, 'rep', None),
    )


def _cmd_hbonds(args: argparse.Namespace) -> None:
    from .h_bonds.h_bonds_ip6 import main as hb_main

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    hb_main(
        project,
        viz=args.viz,
        stage=getattr(args, 'stage', 'final'),
        only=only,
        rep=getattr(args, 'rep', None),
    )


def _cmd_pp_matrix(args: argparse.Namespace) -> None:
    from .h_bonds.make_pp_matrix_ip6 import main as ppm_main

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    root = os.path.join(project, "process")
    ppm_main(root, args.kind, args.mode, args.width_mode, args.width_ref, stage=getattr(args, 'stage', 'final'), only=only)


def _cmd_hbonds_compare(args: argparse.Namespace) -> None:
    from .h_bonds.h_bonds_comparison_ip6_microstates import (
        main as hbonds_compare_main,
    )

    project = os.path.abspath(args.project)
    only = parse_molecule_args(getattr(args, 'molecule', None))
    hbonds_compare_main(args.type, project, stage=getattr(args, 'stage', 'final'), only=only)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ip6-analysis",
        description="Unified CLI for IP6 analysis tools.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # distance
    p_dist = sub.add_parser(
        "distance",
        help="Compute Na–P distances for the IP6 project (per microstate).",
    )
    p_dist.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_dist.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Analyze final (run_final_gmx_simulation), sample (run_gmx_simulation), or remd (run_remd_gmx_simulation) trajectories.",
    )
    p_dist.add_argument(
        "--rep",
        type=int,
        default=None,
        help="Replica index for stage=remd (e.g. 0..N-1).",
    )
    p_dist.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated). Accepts 'IP_…' or bare names.")
    p_dist.set_defaults(func=_cmd_distance)

    # count
    p_count = sub.add_parser(
        "dist-count",
        help=(
            "Run Na–P distances, Na+ counts, buried volume, and excess analysis for the IP6 project."
        ),
    )
    p_count.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_count.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Analyze final (run_final_gmx_simulation), sample (run_gmx_simulation), or remd (run_remd_gmx_simulation) trajectories.",
    )
    p_count.add_argument(
        "--rep",
        type=int,
        default=None,
        help="Replica index for stage=remd (e.g. 0..N-1).",
    )
    p_count.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_count.set_defaults(func=_cmd_count)


    # summarize-nap
    p_sum = sub.add_parser(
        "summarize-nap",
        help=(
            "Summarize Na–P distances/counts/coordination across patterns. Requires an output directory.\n"
            "This command aggregates results from all microstates and produces summary reports.\n"
            "Run 'ip6-analysis dist-count -p <project>' and/or 'ip6-analysis coordination -p <project>'\n"
            "before using this, to generate the necessary input files."
        ),
    )
    p_sum.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_sum.add_argument(
        "-o", "--out", "--output", dest="output_dir", required=True, help="Output directory."
    )

    # NEW: metric selectors (non-mutually-exclusive)
    p_sum.add_argument(
        "--distance", action="store_true",
        help="Generate distance summaries (requires: distances_NaP*.xvg). If none of --distance/--count/--coordination is given, all are run."
    )
    p_sum.add_argument(
        "--count", action="store_true",
        help="Generate ion-count summaries (requires: ion_count_P*.xvg)."
    )
    p_sum.add_argument(
        "--coordination", action="store_true",
        help="Generate coordination summaries (requires: NaP_coordination_bool.npy)."
    )

    # NEW: plot-type toggles
    p_sum.add_argument(
        "--plots", default="means,box,overlaid,stats",
        help="Comma-separated plot types for distance/count. Choices: means,box,overlaid,stats. Default: all."
    )
    p_sum.add_argument(
        "--coord-plots", default="box",
        help="Comma-separated plot types for coordination. Choices: box,stats. Default: box."
    )
    p_sum.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Read from analyze_final_sim, analyze_sample_sim, or analyze_remd_sim.",
    )

    p_sum.set_defaults(func=_cmd_summarize_nap)

    # plot-2d
    p_plot = sub.add_parser(
        "plot-2d", help="Render 2D molecule depictions (SVG) for microstates."
    )
    p_plot.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_plot.add_argument(
        "--subdir",
        default="process",
        help="Relative subdirectory under the project root to traverse (default: process).",
    )
    p_plot.set_defaults(func=_cmd_plot_2d)

    # coordination
    p_coord = sub.add_parser(
        "coordination",
        help=(
            "Analyze Na–IP6 coordination across microstates (counts, RDFs, optional projections)."
        ),
    )
    p_coord.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_coord.add_argument(
        "--subdir",
        default="process",
        help="Relative subdirectory under the project root to traverse (default: process).",
    )
    p_coord.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Analyze final, sample, or REMD trajectories.",
    )
    p_coord.add_argument(
        "--rep",
        type=int,
        default=None,
        help="Replica index for stage=remd (e.g. 0..N-1).",
    )
    p_coord.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_coord.add_argument(
        "--determine-global-y",
        action="store_true",
        help="Scan all microstates to determine a common RDF y-limit.",
    )
    p_coord.add_argument(
        "--rdf-y-max",
        type=float,
        default=None,
        help="Override RDF y-limit with a fixed value (skips scanning).",
    )
    p_coord.add_argument(
        "--plot-projection",
        action="store_true",
        help="Also generate 2D/3D network projection snapshots.",
    )
    p_coord.add_argument(
        "--rdf-data",
        default=None,
        help="Path to an RDF cache file (JSON). Reuse or refresh cached RDF curves.",
    )
    p_coord.add_argument(
        "--no-boxplot",
        dest="boxplot",
        action="store_false",
        help="Skip the per-microstate coordination boxplot.",
    )
    p_coord.set_defaults(func=_cmd_coordination)

    # charges
    p_chg = sub.add_parser(
        "charges",
        help=(
            "Plot charge distributions and histograms for atoms, with optional symmetry handling."
        ),
    )
    p_chg.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_chg.add_argument(
        "--remove-outlier",
        action="store_true",
        help="Apply 1.5*IQR outlier removal before computing averages and plots.",
    )
    p_chg.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_chg.set_defaults(func=_cmd_charges)

    # dihedrals
    p_dih = sub.add_parser(
        "dihedrals",
        help="Analyze phosphate and ring dihedral angles from final GMX runs and generate plots.",
    )
    p_dih.add_argument(
        "-p",
        "--project",
        required=False,
        help="Project root (defaults to $ELECTROFIT_PROJECT_PATH or CWD if omitted).",
    )
    p_dih.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Analyze final, sample, or REMD trajectories.",
    )
    p_dih.add_argument(
        "--rep",
        type=int,
        default=None,
        help="Replica index for stage=remd (e.g. 0..N-1).",
    )
    p_dih.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_dih.set_defaults(func=_cmd_dihedrals)

    # hbonds (run per-species analysis and plots)
    p_hb = sub.add_parser(
        "hbonds",
        help="Run hydrogen-bond analysis for each species and produce plots/reports.",
    )
    p_hb.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_hb.add_argument(
        "--viz",
        "-viz",
        action="store_true",
        help="Generate plots/figures (opt-in).",
    )
    p_hb.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Analyze final, sample, or REMD trajectories.",
    )
    p_hb.add_argument(
        "--rep",
        type=int,
        default=None,
        help="Replica index for stage=remd (e.g. 0..N-1).",
    )
    p_hb.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_hb.set_defaults(func=_cmd_hbonds)

    # pp-matrix (P-to-P matrix from H-bond XPM/log)
    p_ppm = sub.add_parser(
        "pp-matrix",
        help=(
            "Compute P-to-P interaction matrices from H-bond data and draw directed diagrams.\n"
            "Run 'ip6-analysis hbonds -p <project>' before using."
        ),
    )
    p_ppm.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_ppm.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Read HBonds outputs from analyze_final_sim, analyze_sample_sim, or analyze_remd_sim.",
    )
    p_ppm.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_ppm.add_argument(
        "-k",
        "--kind",
        choices=["intra", "inter"],
        default="intra",
        help="Choose which H-bond matrices to use (default: intra).",
    )
    p_ppm.add_argument(
        "--mode",
        choices=["union", "sum", "both"],
        default="union",
        help="How to aggregate O–O bonds into P→P fractions (default: union).",
    )
    p_ppm.add_argument(
        "--width-mode",
        choices=["continuous", "coarse"],
        default="continuous",
        help="Edge width scaling mode for diagrams (default: continuous).",
    )
    p_ppm.add_argument(
        "--width-ref",
        choices=["auto", "absolute"],
        default="auto",
        help="Scale relative to graph max (auto) or absolute fraction (absolute).",
    )
    p_ppm.set_defaults(func=_cmd_pp_matrix)

    # hbonds-compare (combined microstate summaries)
    p_hbc = sub.add_parser(
        "hbonds-compare",
        help=(
            "Aggregate H-bond data across microstates: 3x2 figure, summary plot, and textual tables."
        ),
    )
    p_hbc.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_hbc.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Read HBonds outputs from analyze_final_sim, analyze_sample_sim, or analyze_remd_sim.",
    )
    p_hbc.add_argument(
        "-t",
        "--type",
        default="intra",
        help="H-bond type prefix used in filenames (default: intra).",
    )
    p_hbc.add_argument("-m", "--molecule", action="append", help="Molecule(s) to include (repeat or comma-separated).")
    p_hbc.set_defaults(func=_cmd_hbonds_compare)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

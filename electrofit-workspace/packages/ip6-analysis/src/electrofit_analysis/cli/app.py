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
python -m electrofit_analysis.cli.app count -p /path/to/project
python -m electrofit_analysis.cli.app summarize-nap -p /path/to/project -o /path/to/out
python -m electrofit_analysis.cli.app plot-2d -p /path/to/project [--subdir process]
python -m electrofit_analysis.cli.app coordination -p /path/to/project [--subdir process] [--determine-global-y] [--rdf-y-max 1800] [--plot-projection]
"""

from __future__ import annotations

import argparse
import os


def _cmd_distance(args: argparse.Namespace) -> None:
    from .coordination.na_p_distance_ip6 import main as distance_main

    project = os.path.abspath(args.project)
    distance_main(project)


def _cmd_count(args: argparse.Namespace) -> None:
    from .coordination.na_p_dist_count_ip6 import main as count_main

    project = os.path.abspath(args.project)
    count_main(project)


def _cmd_summarize_nap(args: argparse.Namespace) -> None:
    from .coordination.summerize_nap_dist_count_ip6 import main as summarize_main

    project = os.path.abspath(args.project)
    root_dir = os.path.join(project, "process")
    out_dir = os.path.abspath(args.output_dir)
    summarize_main(root_dir, out_dir)


def _cmd_plot_2d(args: argparse.Namespace) -> None:
    from .plot_molecule_2d import main as plot2d_main

    project = os.path.abspath(args.project)
    plot2d_main(project, args.subdir)


def _cmd_coordination(args: argparse.Namespace) -> None:
    from .coordination.Na_IP6_coordination import main as coord_main

    project = os.path.abspath(args.project)
    coord_main(
        project_dir=project,
        subdir=args.subdir,
        determine_global_y=args.determine_global_y,
        rdf_y_max=args.rdf_y_max,
        plot_projection=args.plot_projection,
    )


def _cmd_charges(args: argparse.Namespace) -> None:
    from .charges.plot_charges_ip6 import main as charges_main

    project = os.path.abspath(args.project)
    charges_main(project, remove_outlier=args.remove_outlier)


def _cmd_dihedrals(args: argparse.Namespace) -> None:
    from .dihedrals.dihedrals_ip6 import main as dihedrals_main

    project = os.path.abspath(args.project) if args.project else None
    dihedrals_main(project)


def _cmd_hbonds(args: argparse.Namespace) -> None:
    from .h_bonds.h_bonds_ip6 import main as hb_main

    project = os.path.abspath(args.project)
    hb_main(project, viz=args.viz)


def _cmd_pp_matrix(args: argparse.Namespace) -> None:
    from .h_bonds.make_pp_matrix_ip6 import main as ppm_main

    project = os.path.abspath(args.project)
    root = os.path.join(project, "process")
    ppm_main(root, args.kind, args.mode, args.width_mode, args.width_ref)


def _cmd_hbonds_compare(args: argparse.Namespace) -> None:
    from .h_bonds.h_bonds_comparison_ip6_microstates import (
        main as hbonds_compare_main,
    )

    project = os.path.abspath(args.project)
    hbonds_compare_main(args.type, project)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ip6-analysis",
        description="Unified CLI for IP6 analysis tools (formerly electrofit-analysis).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # distance
    p_dist = sub.add_parser(
        "distance",
        help="Compute Na–P distances for an IP6 project (per microstate).",
    )
    p_dist.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_dist.set_defaults(func=_cmd_distance)

    # count
    p_count = sub.add_parser(
        "count",
        help=(
            "Run Na–P distances, Na+ counts, buried volume, and excess analysis for an IP6 project."
        ),
    )
    p_count.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_count.set_defaults(func=_cmd_count)

    # summarize-nap
    p_sum = sub.add_parser(
        "summarize-nap",
        help=(
            "Summarize Na–P distances/counts/coordination across patterns. Requires an output directory."
        ),
    )
    p_sum.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_sum.add_argument(
        "-o", "--out", "--output", dest="output_dir", required=True, help="Output directory."
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
    p_chg.set_defaults(func=_cmd_charges)

    # dihedrals
    p_dih = sub.add_parser(
        "dihedrals",
        help="Analyze dihedral angles from final GMX runs and generate plots.",
    )
    p_dih.add_argument(
        "-p",
        "--project",
        required=False,
        help="Project root (defaults to $ELECTROFIT_PROJECT_PATH or CWD if omitted).",
    )
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
    p_hb.set_defaults(func=_cmd_hbonds)

    # pp-matrix (P-to-P matrix from H-bond XPM/log)
    p_ppm = sub.add_parser(
        "pp-matrix",
        help=(
            "Compute P-to-P interaction matrices from H-bond data and draw directed diagrams."
        ),
    )
    p_ppm.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_ppm.add_argument(
        "-k",
        "--kind",
        choices=["intra", "inter"],
        default="intra",
        help="Choose which H-bond matrices to use (default: intra).",
    )
    p_ppm.add_argument(
        "-m",
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
        "-t",
        "--type",
        default="intra",
        help="H-bond type prefix used in filenames (default: intra).",
    )
    p_hbc.set_defaults(func=_cmd_hbonds_compare)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

import logging
import argparse
import os
import sys

import seaborn as sns
from electrofit_analysis.cli.common import resolve_stage
from electrofit.infra.logging import setup_logging
from electrofit_analysis.structure.util.common_gmx import (
    create_index_file,
    plot_all_distances_subplots as _plot_all_distances_subplots,
    run_pairdist_commands,
)

# Set the style for seaborn
sns.set_context("talk")

def plot_all_distances_subplots_png(output_prefix, num_groups, plot_filename="all_distances_subplots.png"):
    """Thin wrapper to call shared utility with a PNG default and title."""
    _plot_all_distances_subplots(
        output_prefix=output_prefix,
        num_groups=num_groups,
        plot_filename=plot_filename,
        suptitle="Minimum Distances Between Phosphorus Groups and Na-Ions",
    )


# ----------------------------
# Main Execution Flow
# ----------------------------


def main(project_dir: str, stage: str = 'final') -> None:
    """Run Na–P distance analysis for all molecules under project_dir.

    Args:
        project_dir: Path to the project root (must contain a 'process' subdirectory).
    """
    # Define the base process directory
    process_dir = os.path.join(project_dir, "process")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)

        # Check if it's a directory
        if os.path.isdir(folder_path):
            run_dir_name, analyze_base = resolve_stage(stage)
            run_sim_dir = os.path.join(folder_path, run_dir_name)

            # Check if run dir exists
            if os.path.isdir(run_sim_dir):
                # Define the destination directory analyze_final_sim/NaP_dist_count
                dest_dir = os.path.join(folder_path, analyze_base, "NaP_dist_count")
                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)

                # Define paths (adjust these as necessary)
                structure_file = os.path.join(run_sim_dir, "md.gro")
                trajectory_file = os.path.join(run_sim_dir, "md_center.xtc")
                # place index file inside the destination directory for consistency
                index_file = os.path.join(dest_dir, "NA_P_index.ndx")
                topology_file = os.path.join(run_sim_dir, "md.tpr")
                selection_group = "NA"
                output_prefix = (
                    "distances_NaP"  # Generates distances_P1.xvg to distances_P6.xvg
                )
                log_file = "distances_NaP_gmx.log"

                # List of phosphorus groups
                p_groups = ["P", "P1", "P2", "P3", "P4", "P5"]

                # Set up logging
                setup_logging(log_file)
                logging.info("Logging is set up.")

                # Check if input files exist
                input_files = [structure_file, trajectory_file, topology_file]
                for file in input_files:
                    if not os.path.isfile(file):
                        logging.error(
                            f"Required input file '{file}' not found. Exiting."
                        )
                        sys.exit(1)
                logging.info("All required input files are present.")

                # Step 1: Create the index file
                create_index_file(structure_file, index_file)

                # Step 2: Run gmx pairdist for each phosphorus group
                run_pairdist_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    selection_group=selection_group,
                    output_prefix=output_prefix,
                )

                # Optional Step 3: Plot the distance data
                plot_all_distances_subplots_png(
                    output_prefix=output_prefix, num_groups=len(p_groups)
                )

                logging.info("All tasks completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compute Na–P distances for the IP6 project. Provide the project root via --project."
        )
    )
    parser.add_argument(
        "--project",
        "-p",
        required=True,
        help="Path to the project root directory (contains 'process').",
    )
    args = parser.parse_args()
    main(os.path.abspath(args.project))

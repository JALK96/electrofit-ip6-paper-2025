import logging
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from electrofit.infra.logging import setup_logging
from electrofit_analysis.structure.util.common_gmx import (
    create_index_file,
    plot_all_distances_subplots,
    run_pairdist_commands,
    plot_timeseries_grid,
    run_ion_count_commands,
    run_buried_volume_commands,
    run_ion_count_whole_molecule
)

# Set the style for seaborn
sns.set_context("talk")


# keep only the shared utility; local duplicate removed to satisfy linter

def plot_ion_counts_subplots(
    output_prefix, groups, plot_filename="ion_counts_subplots.pdf"
):
    """Plot number of ions within cutoff for each group using generic grid plotter."""
    files = [f"{output_prefix}{g}.xvg" for g in groups]
    labels = groups
    plot_timeseries_grid(
        files=files,
        labels=labels,
        ylabel="Na+ ions",
        outfile=plot_filename,
        suptitle=None,  # was commented out zuvor
        unit_scale=1000.0,  # ps → ns
        nrows=2,
        ncols=3,
        line_color="darkblue",
        mean_line=True,
        mean_color="red",
        mean_fmt="{mean:.1f}",
        xlabel="Time (ns)",
    )

def plot_whole_molecule_ion_count(xvg_file, plot_filename="ion_count_whole_mol.pdf"):
    """
    Plots a single line: # of NA ions within cutoff for the entire molecule vs time.
    """
    logging.info("Plotting ion count for the entire molecule...")
    try:
        data = np.loadtxt(xvg_file, comments=("#", "@"))
        time = data[:, 0] / 1000.0
        ion_count = data[:, 1]
        mean_count = np.mean(ion_count)

        plt.figure(figsize=(6, 4))
        plt.plot(time, ion_count, color="darkblue", linestyle="-", linewidth=0.5)
        plt.axhline(y=mean_count, color="red", linestyle="--", linewidth=1.5)
        plt.text(
            0.95,
            0.1,
            f"Mean: {mean_count:.1f}",
            ha="right",
            va="top",
            transform=plt.gca().transAxes,
            color="red",
            fontsize=12,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )
        # plt.title('Na+ Ion Count Near Entire Molecule')
        plt.xlabel("Time (ns)")
        plt.ylabel("# of Na+ ions")
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=300)
        # plt.show()
        logging.info("Ion count for entire molecule plotted successfully.")
    except Exception as e:
        logging.error(f"Error plotting whole molecule ion count: {e}")


# --------------------------------------------------------------------
#  Third-order (frame-specific) excess definition
#
#  N_excess(t) = N_in_sphere(t) / [ c · ( 4/3·π·r³  −  V_IP6(t) ) ]
#
#  c           = N_total_Na / V_box  (bulk concentration)
#  V_IP6(t)    = n_IP6_atoms_in_sphere(t) · v_atom_IP6  (per-atom excluded volume)
#  r           = cutoff_radius used in the gmx select commands.
#
#  The code logs every piece of the formula so you can reconstruct how
#  each number was obtained when reading the log file.


def _parse_box_volume_from_gro(gro_file):
    """
    Returns simulation box volume (nm³) from the last line of a .gro file.
    For orthorhombic boxes the last line contains three floats lx ly lz.
    For triclinic boxes the first three floats are the box vector lengths
    (off‑diagonal terms are ignored → gives upper bound to actual volume).
    """
    with open(gro_file, "r") as fh:
        lines = fh.readlines()
        # last non‑empty line should contain box vectors
        box_line = lines[-1].strip().split()
        if len(box_line) < 3:
            raise ValueError("Could not read box vectors from .gro file.")
        lx, ly, lz = map(float, box_line[:3])
        return lx * ly * lz  # nm³


def _count_na_atoms_in_gro(gro_file):
    """
    Counts how many atoms have the string 'NA' in columns 11‑15 (atom name)
    of a .gro file. Works for standard GROMACS naming.
    """
    count = 0
    with open(gro_file, "r") as fh:
        lines = fh.readlines()
    # skip first 2 title/atom‑number lines, last line is box
    for line in lines[2:-1]:
        atom_name = line[10:15].strip()
        if atom_name.startswith("NA"):
            count += 1
    if count == 0:
        raise ValueError("No Na⁺ atoms found in structure file; check naming.")
    return count


def compute_excess_ion_counts(
    ion_count_prefix,
    buried_prefix,
    groups,
    structure_file,
    cutoff_radius,
    output_prefix="excess_",
):
    """
    Excess‑ion factor using *true* solvent‑accessible volume per frame:

        N_excess(t) = N_in_sphere(t) /
                      [ c_Na · (4/3·π·r³ − V_buried(t)) ]

    where V_buried(t) is obtained from gmx sasa (-tv) and already has units nm³.
    """
    import math

    # Bulk Na⁺ concentration
    n_na_total = _count_na_atoms_in_gro(structure_file)
    v_box = _parse_box_volume_from_gro(structure_file)
    c_na = n_na_total / v_box

    v_sphere = (4.0 / 3.0) * math.pi * cutoff_radius**3
    logging.info("============================================================")
    logging.info("EXCESS‑ION FORMULA (grid/SAV version):")
    logging.info("  N_excess(t) = N_in / [ c_Na · (4/3·π·r³ − V_buried(t)) ]")
    logging.info(f"  r = {cutoff_radius:.3f} nm ; 4/3·π·r³ = {v_sphere:.6f} nm³")
    logging.info(
        f"  c_Na = N_Na_tot / V_Box = {n_na_total}/{v_box:.6f} = {c_na:.6f} nm⁻³"
    )
    logging.info("============================================================")

    for group in groups:
        na_file = f"{ion_count_prefix}{group}.xvg"
        buried_file = f"{buried_prefix}{group}.xvg"
        out_file = f"{output_prefix}{group}.xvg"
        try:
            na_data = np.loadtxt(na_file, comments=("#", "@"))
            buried_data = np.loadtxt(buried_file, comments=("#", "@"))
            if na_data.shape[0] != buried_data.shape[0]:
                raise ValueError("Frame mismatch between Na and buried‑volume files")

            time_ps = na_data[:, 0]
            n_in = na_data[:, 1]
            v_buried = buried_data[:, 1]
            v_eff = v_sphere - v_buried
            # ensure positive effective volume
            v_eff[v_eff <= 1e-9] = 1e-9

            n_solution = c_na * v_eff
            n_excess = n_in / n_solution

            header = (
                '@    title "Excess Na+ ions (grid SAV)"\n'
                '@    xaxis  label "Time (ps)"\n'
                '@    yaxis  label "N_excess"\n'
            )
            np.savetxt(
                out_file,
                np.column_stack((time_ps, n_excess)),
                header=header,
                comments="",
            )

            # log means
            logging.info(
                f"[{group}] mean N_in = {np.mean(n_in):.4f}, "
                f"mean V_buried = {np.mean(v_buried):.4f} nm³, "
                f"mean V_eff = {np.mean(v_eff):.4f} nm³, "
                f"mean N_excess = {np.mean(n_excess):.4f}"
            )
        except Exception as e:
            logging.error(f"Failed for {group}: {e}")


def plot_excess_ion_counts_subplots(
    output_prefix, groups, plot_filename="excess_ion_counts_subplots.pdf"
):
    """Plot N_excess vs time for each group using generic grid plotter."""
    files = [f"{output_prefix}{g}.xvg" for g in groups]
    labels = groups
    plot_timeseries_grid(
        files=files,
        labels=labels,
        ylabel="N_excess",
        outfile=plot_filename,
        suptitle=None,
        unit_scale=1000.0,
        nrows=2,
        ncols=3,
        line_color=None,  # vorher Standardfarbe
        mean_line=True,
        mean_color="black",  # vorher keine spezifische Farbe → neutral
        mean_fmt="{mean:.2f}",
        xlabel="Time (ns)",
    )


# ----------------------------
# Main Execution Flow
# ----------------------------

def main(project_dir: str) -> None:
    """Run Na–P distance and ion-count analysis for all molecules under project_dir.

    Args:
        project_dir: Path to the project root (must contain a 'process' subdirectory).
    """
    # Define the base process directory
    process_dir = os.path.join(project_dir, "process")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)

        if os.path.isdir(folder_path):
            # Define the 'run_final_gmx_simulation' directory within this folder
            run_final_sim_dir = os.path.join(folder_path, "run_final_gmx_simulation")

            if os.path.isdir(run_final_sim_dir):
                # Define the destination directory 'analyze_final_sim'
                dest_dir = os.path.join(folder_path, "analyze_final_sim")
                dest_dir = os.path.join(dest_dir, "NaP_dist_count")

                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)

                # Define paths (adjust these as necessary)
                structure_file = os.path.join(run_final_sim_dir, "md.gro")
                trajectory_file = os.path.join(run_final_sim_dir, "md_center.xtc")
                topology_file = os.path.join(run_final_sim_dir, "md.tpr")
                index_file = os.path.join(dest_dir, "NA_P_index.ndx")
                selection_group = "NA"
                # We'll produce distances like distances_NaP1.xvg -> distances_NaP6.xvg
                output_prefix = "distances_NaP"
                log_file = "distances_NaP_gmx.log"

                # List of phosphorus groups
                p_groups = ["P", "P1", "P2", "P3", "P4", "P5"]

                # Setup logging
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

                # Distance stage: only run if outputs are missing
                distance_files = [f"{output_prefix}{i}.xvg" for i in range(1, len(p_groups) + 1)]
                distances_exist = all(os.path.isfile(f) for f in distance_files)

                if not distances_exist:
                    # Ensure index file exists before running pairdist
                    if not os.path.isfile(index_file):
                        logging.info("Index file not found → creating it before pairdist…")
                        create_index_file(structure_file, index_file)

                    logging.info("Distance .xvg files missing → running gmx pairdist stage…")
                    run_pairdist_commands(
                        trajectory=trajectory_file,
                        topology=topology_file,
                        index_file=index_file,
                        groups=p_groups,
                        selection_group=selection_group,
                        output_prefix=output_prefix,
                    )
                else:
                    logging.info("Found all distance files → skipping pairdist stage.")
                    # Later steps still need an index → ensure it exists
                    if not os.path.isfile(index_file):
                        logging.info("Index file missing, creating it for subsequent steps…")
                        create_index_file(structure_file, index_file)

                # Step 3: Plot the distance data
                plot_all_distances_subplots(
                    output_prefix=output_prefix,
                    num_groups=len(p_groups),
                    plot_filename="all_distances_subplots.pdf",
                )

                # NEW STEPS: Count how many ions are close to each phosphate group
                # Adjust cutoff as desired, e.g., 0.5 nm
                cutoff_distance = 0.5
                ion_count_prefix = "ion_count_"

                run_ion_count_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    selection_group=selection_group,
                    cutoff=cutoff_distance,
                    output_prefix=ion_count_prefix,
                )

                # Plot the ion count data for each group
                plot_ion_counts_subplots(
                    output_prefix=ion_count_prefix,
                    groups=p_groups,
                    plot_filename="ion_counts_subplots.pdf",
                )

                # ----------------------------------------------------------------
                # NEW FEATURE 2: compute & plot excess Na⁺ counts
                # ----------------------------------------------------------------
                excess_prefix = "excess_NaP_"
                buried_prefix = "buriedVol_"
                run_buried_volume_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    ip6_group="Other",
                    cutoff=cutoff_distance,
                    output_prefix=buried_prefix,
                )
                compute_excess_ion_counts(
                    ion_count_prefix=ion_count_prefix,
                    buried_prefix=buried_prefix,
                    groups=p_groups,
                    structure_file=structure_file,
                    cutoff_radius=cutoff_distance,
                    output_prefix=excess_prefix,
                )
                plot_excess_ion_counts_subplots(
                    output_prefix=excess_prefix,
                    groups=p_groups,
                    plot_filename="excess_ion_counts_subplots.pdf",
                )

                # OPTIONAL: If you have an index group for the entire molecule
                # (e.g., "MOL") in your index file, you could do:
                run_ion_count_whole_molecule(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    whole_group="Other",  # must exist in your index
                    selection_group="NA",
                    cutoff=cutoff_distance,
                    output_xvg="ion_count_MOL.xvg",
                )

                plot_whole_molecule_ion_count(
                    xvg_file="ion_count_MOL.xvg",
                    plot_filename="ion_count_whole_mol.pdf",
                )

                logging.info(
                    "All tasks (distance + ion count) completed successfully.\n"
                )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compute Na–P distances, Na+ counts, and excess metrics for an IP6 project."
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

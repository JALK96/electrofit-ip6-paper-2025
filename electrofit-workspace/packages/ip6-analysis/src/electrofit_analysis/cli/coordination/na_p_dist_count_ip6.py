import logging
import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from electrofit_analysis.cli.common import (
    resolve_stage,
    normalize_micro_name,
    resolve_run_and_analyze_dirs,
)
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

def _format_ion_label(ion_name: str) -> str:
    ion_upper = (ion_name or "").strip().upper()
    mapping = {
        "NA": r"Na$^+$",
        "K": r"K$^+$",
        "MG": r"Mg$^{2+}$",
        "CA": r"Ca$^{2+}$",
    }
    return mapping.get(ion_upper, ion_upper)


def plot_ion_counts_subplots(
    output_prefix, groups, plot_filename="ion_counts_subplots.pdf", ion_label: str = "Na+"
):
    """Plot number of ions within cutoff for each group using generic grid plotter."""
    files = [f"{output_prefix}{g}.xvg" for g in groups]
    labels = groups
    plot_timeseries_grid(
        files=files,
        labels=labels,
        ylabel=f"{ion_label} ions",
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

def plot_whole_molecule_ion_count(
    xvg_file,
    plot_filename="ion_count_whole_mol.pdf",
    ion_label: str = "Na+",
):
    """
    Plots a single line: number of ions within cutoff for the entire molecule vs time.
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
        # plt.title('Ion Count Near Entire Molecule')
        plt.xlabel("Time (ns)")
        plt.ylabel(f"# of {ion_label} ions")
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
#  c           = N_total_ion / V_box  (bulk concentration)
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


def _count_ion_atoms_in_gro(gro_file, ion_name: str):
    """
    Counts how many atoms match ion_name in columns 11‑15 (atom name)
    of a .gro file. Works for standard GROMACS naming.
    """
    ion_upper = ion_name.strip().upper()
    count = 0
    with open(gro_file, "r") as fh:
        lines = fh.readlines()
    # skip first 2 title/atom‑number lines, last line is box
    for line in lines[2:-1]:
        atom_name = line[10:15].strip()
        if atom_name.upper().startswith(ion_upper):
            count += 1
    if count == 0:
        raise ValueError(f"No {ion_upper} ions found in structure file; check naming.")
    return count


def compute_excess_ion_counts(
    ion_count_prefix,
    buried_prefix,
    groups,
    structure_file,
    cutoff_radius,
    output_prefix="excess_",
    ion_name: str = "NA",
):
    """
    Excess‑ion factor using *true* solvent‑accessible volume per frame:

        N_excess(t) = N_in_sphere(t) /
                      [ c_ion · (4/3·π·r³ − V_buried(t)) ]

    where V_buried(t) is obtained from gmx sasa (-tv) and already has units nm³.
    """
    import math

    ion_upper = ion_name.strip().upper()
    # Bulk ion concentration
    n_ion_total = _count_ion_atoms_in_gro(structure_file, ion_upper)
    v_box = _parse_box_volume_from_gro(structure_file)
    c_ion = n_ion_total / v_box

    v_sphere = (4.0 / 3.0) * math.pi * cutoff_radius**3
    logging.info("============================================================")
    logging.info("EXCESS‑ION FORMULA (grid/SAV version):")
    logging.info("  N_excess(t) = N_in / [ c_ion · (4/3·π·r³ − V_buried(t)) ]")
    logging.info(f"  r = {cutoff_radius:.3f} nm ; 4/3·π·r³ = {v_sphere:.6f} nm³")
    logging.info(
        f"  c_ion = N_ion_tot / V_Box = {n_ion_total}/{v_box:.6f} = {c_ion:.6f} nm⁻³"
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
                raise ValueError("Frame mismatch between ion-count and buried‑volume files")

            time_ps = na_data[:, 0]
            n_in = na_data[:, 1]
            v_buried = buried_data[:, 1]
            v_eff = v_sphere - v_buried
            # ensure positive effective volume
            v_eff[v_eff <= 1e-9] = 1e-9

            n_solution = c_ion * v_eff
            n_excess = n_in / n_solution

            header = (
                f'@    title "Excess {ion_upper} ions (grid SAV)"\n'
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

def main(
    project_dir: str,
    stage: str = 'final',
    only=None,
    rep: int | None = None,
    ion_name: str = "NA",
    cutoff_nm: float = 0.5,
) -> None:
    """Run ion–P distance and ion-count analysis for all molecules under project_dir.

    Args:
        project_dir: Path to the project root (must contain a 'process' subdirectory).
    """
    ion_name = ion_name.strip().upper()
    ion_label = _format_ion_label(ion_name)
    if cutoff_nm <= 0:
        raise ValueError(f"cutoff_nm must be > 0, got {cutoff_nm}")

    # Define the base process directory
    process_dir = os.path.join(project_dir, "process")

    # Loop through each subdirectory in the process directory
    only_norm = {normalize_micro_name(x) for x in only} if only else None
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)

        if os.path.isdir(folder_path):
            if only_norm and folder_name not in only_norm:
                continue
            run_dir_name, analyze_base = resolve_stage(stage)
            run_sim_dir, analyze_base_dir = resolve_run_and_analyze_dirs(
                Path(folder_path), stage, run_dir_name, analyze_base, rep
            )

            if os.path.isdir(run_sim_dir):
                # Define the destination directory based on stage
                dest_dir = os.path.join(analyze_base_dir, "IonP_dist_count")

                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)

                # Define paths (adjust these as necessary)
                if stage.strip().lower() == "remd":
                    structure_file = os.path.join(run_sim_dir, "remd.gro")
                    trajectory_file = os.path.join(run_sim_dir, "remd_center.xtc")
                    topology_file = os.path.join(run_sim_dir, "remd.tpr")
                else:
                    structure_file = os.path.join(run_sim_dir, "md.gro")
                    trajectory_file = os.path.join(run_sim_dir, "md_center.xtc")
                    topology_file = os.path.join(run_sim_dir, "md.tpr")
                index_file = os.path.join(dest_dir, "Ion_P_index.ndx")
                selection_group = ion_name
                # We'll produce distances like distances_IonP1.xvg -> distances_IonP6.xvg
                output_prefix = "distances_IonP"
                log_file = "distances_IonP_gmx.log"

                # List of phosphorus groups
                p_groups = ["P", "P1", "P2", "P3", "P4", "P5"]

                # Setup logging
                setup_logging(log_file)
                logging.info("Logging is set up.")
                logging.info(
                    "dist-count parameters: ion_name=%s, cutoff=%.3f nm, stage=%s",
                    ion_name,
                    cutoff_nm,
                    stage,
                )

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
                    ion_label=ion_label,
                )

                # NEW STEPS: Count how many ions are close to each phosphate group
                cutoff_distance = cutoff_nm
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
                    ion_label=ion_label,
                )

                # ----------------------------------------------------------------
                # compute & plot excess-ion counts
                # ----------------------------------------------------------------
                excess_prefix = "excess_IonP_"
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
                    ion_name=ion_name,
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
                    selection_group=ion_name,
                    cutoff=cutoff_distance,
                    output_xvg="ion_count_MOL.xvg",
                )

                plot_whole_molecule_ion_count(
                    xvg_file="ion_count_MOL.xvg",
                    plot_filename="ion_count_whole_mol.pdf",
                    ion_label=ion_label,
                )

                logging.info(
                    "All tasks (distance + ion count) completed successfully.\n"
                )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compute ion–P distances, ion counts, buried volume, and excess metrics for the IP6 project."
        )
    )
    parser.add_argument(
        "--project",
        "-p",
        required=True,
        help="Path to the project root directory (contains 'process').",
    )
    parser.add_argument(
        "--stage",
        choices=["final", "sample", "remd"],
        default="final",
        help="Analyze final, sample, or remd trajectories (default: final).",
    )
    parser.add_argument(
        "--ion-name",
        default="NA",
        help="Ion atom name (e.g., NA, K, MG, CA).",
    )
    parser.add_argument(
        "--cutoff-nm",
        type=float,
        default=0.5,
        help="Cutoff in nm for ion-count and buried-volume calculations (default: 0.5).",
    )
    args = parser.parse_args()
    main(
        os.path.abspath(args.project),
        stage=args.stage,
        ion_name=args.ion_name,
        cutoff_nm=args.cutoff_nm,
    )

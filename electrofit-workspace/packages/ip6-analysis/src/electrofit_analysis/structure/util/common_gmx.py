"""Common GMX utilities used across structure analysis scripts.

This module centralizes small helpers that were duplicated in multiple
scripts, such as building index files, running pair distance calculations,
and plotting the resulting time series in subplots.

All functions use logging instead of printing and are side‑effect free
except for writing files and calling external GMX commands via
electrofit.cli.run_commands.run_command.
"""

from __future__ import annotations

import logging
from typing import Iterable
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from electrofit.cli.run_commands import run_command

# -----------------------------
# Shared constants
# -----------------------------
DEFAULT_P_GROUPS = ["P", "P1", "P2", "P3", "P4", "P5"]

def create_index_file(structure_file: str, index_file: str) -> None:
    """Create a GMX index file with phosphate groups P, P1..P5.

    Parameters
    ----------
    structure_file: Path to the structure (.gro)
    index_file: Output .ndx file
    """
    logging.info("Starting index file creation with gmx make_ndx…")
    make_ndx_commands = "a P\na P1\na P2\na P3\na P4\na P5\nq\n"
    cmd = f'printf "{make_ndx_commands}" | gmx make_ndx -f {structure_file} -o {index_file}'
    run_command(cmd)
    logging.info("Index file created successfully.")

def run_pairdist_commands(
    trajectory: str,
    topology: str,
    index_file: str,
    groups: Iterable[str],
    selection_group: str,
    output_prefix: str,
) -> None:
    """Run gmx pairdist for each group against a selection group.

    Produces <output_prefix><i>.xvg files (min distance vs time).
    """
    logging.info("Starting pair distance calculations with gmx pairdist…")
    for i, group in enumerate(groups, start=1):
        output_file = f"{output_prefix}{i}.xvg"
        cmd = (
            f"gmx pairdist -f {trajectory} -s {topology} -n {index_file} "
            f'-ref "{group}" -sel "{selection_group}" -pbc no -rmpbc yes -o {output_file}'
        )
        logging.info(f"Running gmx pairdist for group '{group}' → '{output_file}'…")
        run_command(cmd)
        logging.info(f"Distance calculation for group '{group}' completed.")
    logging.info("All gmx pairdist commands executed successfully.")

def plot_all_distances_subplots(
    output_prefix: str,
    num_groups: int,
    plot_filename: str = "all_distances_subplots.pdf",
    suptitle: str | None = None,
    ion_label: str = "Na+",
) -> None:
    """Plot many distance time series as a grid of subplots.

    Parameters
    ----------
    output_prefix: Prefix of input .xvg files (e.g. 'distances_P')
    num_groups: How many sequential files to plot (…1.xvg, …2.xvg, …)
    plot_filename: Output figure filename
    suptitle: Optional figure title (adds spacing above with y=1.05)
    """
    logging.info("Starting to plot all distance data as subplots…")
    # keep console clean while reading many files
    logging.disable(logging.CRITICAL)

    nrows, ncols = 2, 3  # tuned for 6 groups
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i in range(1, num_groups + 1):
        filename = f"{output_prefix}{i}.xvg"
        ax = axes[i - 1]
        try:
            data = np.loadtxt(filename, comments=("#", "@"))
            time = data[:, 0] / 1000.0  # ps → ns
            distance = data[:, 1]
            mean_distance = float(np.mean(distance))
            group_label = f"P{i}"

            ax.plot(time, distance, color="black", linestyle="-", linewidth=0.5)
            ax.axhline(y=mean_distance, color="red", linestyle="--", linewidth=1.5)
            ax.text(
                0.95,
                0.95,
                f"Mean: {mean_distance:.3f} nm",
                ha="right",
                va="top",
                transform=ax.transAxes,
                color="red",
                fontsize=12,
                bbox=dict(facecolor="white", alpha=0, edgecolor="none"),
            )
            ax.set_title(f"{group_label} - {ion_label}", fontsize=14)
            ax.set_xlabel("Time (ns)")
            if i % ncols == 1:
                ax.set_ylabel("Distance (nm)")
            ax.grid(False)
        except Exception:
            ax.text(
                0.5,
                0.5,
                "Error loading data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="red",
            )
            ax.set_title(f"Distance: P{i} - NA (Error)", fontsize=14)

    plt.tight_layout()
    if suptitle:
        fig.suptitle(suptitle, fontsize=18, y=1.05)
    plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
    logging.disable(logging.NOTSET)
    logging.info("Subplot plotting completed successfully.")

def plot_timeseries_grid(
    files: list[str],
    labels: list[str],
    ylabel: str,
    outfile: str,
    suptitle: str | None = None,
    unit_scale: float = 1000.0,
    nrows: int = 2,
    ncols: int = 3,
    line_color: str | None = None,
    mean_line: bool = True,
    mean_color: str = "red",
    mean_fmt: str = "{mean:.3f}",
    xlabel: str = "Time (ns)",
) -> None:
    """Generic helper to plot multiple .xvg time series as a subplot grid.

    Parameters
    ----------
    files : list[str]
        Paths to .xvg files, one per subplot (order matters).
    labels : list[str]
        Per-subplot labels used as small titles.
    ylabel : str
        Y-axis label shown on the first column.
    outfile : str
        Output figure path.
    suptitle : str | None
        Optional figure suptitle.
    unit_scale : float
        Divide time column by this (e.g. 1000.0 for ps→ns).
    nrows, ncols : int
        Grid layout. Extra axes will be removed if files < nrows*ncols.
    line_color : Optional[str]
        Matplotlib color for the series; if None, uses default cycle.
    mean_line : bool
        Whether to draw a dashed horizontal mean line.
    mean_color : str
        Color of the mean line and mean text.
    mean_fmt : str
        Format string for mean annotation. Must include '{mean}'.
    xlabel : str
        X-axis label for all subplots.
    """
    assert len(files) == len(labels), "files and labels length must match"

    logging.info("Plotting time series grid…")
    logging.disable(logging.CRITICAL)

    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, (path, lab) in enumerate(zip(files, labels)):
        ax = axes[i]
        try:
            data = np.loadtxt(path, comments=("#", "@"))
            time = data[:, 0] / unit_scale if unit_scale else data[:, 0]
            series = data[:, 1]
            mean_val = float(np.mean(series))

            if line_color:
                ax.plot(time, series, linestyle="-", linewidth=0.5, color=line_color)
            else:
                ax.plot(time, series, linestyle="-", linewidth=0.5)

            if mean_line:
                ax.axhline(y=mean_val, color=mean_color, linestyle="--", linewidth=1.5)
                ax.text(
                    0.95,
                    0.95,
                    f"Mean: {mean_fmt.format(mean=mean_val)}",
                    ha="right",
                    va="top",
                    transform=ax.transAxes,
                    color=mean_color,
                    fontsize=12,
                    bbox=dict(facecolor="white", alpha=0.7 if line_color else 0, edgecolor="none"),
                )

            ax.set_title(str(lab), fontsize=14)
            ax.set_xlabel(xlabel)
            # y-label only in first column
            if (i % ncols) == 0:
                ax.set_ylabel(ylabel)
            ax.grid(False)
        except Exception:
            ax.text(
                0.5,
                0.5,
                "Error loading data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="red",
            )
            ax.set_title(f"{lab} (Error)", fontsize=14)

    # remove unused axes
    total = nrows * ncols
    for j in range(len(files), total):
        fig.delaxes(axes[j])

    plt.tight_layout()
    if suptitle:
        fig.suptitle(suptitle, fontsize=18, y=1.05)
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    logging.disable(logging.NOTSET)
    logging.info("Finished plotting time series grid.")

def run_ion_count_commands(
    trajectory,
    topology,
    index_file,
    groups,
    selection_group,
    cutoff=0.5,
    output_prefix="ion_count_",
):
    """
    Uses gmx select to count how many ions (selection_group) are within a given
    distance cutoff of each phosphate group. Produces .xvg files with counts over time.
    """
    logging.info(f"Starting ion count within {cutoff} nm for each group in {groups}...")

    for group in groups:
        # We'll name each output file based on the group, e.g. ion_count_P.xvg, ion_count_P1.xvg, etc.
        output_file = f"{output_prefix}{group}.xvg"
        # Construct the 'gmx select' expression
        select_expr = (
            f'group "{selection_group}" and within {cutoff} of group "{group}"'
        )
        command = (
            f"gmx select -f {trajectory} -s {topology} -n {index_file} "
            f"-select '{select_expr}' "
            f"-os {output_file} "  # .xvg with # of selected atoms vs time
            f"-on dummy_index.ndx "  # optional index file (not used further, but needed by gmx)
            f"2>&1"
        )
        logging.info(
            f"Running gmx select for group '{group}', output -> '{output_file}'..."
        )
        run_command(command)

    logging.info("All ion count commands executed successfully.")

# Optional: Count how many ions are near the *entire molecule* (e.g., 'MOL' group)
def run_ion_count_whole_molecule(
    trajectory,
    topology,
    index_file,
    whole_group="Other",
    selection_group="NA",
    cutoff=0.5,
    output_xvg="ion_count_MOL.xvg",
):
    """
    If you have an index group for the entire molecule (named e.g. 'MOL'),
    this function will count how many 'selection_group' atoms are within
    `cutoff` nm of that entire molecule.
    """
    logging.info(f"Starting ion count for entire molecule group '{whole_group}'...")
    select_expr = (
        f'group "{selection_group}" and within {cutoff} of group "{whole_group}"'
    )
    command = (
        f"gmx select -f {trajectory} -s {topology} -n {index_file} "
        f"-select '{select_expr}' "
        f"-os {output_xvg} "
        f"-on dummy_mol_index.ndx "
        f"2>&1"
    )
    run_command(command)
    logging.info("Ion count for entire molecule completed successfully.")

def run_ip6_exclusion_count_commands(
    trajectory,
    topology,
    index_file,
    groups,
    ip6_group="Other",
    cutoff=0.5,
    output_prefix="ip6_count_",
):
    """
    Uses gmx select to count how many IP6 atoms (or whatever group name is
    passed via `ip6_group`) are found within `cutoff` nm of each phosphate
    group *for every frame* of the trajectory.

    The output is written to <output_prefix><group>.xvg files in exactly the
    same format that `run_ion_count_commands` produces, enabling frame‑wise
    combination with the Na⁺ counts.

    Parameters
    ----------
    trajectory : str
        Path to the .xtc/.trr trajectory file.
    topology : str
        Path to the .tpr (or any structure readable by GROMACS).
    index_file : str
        .ndx file containing both the phosphate groups and the IP6 group.
    groups : list[str]
        List of phosphate group names (e.g. ["P", "P1", ...]).
    ip6_group : str, optional
        Name of the index group that represents all IP6 atoms. Default "Other".
    cutoff : float, optional
        Cut‑off radius in nm (must be identical to the one used for the Na⁺
        counts). Default 0.5 nm.
    output_prefix : str, optional
        Prefix for the generated .xvg files. Default "ip6_count_".
    """
    logging.info(
        f"Starting IP6 exclusion count within {cutoff} nm for groups: {groups}"
    )
    for group in groups:
        out_file = f"{output_prefix}{group}.xvg"
        select_expr = f'group "{ip6_group}" and within {cutoff} of group "{group}"'
        cmd = (
            f"gmx select -f {trajectory} -s {topology} -n {index_file} "
            f"-select '{select_expr}' "
            f"-os {out_file} "  # writes time‑series of #selected atoms
            f"-on dummy_ip6_index.ndx "  # dummy index output (ignored later)
            f"2>&1"
        )
        logging.info(f'gmx select (IP6) for group "{group}" → "{out_file}"')
        run_command(cmd)
    logging.info("IP6 exclusion counting completed for all groups.")

def run_buried_volume_commands(
    trajectory,
    topology,
    index_file,
    groups,
    ip6_group="Other",
    cutoff=0.5,
    output_prefix="buriedVol_",
    ndots=384,
):
    traj = str(Path(trajectory))
    top  = str(Path(topology))
    ndx  = str(Path(index_file))

    # quick sanity checks
    for p in (traj, top, ndx):
        if not Path(p).is_file():
            raise FileNotFoundError(f"Required input not found: {p}")
    if not groups:
        raise ValueError("No groups provided for buried-volume calculation.")
    if ndots <= 0:
        raise ValueError("ndots must be > 0")

    cmd_common = [
        "gmx", "sasa",
        "-f", traj,
        "-s", top,
        "-n", ndx,
        "-probe", "0",
        "-ndots", str(ndots),
        "-quiet",
    ]

    logging.info("Starting buried-volume calculations…")
    for grp in groups:
        outfile = f"{output_prefix}{grp}.xvg"
        sel = f'group "{ip6_group}" and within {cutoff:g} of group "{grp}"'
        cmd = cmd_common + ["-surface", sel, "-tv", outfile]
        logging.info("  %s  ←  %s", outfile, sel)
        run_command(cmd)
    logging.info("Finished.")

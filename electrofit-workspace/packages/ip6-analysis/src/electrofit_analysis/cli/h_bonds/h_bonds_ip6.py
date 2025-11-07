import math
import argparse
import os
import re
import subprocess
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from electrofit.cli.run_commands import run_command
from electrofit.io.files import find_file_with_extension

sns.set_context("talk")


def load_hb_num_xvg(filename):
    """
    Load data from an .xvg file, skipping lines that start with @ or #.

    Parameters:
        filename (str): Path to the .xvg file.

    Returns:
        np.ndarray: 2D array with the data.
    """
    data = []
    with open(filename, "r") as file:
        for line in file:
            if line.startswith(("@", "#")):
                continue  # Skip header lines
            parts = line.strip().split()
            if len(parts) >= 2:  # Ensure there are at least two columns
                try:
                    # Convert all parts to float
                    floats = [float(part) for part in parts]
                    data.append(floats)
                except ValueError:
                    # Skip lines that don't contain valid floats
                    continue
    return np.array(data)


def refine_atom_name(atom):
    """
    Refines the atom name based on specified rules:
    - If the atom is named 'O', change it to 'O1'.
    - If the atom is named 'O<number>', increment the number by 1 (e.g., 'O10' to 'O11').

    Parameters:
    - atom (str): Original atom name.

    Returns:
    - str: Refined atom name.
    """
    # Match atoms like 'O' (no number)
    if re.fullmatch(r"[A-Za-z]+", atom):
        return atom + "1"

    # Match atoms like 'O10', 'O2', etc.
    match = re.fullmatch(r"([A-Za-z]+)(\d+)", atom)
    if match:
        name = match.group(1)
        number = int(match.group(2)) + 1  # Increment the number by 1
        return f"{name}{number}"

    # If atom name doesn't match expected patterns, return as is
    return atom


def parse_xpm(file_path):
    """
    Parses an XPM file and converts it to a binary NumPy array.

    Parameters:
    - file_path: str, path to the XPM file.

    Returns:
    - data_matrix: np.ndarray, binary matrix representing hydrogen bonds.
    - metadata: dict, contains title, labels, etc.
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Initialize variables
    metadata = {}
    data_lines = []
    color_map = {}
    header_found = False
    data_started = False

    for line in lines:
        line = line.strip()

        # Extract metadata
        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ":" in comment:
                key, value = comment.split(":", 1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        # Identify the start of data
        if line.startswith("static char"):
            continue  # Skip the static declaration line

        if line.startswith('"') and not header_found:
            # Header line containing dimensions and color info
            header = line.strip('",')
            tokens = header.split()
            if len(tokens) >= 4:
                width = int(tokens[0])
                height = int(tokens[1])
                num_colors = int(tokens[2])
                chars_per_pixel = int(tokens[3])
                header_found = True
            continue

        if header_found and not data_started:
            # Read color definitions
            color_def = line.strip('",')
            # Example: "   c #FFFFFF " or "o  c #FF0000 "
            match = re.match(r"(.{%d})\s+c\s+(\S+)" % chars_per_pixel, color_def)
            if match:
                symbol = match.group(1)
                color = match.group(2)
                color_map[symbol] = color
            if len(color_map) == num_colors:
                data_started = True
            continue

        if data_started:
            # Read data lines
            if line.startswith('"'):
                data_line = line.strip('",')
                data_lines.append(data_line)

    # Convert data lines to binary matrix
    data_matrix = np.zeros((height, width), dtype=int)

    for y, line in enumerate(data_lines):
        for x, char in enumerate(line):
            if char in color_map:
                if color_map[char] == "#FF0000":  # Present bond
                    data_matrix[y, x] = 1
                else:
                    data_matrix[y, x] = 0
            else:
                data_matrix[y, x] = 0  # Default to 0 if unknown

    return data_matrix, metadata


def analyze_hydrogen_bonds(data_matrix, metadata):
    """
    Analyzes the hydrogen bond data.

    Parameters:
    - data_matrix: np.ndarray, binary matrix representing hydrogen bonds.
    - metadata: dict, contains title, labels, etc.

    Returns:
    - analysis_results: dict, contains various analysis metrics.
    """
    analysis_results = {}

    # Total number of hydrogen bonds over time
    hbonds_over_time = np.sum(data_matrix, axis=0)

    # Total occurrence of each hydrogen bond
    hbonds_per_index = np.sum(data_matrix, axis=1)

    # Lifetime distribution (how long each bond persists)
    # This is more complex and would require tracking continuous sequences of 1s
    lifetimes = []
    for bond in data_matrix:
        current_lifetime = 0
        bond_lifetimes = []
        for state in bond:
            if state == 1:
                current_lifetime += 1
            else:
                if current_lifetime > 0:
                    bond_lifetimes.append(current_lifetime)
                    current_lifetime = 0
        if current_lifetime > 0:
            bond_lifetimes.append(current_lifetime)
        lifetimes.append(bond_lifetimes)
    analysis_results["hbonds_over_time"] = hbonds_over_time
    analysis_results["hbonds_per_index"] = hbonds_per_index
    analysis_results["lifetimes"] = lifetimes

    return analysis_results


def visualize_data_donor_accpetor_pair(
    xpm_file="intra_hb_matrix.xpm",
    hbond_df=None,
    output_prefix="intra",
    bin_size_ns=0.1,
    time_per_frame_ns=0.01,
):
    """
    Visualizes the hydrogen bond data and analysis with time binning.

    Parameters:
    - xpm_file: str, path to the XPM file.
    - hbond_df: pd.DataFrame, DataFrame containing donor and acceptor labels.
    - output_prefix: str, prefix for the output file names.
    - bin_size_ns: float, size of each time bin in nanoseconds.
    - time_per_frame_ns: float, time duration of each frame in nanoseconds.

    Returns:
    - None
    """
    data_matrix, metadata = parse_xpm(xpm_file)
    analysis_results = analyze_hydrogen_bonds(data_matrix, metadata)

    matrix_shape = np.shape(data_matrix)
    print(f"Data matrix shape: {matrix_shape}")

    # Extract donor-acceptor labels from hbond_df
    if hbond_df is not None:
        hbond_labels = hbond_df["donor"] + "-" + hbond_df["acceptor"]
    else:
        hbond_labels = [str(i) for i in range(matrix_shape[0])]

    # Verify alignment
    if len(hbond_labels) != data_matrix.shape[0]:
        print(
            f"Number of labels: {len(hbond_labels)}, Number of rows in data_matrix: {data_matrix.shape[0]}"
        )
        raise ValueError(
            "Mismatch between number of hydrogen bonds in hbond_df and data_matrix."
        )

    # Calculate the number of frames per bin
    frames_per_bin = int(bin_size_ns / time_per_frame_ns)
    if frames_per_bin <= 0:
        raise ValueError("Bin size must be larger than the time per frame.")

    # Number of bins
    num_bins = data_matrix.shape[1] // frames_per_bin
    if data_matrix.shape[1] % frames_per_bin != 0:
        num_bins += 1  # Include partial bin

    # Initialize binned data matrix
    binned_matrix = np.zeros((data_matrix.shape[0], num_bins))

    # Aggregate data into bins
    for i in range(num_bins):
        start_idx = i * frames_per_bin
        end_idx = min((i + 1) * frames_per_bin, data_matrix.shape[1])
        bin_data = data_matrix[:, start_idx:end_idx]
        # Calculate the fraction of time the bond is present in the bin
        binned_matrix[:, i] = np.mean(bin_data, axis=1)

    # Adjust figure size based on the number of hydrogen bonds
    fig_height = max(6, 0.3 * len(hbond_labels))
    plt.figure(figsize=(8, fig_height))

    # Heatmap of hydrogen bonds
    plt.imshow(
        binned_matrix, aspect="auto", cmap="Reds", origin="lower", vmin=0, vmax=1
    )
    plt.title(f"{metadata.get('title', 'Hydrogen Bond Existence Map')} (Binned)")
    plt.xlabel("Time (ns)")
    plt.ylabel("Donor-Acceptor Pairs")

    # Adjust time axis labels
    bin_times_ns = np.arange(num_bins) * bin_size_ns
    # Determine the number of ticks you want to display
    num_ticks = 5

    # Generate tick positions evenly spaced across the bins
    tick_positions = np.linspace(0, num_bins - 1, num_ticks, dtype=int)

    # Generate tick labels corresponding to the bin times
    tick_labels = [f"{bin_times_ns[pos]:.1f}" for pos in tick_positions]

    # Set the x-ticks with the selected positions and labels
    plt.xticks(tick_positions, labels=tick_labels)

    # Set y-ticks with donor-acceptor labels
    plt.yticks(np.arange(len(hbond_labels)), hbond_labels)

    # Adjust y-axis tick label font size if necessary
    plt.tick_params(axis="y", which="major", labelsize=8)

    # Color bar representing fraction of bond presence
    cbar = plt.colorbar(label="Fraction of Bond Presence")
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(["0%", "25%", "50%", "75%", "100%"])

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_existence_map_binned.pdf")
    plt.close()

    print("Hydrogen bonds per index:", analysis_results["hbonds_per_index"])

    # Histogram of Hydrogen Bonds per Index (Horizontal Bar Plot)
    plt.figure(figsize=(8, fig_height))
    plt.barh(
        range(len(analysis_results["hbonds_per_index"])),
        analysis_results["hbonds_per_index"],
        color="darkred",
    )
    plt.title("Total Occurrence of Each Hydrogen Bond")
    plt.xlabel("Occurrences")
    plt.ylabel("Donor-Acceptor Pairs")
    plt.yticks(np.arange(len(hbond_labels)), hbond_labels)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_occurrences.pdf")
    plt.close()

    # Lifetime Distribution
    # Suppose each frame is 10 ps = 0.01 ns
    time_per_frame = 0.01  # ns per frame

    # Lifetime Distribution
    # Flatten all lifetimes (still in frames)
    all_lifetimes_frames = [
        lifetime
        for bond_lifetimes in analysis_results["lifetimes"]
        for lifetime in bond_lifetimes
    ]

    # Convert frames --> ns
    all_lifetimes_ns = [lifetime * time_per_frame for lifetime in all_lifetimes_frames]

    # Define histogram bins in ns.
    # Bin every 0.01 ns:

    max_lifetime = max(all_lifetimes_ns)
    bin_width = 0.01  # 0.01 ns bin width (adjust to taste)
    num_bins = int(max_lifetime / bin_width) + 2  # +2 for a little padding
    bins = np.arange(0, num_bins * bin_width, bin_width)

    plt.figure(figsize=(8, 6))
    # Plot the histogram of lifetimes in ns
    plt.hist(all_lifetimes_ns, bins=bins, edgecolor=None, color="darkred")

    plt.title("Hydrogen Bond Lifetime Distribution")
    plt.xlabel("Lifetime (ns)")  # Note the change here!
    plt.ylabel("Frequency")
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_lifetime.pdf")
    plt.close()

# --- occurences in ns + total number of H-bonds (event based)

def visualize_data_donor_acceptor(
    xpm_file="inter_hb_matrix.xpm",
    hbond_df=None,
    output_prefix="inter",
    bin_size_ns=0.1,
    time_per_frame_ns=0.01,
):
    """
    Visualizes the hydrogen bond data and analysis with time binning.
    Aggregates data per oxygen atom in the molecule (donor or acceptor), showing both the total time and the number of distinct hydrogen bond occurrences
    associated with that oxygen atom, regardless of the water molecules they are interacting with.

    Parameters:
    - xpm_file: str, path to the XPM file.
    - hbond_df: pd.DataFrame, DataFrame containing donor and acceptor labels.
    - output_prefix: str, prefix for the output file names.
    - bin_size_ns: float, size of each time bin in nanoseconds.
    - time_per_frame_ns: float, time duration of each frame in nanoseconds.

    Returns:
    - None
    """
    # Ensure necessary functions are defined or imported
    # Assuming parse_xpm and analyze_hydrogen_bonds are defined elsewhere
    data_matrix, metadata = parse_xpm(xpm_file)
    analyze_hydrogen_bonds(data_matrix, metadata)

    matrix_shape = np.shape(data_matrix)
    print(f"Data matrix shape: {matrix_shape}")

    # Verify that hbond_df is provided
    if hbond_df is None:
        raise ValueError("hbond_df must be provided.")

    # Ensure that 'idx' in hbond_df matches the rows in data_matrix
    hbond_df = hbond_df.reset_index(drop=True)
    hbond_df["idx"] = hbond_df.index

    # Identify oxygen atoms in the molecule
    def is_molecule_atom(atom_name):
        # Assuming molecule atoms start with 'O' followed by digits
        return bool(re.match(r"^O\d+$", atom_name))

    # Filter donors and acceptors to molecule oxygens
    molecule_donors = hbond_df["donor"][hbond_df["donor"].apply(is_molecule_atom)]
    molecule_acceptors = hbond_df["acceptor"][
        hbond_df["acceptor"].apply(is_molecule_atom)
    ]

    # Get unique oxygen atoms in molecule
    molecule_oxygens = set(molecule_donors).union(set(molecule_acceptors))

    # Create mapping from oxygen to indices
    oxygen_to_indices = {}
    for oxygen in molecule_oxygens:
        indices = hbond_df[
            (hbond_df["donor"] == oxygen) | (hbond_df["acceptor"] == oxygen)
        ]["idx"].tolist()
        oxygen_to_indices[oxygen] = indices

    # Calculate the number of frames per bin
    frames_per_bin = int(bin_size_ns / time_per_frame_ns)
    if frames_per_bin <= 0:
        raise ValueError("Bin size must be larger than the time per frame.")

    # Number of bins
    num_bins = data_matrix.shape[1] // frames_per_bin
    if data_matrix.shape[1] % frames_per_bin != 0:
        num_bins += 1  # Include partial bin

    # Initialize binned data matrix
    binned_matrix = np.zeros((data_matrix.shape[0], num_bins))

    # Aggregate data into bins
    for i in range(num_bins):
        start_idx = i * frames_per_bin
        end_idx = min((i + 1) * frames_per_bin, data_matrix.shape[1])
        bin_data = data_matrix[:, start_idx:end_idx]
        # Calculate the fraction of time the bond is present in the bin
        binned_matrix[:, i] = np.mean(bin_data, axis=1)

    # Now, aggregate data per oxygen atom for the heatmap
    num_oxygens = len(oxygen_to_indices)
    aggregated_binned_matrix = np.zeros((num_oxygens, num_bins))

    oxygen_list = list(oxygen_to_indices.keys())
    for i, oxygen in enumerate(oxygen_list):
        indices = oxygen_to_indices[oxygen]
        oxygen_data = binned_matrix[
            indices, :
        ]  # shape: (num_hbonds_for_oxygen, num_bins)
        # Aggregate over the rows
        aggregated_data = np.max(
            oxygen_data, axis=0
        )  # For binary data, max indicates presence
        aggregated_binned_matrix[i, :] = aggregated_data

    # Prepare roles for coloring y-tick labels
    oxygen_roles = {}
    for oxygen in oxygen_list:
        is_donor = any(hbond_df["donor"] == oxygen)
        is_acceptor = any(hbond_df["acceptor"] == oxygen)
        if is_donor and is_acceptor:
            role = "both"
        elif is_donor:
            role = "donor"
        elif is_acceptor:
            role = "acceptor"
        else:
            role = "unknown"  # This shouldn't happen
        oxygen_roles[oxygen] = role

    # Adjust figure size based on the number of oxygens
    fig_height = max(6, 0.3 * len(oxygen_list))
    plt.figure(figsize=(8, fig_height))

    # Heatmap of hydrogen bonds per oxygen atom
    plt.imshow(
        aggregated_binned_matrix,
        aspect="auto",
        cmap="Reds",
        origin="lower",
        vmin=0,
        vmax=1,
    )
    plt.title(f"{metadata.get('title', 'Hydrogen Bond Existence Map')} (Binned)")
    plt.xlabel("Time (ns)")
    plt.ylabel("Oxygen Atoms")

    # Adjust time axis labels
    bin_times_ns = np.arange(num_bins) * bin_size_ns
    num_ticks = 5  # Adjust as needed
    tick_positions = np.linspace(0, num_bins - 1, num_ticks, dtype=int)
    tick_labels = [f"{bin_times_ns[pos]:.1f}" for pos in tick_positions]
    plt.xticks(tick_positions, labels=tick_labels)

    # Set y-ticks with oxygen labels
    plt.yticks(np.arange(len(oxygen_list)), oxygen_list)

    # Adjust y-axis tick label font size if necessary
    plt.tick_params(axis="y", which="major", labelsize=8)

    # Color y-tick labels based on role
    ax = plt.gca()
    yticks = ax.get_yticklabels()
    for tick_label in yticks:
        text = tick_label.get_text()
        role = oxygen_roles.get(text, "unknown")
        if role == "donor":
            tick_label.set_color("blue")
        elif role == "acceptor":
            tick_label.set_color("green")
        elif role == "both":
            tick_label.set_color("purple")

    # Color bar representing fraction of bond presence
    cbar = plt.colorbar(label="Fraction of Bond Presence")
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(["0%", "25%", "50%", "75%", "100%"])

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_existence_map_binned.pdf")
    plt.close()

    # ===== Occurrence Counting =====

    # --- Time-Resolved Occurrence ---
    # Calculate the total time each hydrogen bond is present
    hbonds_time_resolved = np.sum(binned_matrix, axis=1) * bin_size_ns  # in ns

    # --- Event-Based Occurrence ---
    # Define an occurrence as a transition from 0 to 1
    # Vectorized approach for efficiency
    # Prepend a column of zeros to detect transitions at the first frame
    bond_transitions = np.diff(data_matrix, prepend=0, axis=1)
    hbonds_event_based = np.sum(
        bond_transitions == 1, axis=1
    )  # Number of events per hbond

    # Assign both counts to the DataFrame
    hbond_df["count_time_ns"] = hbonds_time_resolved
    hbond_df["count_events"] = hbonds_event_based

    # --- Aggregating Counts per Oxygen Atom ---

    # Melt the DataFrame for time-resolved counts
    melted_time = hbond_df.melt(
        id_vars=["idx", "count_time_ns"],
        value_vars=["donor", "acceptor"],
        value_name="oxygen",
    )
    melted_time = melted_time.dropna(subset=["oxygen"])
    oxygen_time_occurrences = melted_time.groupby("oxygen")["count_time_ns"].sum()
    oxygen_time_occurrences = oxygen_time_occurrences.reindex(
        molecule_oxygens, fill_value=0
    )

    # Melt the DataFrame for event-based counts
    melted_events = hbond_df.melt(
        id_vars=["idx", "count_events"],
        value_vars=["donor", "acceptor"],
        value_name="oxygen",
    )
    melted_events = melted_events.dropna(subset=["oxygen"])
    oxygen_event_occurrences = melted_events.groupby("oxygen")["count_events"].sum()
    oxygen_event_occurrences = oxygen_event_occurrences.reindex(
        molecule_oxygens, fill_value=0
    )

    # --- Plotting Histograms ---

    # Define a combined figure with two subplots
    fig, axes = plt.subplots(
        nrows=2, ncols=1, figsize=(10, 2 * fig_height), constrained_layout=True
    )

    # --- Histogram 1: Total Time in ns ---
    axes[0].barh(
        range(len(oxygen_time_occurrences)),
        oxygen_time_occurrences.values,
        color="darkred",
    )
    axes[0].set_title("Total Time of Hydrogen Bonds per Oxygen Atom")
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("Oxygen Atoms")
    axes[0].set_yticks(np.arange(len(oxygen_list)))
    axes[0].set_yticklabels(oxygen_list)
    axes[0].grid(False)

    # Color y-tick labels based on role
    yticks = axes[0].get_yticklabels()
    for tick_label in yticks:
        text = tick_label.get_text()
        role = oxygen_roles.get(text, "unknown")
        if role == "donor":
            tick_label.set_color("blue")
        elif role == "acceptor":
            tick_label.set_color("green")
        elif role == "both":
            tick_label.set_color("purple")

    # --- Histogram 2: Total Number of Occurrences ---
    axes[1].barh(
        range(len(oxygen_event_occurrences)),
        oxygen_event_occurrences.values,
        color="darkgreen",
    )
    axes[1].set_title("Total Number of Hydrogen Bond Occurrences per Oxygen Atom")
    axes[1].set_xlabel("Number of Occurrences")
    axes[1].set_ylabel("Oxygen Atoms")
    axes[1].set_yticks(np.arange(len(oxygen_list)))
    axes[1].set_yticklabels(oxygen_list)
    axes[1].grid(False)

    # Color y-tick labels based on role
    yticks = axes[1].get_yticklabels()
    for tick_label in yticks:
        text = tick_label.get_text()
        role = oxygen_roles.get(text, "unknown")
        if role == "donor":
            tick_label.set_color("blue")
        elif role == "acceptor":
            tick_label.set_color("green")
        elif role == "both":
            tick_label.set_color("purple")

    # Save the combined histograms
    plt.savefig(f"{output_prefix}_hb_occurrences.pdf")
    plt.close()

    print(f"Visualization completed. Files saved with prefix '{output_prefix}_hb_'.")


# ---- stop occurences


def plot_hb_dist_xvg(file_path, plot_file_name="hb_distribution.pdf"):
    """
    Reads an XVG file and plots the data using matplotlib.

    Parameters:
    - file_path: str, path to the XVG file.

    Returns:
    - None
    """
    title = "XVG Plot"
    x_label = "X-axis"
    y_label = "Y-axis"
    data = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Handle metadata lines
            if line.startswith("#"):
                continue  # Ignore comments
            elif line.startswith("@"):
                tokens = line.split()
                if len(tokens) >= 3:
                    if tokens[1] == "title":
                        title = " ".join(tokens[2:]).strip('"')
                    elif tokens[1] == "xaxis":
                        if tokens[2] == "label":
                            x_label = " ".join(tokens[3:]).strip('"')
                    elif tokens[1] == "yaxis":
                        if tokens[2] == "label":
                            y_label = " ".join(tokens[3:]).strip('"')
            else:
                # Data lines
                try:
                    x, y = map(float, line.split())
                    data.append((x, y))
                except ValueError:
                    # Handle lines that do not have exactly two float numbers
                    continue

    if not data:
        print("No data found in the XVG file.")
        return

    # Convert data to NumPy arrays for easier handling
    data = np.array(data)
    x = data[:, 0]
    y = data[:, 1]

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker="o", linestyle="-", color="darkblue")
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{plot_file_name}.pdf", format="pdf")


def add_dashed_hbonds(
    svg_content,
    hbond_pairs,
    atom_coords,
    dasharray="5,5",
    stroke_width="2",
    angle_threshold=5,
    offset_step=5,
):
    """
    Adds dashed lines between specified atom pairs in the SVG content, adjusting overlapping lines to prevent them from overlaying.
    The lines have gradients to indicate direction from black to light grey.

    Parameters:
        svg_content (str): The original SVG content.
        hbond_pairs (list of tuples): List of atom index pairs representing H-bonds (donor_idx, acceptor_idx).
        atom_coords (dict): Dictionary mapping atom indices to (x, y) drawing coordinates.
        dasharray (str): Dash pattern for dashed lines.
        stroke_width (str): Width of the dashed lines.
        angle_threshold (float): Threshold angle in degrees to consider lines as overlapping.
        offset_step (float): Step size for offsetting overlapping lines.

    Returns:
        str: Modified SVG content with adjusted dashed H-bond lines inserted after </rect>.
    """
    import uuid

    # Collect all lines with their data
    lines_data = []

    # Prepare gradient definitions
    gradient_defs = []

    for idx, (donor_idx, acceptor_idx) in enumerate(hbond_pairs):
        if donor_idx not in atom_coords or acceptor_idx not in atom_coords:
            print(
                f"Warning: Atom indices {donor_idx} or {acceptor_idx} not found in molecule."
            )
            continue

        x1, y1 = atom_coords[donor_idx]
        x2, y2 = atom_coords[acceptor_idx]

        # Calculate direction vector
        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy)
        if length == 0:
            print(
                f"Warning: Zero length line between atom {donor_idx} and {acceptor_idx}. Line skipped."
            )
            continue
        dir_vector = (dx / length, dy / length)

        # Prepare data for this line
        line_data = {
            "donor_idx": donor_idx,
            "acceptor_idx": acceptor_idx,
            "x1": x1,
            "y1": y1,
            "x2": x2,
            "y2": y2,
            "dir_vector": dir_vector,
            "offset": 0,  # Will be updated if overlap detected
            "gradient_id": f"grad_{uuid.uuid4().hex}",  # Unique gradient ID
        }
        lines_data.append(line_data)

    # Detect overlapping lines and apply offsets
    for line1, line2 in combinations(lines_data, 2):
        # Check if lines share a common atom
        common_atoms = {line1["donor_idx"], line1["acceptor_idx"]}.intersection(
            {line2["donor_idx"], line2["acceptor_idx"]}
        )
        if common_atoms:
            # Calculate angle between direction vectors
            dot_product = (
                line1["dir_vector"][0] * line2["dir_vector"][0]
                + line1["dir_vector"][1] * line2["dir_vector"][1]
            )
            angle_rad = math.acos(
                min(max(dot_product, -1.0), 1.0)
            )  # Clamp value to valid range
            angle_deg = math.degrees(angle_rad)
            # Check if lines are nearly colinear
            if (
                abs(angle_deg) < angle_threshold
                or abs(angle_deg - 180) < angle_threshold
            ):
                # Apply offsets to both lines
                line1["offset"] += offset_step
                line2["offset"] -= offset_step  # Alternate offset direction

    # Create SVG elements for lines and gradient definitions
    dashed_lines = []
    for line_data in lines_data:
        x1, y1 = line_data["x1"], line_data["y1"]
        x2, y2 = line_data["x2"], line_data["y2"]
        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy)
        nx = -dy / length
        ny = dx / length
        # Apply offset
        offset_distance = line_data["offset"]
        x1_adj = x1 + nx * offset_distance
        y1_adj = y1 + ny * offset_distance
        x2_adj = x2 + nx * offset_distance
        y2_adj = y2 + ny * offset_distance

        # Use black and light grey for gradient
        donor_color_hex = "#000000"  # Black
        acceptor_color_hex = "#D3D3D3"  # Light grey

        # Define gradient for this line
        gradient_def = f'''
        <linearGradient id="{line_data["gradient_id"]}" gradientUnits="userSpaceOnUse" x1="{x1_adj}" y1="{y1_adj}" x2="{x2_adj}" y2="{y2_adj}">
            <stop offset="0%" stop-color="{donor_color_hex}" />
            <stop offset="100%" stop-color="{acceptor_color_hex}" />
        </linearGradient>
        '''
        gradient_defs.append(gradient_def)

        # Create the dashed line SVG element with gradient stroke
        line = f'<line x1="{x1_adj}" y1="{y1_adj}" x2="{x2_adj}" y2="{y2_adj}" stroke="url(#{line_data["gradient_id"]})" stroke-width="{stroke_width}" stroke-dasharray="{dasharray}" />'
        dashed_lines.append(line)

    if dashed_lines:
        # Join all gradient definitions
        gradients_svg = "\n".join(gradient_defs)

        # Insert gradient definitions into the SVG before <!-- END OF HEADER --> comment
        header_end_pos = svg_content.find("<!-- END OF HEADER -->")
        if header_end_pos != -1:
            svg_content = (
                svg_content[:header_end_pos]
                + f"\n<defs>\n{gradients_svg}\n</defs>\n"
                + svg_content[header_end_pos:]
            )
        else:
            # Fallback if <!-- END OF HEADER --> not found
            insert_pos = svg_content.find(">") + 1
            svg_content = (
                svg_content[:insert_pos]
                + f"\n<defs>\n{gradients_svg}\n</defs>\n"
                + svg_content[insert_pos:]
            )
            print(
                "Warning: <!-- END OF HEADER --> comment not found. Inserting <defs> after <svg> tag."
            )

        # Join all dashed line elements
        dashed_svg_content = "\n".join(dashed_lines)

        # Find the position after the closing </rect> tag
        insert_pos = svg_content.find("</rect>") + len("</rect>")
        if insert_pos == -1 + len("</rect>"):
            # </rect> not found, default to inserting after </defs>
            insert_pos = svg_content.find("</defs>") + len("</defs>")
            if insert_pos == -1 + len("</defs>"):
                # Neither </rect> nor </defs> found, insert after opening <svg> tag
                insert_pos = svg_content.find(">") + 1
                print(
                    "Warning: Neither </rect> nor </defs> tag found. Inserting dashed lines after <svg> tag."
                )
        svg_content = (
            svg_content[:insert_pos]
            + f"\n{dashed_svg_content}\n"
            + svg_content[insert_pos:]
        )
        print(f"Added {len(dashed_lines)} dashed H-bond line(s) with gradients to SVG.")
    else:
        print("No dashed lines were added to the SVG content.")

    return svg_content


# Legend
def create_legend(x, y, legend_entries, font_size=12):
    """
    Creates SVG elements for the legend.

    Parameters:
        x (float): The x-coordinate of the top-left corner of the legend.
        y (float): The y-coordinate of the top-left corner of the legend.
        legend_entries (list of dict): List of legend entries to include.
        font_size (int): Font size for the legend text.

    Returns:
        str: SVG content representing the legend.
    """
    legend_svg = []
    entry_height = font_size + 10
    for i, entry in enumerate(legend_entries):
        entry_y = y + i * entry_height

        if "dasharray" in entry:
            # Draw a dashed line for the hydrogen bond legend entry
            line = f'<line x1="{x}" y1="{entry_y + font_size / 2}" x2="{x + 20}" y2="{entry_y + font_size / 2}" stroke="{entry["color"]}" stroke-width="2" stroke-dasharray="{entry["dasharray"]}" />'
            legend_svg.append(line)
        else:
            # Draw a colored rectangle for the atom legend entries
            rect = f'<rect x="{x}" y="{entry_y}" width="20" height="{font_size}" fill="{entry["color"]}" stroke="black" />'
            legend_svg.append(rect)

        # Add the text label
        text_x = x + 25
        text_y = entry_y + font_size
        text = f'<text x="{text_x}" y="{text_y}" font-size="{font_size}" font-family="Arial">{entry["label"]}</text>'
        legend_svg.append(text)

    return "\n".join(legend_svg)


def add_legend_to_svg(svg_content, legend_svg_content):
    """
    Inserts the legend SVG content into the main SVG content.

    Parameters:
        svg_content (str): The original SVG content.
        legend_svg_content (str): The SVG content for the legend.

    Returns:
        str: Modified SVG content with the legend included.
    """
    # Find the position before the closing </svg> tag
    insert_pos = svg_content.rfind("</svg>")
    if insert_pos == -1:
        print("Warning: </svg> tag not found. Legend will not be added.")
        return svg_content

    # Insert the legend before the </svg> tag
    svg_content = (
        svg_content[:insert_pos]
        + f"\n{legend_svg_content}\n"
        + svg_content[insert_pos:]
    )
    print("Legend added to SVG.")
    return svg_content


def parse_hbond_log_to_dataframe(file_path):
    """
    Parses a GROMACS .log file to extract donor-acceptor pairs with indices and refined atom names.

    Parameters:
    - file_path (str): Path to the .log file.

    Returns:
    - pd.DataFrame: DataFrame containing idx, donor, and acceptor columns.
    """
    hbond_pairs = []

    # Regular expression patterns
    # Pattern to match lines with donor and acceptor information
    # Example line: I211O10              -      I211O2
    line_pattern = re.compile(r"^\s*(\S+)\s+-\s+(\S+)\s*$")

    # Pattern to extract atom name and index from a string like 'I211O10'
    # Assuming residue identifier is 'I211' and atom name is 'O10'
    atom_pattern = re.compile(r"^[A-Za-z]+\d+([A-Za-z]+\d*)$")

    with open(file_path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            # Skip empty lines and header lines
            if (
                not line
                or line.startswith("#")
                or line.startswith('"""')
                or line.startswith("*")
            ):
                continue

            # Match the line with donor and acceptor
            match = line_pattern.match(line)
            if match:
                donor_full = match.group(1)  # e.g., 'I211O10'
                acceptor_full = match.group(2)  # e.g., 'I211O2'

                # Extract atom names using regex
                donor_match = atom_pattern.match(donor_full)
                acceptor_match = atom_pattern.match(acceptor_full)

                if donor_match and acceptor_match:
                    donor_atom = donor_match.group(1)  # e.g., 'O10'
                    acceptor_atom = acceptor_match.group(1)  # e.g., 'O2'

                    # Refine atom names
                    refined_donor = refine_atom_name(donor_atom)
                    refined_acceptor = refine_atom_name(acceptor_atom)

                    # Combine donor and acceptor with colon
                    pair = {"donor": refined_donor, "acceptor": refined_acceptor}
                    hbond_pairs.append(pair)
                else:
                    print(
                        f"Warning (Line {line_number}): Couldn't parse atoms in line: {line}"
                    )
            else:
                print(
                    f"Warning (Line {line_number}): Line didn't match expected format and was skipped: {line}"
                )

    # Create DataFrame with index
    df = pd.DataFrame(hbond_pairs)
    df.reset_index(inplace=True)
    df.rename(columns={"index": "idx"}, inplace=True)

    # Adjust idx to start from 0
    df["idx"] = df.index

    return df


def _determine_hbond_command() -> str:
    """Determine which gmx hbond variant supports the expected flags."""
    try:
        result = subprocess.run(
            ["gmx", "--version"],
            capture_output=True,
            check=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "gmx hbond"

    match = re.search(r"GROMACS version:\s*(\d+)\.(\d+)", result.stdout)
    if not match:
        return "gmx hbond"

    major, _minor = int(match.group(1)), int(match.group(2))
    if major > 2023:
        return "gmx hbond-legacy"
    return "gmx hbond"


from electrofit_analysis.cli.common import resolve_stage

def main(project_dir: str, viz: bool = False, stage: str = 'final') -> None:
    """Run hydrogen bond analysis for all molecules in the project's process directory.

    Args:
        project_dir: Path to the project root (must contain a 'process' subdirectory).
    """
    preferred_hbond_command = _determine_hbond_command()
    process_dir = os.path.join(project_dir, "process")

    if not os.path.isdir(process_dir):
        print(f"No 'process' directory found under: {project_dir}")
        return

    # Loop through each subdirectory in the process directory
    run_dir_name, analyze_base = resolve_stage(stage)
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)
        print(f"Running hydrogen bond analysis for: {folder_name}")

        # Check if it's a directory
        if os.path.isdir(folder_path):
            
            run_final_sim_dir = os.path.join(folder_path, run_dir_name)
            run_gau_create_gmx_in_dir = os.path.join(folder_path, "run_gau_create_gmx_in")

            # Check if 'run_final_gmx_simulation' exists
            if os.path.isdir(run_final_sim_dir):
                dest_dir = os.path.join(folder_path, analyze_base)
                dest_dir = os.path.join(dest_dir, "h_bonds")
                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)

                gro_file_path = os.path.join(run_final_sim_dir, "md.tpr")
                xtc_file_path = os.path.join(run_final_sim_dir, "md_center.xtc")

                print(f"reading gro file: {gro_file_path}")
                print(f"reading xtc file: {xtc_file_path}")
                hbond_candidates = [preferred_hbond_command]
                if preferred_hbond_command != "gmx hbond-legacy":
                    hbond_candidates.append("gmx hbond-legacy")

                selected_hbond_command = None
                last_error = None

                for candidate in hbond_candidates:
                    try:
                        run_command(
                            f'echo "2\n2\n" | {candidate} -s {gro_file_path} -f {xtc_file_path} -hbn intra_hb_idx.ndx -num intra_hb_num.xvg -dist intra_hb_dist.xvg -g intra_hb.log -hbm intra_hb_matrix.xpm',
                            cwd=dest_dir,
                        )
                    except subprocess.CalledProcessError as err:
                        last_error = err
                        if candidate != hbond_candidates[-1]:
                            print(
                                "gmx hbond failed with legacy flags; retrying with gmx hbond-legacy."
                            )
                        continue
                    selected_hbond_command = candidate
                    break

                if selected_hbond_command is None:
                    raise last_error  # Propagate the last failure if everything failed.

                if selected_hbond_command.endswith("hbond-legacy"):
                    print("Using gmx hbond-legacy for compatibility with legacy flags.")

                run_command(
                    f'echo "2\n5\n" | {selected_hbond_command} -s {gro_file_path} -f {xtc_file_path} -hbn inter_hb_idx.ndx -num inter_hb_num.xvg -dist inter_hb_dist.xvg -g inter_hb.log -hbm inter_hb_matrix.xpm',
                    cwd=dest_dir,
                )

                # ---- Plot inter (optional) ----
                if viz:
                    # Load the hb_num.xvg file
                    xvg_filename = "inter_hb_num.xvg"
                    data = load_hb_num_xvg(xvg_filename)

                    # Check if data was loaded correctly
                    if data.size == 0:
                        raise ValueError(f"No data found in {xvg_filename}.")

                    # Extract columns
                    time = data[:, 0]  # Time in ps
                    time = time / 10**3
                    num_hbonds = data[:, 1]  # Number of Hydrogen Bonds (s0)
                    pairs_within_0_35_nm = (
                        data[:, 2] if data.shape[1] > 2 else None
                    )  # Pairs within 0.35 nm (s1)

                    # Plotting
                    plt.figure(figsize=(8, 6))
                    plt.plot(
                        time, num_hbonds, label="Hydrogen Bonds", color="darkblue", linewidth=2
                    )
                    if pairs_within_0_35_nm is not None:
                        plt.plot(
                            time,
                            pairs_within_0_35_nm,
                            label="Pairs within 0.35 nm",
                            color="darkred",
                            linewidth=2,
                        )

                    plt.title("Hydrogen Bonds Over Time")
                    plt.xlabel("Time (ns)")
                    plt.ylabel("Number")
                    plt.legend()
                    plt.grid(False)

                    # Save the plot as PDF
                    output_svg = "inter_hb_num_over_time.pdf"
                    plt.savefig(output_svg, format="pdf")
                    plt.close()

                    print(f"Plot saved as {output_svg}")

                    # Plot and save hb donor acceptor distance distribution:
                    plot_hb_dist_xvg(
                        "inter_hb_dist.xvg", plot_file_name="inter_hb_distriution"
                    )

                # ------- work in progress

                inter_log_file = "inter_hb.log"

                # Parse the log file and get the DataFrame (always available for downstream use)
                inter_hbond_df = parse_hbond_log_to_dataframe(inter_log_file)

                print(inter_hbond_df)
                if viz:
                    visualize_data_donor_acceptor(
                        xpm_file="inter_hb_matrix.xpm",
                        hbond_df=inter_hbond_df,
                        output_prefix="inter",
                    )

                # --------------- Plot intra 2 ---------------
                log_file = "intra_hb.log"
                output_csv = "intra_hbond_pairs.csv"

                # Parse the log file and get the DataFrame
                hbond_df = parse_hbond_log_to_dataframe(log_file)

                # Display the DataFrame
                print("Extracted Donor-Acceptor Pairs DataFrame:")
                print(hbond_df)

                # Save the DataFrame to a CSV file
                hbond_df.to_csv(output_csv, index=False)
                print(f"\nDonor-Acceptor pairs have been saved to '{output_csv}'.")

                if viz:
                    visualize_data_donor_acceptor(
                        xpm_file="intra_hb_matrix.xpm",
                        hbond_df=hbond_df,
                        output_prefix="intra2",
                    )

                    data_matrix, meta_data = parse_xpm("intra_hb_matrix.xpm")

                    analysis_results = analyze_hydrogen_bonds(
                        data_matrix=data_matrix, metadata=meta_data
                    )

                    hbonds_per_index = analysis_results["hbonds_per_index"]
                    # Add 'count' column to hbond_df
                    hbond_df["count"] = hbonds_per_index

                    # For each donor, get the two most abundant hydrogen bonds
                    top_hbonds_df = (
                        hbond_df.groupby("donor")
                        .apply(lambda x: x.nlargest(2, "count"))
                        .reset_index(drop=True)
                    )

                    print("Top hydrogen bonds:")
                    print(top_hbonds_df)

                    # Path to your MOL2 file
                    os.chdir(run_gau_create_gmx_in_dir)
                    mol2_file_name = find_file_with_extension("mol2")
                    os.chdir(dest_dir)
                    mol2_file = os.path.join(run_gau_create_gmx_in_dir, mol2_file_name)

                    # Load the molecule
                    molecule = Chem.MolFromMol2File(mol2_file, removeHs=False)

                    if molecule is None:
                        raise ValueError(f"Failed to load molecule from {mol2_file}")

                    print("Molecule loaded successfully!")

                    # Set preference to use CoordGen
                    rdDepictor.SetPreferCoordGen(True)

                    # Generate 2D coordinates
                    rdDepictor.Compute2DCoords(molecule)

                    # Create custom atom labels
                    atom_counters = {}
                    atom_labels = {}

                    for atom in molecule.GetAtoms():
                        symbol = atom.GetSymbol()
                        idx = atom.GetIdx()

                        # Update the counter for this atom type
                        if symbol not in atom_counters:
                            atom_counters[symbol] = 1
                        else:
                            atom_counters[symbol] += 1

                        # Assign the custom label
                        atom_labels[idx] = f"{symbol}{atom_counters[symbol]}"

                    # Build a mapping from labels to atom indices
                    label_to_atom_idx = {label: idx for idx, label in atom_labels.items()}

                    # Now, get the donor and acceptor atom indices for highlighting
                    donor_atom_indices = []
                    acceptor_atom_indices = []

                    for _, row in hbond_df.iterrows():
                        donor_label = row["donor"]
                        acceptor_label = row["acceptor"]
                        donor_idx = label_to_atom_idx.get(donor_label)
                        acceptor_idx = label_to_atom_idx.get(acceptor_label)
                        if donor_idx is not None:
                            donor_atom_indices.append(donor_idx)
                        if acceptor_idx is not None:
                            acceptor_atom_indices.append(acceptor_idx)

                    # Remove duplicates
                    donor_atom_indices = list(set(donor_atom_indices))
                    acceptor_atom_indices = list(set(acceptor_atom_indices))

                    # Identify atoms that are both donors and acceptors
                    both_donor_acceptor_indices = set(donor_atom_indices) & set(
                        acceptor_atom_indices
                    )

                    # Update donor and acceptor lists to exclude atoms that are both
                    donor_only_indices = set(donor_atom_indices) - both_donor_acceptor_indices
                    acceptor_only_indices = (
                        set(acceptor_atom_indices) - both_donor_acceptor_indices
                    )

                    print("Donor only atom indices:", donor_only_indices)
                    print("Acceptor only atom indices:", acceptor_only_indices)
                    print(
                        "Both donor and acceptor atom indices:", both_donor_acceptor_indices
                    )

                    # Collect hbond_pairs from top_hbonds_df for drawing dashed lines
                    hbond_pairs = []

                    for _, row in top_hbonds_df.iterrows():
                        donor_label = row["donor"]
                        acceptor_label = row["acceptor"]
                        donor_idx = label_to_atom_idx.get(donor_label)
                        acceptor_idx = label_to_atom_idx.get(acceptor_label)
                        if donor_idx is not None and acceptor_idx is not None:
                            hbond_pairs.append((donor_idx, acceptor_idx))
                        else:
                            print(
                                f"Warning: Atom label {donor_label} or {acceptor_label} not found in molecule."
                            )

                    print("Hydrogen bond pairs for dashed lines:", hbond_pairs)

                    # Combine all atom indices for highlighting
                    highlight_atoms = donor_only_indices.union(acceptor_only_indices).union(
                        both_donor_acceptor_indices
                    )

                    # Define colors
                    # Donor only atoms: light blue
                    # Acceptor only atoms: light green
                    # Both donor and acceptor atoms: teal

                    highlightAtomColors = {}

                    # Donor only atoms (light blue)
                    for idx in donor_only_indices:
                        highlightAtomColors[idx] = (0.6, 0.8, 1)  # Light blue

                    # Acceptor only atoms (light green)
                    for idx in acceptor_only_indices:
                        highlightAtomColors[idx] = (0.6, 1, 0.6)  # Light green

                    # Both donor and acceptor atoms (teal)
                    for idx in both_donor_acceptor_indices:
                        highlightAtomColors[idx] = (0.6, 0.9, 0.8)  # Teal color

                    # Set SVG size
                    svg_size = 500

                    # Draw the molecule with highlighted atoms
                    drawer = rdMolDraw2D.MolDraw2DSVG(svg_size, svg_size)
                    opts = drawer.drawOptions()
                    opts.addAtomIndices = False
                    opts.addBondIndices = False
                    # Set font size
                    opts.baseFontSize = 0.3

                    # Assign custom atom labels individually
                    for idx, label in atom_labels.items():
                        opts.atomLabels[idx] = label

                    # Prepare the molecule for drawing
                    rdMolDraw2D.PrepareMolForDrawing(molecule)

                    # Draw the molecule
                    drawer.DrawMolecule(
                        molecule,
                        highlightAtoms=highlight_atoms,
                        highlightAtomColors=highlightAtomColors,
                    )
                    drawer.FinishDrawing()
                    svg = drawer.GetDrawingText()

                    # Extract drawing coordinates for each atom
                    draw_coords = {}
                    for atom_idx in range(molecule.GetNumAtoms()):
                        # Get the drawing coordinates (in pixels)
                        point = drawer.GetDrawCoords(atom_idx)
                        draw_coords[atom_idx] = (point.x, point.y)

                    # Add dashed lines for H-bonds
                    modified_svg = add_dashed_hbonds(svg, hbond_pairs, draw_coords)

                    # Save the modified SVG
                    with open("molecule_dashed_hbonds.svg", "w") as f:
                        f.write(modified_svg)

                    print(
                        "Modified SVG with dashed H-bonds saved as 'molecule_dashed_hbonds.svg'"
                    )

                    # Create legend SVG content
                    # Convert RGB tuples to hex codes
                    def rgb_to_hex(rgb_tuple):
                        return "#" + "".join(f"{int(255 * x):02X}" for x in rgb_tuple)

                    # Prepare the legend entries based on the data
                    legend_entries = []

                    # Include 'Donor' entry if there are donor-only atoms
                    if donor_only_indices:
                        legend_entries.append(
                            {"label": "Donor (D)", "color": rgb_to_hex((0.6, 0.8, 1))}
                        )  # Light blue

                    # Include 'Acceptor' entry if there are acceptor-only atoms
                    if acceptor_only_indices:
                        legend_entries.append(
                            {"label": "Acceptor (A)", "color": rgb_to_hex((0.6, 1, 0.6))}
                        )  # Light green

                    # Include 'D/A' entry if there are atoms that are both donors and acceptors
                    if both_donor_acceptor_indices:
                        legend_entries.append(
                            {"label": "D/A", "color": rgb_to_hex((0.6, 0.9, 0.8))}
                        )  # Teal

                    # Include 'H-Bond' entry if there are any hydrogen bonds (dashed lines)
                    if hbond_pairs:
                        legend_entries.append(
                            {"label": "H-Bond", "color": "lightgray", "dasharray": "5,5"}
                        )

                    legend_x = 20  # X position of the legend
                    legend_y = 20  # Y position of the legend
                    legend_svg_content = create_legend(
                        legend_x, legend_y, legend_entries, font_size=14
                    )

                    # Add legend to SVG content
                    modified_svg_with_legend = add_legend_to_svg(
                        modified_svg, legend_svg_content
                    )

                    # Save the modified SVG with legend
                    with open("molecule_dashed_hbonds_with_legend.svg", "w") as f:
                        f.write(modified_svg_with_legend)

                    print(
                        "Modified SVG with dashed H-bonds and legend saved as 'molecule_dashed_hbonds_with_legend.svg'"
                    )

                    # ---- Plot intra ----
                    visualize_data_donor_accpetor_pair(
                        xpm_file="intra_hb_matrix.xpm", hbond_df=hbond_df
                    )
                    # Load the hb_num.xvg file
                    xvg_filename = "intra_hb_num.xvg"
                    data = load_hb_num_xvg(xvg_filename)

                    # Check if data was loaded correctly
                    if data.size == 0:
                        raise ValueError(f"No data found in {xvg_filename}.")

                    # Extract columns
                    time = data[:, 0]  # Time in ps
                    num_hbonds = data[:, 1]  # Number of Hydrogen Bonds (s0)
                    pairs_within_0_35_nm = (
                        data[:, 2] if data.shape[1] > 2 else None
                    )  # Pairs within 0.35 nm (s1)

                    # Plotting
                    plt.figure(figsize=(8, 6))
                    plt.plot(
                        time, num_hbonds, label="Hydrogen Bonds", color="darkblue", linewidth=2
                    )
                    if pairs_within_0_35_nm is not None:
                        plt.plot(
                            time,
                            pairs_within_0_35_nm,
                            label="Pairs within 0.35 nm",
                            color="darkred",
                            linewidth=2,
                        )

                    plt.title("Hydrogen Bonds Over Time")
                    plt.xlabel("Time (ps)")
                    plt.ylabel("Number")
                    plt.legend()
                    plt.grid(False)

                    # Save the plot as PDF
                    output_svg = "intra_hb_num_over_time.pdf"
                    plt.savefig(output_svg, format="pdf")
                    plt.close()

                    print(f"Plot saved as {output_svg}")

                    # Plot and save hb donor acceptor distance distribution:
                    plot_hb_dist_xvg(
                        "intra_hb_dist.xvg", plot_file_name="intra_hb_distriution"
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Analyze hydrogen bonds for IP6 project. Provide the project root via --project."
        )
    )
    parser.add_argument(
        "--project",
        "-p",
        required=True,
        help="Path to the project root directory (contains 'process').",
    )
    parser.add_argument(
        "--viz",
        "-viz",
        action="store_true",
        help="Generate and save plots/figures (default: off).",
    )
    args = parser.parse_args()
    main(os.path.abspath(args.project), viz=args.viz)

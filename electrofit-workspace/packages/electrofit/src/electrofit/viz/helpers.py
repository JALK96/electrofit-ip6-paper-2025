"""Visualisierungs-Hilfsfunktionen (umgezogen aus ``electrofit.plotting.helpers``).

Diese Datei konsolidiert die zuvor im veralteten Namespace ``electrofit.plotting``
liegenden Funktionen. Importiere künftig direkt über ``electrofit.viz`` bzw.
``electrofit.viz.helpers``.
"""

import os
import re
import logging
import matplotlib.pyplot as plt  # Hard dependency: removed lazy/skip logic
import numpy as np

from .curly_brace import draw_curly_brace

logger = logging.getLogger(__name__)

__all__ = [
    "plot_charges_by_atom",
    "plot_charges_by_atom_sym",
    "create_atom_color_mapping",
    "plot_charges_by_symmetry",
]


ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)(\d*)")


## Removed _import_matplotlib: we now require matplotlib at runtime. Any ImportError should surface.


def _write_avg_charges(base_dir: str, filename: str, avg1, avg2=None) -> None:
    """Hilfsfunktion (Task 1) zum Schreiben der durchschnittlichen Ladungen."""
    path = os.path.join(base_dir, filename)
    with open(path, "w") as f:
        for v in avg1:
            f.write(f"{round(v,4)}\n")
        if avg2 is not None:
            f.write("\n")
            for v in avg2:
                f.write(f"{round(v,4)}\n")


def plot_charges_by_atom(atoms_dict1, initial_charges_dict, base_dir, atoms_dict2=None):
    """Plot charge distributions (Dataset 1 vs optional Dataset 2) – best effort.

    Fällt der Matplotlib Import, werden nur Durchschnittsladungen persistiert.
    """
    # Hard import already executed at module import; if unavailable this function will not run.
    atom_names = list(atoms_dict1.keys())
    charges_data1 = [atoms_dict1[atom]["charges"] for atom in atom_names]
    avg_charges1 = [atoms_dict1[atom]["average_charge"] for atom in atom_names]
    # Initial charges dict holds per-conformer charge lists; for the scatter we want a single scalar.
    # Verwende den Mittelwert der initialen Serie (Legacy-Parität: einzelne Werte pro Atom).
    def _mean(lst):
        return sum(lst) / len(lst) if lst else 0.0
    init_charges = [_mean(initial_charges_dict[atom]["charges"]) for atom in atom_names]
    avg_sum1 = sum(avg_charges1)

    if atoms_dict2 is not None:
        charges_data2 = [atoms_dict2[atom]["charges"] for atom in atom_names]
        avg_charges2 = [atoms_dict2[atom]["average_charge"] for atom in atom_names]
        avg_sum2 = sum(avg_charges2)
    else:
        charges_data2 = None
        avg_charges2 = None
        avg_sum2 = None

    color_map = {"H": "royalblue", "O": "darkred", "P": "darkorange", "C": "darkblue"}
    fig, ax = plt.subplots(figsize=(8, 6))
    positions = np.arange(len(atom_names))
    width = 0.35
    positions1 = positions - width / 2
    positions2 = positions + width / 2
    vp1 = ax.violinplot(
        charges_data1,
        positions=positions1,
        widths=width,
        showmeans=True,
        showmedians=False,
        showextrema=True,
    )
    if charges_data2 is not None:
        vp2 = ax.violinplot(
            charges_data2,
            positions=positions2,
            widths=width,
            showmeans=True,
            showmedians=False,
            showextrema=True,
        )

    # Early guard: skip plotting if any dataset contains only empty charge series
    if not atom_names or any(len(c) == 0 for c in charges_data1):
        logger.warning(
            "Skipping violin plot (plot_charges_by_atom): empty charge series encountered."
        )
        # Still write average charges file for downstream expectations
        charges_name = (
            "average_charges.chg"
            if atoms_dict2 is None
            else "average_charges_comparison.chg"
        )
        _write_avg_charges(base_dir, charges_name, avg_charges1, None)
        return

    def set_violin_colors(vp, alpha_value):
        for i, body in enumerate(vp["bodies"]):
            atom_name = atom_names[i]
            match = ELEMENT_RE.match(atom_name)
            element = match.group(1) if match else "Unknown"
            color = color_map.get(element, "black")
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha_value)
        for partname in ("cmeans", "cmins", "cmaxes", "cbars"):
            if partname in vp:
                items = vp[partname]
                if hasattr(items, "set_color") and not isinstance(items, list):
                    line_colors = []
                    for atom_name in atom_names:
                        match = ELEMENT_RE.match(atom_name)
                        element = match.group(1) if match else "Unknown"
                        line_colors.append(color_map.get(element, "black"))
                    items.set_color(line_colors)
                    items.set_linewidth(1.5)
                    items.set_alpha(alpha_value)

    set_violin_colors(vp1, 0.5)
    if charges_data2 is not None:
        set_violin_colors(vp2, 0.2)

    ax.scatter(positions1, init_charges, color="black", marker="o", label="Initial Charge", zorder=5, s=5)
    ax.set_xticks(positions)
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel("Atom Names")
    ax.set_ylabel("Atomic Partial Charge (e)")
    if avg_sum2 is not None:
        ax.set_title(
            f"Charge Distribution (Σ Avg 1: {round(avg_sum1,4)}, Σ Avg 2: {round(avg_sum2,4)})"
        )
    else:
        ax.set_title(f"Charge Distribution (Σ Avg: {round(avg_sum1,4)})")
    handles = [
        plt.Line2D(
            [],
            [],
            color="royalblue",
            marker="s",
            linestyle="None",
            markersize=10,
            label="Hydrogen",
        ),
        plt.Line2D(
            [],
            [],
            color="darkred",
            marker="s",
            linestyle="None",
            markersize=10,
            label="Oxygen",
        ),
        plt.Line2D(
            [],
            [],
            color="darkorange",
            marker="s",
            linestyle="None",
            markersize=10,
            label="Phosphorus",
        ),
        plt.Line2D(
            [],
            [],
            color="darkblue",
            marker="s",
            linestyle="None",
            markersize=10,
            label="Carbon",
        ),
        plt.Line2D(
            [],
            [],
            color="black",
            marker="o",
            linestyle="None",
            markersize=5,
            label="Initial Charge",
        ),
    ]
    if charges_data2 is not None:
        handles.extend(
            [
                plt.Line2D(
                    [],
                    [],
                    color="grey",
                    alpha=0.7,
                    marker="s",
                    linestyle="None",
                    markersize=10,
                    label="Dataset 1",
                ),
                plt.Line2D(
                    [],
                    [],
                    color="grey",
                    alpha=0.4,
                    marker="s",
                    linestyle="None",
                    markersize=10,
                    label="Dataset 2",
                ),
            ]
        )
    ax.legend(handles=handles, frameon=False, title="Legend")
    plt.tight_layout()
    figure_name = "charges.pdf" if atoms_dict2 is None else "charges_comparison.pdf"
    plt.savefig(os.path.join(base_dir, figure_name))
    plt.close(fig)
    charges_name = (
        "average_charges.chg"
        if atoms_dict2 is None
        else "average_charges_comparison.chg"
    )
    _write_avg_charges(base_dir, charges_name, avg_charges1, avg_charges2)

# Note: This plot only works for IP6 Configs
def plot_charges_by_atom_sym(
    atoms_dict1,
    initial_charges_dict,
    base_dir,
    atoms_dict2=None,
    equivalent_groups=None,
):
    """
    Plot the charge distributions of atoms along with their average and initial charges for IP6 microstates,
    coloring the violin plots based on atom groups. If a second charges dictionary is provided,
    compare the two datasets in the same plot.

    Parameters:
    - atoms_dict1: Dictionary of atoms and their charges (first dataset).
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    - atoms_dict2: (Optional) Dictionary of atoms and their charges (second dataset for comparison).
    - equivalent_groups: (Optional) Dictionary of symmetric atom groups.
    """

    atom_names = list(atoms_dict1.keys())
    charges_data1 = [atoms_dict1[atom]["charges"] for atom in atom_names]
    avg_charges1 = [atoms_dict1[atom]["average_charge"] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]["charges"] for atom in atom_names]
    # (avg_charges1 Summe nicht benötigt -> entfernt)

    # Build set of symmetric atoms
    symmetric_atoms = set()
    if equivalent_groups is not None:
        for representative, group in equivalent_groups.items():
            symmetric_atoms.add(representative)
            symmetric_atoms.update(group)
    else:
        symmetric_atoms = set()  # Empty set if no equivalent_groups provided

    # Build indices for symmetric and non-symmetric atoms
    symmetric_atom_indices = []
    non_symmetric_atom_indices = []

    for i, atom in enumerate(atom_names):
        if atom in symmetric_atoms:
            symmetric_atom_indices.append(i)
        else:
            non_symmetric_atom_indices.append(i)

    if atoms_dict2 is not None:
        # For symmetric atoms only
        atom_names_symmetric = [atom_names[i] for i in symmetric_atom_indices]
        charges_data2_symmetric = [
            atoms_dict2[atom]["charges"] for atom in atom_names_symmetric
        ]
        avg_charges2_symmetric = [
            atoms_dict2[atom]["average_charge"] for atom in atom_names_symmetric
        ]
    # (avg_charges2_symmetric Summe nicht benötigt -> entfernt)
    else:
        charges_data2_symmetric = None
        avg_charges2_symmetric = None

    # Define groups (same as your existing code)
    groups = {
        "C1-C6": [],
        "O1-O6": [],
        "O7-O24": [],
        "P1-P6": [],
        "H1-H6": [],
        # 'H7-Hx' group will be added dynamically
    }

    # Map each atom to its group
    atom_to_group = {}

    max_H_number = None  # Will store the maximum number of H atoms >=7
    H7_indices = []  # To store indices of H atoms with number >=7

    for i, atom_name in enumerate(atom_names):
        match = re.match(r"^([A-Z][a-z]?)(\d+)$", atom_name)
        if match:
            element = match.group(1)
            num = int(match.group(2))
            group_name = None
            if element == "C" and 1 <= num <= 6:
                groups["C1-C6"].append(i)
                group_name = "C1-C6"
            elif element == "O":
                if 1 <= num <= 6:
                    groups["O1-O6"].append(i)
                    group_name = "O1-O6"
                elif 7 <= num <= 24:
                    groups["O7-O24"].append(i)
                    group_name = "O7-O24"
            elif element == "P" and 1 <= num <= 6:
                groups["P1-P6"].append(i)
                group_name = "P1-P6"
            elif element == "H":
                if 1 <= num <= 6:
                    groups["H1-H6"].append(i)
                    group_name = "H1-H6"
                elif num >= 7:
                    H7_indices.append(i)
                    if max_H_number is None or num > max_H_number:
                        max_H_number = num
                    # group_name will be assigned later
            if group_name:
                atom_to_group[atom_name] = group_name
        else:
            # Handle non-matching atom names if necessary
            pass

    # After the loop, handle the H7-Hx group
    if H7_indices:
        if max_H_number == 7:
            H7_label = "H7"
        else:
            H7_label = f"H7-H{max_H_number}"
        # Add this group to groups
        groups[H7_label] = H7_indices
        # Assign group names to atom_to_group
        for idx in H7_indices:
            atom_name = atom_names[idx]
            atom_to_group[atom_name] = H7_label
    else:
        H7_label = None  # No H7 or higher hydrogens present

    # Define a color mapping for each group
    color_map = {
        "C1-C6": "black",
        "O1-O6": "lightcoral",
        "O7-O24": "red",
        "P1-P6": "darkorange",
        "H1-H6": "lightgrey",
    }

    # If H7 group exists, add it to the color map
    if H7_label:
        color_map[H7_label] = "dimgrey"

    # Hard import path (no graceful degradation)

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 6))

    positions = np.arange(len(atom_names))
    width = 0.35  # Width of each violin

    positions1 = positions - width / 2
    positions2 = positions + width / 2

    # Plot violin plots for Dataset 1 (all atoms)
    vp1 = ax.violinplot(
        charges_data1,
        positions=positions1,
        widths=width,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )

    # Plot violin plots for Dataset 2 (symmetric atoms only)
    if charges_data2_symmetric is not None:
        positions2_symmetric = positions2[symmetric_atom_indices]
        vp2 = ax.violinplot(
            charges_data2_symmetric,
            positions=positions2_symmetric,
            widths=width,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

    # Function to set colors based on atom group
    def set_violin_colors(vp, atom_names_list, alpha_value, red_mean=False):
        for i, body in enumerate(vp["bodies"]):
            atom_name = atom_names_list[i]
            group_name = atom_to_group.get(atom_name, "Unknown")
            color = color_map.get(
                group_name, "black"
            )  # Default to black if group not in color_map
            # Set the face and edge color of the violin body
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha_value)

    # Set colors for Dataset 1
    set_violin_colors(vp1, atom_names, alpha_value=0.4)

    # Set colors for Dataset 2
    if charges_data2_symmetric is not None:
        set_violin_colors(vp2, atom_names_symmetric, alpha_value=0.2)

    # Add scatter points for initial charges
    ax.scatter(
        positions1,
        init_charges,
        color="black",
        marker="o",
        label="Initial Charge",
        zorder=5,
        s=7,
    )

    # For non-symmetric atoms, plot their mean charges at positions1 with black rectangles
    non_symmetric_positions1 = positions1[non_symmetric_atom_indices]
    non_symmetric_avg_charges1 = [avg_charges1[i] for i in non_symmetric_atom_indices]
    ax.scatter(
        non_symmetric_positions1,
        non_symmetric_avg_charges1,
        color="black",
        marker="s",
        label="Mean (Non-Symmetric)",
        zorder=5,
        s=7,
    )

    # For symmetric atoms, plot their average charges from Dataset 2 at positions2_symmetric with red rectangles
    if charges_data2_symmetric is not None:
        ax.scatter(
            positions2_symmetric,
            avg_charges2_symmetric,
            color="red",
            marker="s",
            label="Combined Mean",
            zorder=5,
            s=7,
        )

    # Remove x-tick labels
    ax.set_xticks([])

    # Adjust y-limits to make space for braces
    y_min, y_max = ax.get_ylim()
    ax.set_ylim(bottom=y_min - 0.5)  # Adjust as needed

    # Define font dictionary for labels
    label_font = {
        "family": "serif",
        "color": "black",
        "weight": "normal",
        "style": "normal",
        "size": 12,
    }

    # Bracket coefficient (controls curvature)
    # (Your existing code for curly braces)

    # Draw curly braces and labels using curlyBrace
    for group_label, indices in groups.items():
        if not indices:
            continue  # Skip empty groups
        x_start = positions[min(indices)] - width
        x_end = positions[max(indices)] + width
        y = y_min - 0.1  # Position below the plot; adjust as needed

        # Define the points for the brace
        p1 = (x_start, y)
        p2 = (x_end, y)

        k_r = 0.09 / (x_end - x_start)

        # Choose color based on group (optional)
        brace_color = color_map.get(group_label, "black")

        # Draw the curly brace
        draw_curly_brace(
            fig,
            ax,
            p2,
            p1,
            k_r=k_r,
            auto=False,
            text=group_label,
            color=brace_color,
            lw=1.5,
            text_offset_lines=2,
            fontdict=label_font,
        )

    # Labeling
    ax.set_ylabel("Atomic Partial Charge (e)")

    # Create custom legend
    handles = []
    handles.append(
        plt.Line2D(
            [],
            [],
            color="black",
            marker="o",
            linestyle="None",
            markersize=5,
            label="Initial Charge",
        )
    )
    handles.append(
        plt.Line2D(
            [],
            [],
            color="black",
            marker="s",
            linestyle="None",
            markersize=5,
            label="Mean (Non-Sym)",
        )
    )
    if charges_data2_symmetric is not None:
        handles.append(
            plt.Line2D(
                [],
                [],
                color="red",
                marker="s",
                linestyle="None",
                markersize=5,
                label="Combined Mean",
            )
        )
        handles.extend(
            [
                plt.Line2D(
                    [],
                    [],
                    color="grey",
                    alpha=0.7,
                    marker="s",
                    linestyle="None",
                    markersize=10,
                    label="Sampled Charges",
                ),
                plt.Line2D(
                    [],
                    [],
                    color="grey",
                    alpha=0.4,
                    marker="s",
                    linestyle="None",
                    markersize=10,
                    label="Combined Charges",
                ),
            ]
        )
    else:
        handles.append(
            plt.Line2D(
                [],
                [],
                color="grey",
                alpha=0.7,
                marker="s",
                linestyle="None",
                markersize=10,
                label="Dataset 1",
            )
        )

    # Legende (ohne gruppenspezifische Farbfelder)
    ax.legend(handles=handles, title="Legend", frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_name = "charges_ip6.pdf" if atoms_dict2 is None else "charges_comparison_ip6.pdf"
    figure_path = os.path.join(base_dir, figure_name)
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_name = (
        "average_charges.chg"
        if atoms_dict2 is None
        else "average_charges_comparison.chg"
    )
    _write_avg_charges(base_dir, charges_name, avg_charges1, avg_charges2_symmetric)


def create_atom_color_mapping(atom_names, symmetry_groups):
    # List of colors to use for the groups
    group_colors_list = [
        "darkred",
        "darkgreen",
        "darkorange",
        "purple",
        "royalblue",
        "lightcoral",
        "deepskyblue",
        "mediumvioletred",
        "orange",
        "olive",
        "teal",
        "dodgerblue",
        "darkkhaki",
        "salmon",
        "firebrick",
        "olivedrab",
        "palevioletred",
    ]
    group_colors = {}

    # Map each group to a color
    for i, (group_representative, group_atoms) in enumerate(symmetry_groups.items()):
        color = group_colors_list[
            i % len(group_colors_list)
        ]  # Cycle through colors if needed
        # Include the representative atom in the group
        group = [group_representative] + group_atoms
        for atom in group:
            group_colors[atom] = color

    # Assign 'darkblue' to atoms not in any group
    atom_to_color = {}
    for atom in atom_names:
        color = group_colors.get(atom, "darkblue")
        atom_to_color[atom] = color

    return atom_to_color


def plot_charges_by_symmetry(
    atoms_dict1, initial_charges_dict, base_dir, symmetry_groups, atoms_dict2=None
):
    """
    Plot the charge distributions of atoms, coloring the violin plots based on symmetry groups.
    If a second charges dictionary is provided, compare the two datasets in the same plot.

    Parameters:
    - atoms_dict1: Dictionary of atoms and their charges (first dataset).
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    - symmetry_groups: Dictionary mapping representative atoms to lists of equivalent atoms.
    - atoms_dict2: (Optional) Dictionary of atoms and their charges (second dataset for comparison).
    """
    atom_names = list(atoms_dict1.keys())
    charges_data1 = [atoms_dict1[atom]["charges"] for atom in atom_names]
    avg_charges1 = [atoms_dict1[atom]["average_charge"] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]["charges"] for atom in atom_names]
    avg_sum1 = sum(avg_charges1)

    if atoms_dict2 is not None:
        charges_data2 = [atoms_dict2[atom]["charges"] for atom in atom_names]
        avg_charges2 = [atoms_dict2[atom]["average_charge"] for atom in atom_names]
        avg_sum2 = sum(avg_charges2)
    else:
        charges_data2 = None
        avg_charges2 = None
        avg_sum2 = None

    # Guard against empty data (avoid matplotlib reductions on empty arrays)
    if (not atom_names) or any(len(c) == 0 for c in charges_data1):
        logger.warning(
            "Skipping violin plot (plot_charges_by_symmetry): empty charge series encountered."
        )
        charges_name = (
            "average_charges.chg"
            if atoms_dict2 is None
            else "average_charges_comparison.chg"
        )
        _write_avg_charges(base_dir, charges_name, avg_charges1, avg_charges2)
        return

    # Create atom-to-color mapping
    atom_to_color = create_atom_color_mapping(atom_names, symmetry_groups)

    # Hard import path (no graceful degradation)

    fig, ax = plt.subplots(figsize=(8, 6))

    positions = np.arange(len(atom_names))
    width = 0.35  # Width of each violin

    positions1 = positions - width / 2
    positions2 = positions + width / 2

    # Plot violin plots for dataset 1
    vp1 = ax.violinplot(
        charges_data1,
        positions=positions1,
        widths=width,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )

    # Plot violin plots for dataset 2 if provided
    if charges_data2 is not None:
        vp2 = ax.violinplot(
            charges_data2,
            positions=positions2,
            widths=width,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

    # Function to set colors based on symmetry groups
    def set_violin_colors(vp, alpha_value):
        for i, body in enumerate(vp["bodies"]):
            atom_name = atom_names[i]
            color = atom_to_color[atom_name]
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha_value)

        # Customize the mean lines and other parts
        for partname in ("cmeans", "cmins", "cmaxes", "cbars"):
            if partname in vp:
                items = vp[partname]
                colors = [atom_to_color[atom_names[i]] for i in range(len(atom_names))]
                items.set_color(colors)
                items.set_linewidth(1.5)
                items.set_alpha(alpha_value)

    # Set colors for dataset 1
    set_violin_colors(vp1, alpha_value=0.5)

    # Set colors for dataset 2
    if charges_data2 is not None:
        set_violin_colors(vp2, alpha_value=0.2)

    # Add scatter points for initial charges
    ax.scatter(
        positions1,
        init_charges,
        color="black",
        marker="o",
        label="Initial Charge",
        zorder=5,
        alpha=1,
        s=5,
    )

    # Add scatter points for average charges
    ax.scatter(
        positions1,
        avg_charges1,
        color="black",
        marker="^",
        label="Average Charge",
        zorder=5,
        alpha=1,
        s=5,
    )

    if charges_data2 is not None:
        # Add scatter points for average charges
        ax.scatter(
            positions2,
            avg_charges2,
            color="black",
            marker="^",
            label="Initial Charge",
            zorder=5,
            alpha=0.5,
            s=5,
        )

    # Labeling
    ax.set_xticks(positions)
    ax.set_xticklabels(atom_names, rotation=60)
    ax.set_xlabel("Atom Names")
    ax.set_ylabel("Atomic Partial Charge (e)")
    if avg_sum2 is not None:
        ax.set_title(
            f"Charge Distribution for Each Atom\n(Average Sum Dataset 1: {round(avg_sum1, 2)}, Dataset 2: {round(avg_sum2, 2)})"
        )
    else:
        ax.set_title(
            f"Charge Distribution for Each Atom (Average Sum: {round(avg_sum1, 2)})"
        )

    # Create custom legend
    list(set(atom_to_color.values()))
    handles = []
    # for color in group_colors:
    #    handles.append(plt.Line2D([], [], color=color, marker='s', linestyle='None', markersize=10, label=color))

    handles.append(
        plt.Line2D(
            [],
            [],
            color="black",
            marker="o",
            linestyle="None",
            markersize=5,
            label="Initial Charge",
        )
    )
    handles.append(
        plt.Line2D(
            [],
            [],
            color="black",
            marker="^",
            linestyle="None",
            markersize=5,
            label="Average Charge",
        )
    )

    # Add dataset labels to legend if comparing
    if charges_data2 is not None:
        handles.extend(
            [
                plt.Line2D(
                    [],
                    [],
                    color="darkblue",
                    alpha=0.5,
                    marker="s",
                    linestyle="None",
                    markersize=10,
                    label="Dataset 1",
                ),
                plt.Line2D(
                    [],
                    [],
                    color="darkblue",
                    alpha=0.2,
                    marker="s",
                    linestyle="None",
                    markersize=10,
                    label="Dataset 2",
                ),
            ]
        )

    ax.legend(handles=handles, frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_name = (
        "charges_by_symmetry.pdf"
        if atoms_dict2 is None
        else "charges_by_symmetry_comparison.pdf"
    )
    figure_path = os.path.join(base_dir, figure_name)
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_name = (
        "average_charges.chg"
        if atoms_dict2 is None
        else "average_charges_comparison.chg"
    )
    _write_avg_charges(base_dir, charges_name, avg_charges1, avg_charges2)
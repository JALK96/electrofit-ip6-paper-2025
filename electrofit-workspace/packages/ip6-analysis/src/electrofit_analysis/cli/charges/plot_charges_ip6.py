import fnmatch
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import necessary modules from the project
from electrofit.config.loader import load_config, resolve_symmetry_flags
from electrofit.config.legacy import ConfigParser
from electrofit.io.files import (
    adjust_atom_names,
    extract_charges_from_subdirectories,
    find_file_with_extension,
    parse_charges_from_mol2,
)
from electrofit.io.symmetry import load_symmetry_groups
from electrofit.viz.helpers import (
    plot_charges_by_atom,
    plot_charges_by_atom_sym,
    plot_charges_by_symmetry,
)
from electrofit_analysis.cli.common import parse_molecule_args, normalize_micro_name

def calculate_symmetric_group_averages(charges_dict, equivalent_groups):
    """
    Calculate the mean of the average charges for symmetric atoms.

    Parameters:
    - charges_dict: Dictionary containing charges data for atoms.
    - equivalent_groups: Dictionary containing symmetric atom groups.

    Returns:
    - updated_charges_dict: Dictionary with updated average charges for symmetric groups.
    """
    updated_charges_dict = charges_dict.copy()  # Create a copy to modify

    # Iterate over each group in the equivalent groups
    for representative, group in equivalent_groups.items():
        # Add the representative atom to the group
        full_group = [representative] + group

        # Collect the charges for all atoms in the group
        group_charges = []
        for atom in full_group:
            if atom in charges_dict:
                group_charges.extend(charges_dict[atom]['charges'])

        # Calculate the mean of the charges
        if group_charges:
            group_average = sum(group_charges) / len(group_charges)

            # Update the average charge for all atoms in the group
            for atom in full_group:
                if atom in updated_charges_dict:
                    updated_charges_dict[atom]['average_charge'] = group_average

    return updated_charges_dict

def combine_and_calculate_symmetric_group_averages(charges_dict, equivalent_groups):
    """
    Calculate the mean of the average charges for symmetric atoms, update the average charges,
    and create a combined dataset with expanded charges lists for symmetric atoms.

    Parameters:
    - charges_dict (dict): Dictionary containing charges data for atoms.
      Expected format:
          {
              'C1': {'charges': [...], 'average_charge': ...},
              'C2': {'charges': [...], 'average_charge': ...},
              ...
          }

    - equivalent_groups (dict): Dictionary containing symmetric atom groups.
      Expected format:
          {
              'C1': ['C2', 'C3'],
              'O1': ['O2', 'O3'],
              ...
          }

    Returns:
    - updated_charges_dict (dict): Dictionary with updated average charges for symmetric groups.
    - combined_charges_dict (dict): Dictionary with updated average charges and expanded charges lists of combined datasets for symmetric atom entries.
    """
    import copy
    # Deep copy to avoid modifying the original dictionaries
    updated_charges_dict = copy.deepcopy(charges_dict)
    combined_charges_dict = copy.deepcopy(charges_dict)

    # Iterate over each group in equivalent_groups
    for representative, group in equivalent_groups.items():
        # Form the complete group including the representative
        full_group = [representative] + group

        # Collect charges from each atom in the group
        group_charges = []
        for atom in full_group:
            if atom in charges_dict:
                group_charges.extend(charges_dict[atom]['charges'])
            else:
                print(f"Warning: Atom '{atom}' not found in charges_dict.")

        # Calculate the mean charge for the group
        if group_charges:
            group_average = sum(group_charges) / len(group_charges)

            # Update the average_charge for each atom in updated_charges_dict
            for atom in full_group:
                if atom in updated_charges_dict:
                    updated_charges_dict[atom]['average_charge'] = group_average
                else:
                    print(f"Warning: Atom '{atom}' not found in updated_charges_dict.")

            # Update the average_charge and charges list for each atom in combined_charges_dict
            for atom in full_group:
                if atom in combined_charges_dict:
                    combined_charges_dict[atom]['average_charge'] = group_average
                    combined_charges_dict[atom]['charges'] = group_charges.copy()
                else:
                    print(f"Warning: Atom '{atom}' not found in combined_charges_dict.")
        else:
            print(f"Warning: No charges found for group '{representative}'.")

    # Atoms not in symmetric groups remain unchanged in both dictionaries
    return updated_charges_dict, combined_charges_dict

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

def plot_histograms(
    df_original,
    df_adjusted,
    title,
    filename,
    adjusted_average_charges=None,
    symmetric_atoms=None,
    atom_to_color=None,
    symmetric_groups=None,
    combine_original_data=False,
    remove_outlier=False,
):
    """
    Plot histograms of the DataFrame and save the plot.

    Parameters:
        df_original (pd.DataFrame): Original charges data.
        df_adjusted (pd.DataFrame): Adjusted charges data.
        title (str): Title of the plot.
        filename (str): Filename to save the plot.
        adjusted_average_charges (dict, optional): Dictionary of adjusted average charges.
        symmetric_atoms (set, optional): Set of atom names that are symmetric.
        atom_to_color (dict, optional): Mapping of atom names to colors.
        symmetric_groups (dict, optional): Dictionary of symmetry groups.
    """
    num_atoms = len(df_original.columns)
    num_cols = min(4, num_atoms)
    num_rows = (num_atoms + num_cols - 1) // num_cols  # Calculate rows needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(4 * num_cols, 3 * num_rows))
    axes = axes.flatten()

    for i, col in enumerate(df_original.columns):
        ax = axes[i]
        data_original = df_original[col].dropna()  # Original data
        data_adjusted = df_adjusted[col].dropna()  # Adjusted data

        # Determine if the atom is symmetric
        is_symmetric = symmetric_atoms is not None and col in symmetric_atoms
        # Get the color for the atom
        if atom_to_color is not None:
            color = atom_to_color.get(col, "darkblue")
        else:
            color = "darkblue"

        # Calculate common bins for both datasets
        if is_symmetric and symmetric_groups is not None:
            # For symmetric atoms, combine data from the symmetry group
            group_found = False
            for rep, group in symmetric_groups.items():
                group_atoms = [rep] + group
                if col in group_atoms:
                    # Combine charges for the group
                    combined_data_adjusted = pd.Series(dtype=float)
                    combined_data_original = pd.Series(dtype=float)

                    for atom in group_atoms:
                        combined_data_adjusted = pd.concat(
                            [combined_data_adjusted, df_adjusted[atom].dropna()],
                            ignore_index=True,
                        )
                        combined_data_original = pd.concat(
                            [combined_data_original, df_original[atom].dropna()],
                            ignore_index=True,
                        )

                    if not combine_original_data:
                        combined_data_original = data_original
                        # Use combined data to calculate bins
                        data_combined = pd.concat(
                            [combined_data_original, combined_data_adjusted],
                            ignore_index=True,
                        )
                        bins = np.histogram_bin_edges(data_combined, bins=20)
                        # Plot original histogram in background
                        ax.hist(
                            combined_data_original,
                            bins=bins,
                            color=color,
                            alpha=0.9,
                            edgecolor="black",
                            label="Original",
                        )
                        # Plot combined adjusted histogram
                        ax.hist(
                            combined_data_adjusted,
                            bins=bins,
                            color=color,
                            alpha=0.5,
                            edgecolor="black",
                            label="Combined",
                        )
                        group_found = True
                        break

                    else:
                        # Use combined data to calculate bins
                        data_combined = pd.concat(
                            [combined_data_original, combined_data_adjusted],
                            ignore_index=True,
                        )
                        bins = np.histogram_bin_edges(data_combined, bins=20)
                        # Plot original histogram in background
                        ax.hist(
                            combined_data_original,
                            bins=bins,
                            color=color,
                            alpha=0.5,
                            edgecolor="black",
                            label="Comb. Original",
                        )
                        # Plot combined adjusted histogram
                        ax.hist(
                            combined_data_adjusted,
                            bins=bins,
                            color=color,
                            alpha=0.9,
                            edgecolor="black",
                            label="Comb. Clipped",
                        )
                        group_found = True
                        break
            if not group_found:
                # Atom is symmetric but not found in groups (shouldn't happen)
                # Use data from the individual atom
                data_combined = pd.concat(
                    [data_original, data_adjusted], ignore_index=True
                )
                bins = np.histogram_bin_edges(data_combined, bins=20)
                # Plot original and adjusted histograms
                ax.hist(
                    data_original,
                    bins=bins,
                    color=color,
                    alpha=0.5,
                    edgecolor="black",
                    label="Original",
                )
                if remove_outlier:
                    ax.hist(
                        data_adjusted,
                        bins=bins,
                        color=color,
                        alpha=0.9,
                        edgecolor="black",
                        label="Clipped",
                    )
        else:
            # Non-symmetric atom
            # Use data from the individual atom
            data_combined = pd.concat([data_original, data_adjusted], ignore_index=True)
            bins = np.histogram_bin_edges(data_combined, bins=20)
            # Plot original and adjusted histograms
            ax.hist(
                data_original,
                bins=bins,
                color=color,
                alpha=0.5,
                edgecolor="black",
                label="Original",
            )
            if remove_outlier:
                ax.hist(
                    data_adjusted,
                    bins=bins,
                    color=color,
                    alpha=0.9,
                    edgecolor="black",
                    label="Clipped",
                )

        ax.set_title(col)
        ax.set_xlabel("Charge")
        ax.set_ylabel("Frequency")

        if adjusted_average_charges is not None:
            # For symmetric atoms, use the group average charge
            if is_symmetric:
                mean_value = adjusted_average_charges.get(col, None)
            else:
                # For non-symmetric atoms, use their individual mean
                mean_value = data_adjusted.mean()
        else:
            # Calculate the mean of the adjusted data
            mean_value = data_adjusted.mean()

        if mean_value is not None:
            # Set line color based on symmetry
            line_color = "red" if is_symmetric else "black"
            # Plot a vertical dashed line at the mean
            ax.axvline(
                mean_value,
                color=line_color,
                linestyle="dashed",
                linewidth=2,
                label=f"{mean_value:.2f}",
            )
            # Add a text label with the mean value
            # ax.text(
            #    0.95, 0.95,
            #    f'{mean_value:.2f}',
            #    color=line_color,
            #    fontsize=10,
            #    ha='right', va='top',
            #    transform=ax.transAxes
            # )

        # Add legend
        ax.legend()

    # Remove any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename)
    plt.close()

def main(project_root: str, remove_outlier: bool = False, only=None) -> None:
    """Run IP6 charge plots for all molecules under <project_root>/process.

    Parameters
    ----------
    project_root : str
        Electrofit project root directory.
    remove_outlier : bool, optional
        Apply 1.5*IQR outlier removal before computing averages and plots.
    """
    process_dir = os.path.join(project_root, "process")
    if not os.path.isdir(process_dir):
        print(f"process directory not found: {process_dir}")
        return

    # normalize filter set if provided
    only_norm = None
    if only:
        only_norm = {normalize_micro_name(x) for x in only}

    for sub_dir in sorted(os.listdir(process_dir)):
        sub_dir_path = os.path.join(process_dir, sub_dir)
        if not os.path.isdir(sub_dir_path):
            print(f"Skipping '{sub_dir_path}' as it is not a directory.")
            continue
        if only_norm and sub_dir not in only_norm:
            continue

        base_dir = os.path.join(sub_dir_path, "extracted_conforms")
        pis_dir = os.path.join(sub_dir_path, "run_gau_create_gmx_in")

        ac_files = glob.glob(os.path.join(pis_dir, "*.acpype"))
        if not ac_files:
            print(
                f"No '.acpype' directories/files found in '{pis_dir}'. Skipping '{sub_dir_path}'."
            )
            continue
        print(f"Acpype files found in {pis_dir}: {ac_files}")
        ac_dir = ac_files[0]

        results_dir = os.path.join(sub_dir_path, "results")
        os.makedirs(results_dir, exist_ok=True)
        plots_dir = os.path.join(results_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        # Prefer TOML in process/<mol>/results; fallback to legacy .ef in base_dir
        toml_path = os.path.join(results_dir, "electrofit.toml")
        if os.path.isfile(toml_path):
            cfg = load_config(project_root=project_root, context_dir=sub_dir_path)
            molecule_name = cfg.project.molecule_name or sub_dir
            print(f"Processing: {molecule_name}")
            charge = cfg.project.charge
            print(f"Charge set to: {charge}")
            atom_type = cfg.project.atom_type
            print(f"AtomType set to: {atom_type}")
            adjust_sym, ignore_sym = resolve_symmetry_flags(cfg, "ensemble")
            cfg.project.adjust_symmetry = adjust_sym  # type: ignore[attr-defined]
            cfg.project.ignore_symmetry = ignore_sym  # type: ignore[attr-defined]
            print("AdjustSymmetry set to:", adjust_sym)
            print("IgnoreSymmetry set to:", ignore_sym)
            calc_group_average = bool(cfg.project.calculate_group_average)
            print("CalculateGroupAverage set to:", calc_group_average)
        else:
            # Legacy .ef fallback
            cwd0 = os.getcwd()
            try:
                os.chdir(base_dir)
                ef_name = find_file_with_extension("ef")
            finally:
                os.chdir(cwd0)
            if not ef_name:
                print(f"No TOML at {toml_path} and no .ef in {base_dir}. Skipping.")
                continue
            config_file_path = os.path.join(base_dir, ef_name)
            config = ConfigParser(config_file_path)
            molecule_name = config.MoleculeName
            print(f"Processing: {molecule_name}")
            charge = config.Charge
            print(f"Charge set to: {charge}")
            atom_type = config.AtomType
            print(f"AtomType set to: {atom_type}")
            adjust_sym = config.AdjustSymmetry
            print("AdjustSymmetry set to:", adjust_sym)
            ignore_sym = config.IgnoreSymmetry
            print("IgnoreSymmetry set to:", ignore_sym)
            calc_group_average = config.CalculateGroupAverage
            print("CalculateGroupAverage set to:", calc_group_average)

        # Determine the mol2 file to use
        mol2_source_file_path = None
        ac_entries = os.listdir(ac_dir)
        if atom_type:
            mol2_file_pattern = f"*{atom_type}.mol2"
            for file_name in ac_entries:
                if fnmatch.fnmatch(file_name, mol2_file_pattern):
                    mol2_source_file_path = os.path.join(ac_dir, file_name)
                    break
        if mol2_source_file_path is None:
            for cand in ac_entries:
                if cand.endswith("_bcc_gaff2.mol2"):
                    mol2_source_file_path = os.path.join(ac_dir, cand)
                    break
        if mol2_source_file_path is None:
            expected = f"{sub_dir}.mol2"
            if expected in ac_entries:
                mol2_source_file_path = os.path.join(ac_dir, expected)
        if mol2_source_file_path is None:
            mol2s = [f for f in ac_entries if f.endswith(".mol2")]
            if mol2s:
                mol2_source_file_path = os.path.join(ac_dir, mol2s[0])
        if mol2_source_file_path is None:
            print(f"No suitable .mol2 file found in '{ac_dir}'. Skipping.")
            continue

        initial_charges_dict = adjust_atom_names(
            parse_charges_from_mol2(mol2_source_file_path)
        )

        # Extract charges across conformers
        atoms_dict = extract_charges_from_subdirectories(base_dir, results_dir)

        # Plotting charges
        equiv_group = None
        if adjust_sym:
            cwd0 = os.getcwd()
            try:
                os.chdir(pis_dir)
                eq_file = find_file_with_extension("json")
            finally:
                os.chdir(cwd0)
            if eq_file:
                equiv_group = load_symmetry_groups(os.path.join(pis_dir, eq_file))
                plot_charges_by_symmetry(
                    atoms_dict, initial_charges_dict, plots_dir, equiv_group
                )
                plot_charges_by_atom_sym(atoms_dict, initial_charges_dict, plots_dir)
            else:
                print(f"No symmetry JSON in {pis_dir}; plotting without symmetry.")
                plot_charges_by_atom(atoms_dict, initial_charges_dict, plots_dir)
                adjust_sym = False  # fall back
        else:
            plot_charges_by_atom(atoms_dict, initial_charges_dict, plots_dir)

        # Build DataFrame for histograms
        charges_data = {k: v["charges"] for k, v in atoms_dict.items()}
        df = pd.DataFrame(charges_data)

        # Color mapping
        if adjust_sym and equiv_group is not None:
            symmetric_groups = equiv_group
            symmetric_atoms = set()
            for group in symmetric_groups.values():
                symmetric_atoms.update(group)
            symmetric_atoms.update(symmetric_groups.keys())
            atom_to_color = create_atom_color_mapping(df.columns.tolist(), symmetric_groups)
        else:
            symmetric_groups = {}
            symmetric_atoms = set()
            atom_to_color = {atom: "darkblue" for atom in df.columns}

        # Base histograms
        plot_histograms(
            df_original=df,
            df_adjusted=df,
            title="Charge Distribution",
            filename=os.path.join(plots_dir, "hist.pdf"),
            atom_to_color=atom_to_color,
            symmetric_atoms=symmetric_atoms,
        )

        # Outlier handling and additional plots
        if remove_outlier:
            print("Remove Outliers ...")

            def get_outlier_mask(series: pd.Series) -> pd.Series:
                Q1 = series.quantile(0.25)
                Q3 = series.quantile(0.75)
                IQR = Q3 - Q1
                return (series < (Q1 - 1.5 * IQR)) | (series > (Q3 + 1.5 * IQR))

            outlier_mask = pd.Series(False, index=df.index)
            for column in df.columns:
                outlier_mask |= get_outlier_mask(df[column])

            df_no_outliers = df[~outlier_mask].reset_index(drop=True)
            adjusted_average_charges = df_no_outliers.mean().to_dict()

            if calc_group_average and adjust_sym and equiv_group is not None:
                cleaned_atoms_dict = {
                    atom: {
                        "charges": df_no_outliers[atom].tolist(),
                        "average_charge": adjusted_average_charges[atom],
                    }
                    for atom in df_no_outliers.columns
                }
                updated_charges_dict = calculate_symmetric_group_averages(
                    cleaned_atoms_dict, symmetric_groups
                )
                adjusted_average_charges = {
                    atom: data["average_charge"] for atom, data in updated_charges_dict.items()
                }
                plot_histograms(
                    df_original=df,
                    df_adjusted=df_no_outliers,
                    title="Charge Distribution with Group Average Charges",
                    filename=os.path.join(plots_dir, "hist_group_average_clipped_charges.pdf"),
                    adjusted_average_charges=adjusted_average_charges,
                    symmetric_atoms=symmetric_atoms,
                    atom_to_color=atom_to_color,
                    symmetric_groups=symmetric_groups,
                    combine_original_data=True,
                    remove_outlier=True,
                )
                plot_charges_by_symmetry(
                    updated_charges_dict, initial_charges_dict, plots_dir, equiv_group
                )
            elif adjust_sym and equiv_group is not None:
                plot_histograms(
                    df_original=df,
                    df_adjusted=df_no_outliers,
                    title="Clipped Charge Distribution with Symmetry and Clipped Average Charges",
                    filename=os.path.join(plots_dir, "hist_average_clipped_charges.pdf"),
                    adjusted_average_charges=adjusted_average_charges,
                    atom_to_color=atom_to_color,
                    symmetric_atoms=symmetric_atoms,
                    remove_outlier=True,
                )
            else:
                plot_histograms(
                    df_original=df,
                    df_adjusted=df_no_outliers,
                    title="Clipped Charge Distribution with no Symmetry and Clipped Average Charges",
                    filename=os.path.join(plots_dir, "hist_average_clipped_charges.pdf"),
                    adjusted_average_charges=adjusted_average_charges,
                    atom_to_color=atom_to_color,
                    remove_outlier=True,
                )
                cleaned_atoms_dict = {
                    atom: {
                        "charges": df_no_outliers[atom].tolist(),
                        "average_charge": adjusted_average_charges[atom],
                    }
                    for atom in df_no_outliers.columns
                }
                plot_charges_by_atom(cleaned_atoms_dict, initial_charges_dict, plots_dir)
        else:
            if calc_group_average and adjust_sym and equiv_group is not None:
                updated_charges_dict, combind_charges_dict = (
                    combine_and_calculate_symmetric_group_averages(
                        atoms_dict, symmetric_groups
                    )
                )
                group_average_charges = {
                    atom: data["average_charge"] for atom, data in updated_charges_dict.items()
                }
                plot_histograms(
                    df_original=df,
                    df_adjusted=df,
                    title="Charge Distribution with Group Average Charges",
                    filename=os.path.join(plots_dir, "hist_group_average_charges.pdf"),
                    adjusted_average_charges=group_average_charges,
                    symmetric_atoms=symmetric_atoms,
                    atom_to_color=atom_to_color,
                    symmetric_groups=symmetric_groups,
                    combine_original_data=False,
                )
                plot_charges_by_symmetry(
                    updated_charges_dict, initial_charges_dict, plots_dir, equiv_group
                )
                plot_charges_by_atom_sym(
                    atoms_dict,
                    initial_charges_dict,
                    plots_dir,
                    combind_charges_dict,
                    equivalent_groups=symmetric_groups,
                )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot charge distributions for IP6 project directories."
    )
    parser.add_argument(
        "--project",
        type=str,
        required=True,
        help="Path to the electrofit project root (contains 'process/').",
    )
    parser.add_argument(
        "--remove-outlier",
        action="store_true",
        help="Apply 1.5*IQR outlier removal before plotting.",
    )

    args = parser.parse_args()
    main(args.project, remove_outlier=args.remove_outlier)

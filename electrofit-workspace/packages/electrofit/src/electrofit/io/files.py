
import os
import re
import shutil

import glob
import logging
import numpy as np

from electrofit.io.mol2_ops import parse_charges_from_mol2

def copy_and_rename_folders(source: str, destination: str, nested_folder: str = "run_gau_create_gmx_in") -> None:
    """
    Copy folders from the source directory to the destination with an additional nested folder.
    Rename .mol2 files within each folder according to a binary-to-IP mapping.
    
    Parameters:
    - source (str): Path to the source directory.
    - destination (str): Path to the destination directory.
    - nested_folder (str): Name of the intermediate folder to nest files (default is "run_gau_create_gmx_in").
    """

    # Ensure the source directory exists
    if not os.path.isdir(source):
        print(f"Source directory '{source}' does not exist.")
        return

    # Create the destination directory if it does not exist
    os.makedirs(destination, exist_ok=True)

    # Walk through each folder in the source directory
    for folder_name in os.listdir(source):
        folder_path = os.path.join(source, folder_name)

        # Only process directories
        if os.path.isdir(folder_path):

            # Define the new nested path in the destination
            nested_dest_path = os.path.join(destination, folder_name, nested_folder)

            # Create the nested directory structure
            os.makedirs(nested_dest_path, exist_ok=True)

            # Copy all files from the source folder to the nested destination folder
            for item in os.listdir(folder_path):
                item_source_path = os.path.join(folder_path, item)
                item_dest_path = os.path.join(nested_dest_path, item)

                # Copy file or directory
                if os.path.isfile(item_source_path):
                    shutil.copy2(item_source_path, item_dest_path)
                elif os.path.isdir(item_source_path):
                    shutil.copytree(item_source_path, item_dest_path, dirs_exist_ok=True)

            print(f"Copied '{folder_path}' to '{nested_dest_path}'")

def replace_charge_in_ac_file(file_path, new_charge_float, cwd=None):
    """
    Replaces all CHARGE lines in the file with the new charge value.
    """
    if cwd:
        home = os.getcwd()
        os.chdir(cwd)

    new_charge_int = int(round(new_charge_float))
    new_charge_line = f"CHARGE     {new_charge_float:.2f} ( {new_charge_int} )\n"

    with open(file_path, "r") as file:
        lines = file.readlines()

    charge_replaced = False
    updated_lines = []
    for line in lines:
        if line.startswith("CHARGE"):
            updated_lines.append(new_charge_line)
            charge_replaced = True
            print(
                f"Replaced line:\nOld: {line.strip()}\nNew: {new_charge_line.strip()}"
            )
        else:
            updated_lines.append(line)

    if not charge_replaced:
        print("No CHARGE line found in the file.")
        return

    with open(file_path, "w") as file:
        file.writelines(updated_lines)

    if cwd:
        os.chdir(home)

    print(
        f"All CHARGE entries have been updated to {new_charge_float:.2f} ({new_charge_int})."
    )

def modify_gaussian_input(file_path):
    """
    Modifies the Gaussian input file to remove 'opt' and set charge/multiplicity.

    Parameters:
    - file_path (str): Path to the Gaussian input file.
    """
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        with open(file_path, "w") as file:
            for line in lines:
                # Remove 'opt' from the route section
                if " opt" in line:
                    line = line.replace(" opt", "")
                    logging.info("Removed 'opt' from the route section.")

                file.write(line)

    except Exception as e:
        logging.error(f"Error modifying Gaussian input file: {e}")
        raise

def adjust_atom_names(atoms_dict):
    """
    Adjusts atom names in the provided dictionary by appending a count to each unique element symbol.

    Parameters:
        atoms_dict (dict): Dictionary where keys are atom names and values are properties associated with each atom.

    Returns:
        dict: New dictionary with adjusted atom names.
    """

    # Get the list of atom names from the dictionary keys
    atom_names = list(atoms_dict.keys())

    # Initialize a dictionary to keep track of counts for each element
    counts = {}
    adjusted_atom_names = []

    for name in atom_names:
        # Extract the element symbol (handles one or two-letter symbols)
        match = re.match(r"^([A-Z][a-z]?)(\d*)", name)
        if match:
            element = match.group(1)
        else:
            # If the name doesn't match the pattern, keep it unchanged
            adjusted_atom_names.append(name)
            continue

        # Update the count for the element
        counts.setdefault(element, 0)
        counts[element] += 1

        # Adjust the name with the element symbol and its count
        adjusted_name = f"{element}{counts[element]}"
        adjusted_atom_names.append(adjusted_name)

    # Create a mapping from old names to adjusted names
    name_mapping = dict(zip(atom_names, adjusted_atom_names))

    # Update atoms_dict keys with adjusted names
    adjusted_atoms_dict = {
        new_name: atoms_dict[old_name] for old_name, new_name in name_mapping.items()
    }

    return adjusted_atoms_dict

def extract_charges_from_subdirectories(base_dir, results_dir):
    """
    Walk through subdirectories of the base_dir and extract charges from mol2 files.

    Parameters:
    - base_dir: Path to the directory containing subdirectories with mol2 files.

    Returns:
    - adjusted_atoms_dict: Dictionary of atoms with adjusted names and charges collected from all subdirectories.
    """
    atoms_dict = {}
    # Traverse each subdirectory in 'base_dir'
    subdirs = [
        f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))
    ]

    for subdir in subdirs:
        subdir_path = os.path.join(base_dir, subdir)
        entries = os.listdir(subdir_path)
        resp_files = [f for f in entries if f.endswith("_resp.mol2")]
        json_charge = os.path.join(subdir_path, "charges.json")
        if resp_files:
            for file_name in resp_files:
                mol2_path = os.path.join(subdir_path, file_name)
                atoms = parse_charges_from_mol2(mol2_path)
                if not atoms_dict:
                    for atom_name in atoms.keys():
                        atoms_dict[atom_name] = {"charges": []}
                for atom_name, atom_data in atoms.items():
                    atoms_dict[atom_name]["charges"].extend(atom_data["charges"])
        elif os.path.isfile(json_charge):  # test / synthetic path: simple JSON mapping atom->charge
            try:
                import json as _json
                with open(json_charge, "r") as fh:
                    simple = _json.load(fh)
                if not atoms_dict:
                    for atom_name in simple.keys():
                        atoms_dict[atom_name] = {"charges": []}
                for atom_name, charge_val in simple.items():
                    atoms_dict.setdefault(atom_name, {"charges": []})
                    atoms_dict[atom_name]["charges"].append(charge_val)
            except Exception as e:  # pragma: no cover
                print(f"[extract] warn: failed reading synthetic charges.json in {subdir_path}: {e}")

    # Calculate the mean charge for each atom
    for atom_name, atom_data in atoms_dict.items():
        charges = atom_data["charges"]
        atom_data["average_charge"] = np.mean(charges) if charges else 0

    adjusted_atoms_dict = adjust_atom_names(atoms_dict)

    # Write the average charges to the output file
    output_file = os.path.join(results_dir, "average_charges.txt")
    try:
        with open(output_file, "w") as f:
            f.write("#Atom_Name\tAverage_Charge\n")
            for atom_name, atom_data in adjusted_atoms_dict.items():
                f.write(f"{atom_name}\t{atom_data['average_charge']:.4f}\n")
        print(f"Average charges successfully written to {output_file}")
    except Exception as e:
        print(f"An error occurred while writing to {output_file}: {e}")

    return adjusted_atoms_dict

def find_file_with_extension(extension):
    cwd = os.getcwd()
    search_pattern = os.path.join(cwd, f"*.{extension}")
    files = glob.glob(search_pattern)

    if files:
        file_name = os.path.basename(files[0])  # Take the first found result
        print(f"Found file: {file_name}")
        return file_name
    else:
        print(f"No files with .{extension} extension found in {cwd}")
        return None

def strip_extension(file_name):
    # Split the file name into the name and the extension
    name, extension = os.path.splitext(file_name)

    print(f"File name: {name}")
    return name

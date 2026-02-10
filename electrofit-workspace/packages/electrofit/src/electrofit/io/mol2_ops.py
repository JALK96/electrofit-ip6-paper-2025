import os
import uuid  # for unique temporary filenames in mol2_to_pdb_with_bonds
from openbabel import openbabel
import logging

# Delay heavy/native import until actually needed (for crash isolation)
try:
    from openmol import tripos_mol2 as mol2  # type: ignore
except Exception as _imp_err:  # pragma: no cover - diagnostic path
    mol2 = None  # will trigger fallback/explicit error when used
    logging.warning("Deferred openmol import failed at module import: %s", _imp_err)

class Mol2ChargeError(Exception):
    """Domain-specific exception for MOL2 charge update failures."""


def read_charges(chg_file_path):
    """
    Reads charges from a .chg file and returns a list of floats.
    Assumes charges are separated by spaces and may span multiple lines.
    """
    charges = []
    try:
        with open(chg_file_path, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                try:
                    line_charges = [float(charge) for charge in line.split()]
                except ValueError as ve:
                    raise Mol2ChargeError(
                        f"Non-numeric value in charges file '{chg_file_path}': {ve}"
                    ) from ve
                charges.extend(line_charges)
    except FileNotFoundError as fnf:
        raise Mol2ChargeError(f"Charges file '{chg_file_path}' not found") from fnf
    return charges

def update_mol2_charges(mol2_input_path, chg_file_path, mol2_output_path):
    """
    Updates the charges in a MOL2 file using charges from a .chg file.

    Parameters:
    - mol2_input_path: Path to the input MOL2 file.
    - chg_file_path: Path to the .chg file containing new charges.
    - mol2_output_path: Path where the updated MOL2 file will be saved.
    """
    # Step 1: Check if input MOL2 file exists
    logging.info("[mol2-update] ENTER mol2_input=%s chg=%s out=%s", mol2_input_path, chg_file_path, mol2_output_path)
    if not os.path.isfile(mol2_input_path):
        raise Mol2ChargeError(f"MOL2 input file '{mol2_input_path}' does not exist")

    # Step 2: Read the MOL2 file
    logging.info("[mol2-update] Reading MOL2 file ...")
    if mol2 is None:
        raise Mol2ChargeError("openmol.tripos_mol2 not available (import failed earlier)")
    try:
        p = mol2.read(mol2_input_path)
    except Exception as e:
        raise Mol2ChargeError(f"Failed to read MOL2 file '{mol2_input_path}': {e}") from e

    num_atoms = len(p.atom_name)
    logging.info("[mol2-update] Atoms=%d", num_atoms)

    # Step 3: Read the charges from the .chg file
    logging.info("[mol2-update] Reading charges ...")
    charges = read_charges(chg_file_path)
    num_charges = len(charges)
    logging.info("[mol2-update] Charges read=%d", num_charges)

    # Step 4: Validate the number of charges matches the number of atoms
    if num_charges != num_atoms:
        raise Mol2ChargeError(
            f"Charge/atom count mismatch: charges={num_charges} atoms={num_atoms}"
        )
    else:
        logging.info("[mol2-update] Count match OK")

    # Step 5: Update the charges in the MOL2 structure
    logging.info("[mol2-update] Updating atom charges ...")
    for i in range(num_atoms):
        original_charge = p.atom_q[i]
        p.atom_q[i] = charges[i]
        atom_name = p.atom_name[i]
        logging.debug(
            "[mol2-update] atom=%d name=%s old=%s new=%s",
            i + 1,
            atom_name,
            original_charge,
            charges[i],
        )

    # Step 6: Rebuild the MOL2 structure
    logging.info("[mol2-update] Rebuilding MOL2 structure ...")
    try:
        p = mol2.build(p)
    except Exception as e:
        raise Mol2ChargeError(f"Failed to rebuild MOL2 structure: {e}") from e

    # Step 7: Write the updated MOL2 file
    logging.info("[mol2-update] Writing updated MOL2 file ...")
    try:
        mol2.Writer(p, mol2_output_path).write()
    except Exception as e:
        raise Mol2ChargeError(f"Failed to write MOL2 file '{mol2_output_path}': {e}") from e
    logging.info("[mol2-update] DONE")


def mol2_to_pdb_and_back(input_file, output_file, residue_name, cwd=None):
    if cwd:
        os.chdir(cwd)

    obConversion = openbabel.OBConversion()

    # Convert MOL2 to PDB
    obConversion.SetInAndOutFormats("mol2", "pdb")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)

    # Write to PDB format (in memory)
    pdb_string = obConversion.WriteString(mol)

    # Read PDB string back into OBMol
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol2 = openbabel.OBMol()
    obConversion.ReadString(mol2, pdb_string)

    # Set the residue name
    for residue in openbabel.OBResidueIter(mol2):
        residue.SetName(residue_name)

    # Write back to MOL2 format
    obConversion.WriteFile(mol2, output_file)


def pdb_to_mol2(input_file, output_file, residue_name, cwd=None):
    if cwd:
        os.chdir(cwd)

    obConversion = openbabel.OBConversion()

    # Convert PDB to MOL2
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)
    print("Converting ", input_file, " to ", output_file, "...")

    # Set the residue name
    for residue in openbabel.OBResidueIter(mol):
        residue.SetName(residue_name)
    print("Changed Residue Name to:", residue_name)

    # Write to MOL2 format
    obConversion.WriteFile(mol, output_file)


def mol2_to_pdb_with_bonds(input_file, existing_pdb_file, cwd=None, verbose: bool = False):
    """
    Convert a MOL2 file to a PDB format, extract the bond information,
    and insert it into an existing PDB file. Inserts bond info before ENDMDL if present;
    otherwise, inserts before END. Adds a REMARK line documenting the bond information addition.

    Args:
        input_file (str): Path to the input MOL2 file.
        existing_pdb_file (str): Path to the existing PDB file to append bond information.
        cwd (str, optional): Working directory to change to before running. Defaults to None.
    """
    if cwd:
        os.chdir(cwd)

    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"mol2_to_pdb_with_bonds: source mol2 missing: {input_file}")
    if not os.path.isfile(existing_pdb_file):
        raise FileNotFoundError(f"mol2_to_pdb_with_bonds: target PDB missing: {existing_pdb_file}")

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "pdb")
    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, input_file):
        raise RuntimeError(f"OpenBabel failed to read mol2 {input_file}")

    # Use a unique temp filename to avoid collisions in parallel execution
    temp_pdb_path = f"temp_with_bonds_{uuid.uuid4().hex}.pdb"
    if not obConversion.WriteFile(mol, temp_pdb_path):
        raise RuntimeError(f"OpenBabel failed to write temporary PDB from {input_file}")
    if not os.path.isfile(temp_pdb_path):
        raise FileNotFoundError(f"Expected temporary PDB not created: {temp_pdb_path}")

    bond_info = []
    try:
        with open(temp_pdb_path, "r") as pdb_temp:
            for line in pdb_temp:
                if line.startswith("CONECT"):
                    bond_info.append(line)
    finally:
        try:
            os.remove(temp_pdb_path)
        except OSError:
            pass

    # Create the REMARK line to document the added bond information
    remark_line = f"REMARK   1 EDITED BY electrofit: Added bond information from {os.path.basename(input_file)}\n"

    # Read and modify the existing PDB content
    modified_pdb_content = [remark_line]  # Add the REMARK line at the beginning
    endmdl_found = False
    with open(existing_pdb_file, "r") as existing_pdb:
        for line in existing_pdb:
            # Insert bond info right before ENDMDL or END, based on the first found
            if not endmdl_found and line.strip() == "ENDMDL":
                modified_pdb_content.extend(bond_info)
                endmdl_found = True
            elif line.strip() == "END" and not endmdl_found:
                modified_pdb_content.extend(bond_info)
                endmdl_found = True
            modified_pdb_content.append(line)  # Add current line

    # Write the modified content back to the PDB file
    with open(existing_pdb_file, "w") as existing_pdb:
        existing_pdb.writelines(modified_pdb_content)

    import logging
    msg1 = f"Bond information inserted before ENDMDL or END in {existing_pdb_file}"
    msg2 = f"REMARK added: {remark_line.strip()}"
    if verbose:
        logging.info(msg1)
        logging.info(msg2)


def parse_mol2(mol2_file):
    """
    Parses the MOL2 file to extract atoms and bonds.

    Parameters:
    - mol2_file: Path to the MOL2 file.

    Returns:
    - atoms: Dictionary mapping atom_id to atom properties.
    - bonds: List of bonds with origin_atom_id, target_atom_id, and bond_type.
    """
    atoms = {}
    bonds = []
    with open(mol2_file, "r") as f:
        lines = f.readlines()

    section = None
    for line in lines:
        line = line.strip()
        if line.startswith("@<TRIPOS>ATOM"):
            section = "ATOM"
            continue
        elif line.startswith("@<TRIPOS>BOND"):
            section = "BOND"
            continue
        elif line.startswith("@<TRIPOS>"):
            section = None
            continue

        if section == "ATOM":
            if not line:
                continue
            parts = line.split()
            if len(parts) < 9:
                continue  # Incomplete atom line
            atom_id = int(parts[0])
            atom_name = parts[1]
            atom_type = parts[5]
            atoms[atom_id] = {
                "atom_name": atom_name,
                "atom_type": atom_type,
                "connections": [],
            }
        elif section == "BOND":
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue  # Incomplete bond line
            origin = int(parts[1])
            target = int(parts[2])
            bond_type = parts[3]
            bonds.append((origin, target, bond_type))
    return atoms, bonds


def parse_charges_from_mol2(mol2_file):
    """
    Parses the MOL2 file to extract per-atom charges.

    Parameters:
    - mol2_file: Path to the MOL2 file.

    Returns:
    - atoms: Ordered dictionary mapping *unique* atom keys (in MOL2 atom order)
      to a list of charges.

    Notes
    -----
    MOL2 atom names are not guaranteed to be unique (e.g. many files name atoms
    simply "C", "O", "H"). Downstream Step6 aggregation expects a 1:1 mapping
    between atoms and charge values, so we normalise names into a stable,
    unique scheme based on element symbol + running index in file order
    (e.g. O1..On, C1..Cn, H1..Hn).
    """
    import re

    atoms: dict[str, dict] = {}
    with open(mol2_file, "r") as f:
        lines = f.readlines()

    section = None
    element_counts: dict[str, int] = {}
    for line in lines:
        line = line.strip()
        if line.startswith("@<TRIPOS>ATOM"):
            section = "ATOM"
            continue
        elif line.startswith("@<TRIPOS>"):
            section = None
            continue

        if section == "ATOM":
            if not line:
                continue
            parts = line.split()
            if len(parts) < 9:
                continue  # Incomplete atom line
            atom_name = parts[1]
            atom_charge = float(parts[8])
            # Derive element symbol from atom name (1–2 letters).
            # Fall back to the full atom_name if it doesn't match the pattern.
            elem_match = re.match(r"^([A-Z][a-z]?)", atom_name)
            element = elem_match.group(1) if elem_match else atom_name
            element_counts[element] = element_counts.get(element, 0) + 1
            unique_name = f"{element}{element_counts[element]}"
            atoms[unique_name] = {"charges": [atom_charge]}
        else:
            continue

    return atoms

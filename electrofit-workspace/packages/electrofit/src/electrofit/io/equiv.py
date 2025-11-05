import json
import sys
from typing import Dict, Iterable, List, Tuple

def build_connections(atoms, bonds):
    """
    Establishes connections between atoms based on bonds.

    Parameters:
    - atoms: Dictionary of atoms.
    - bonds: List of bonds.

    Returns:
    - None (modifies atoms in place).
    """
    for origin, target, bond_type in bonds:
        atoms[origin]["connections"].append({"atom_id": target, "bond_type": bond_type})
        atoms[target]["connections"].append({"atom_id": origin, "bond_type": bond_type})


def create_equivalence_groups(atoms):
    """
    Creates equivalence groups for oxygen atoms connected exclusively to the same phosphorus atom.

    Parameters:
    - atoms: Dictionary of atoms with connection information.

    Returns:
    - equiv_groups: Dictionary for equiv_groups.json.
    """
    equiv_groups = {}
    # Collect all oxygen atom_ids
    oxygen_atom_ids = [
        atom_id for atom_id, atom in atoms.items() if atom["atom_type"].startswith("O")
    ]
    if not oxygen_atom_ids:
        print("No oxygen atoms found in the MOL2 file.")
        return {}
    min_oxygen_id = min(oxygen_atom_ids)

    # Create label map: atom_id -> label (e.g., 19 -> "O7" if min_oxygen_id=13)
    label_map = {
        atom_id: f"O{atom_id - min_oxygen_id + 1}" for atom_id in oxygen_atom_ids
    }

    # Debug: Print label map
    print("\nLabel Map:")
    for atom_id in sorted(label_map.keys()):
        print(f"  Atom ID {atom_id}: {label_map[atom_id]}")
    print()

    for atom_id, atom in atoms.items():
        # Identify phosphorus atoms
        if not atom["atom_type"].startswith("P"):
            continue

        # Find connected oxygen atoms via single or double bonds
        connected_oxys = []
        for conn in atom["connections"]:
            conn_atom = atoms.get(conn["atom_id"])
            if conn_atom and conn_atom["atom_type"].startswith("O"):
                # Check if this oxygen is bonded only to this phosphorus atom
                if len(conn_atom["connections"]) == 1:
                    connected_oxys.append(
                        {"atom_id": conn["atom_id"], "bond_type": conn["bond_type"]}
                    )
                elif len(conn_atom["connections"]) == 2:
                    # If oxygen is bonded to phosphorus and one more atom, exclude it
                    bonded_atoms = conn_atom["connections"]
                    bonded_to_p = any(
                        other["atom_id"] == atom_id for other in bonded_atoms
                    )
                    bonded_to_others = any(
                        other["atom_id"] != atom_id for other in bonded_atoms
                    )
                    if bonded_to_p and not bonded_to_others:
                        connected_oxys.append(
                            {"atom_id": conn["atom_id"], "bond_type": conn["bond_type"]}
                        )
                else:
                    # More than two connections, exclude
                    continue

        if len(connected_oxys) < 2:
            continue  # Need at least two oxygens to form equivalence

        # Select central oxygen with the lowest atom_id
        central_oxy = min(connected_oxys, key=lambda oxy: oxy["atom_id"])

        # Generate unique labels based on atom_id
        central_label = label_map[central_oxy["atom_id"]]
        equiv_labels = [
            label_map[oxy["atom_id"]]
            for oxy in connected_oxys
            if oxy["atom_id"] != central_oxy["atom_id"]
        ]

        # Debugging Output
        print(f"Phosphorus Atom ID: {atom_id}")
        print(f"  Central Oxygen: {central_label}")
        print(f"  Equivalent Oxygens: {equiv_labels}\n")

        # Ensure no duplicates and valid labels
        if equiv_labels:
            equiv_groups[central_label] = equiv_labels

    return equiv_groups


def write_equivalence_groups(equiv_groups, output_file):
    """
    Writes the equivalence groups to a JSON file.

    Parameters:
    - equiv_groups: Dictionary of equivalence groups.
    - output_file: Path to the output JSON file.

    Returns:
    - None
    """
    with open(output_file, "w") as f:
        json.dump(equiv_groups, f, indent=4)
    print(f"Equivalence groups successfully written to '{output_file}'.")


def generate_atom_labels(atom_numbers):
    element_symbols = {
        1: "H",
        6: "C",
        7: "N",
        8: "O",
        15: "P",
        16: "S",
        17: "Cl",
        35: "Br",
        53: "I",
    }
    label_counters = {}
    atom_labels = []
    for z in atom_numbers:
        symbol = element_symbols.get(z, f"Z{z}")
        count = label_counters.get(symbol, 0) + 1
        label_counters[symbol] = count
        atom_labels.append(f"{symbol}{count}")
    return atom_labels


def load_equivalence_groups(equiv_groups_file):
    try:
        with open(equiv_groups_file, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading equivalence groups: {e}")
        sys.exit(1)


def extract_atom_labels_from_mol2(mol2_file: str) -> List[str]:
    """Return sequential element labels (e.g. C1, O2) derived from the MOL2 ATOM block."""
    labels: List[str] = []
    counts: Dict[str, int] = {}
    try:
        with open(mol2_file, "r") as fh:
            in_atoms = False
            for raw in fh:
                line = raw.strip()
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if not in_atoms:
                    continue
                if not line or line.startswith("@<TRIPOS>"):
                    break
                parts = line.split()
                if len(parts) < 6:
                    continue
                atom_type = parts[5]
                element = ""
                for ch in atom_type:
                    if ch.isalpha():
                        element += ch
                    else:
                        break
                if not element:
                    raw_name = parts[1]
                    element = "".join(ch for ch in raw_name if ch.isalpha()) or raw_name
                element = element.capitalize()
                counts[element] = counts.get(element, 0) + 1
                labels.append(f"{element}{counts[element]}")
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"MOL2 file '{mol2_file}' not found while extracting labels") from exc
    return labels


def apply_symmetry_to_charges(
    mol2_file: str,
    charges: Iterable[float],
    equivalence_groups: Dict[str, List[str]],
) -> Tuple[List[float], List[str]]:
    """Average charges in-place for each symmetry group.

    Returns the updated charge list together with warnings for any missing atoms.
    """
    labels = extract_atom_labels_from_mol2(mol2_file)
    label_to_index = {label: idx for idx, label in enumerate(labels)}
    charges_list = list(charges)
    warnings: List[str] = []
    for representative, group_atoms in equivalence_groups.items():
        full_group = [representative] + list(group_atoms)
        indices: List[int] = []
        missing: List[str] = []
        for atom_label in full_group:
            idx = label_to_index.get(atom_label)
            if idx is None:
                missing.append(atom_label)
            else:
                indices.append(idx)
        if missing:
            warnings.append(
                "Missing atoms for symmetry averaging: "
                + ", ".join(missing)
                + f" (reference: {representative})"
            )
            continue
        if not indices:
            continue
        mean_charge = sum(charges_list[idx] for idx in indices) / len(indices)
        for idx in indices:
            charges_list[idx] = mean_charge
    return charges_list, warnings

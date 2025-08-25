import sys
import re

from electrofit.io.equiv import load_equivalence_groups, generate_atom_labels


def edit_resp_input(input_file, equiv_groups_file, output_file, ignore_sym=False):
    try:
        with open(input_file, "r") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading RESP input file: {e}")
        sys.exit(1)

    charge_mult_pattern = re.compile(r"^\s*-?\d+\s+(\d+)\s*$")
    atom_data_start = None
    num_atoms = None

    for idx, line in enumerate(lines):
        if "Resp charges for organic molecule" in line and "1.0" in lines[idx - 1]:
            charge_mult_line = lines[idx + 1].strip()
            match = charge_mult_pattern.match(charge_mult_line)
            if match:
                num_atoms = int(match.group(1))
                atom_data_start = idx + 2
            else:
                print(
                    f"Error parsing charge and multiplicity line: '{charge_mult_line}'"
                )
                sys.exit(1)
            break
    else:
        print("Error locating atom data in RESP input file.")
        sys.exit(1)

    atom_numbers = []
    for i in range(num_atoms):
        line = lines[atom_data_start + i].strip()
        parts = line.split()
        try:
            atom_numbers.append(int(parts[0]))
        except (IndexError, ValueError):
            atom_numbers.append(0)

    atom_labels = generate_atom_labels(atom_numbers)
    label_to_index = {label: idx + 1 for idx, label in enumerate(atom_labels)}

    equivalence_groups_labels = load_equivalence_groups(equiv_groups_file)

    if ignore_sym:
        equivalence_groups_labels = {}

    equivalence_groups = {}
    for central_label, equiv_labels in equivalence_groups_labels.items():
        central_idx = label_to_index.get(central_label)
        if not central_idx:
            print(f"Warning: Central atom '{central_label}' not found.")
            continue
        equiv_indices = [
            label_to_index.get(lbl) for lbl in equiv_labels if label_to_index.get(lbl)
        ]
        if equiv_indices:
            equivalence_groups[central_idx] = equiv_indices

    for i in range(num_atoms):
        line_number = atom_data_start + i
        parts = lines[line_number].strip().split()
        try:
            z_int = int(parts[0])
            lines[line_number] = f"  {z_int:>3}    0\n"
        except (IndexError, ValueError):
            continue

    for central_idx, equiv_indices in equivalence_groups.items():
        line_number = atom_data_start + (central_idx - 1)
        try:
            z_int = int(lines[line_number].strip().split()[0])
            lines[line_number] = f"  {z_int:>3}    0\n"
        except (IndexError, ValueError):
            continue
        for equiv_idx in equiv_indices:
            line_number = atom_data_start + (equiv_idx - 1)
            try:
                z_int = int(lines[line_number].strip().split()[0])
                lines[line_number] = f"  {z_int:>3}  {central_idx:>3}\n"
            except (IndexError, ValueError):
                continue

    try:
        with open(output_file, "w") as f:
            f.writelines(lines)
        print(f"Modified RESP input file saved as '{output_file}'.")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)
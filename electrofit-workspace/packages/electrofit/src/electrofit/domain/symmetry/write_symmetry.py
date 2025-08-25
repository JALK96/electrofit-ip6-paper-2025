"""Write symmetry mapping from a RESP input file.

Extracted from legacy electrofit.core.symmetry (domain extraction phase).
"""
from __future__ import annotations
import os
from contextlib import contextmanager
from typing import Dict, List

@contextmanager
def _pushd(path: str | None):
    """Temporarily change the working directory if path is provided."""
    if not path:
        yield
        return
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def write_symmetry(
    input_filename: str,
    output_filename: str,
    cwd: str | None = None,
    atomic_number_mapping: Dict[int, str] | None = None,
) -> str:
    """Parse RESP input file and write compact symmetry table.

    Returns absolute path to written file.
    """
    if atomic_number_mapping is None:
        atomic_number_mapping = {1: "H", 6: "C", 8: "O", 15: "P"}

    with _pushd(cwd):
        if not os.path.isfile(input_filename):
            raise FileNotFoundError(f"Input file '{input_filename}' does not exist.")
        with open(input_filename, "r", encoding="utf-8") as f:
            lines: List[str] = f.readlines()

        start_index = None
        occurrences = 0
        for i, line in enumerate(lines):
            if line.strip() == "Resp charges for organic molecule":
                occurrences += 1
                if occurrences == 2:
                    start_index = i + 2
                    break
        if start_index is None:
            raise ValueError("Data section not found in input file.")

        atoms = []
        element_counts: Dict[str, int] = {}
        line_to_atom = {}
        for idx, line in enumerate(lines[start_index:], start=1):
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            try:
                atomic_number = int(parts[0])
                index = int(parts[1])
            except ValueError:
                continue
            element = atomic_number_mapping.get(atomic_number, f"Element{atomic_number}")
            element_counts[element] = element_counts.get(element, 0) + 1
            atom_label = f"{element}{element_counts[element]}"
            atom_info = {
                "label": atom_label,
                "atomic_number": atomic_number,
                "index": index,
                "line_number": idx,
            }
            atoms.append(atom_info)
            line_to_atom[idx] = atom_info

        output_lines = []
        for atom in atoms:
            label = atom["label"]
            index = atom["index"]
            if index == 0:
                output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label}"
            else:
                target_atom = line_to_atom.get(index)
                if target_atom:
                    target_label = target_atom["label"]
                    output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label} = {target_label}"
                else:
                    output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label} (reference atom not found)"
            output_lines.append(output_line)

        with open(output_filename, "w", encoding="utf-8") as f:
            f.write("\n".join(output_lines))
        return os.path.abspath(output_filename)

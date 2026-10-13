from __future__ import annotations

from typing import Dict, Literal, Tuple


CoordDistanceBasis = Literal["oxygen", "phosphorus"]

PHOS_LABELS: Tuple[str, ...] = ("P", "P1", "P2", "P3", "P4", "P5")
PHOS_LABELS_CORRECTED: Tuple[str, ...] = ("P1", "P2", "P3", "P4", "P5", "P6")

# Mapping: phosphate -> three terminal oxygen atoms (no ester/bridging oxygen).
PHOS_OXYGENS: Dict[str, Tuple[str, ...]] = {
    "P": ("O6", "O7", "O8"),
    "P1": ("O9", "O10", "O11"),
    "P2": ("O12", "O13", "O14"),
    "P3": ("O15", "O16", "O17"),
    "P4": ("O18", "O19", "O20"),
    "P5": ("O21", "O22", "O23"),
}

# Mapping: phosphate -> phosphorus center atom.
PHOSPHORUS_ATOMS: Dict[str, Tuple[str, ...]] = {
    "P": ("P",),
    "P1": ("P1",),
    "P2": ("P2",),
    "P3": ("P3",),
    "P4": ("P4",),
    "P5": ("P5",),
}


def normalize_coord_distance_basis(value: str | None) -> CoordDistanceBasis:
    basis = (value or "oxygen").strip().lower()
    if basis not in {"oxygen", "phosphorus"}:
        raise ValueError(
            f"Unsupported coord distance basis '{value}'. Expected 'oxygen' or 'phosphorus'."
        )
    return basis  # type: ignore[return-value]


def site_atom_names_by_basis(
    basis: CoordDistanceBasis,
) -> Dict[str, Tuple[str, ...]]:
    if basis == "oxygen":
        return PHOS_OXYGENS
    if basis == "phosphorus":
        return PHOSPHORUS_ATOMS
    raise ValueError(f"Unsupported coord distance basis '{basis}'.")


def expected_site_atom_count(basis: CoordDistanceBasis) -> int:
    if basis == "oxygen":
        return 3
    if basis == "phosphorus":
        return 1
    raise ValueError(f"Unsupported coord distance basis '{basis}'.")


def site_mode_tag(basis: CoordDistanceBasis) -> str:
    return normalize_coord_distance_basis(basis)


def rdf_curve_tag(basis: CoordDistanceBasis) -> str:
    return "periphO" if basis == "oxygen" else "P"


def select_phosphate_site_atomgroups(
    u,
    basis: CoordDistanceBasis,
    *,
    resname_pattern: str = "I*",
):
    """
    Return a dict {phosphate_label -> AtomGroup} for the chosen coordination basis.

    For `oxygen`, each group must contain exactly 3 atoms.
    For `phosphorus`, each group must contain exactly 1 atom.
    """
    names_map = site_atom_names_by_basis(basis)
    expected = expected_site_atom_count(basis)
    selected = {
        label: u.select_atoms(
            f"resname {resname_pattern} and name " + " ".join(atom_names)
        )
        for label, atom_names in names_map.items()
    }
    for label, ag in selected.items():
        if len(ag) != expected:
            raise ValueError(
                f"Expected {expected} {basis} site atom(s) for {label}, got {len(ag)}."
            )
    return selected

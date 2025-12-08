
"""Helpers for deriving ion/box settings from configuration."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal, TYPE_CHECKING

AVOGADRO = 6.022_140_76e23  # mol^-1
DEFAULT_ANGLES = (90.0, 90.0, 90.0)

if TYPE_CHECKING:  # pragma: no cover - typing aid only
    from electrofit.config.loader import SimulationBoxSection, SimulationIonsSection, IonTargetSection, IonTargetsSection


class IonConfigError(ValueError):
    """Raised when ion/box settings are inconsistent or unsupported."""


@dataclass
class BoxSetup:
    type: str
    edge_nm: float | None = None
    lengths_nm: tuple[float, float, float] | None = None
    angles_deg: tuple[float, float, float] = DEFAULT_ANGLES
    volume_nm3: float | None = None
    notes: list[str] = field(default_factory=list)


@dataclass
class IonSetup:
    mode: Literal["salt", "target"]
    cation: str
    anion: str
    enforce_neutrality: bool
    salt_concentration: float | None = None
    positive_count: int | None = None
    negative_count: int | None = None
    positive_concentration: float | None = None
    inferred_negative: bool = False
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    effective_molarity: float | None = None

    @property
    def net_charge(self) -> int | None:
        if self.positive_count is None or self.negative_count is None:
            return None
        return self.positive_count - self.negative_count


@dataclass
class DerivedSimulation:
    mode: Literal["salt", "target"]
    box: BoxSetup
    ions: IonSetup
    warnings: list[str] = field(default_factory=list)
    log_lines: list[str] = field(default_factory=list)


def _default_valence_map() -> dict[str, int]:
    """Static valence mapping for common ions.

    Sourced from your force field ions.itp and common MD ions.
    """
    return {
        "IB+": +1,
        "CA": +2,
        "CL": -1,
        "NA": +1,
        "MG": +2,
        "K": +1,
        "RB": +1,
        "CS": +1,
        "LI": +1,
        "ZN": +2,
    }


def _resolve_valences(cation: str, anion: str) -> tuple[int, int, list[str]]:
    """Return (z_cation, z_anion, notes) from the static map with sanity checks."""
    notes: list[str] = []
    zmap = _default_valence_map()
    cat_key = (cation or "").strip().upper()
    an_key = (anion or "").strip().upper()
    z_cat = zmap.get(cat_key, +1)
    z_an = zmap.get(an_key, -1)
    if cat_key not in zmap:
        notes.append(f"Cation valence for '{cation}' not in map -> using {z_cat:+d}.")
    if an_key not in zmap:
        notes.append(f"Anion valence for '{anion}' not in map -> using {z_an:+d}.")
    # Enforce conventional signs
    if z_cat < 0:
        z_cat = abs(z_cat)
        notes.append("Cation valence negative -> coerced to positive.")
    if z_an > 0:
        z_an = -abs(z_an)
        notes.append("Anion valence positive -> coerced to negative.")
    if z_cat == 0:
        z_cat = +1
        notes.append("Cation valence 0 invalid -> coerced to +1.")
    if z_an == 0:
        z_an = -1
        notes.append("Anion valence 0 invalid -> coerced to -1.")
    return z_cat, z_an, notes


def _has_target_payload(target: IonTargetSection | None) -> bool:
    return bool(
        target
        and (
            target.desired_count is not None
            or target.concentration is not None
        )
    )


def _compute_box_lengths(volume_nm3: float, box_type: str) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    if volume_nm3 <= 0:
        raise IonConfigError("Computed volume must be positive.")
    bt = (box_type or "cubic").lower()
    if bt in {"cubic", "cube"}:
        length = volume_nm3 ** (1.0 / 3.0)
        return (length, length, length), DEFAULT_ANGLES
    if bt == "dodecahedron":
        # Volume of rhombic dodecahedron: (16/9) * sqrt(3) * a^3
        a = ((9.0 * volume_nm3) / (16.0 * math.sqrt(3.0))) ** (1.0 / 3.0)
        return (a, a, a), DEFAULT_ANGLES
    raise IonConfigError(
        f"Box type '{box_type}' not supported for automatic sizing."
    )


def _targets_section(ions_section: SimulationIonsSection) -> IonTargetsSection | None:
    targets = getattr(ions_section, "targets", None)
    if targets is None:
        return None
    # If both sub-targets are effectively empty, treat as absent
    if not _has_target_payload(getattr(targets, "positive_ion", None)) and not _has_target_payload(getattr(targets, "negative_ion", None)):
        return None
    return targets


def derive_simulation_settings(
    box_section: SimulationBoxSection,
    ions_section: SimulationIonsSection,
    solute_charge: int | float | None = 0,
) -> DerivedSimulation:
    charge = int(solute_charge or 0)
    targets = _targets_section(ions_section)
    positive_target = getattr(targets, "positive_ion", None) if targets else None
    negative_target = getattr(targets, "negative_ion", None) if targets else None
    target_mode = _has_target_payload(positive_target) or _has_target_payload(negative_target)

    log_lines: list[str] = []
    warnings: list[str] = []

    cation = ions_section.cation
    anion = ions_section.anion
    # Determine ion valences (static map; no file I/O)
    z_cat, z_an, z_notes = _resolve_valences(cation, anion)

    if not target_mode:
        # Salt-concentration mode (legacy behaviour)
        box_notes: list[str] = []
        ion_notes: list[str] = []
        concentration = ions_section.salt_concentration
        if concentration is None:
            concentration = 0.15
        log_lines.append(
            f"[config][derived] Ion mode: salt (target concentration={concentration})"
        )
        log_lines.append(
            f"[config][derived] Box padding edge distance = {box_section.edge_nm} nm for type '{box_section.type}'"
        )
        ion_setup = IonSetup(
            mode="salt",
            cation=cation,
            anion=anion,
            enforce_neutrality=ions_section.enforce_neutrality,
            salt_concentration=concentration,
            notes=ion_notes + z_notes,
        )
        box_setup = BoxSetup(
            type=box_section.type,
            edge_nm=box_section.edge_nm,
            notes=box_notes,
        )
        return DerivedSimulation(
            mode="salt",
            box=box_setup,
            ions=ion_setup,
            warnings=warnings,
            log_lines=log_lines,
        )

    # Target-count mode
    if positive_target is None or positive_target.desired_count is None:
        raise IonConfigError(
            "Target ion specification requires 'positive_ion.desired_count'."
        )
    positive_count = positive_target.desired_count
    if positive_count <= 0:
        raise IonConfigError(
            "'positive_ion.desired_count' must be a positive integer."
        )
    if positive_target.concentration is None:
        raise IonConfigError(
            "Target ion specification requires 'positive_ion.concentration'."
        )
    if positive_target.concentration <= 0:
        raise IonConfigError(
            "'positive_ion.concentration' must be > 0 mol/L."
        )
    if negative_target and negative_target.concentration is not None:
        raise IonConfigError(
            "Negative-ion concentration is not supported yet; remove 'negative_ion.concentration'."
        )

    enforce = ions_section.enforce_neutrality
    inferred_negative = False
    negative_count = negative_target.desired_count if negative_target else None
    ion_warnings: list[str] = []
    ion_notes: list[str] = []
    ion_notes.extend(z_notes)

    if enforce:
        # charge + z_cat*n_pos + z_an*n_neg = 0  ->  n_neg = -(charge + z_cat*n_pos)/z_an
        numerator = -(charge + z_cat * positive_count)
        expected_f = numerator / float(z_an)
        expected_negative = int(round(expected_f))
        if abs(expected_f - expected_negative) > 1e-8:
            raise IonConfigError(
                f"Neutrality requires a non-integer number of {anion} ions ({expected_f:.3f}) with valence {z_an:+d}."
                " Adjust counts or ion species so neutrality is attainable."
            )
        if expected_negative < 0:
            raise IonConfigError(
                f"Neutrality would require {expected_negative} {anion} ions (valence {z_an:+d}) based on solute charge {charge}"
                f" and {positive_count} {cation} (valence {z_cat:+d}), which is impossible."
            )
        if negative_count is None:
            negative_count = expected_negative
            inferred_negative = True
            ion_notes.append(
                f"Negative ion count inferred as {negative_count} from charge {charge} with valences {cation}={z_cat:+d}, {anion}={z_an:+d}."
            )
        elif negative_count != expected_negative:
            raise IonConfigError(
                f"Neutrality requires {expected_negative} {anion} ions (valence {z_an:+d}) given {positive_count} {cation} (valence {z_cat:+d})"
                f" and solute charge {charge}, but {negative_count} were specified."
            )
    else:
        if negative_count is None:
            negative_count = 0
            ion_warnings.append(
                f"Neutrality disabled: inserting {positive_count} {cation} and no {anion}."
            )
        net_charge = charge + z_cat * positive_count + z_an * negative_count
        if net_charge != 0:
            ion_warnings.append(
                f"Resulting net system charge will be {net_charge:+d} e."  # sign display
            )

    if negative_count is None or negative_count < 0:
        raise IonConfigError("Derived negative ion count must be >= 0.")

    box_notes: list[str] = []
    if box_section.edge_nm is not None:
        warnings.append(
            f"Ignoring box.edge_nm={box_section.edge_nm} nm because target ion mode controls the volume."
        )
    if ions_section.salt_concentration is not None:
        warnings.append(
            "Ignoring 'salt_concentration' because explicit ion targets were provided."
        )

    # Compute volume from positive ion concentration
    mol_positive = positive_count / AVOGADRO
    volume_litre = mol_positive / positive_target.concentration
    volume_nm3 = volume_litre * 1e24
    lengths, angles = _compute_box_lengths(volume_nm3, box_section.type)

    log_lines.append(
        f"[config][derived] Ion mode: target (enforce_neutrality={enforce})"
    )
    log_lines.append(
        f"[config][derived] Target counts -> {cation}: {positive_count} (z={z_cat:+d}), {anion}: {negative_count} (z={z_an:+d})"
    )
    log_lines.append(
        f"[config][derived] Target concentration {positive_target.concentration} mol/L -> volume {volume_nm3:.3f} nm^3"
    )
    log_lines.append(
        "[config][derived] Box lengths (nm): "
        + " x ".join(f"{value:.3f}" for value in lengths)
        + f" for type '{box_section.type}'"
    )

    box_setup = BoxSetup(
        type=box_section.type,
        lengths_nm=lengths,
        angles_deg=angles,
        volume_nm3=volume_nm3,
        notes=box_notes,
    )
    ion_setup = IonSetup(
        mode="target",
        cation=cation,
        anion=anion,
        enforce_neutrality=enforce,
        positive_count=positive_count,
        negative_count=negative_count,
        positive_concentration=positive_target.concentration,
        inferred_negative=inferred_negative,
        warnings=ion_warnings,
        notes=ion_notes,
        effective_molarity=positive_target.concentration,
    )

    combined_warnings = warnings + ion_warnings
    return DerivedSimulation(
        mode="target",
        box=box_setup,
        ions=ion_setup,
        warnings=combined_warnings,
        log_lines=log_lines,
    )

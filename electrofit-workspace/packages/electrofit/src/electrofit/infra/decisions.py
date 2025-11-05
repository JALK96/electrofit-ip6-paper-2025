"""Decision model & builders for stepwise logging of protocol & symmetry semantics.

English (project guideline).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional
import logging
import os


@dataclass(slots=True)
class DecisionModel:
    stage: str  # 'initial' | 'conformer' | 'sampling' | 'aggregation'
    protocol: str | None
    charges_origin: str  # 'bcc' | 'resp_initial' | 'resp_ensemble' | 'pending' | 'resp_ensemble_symavg'
    symmetry_requested: bool
    symmetry_ignore_flag: bool
    symmetry_effective: str  # 'applied (antechamber defined)' | 'applied (user defined)' | 'none'
    ensemble_mode: bool
    notes: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    # Optional extra key/values (sampling/aggregation metadata) stored as flat list of tuples to keep model stable
    extra: List[tuple[str, str]] = field(default_factory=list)

    def kv_pairs(self) -> list[tuple[str, str]]:
        # Deterministic ordering
        base = [
            ("stage", self.stage),
            ("protocol", str(self.protocol)),
            ("charges_origin", self.charges_origin),
            ("symmetry.adjustment.requested", str(self.symmetry_requested)),
            ("symmetry.ignore_flag", str(self.symmetry_ignore_flag)),
            ("symmetry.effective", self.symmetry_effective),
            ("ensemble_mode", str(self.ensemble_mode)),
            ("notes", '[' + ','.join(self.notes) + ']'),
            ("warnings", '[' + ','.join(self.warnings) + ']'),
        ]
        # Append any extra metadata (already deterministic order if inserted consistently)
        return base + list(self.extra)

    def log(self, step: str) -> None:
       
        # Human readability: grouped pretty table (same env switch as cfg table)
        if os.environ.get('ELECTROFIT_LOG_TABLE', '1') not in {'0', 'false', 'False'}:
            try:
                # Prepare grouped sections
                base_section = [
                    ('stage', self.stage),
                    ('protocol', self.protocol),
                    ('charges_origin', self.charges_origin),
                    ('ensemble_mode', self.ensemble_mode),
                ]
                sym_section = [
                    ('symmetry.adjustment.requested', self.symmetry_requested),
                    ('symmetry.ignore_flag', self.symmetry_ignore_flag),
                    ('symmetry.effective', self.symmetry_effective),
                ]
                extra_section = list(self.extra) if self.extra else []
                notes_section = [(f'note[{i}]', n) for i, n in enumerate(self.notes)] if self.notes else []
                warnings_section = [(f'warning[{i}]', w) for i, w in enumerate(self.warnings)] if self.warnings else []

                def emit_header(title: str):
                    logging.info(f"[{step}][decisions] ─ {title} ─")

                # Compute max key width across all sections for alignment
                all_rows = base_section + sym_section + extra_section + notes_section + warnings_section
                if all_rows:
                    k_width = min(max(len(k) for k, _ in all_rows), 48)
                    logging.info(f"[{step}][decisions] ── decision summary ──")
                    def emit_rows(rows):
                        for k, v in rows:
                            key = (k[:45] + '...') if len(k) > 48 else k
                            logging.info(f"[{step}][decisions] {key.ljust(k_width)} : {v}")
                    emit_header('base')
                    emit_rows(base_section)
                    emit_header('symmetry')
                    emit_rows(sym_section)
                    if extra_section:
                        emit_header('extra')
                        emit_rows(extra_section)
                    if notes_section:
                        emit_header('notes')
                        emit_rows(notes_section)
                    if warnings_section:
                        emit_header('warnings')
                        emit_rows(warnings_section)
            except Exception:
                logging.debug(f"[{step}][decisions] pretty table logging failed", exc_info=True)
        for w in self.warnings:
            logging.warning(f"[{step}][warn] {w}")


def build_initial_decision(protocol: str | None, adjust_sym: bool, ignore_sym: bool) -> DecisionModel:
    notes: list[str] = []
    warnings: list[str] = []
    if protocol == 'bcc':
        charges_origin = 'AM1-BCC'
        if adjust_sym:
            if ignore_sym:
                symmetry_effective = 'none'
                warnings.append('symmetry modifications suppressed (ignore_symmetry=True)')
            else:
                symmetry_effective = 'applied (user defined)'
                notes.append('bcc pipeline: Mulliken → AM1-BCC with symmetry averaging')
        else:
            symmetry_effective = 'applied (antechamber defined)'
    elif protocol == 'opt':
        charges_origin = 'resp_initial'
        if adjust_sym:
            if ignore_sym:
                symmetry_effective = 'none'
                warnings.append('symmetry modifications suppressed (ignore_symmetry=True)')
            else:
                symmetry_effective = 'applied (user defined)'
        else:
            symmetry_effective = 'applied (antechamber defined)'
    else:
        charges_origin = 'AM1-BCC'  # Conservative fallback
        warnings.append(f'unknown protocol={protocol} -> treating as bcc path')
        symmetry_effective = 'applied (antechamber defined)'
    return DecisionModel(
        stage='initial',
        protocol=protocol,
        charges_origin=charges_origin,
        symmetry_requested=adjust_sym,
        symmetry_ignore_flag=ignore_sym,
        symmetry_effective=symmetry_effective,
        ensemble_mode=False,
        notes=notes,
        warnings=warnings,
    )


def build_conformer_decision(protocol: str | None, adjust_sym: bool, ignore_sym: bool) -> DecisionModel:
    notes: list[str] = []
    warnings: list[str] = []
    # Ensemble RESP always single-point; protocol only influences notes.
    charges_origin = 'resp_ensemble'
    if protocol == 'bcc':
        notes.append('overriding initial bcc charges with ensemble RESP')
    elif protocol == 'opt':
        notes.append('ensemble uses single-point RESP (no per-conformer optimisation)')
    else:
        warnings.append(f'unknown protocol={protocol} -> continuing with ensemble RESP')
    if adjust_sym:
        if ignore_sym:
            symmetry_effective = 'none'
            warnings.append('symmetry modifications suppressed (ignore_symmetry=True)')
        else:
            symmetry_effective = 'applied (user defined)'
    else:
        symmetry_effective = 'applied (antechamber defined)'
    return DecisionModel(
        stage='conformer',
        protocol=protocol,
        charges_origin=charges_origin,
        symmetry_requested=adjust_sym,
        symmetry_ignore_flag=ignore_sym,
        symmetry_effective=symmetry_effective,
        ensemble_mode=True,
        notes=notes,
        warnings=warnings,
    )


def build_sampling_decision(
    *,
    protocol: Optional[str],
    adjust_sym: bool,
    ignore_sym: bool,
    sampling_method: str,
    sample_count: int,
    seed: Optional[int],
    symmetry_json_present: bool,
) -> DecisionModel:
    """Decision for Step4 sampling (no charges yet: pending state).

    Symmetry not applied here; only planning / availability is logged.
    """
    notes: List[str] = [f"sampling.method={sampling_method}", f"sampling.count={sample_count}"]
    if seed is not None:
        notes.append(f"sampling.seed={seed}")
    warnings: List[str] = []
    if adjust_sym:
        if ignore_sym:
            notes.append('symmetry requested+ignored (will suppress modifications in later steps)')
        else:
            notes.append('symmetry planned for RESP stages')
        if not symmetry_json_present and protocol in {'bcc','opt'}:
            warnings.append('symmetry groups JSON missing at sampling time')
    # charges not yet aggregated or RESP ensemble run
    return DecisionModel(
        stage='sampling',
        protocol=protocol,
        charges_origin='pending',
        symmetry_requested=adjust_sym,
        symmetry_ignore_flag=ignore_sym,
        symmetry_effective='applied (antechamber defined)',
        ensemble_mode=False,
        notes=notes,
        warnings=warnings,
        extra=[],
    )


def build_aggregation_decision(
    *,
    protocol: Optional[str],
    adjust_sym: bool,
    ignore_sym: bool,
    calc_group_average: bool,
    group_average_applied: bool,
    symmetry_json_found: bool,
) -> DecisionModel:
    """Decision for Step6 aggregation of RESP ensemble charges.

    group_average_applied indicates whether symmetry group averaging actually produced an alternative charge set.
    """
    notes: List[str] = []
    warnings: List[str] = []
    charges_origin = 'resp_ensemble'
    if group_average_applied:
        charges_origin = 'resp_ensemble_symavg'
        notes.append('group_symmetry_average_applied')
    else:
        if calc_group_average and adjust_sym and not ignore_sym:
            warnings.append('group average requested but not applied (missing data or early exit)')
    if adjust_sym:
        if ignore_sym:
            symmetry_effective = 'none'
            warnings.append('symmetry modifications suppressed (ignore_symmetry=True)')
        else:
            symmetry_effective = 'applied (user defined)' if symmetry_json_found else 'applied (antechamber defined)'
            if symmetry_effective == 'applied (antechamber defined)':
                warnings.append('symmetry JSON missing during aggregation')
    else:
        symmetry_effective = 'applied (antechamber defined)'
    if protocol == 'opt':
        notes.append('aggregation after initial opt + ensemble RESP')
    elif protocol == 'bcc':
        notes.append('aggregation overrides initial bcc charges')
    elif protocol not in {None, 'opt', 'bcc'}:
        warnings.append(f'unknown protocol={protocol} during aggregation')
    # extras
    extras: List[tuple[str,str]] = [
        ('aggregation.group_average_requested', str(calc_group_average)),
        ('aggregation.group_average_applied', str(group_average_applied)),
        ('symmetry.json_found', str(symmetry_json_found)),
    ]
    return DecisionModel(
        stage='aggregation',
        protocol=protocol,
        charges_origin=charges_origin,
        symmetry_requested=adjust_sym,
        symmetry_ignore_flag=ignore_sym,
        symmetry_effective=symmetry_effective,
        ensemble_mode=True,
        notes=notes,
        warnings=warnings,
        extra=extras,
    )

__all__ = [
    'DecisionModel',
    'build_initial_decision',
    'build_conformer_decision',
    'build_sampling_decision',
    'build_aggregation_decision',
]

"""Initial structure preparation & charge derivation orchestration.

First extraction pass from legacy `core.process_initial_structure`.
Goals:
- Provide clearer functional boundaries (Gaussian opt vs RESP vs topology generation).
- Reduce side effects in preparation for further adapter abstraction.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import logging
import os
import subprocess

from electrofit.cli.run_commands import (
    create_gaussian_input,
    gaussian_out_to_prepi,
    run_acpype as _run_acpype_real,
    run_espgen,
    run_gaussian_calculation,
    run_resp,
    run_python,
    run_command,
)
from electrofit.io.files import find_file_with_extension
from electrofit.io.mol2_ops import mol2_to_pdb_and_back, update_mol2_charges
from electrofit.domain.symmetry.write_symmetry import write_symmetry
from electrofit.domain.symmetry.resp_constraints import apply_and_optionally_modify, RespSymmetryConfig
from electrofit.io.resp_edit import edit_resp_input
from electrofit.infra.decisions import build_initial_decision
from electrofit.io.equiv import apply_symmetry_to_charges


@dataclass(slots=True)
class InitialPrepConfig:
    """Configuration for the initial preparation step.

    NOTE: This is *input* configuration, not the result. See :class:`InitialResult`.
    """
    molecule_name: str
    mol2_file: str
    net_charge: int
    residue_name: str
    atom_type: str = "gaff2"
    adjust_sym: bool = False
    ignore_sym: bool = False
    symmetry_initial_mode: str | None = None
    protocol: str = "bcc"  # or 'opt'


@dataclass(slots=True)
class InitialResult:
    """Structured outcome of the initial preparation.

    This avoids forcing downstream code (or tests) to re-discover artifact paths
    via globbing. For backwards compatibility existing callers can ignore the
    return value; a future refactor will thread this object through the pipeline.
    """
    protocol: str
    residue_name: str
    charges_source: str  # 'bcc' or 'resp'
    mol2_final: str  # path inside scratch_dir
    acpype_dir: str | None
    gmx_gro: str | None
    gmx_itp: str | None
    gmx_top: str | None
    posre_itps: list[str]


def normalize_residue_roundtrip(mol2_file: str, residue_name: str):
    """Normalize residue tags in a MOL2 file to a target residue_name.

    Performs an OpenBabel mol2->pdb->mol2 roundtrip (unless disabled) and verifies
    that all ATOM records have the desired residue tag. Logs diagnostic info.
    """
    try:
        # Quick scan for atom count
        atom_count = 0
        with open(mol2_file, "r") as fh:
            in_atoms = False
            for line in fh:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if in_atoms:
                    if line.startswith("@<TRIPOS>") or not line.strip():
                        break
                    atom_count += 1
                    if atom_count > 50:  # early stop
                        break
        if atom_count <= 1:
            logging.debug("[resnorm] skip (atom_count<=1) %s", mol2_file)
            return
        if os.environ.get("ELECTROFIT_DISABLE_OBABEL"):
            logging.info("[resnorm] disabled by env ELECTROFIT_DISABLE_OBABEL for %s", mol2_file)
            return
        def collect_resnames() -> set[str]:
            names: set[str] = set()
            with open(mol2_file, "r") as fh2:
                in_atoms2 = False
                for ln in fh2:
                    if ln.startswith("@<TRIPOS>ATOM"):
                        in_atoms2 = True
                        continue
                    if in_atoms2:
                        if ln.startswith("@<TRIPOS>") or not ln.strip():
                            break
                        parts = ln.split()
                        if len(parts) >= 8:
                            names.add(parts[7])
            return names
        before = collect_resnames()
        mol2_to_pdb_and_back(mol2_file, mol2_file, residue_name)
        after = collect_resnames()
        if after == {residue_name}:
            # Exact match
            if before != after:
                logging.info("[resnorm] normalized residues %s -> %s (%s)", before, after, os.path.basename(mol2_file))
            else:
                logging.debug("[resnorm] residues already normalized (%s)", residue_name)
        else:
            # Accept variants like IP61, IP62 ... if all share prefix
            if after and all(a.startswith(residue_name) for a in after):
                if before != after:
                    logging.info("[resnorm] prefix-normalized residues %s -> %s (prefix=%s)", before, after, residue_name)
                else:
                    logging.debug("[resnorm] residues already prefix-aligned (%s)", residue_name)
            else:
                logging.warning("[resnorm] mismatch after normalization: wanted prefix '%s' got %s in %s", residue_name, after, mol2_file)
    except Exception as e:  # pragma: no cover - defensive
        logging.warning("[resnorm] OpenBabel roundtrip suppressed after error: %s", e)


# RESP symmetry handling unified via domain.symmetry.resp_constraints.apply_and_optionally_modify


def _run_antechamber_mul(mol2_file: str, ac_output: str, net_charge: int, scratch_dir: str) -> None:
    """Generate antechamber AC file with Mulliken charges."""
    command = [
        "antechamber",
        "-i",
        mol2_file,
        "-fi",
        "mol2",
        "-fo",
        "ac",
        "-o",
        ac_output,
        "-c",
        "mul",
        "-nc",
        str(net_charge),
    ]
    run_command(command, cwd=scratch_dir)


def _run_am1bcc(ac_input: str, ac_output: str, scratch_dir: str) -> None:
    """Apply AM1-BCC corrections while staying in AC format."""
    command = [
        "am1bcc",
        "-i",
        ac_input,
        "-o",
        ac_output,
        "-j",
        "5",
    ]
    run_command(command, cwd=scratch_dir)


def _extract_charges_from_ac(ac_path: str) -> list[float]:
    """Parse charges from an antechamber AC file."""
    charges: list[float] = []
    with open(ac_path, "r") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                charges.append(float(parts[8]))
            except ValueError:
                continue
    return charges


def _write_charges_file(path: str, charges: list[float]) -> None:
    """Write charges as one value per line with fixed precision."""
    with open(path, "w") as fh:
        for value in charges:
            fh.write(f"{value:.6f}\n")


def process_initial(cfg: InitialPrepConfig, scratch_dir: str, original_dir: str, input_files: list[str], run_acpype=_run_acpype_real) -> InitialResult:
    # Ensure residue name normalization (legacy behaviour) before any downstream tool consumption.
    # Previously this happened in the legacy shim unconditionally; during refactor it was only kept there.
    # Re-introduce here so new callers using domain API directly still obtain consistent residue naming.
    # Ensure both scratch copy and (optionally) original are normalized so downstream
    # ACPYPE/GROMACS topology carries the intended residue tag. Previously normalization
    # happened pre-scratch in the legacy shim; after refactor it ran only on the original
    # file, leaving the scratch copy stale. We now normalize the scratch file explicitly
    # (and original as best-effort) to restore prior guarantees.
    try:
        scratch_mol2 = os.path.join(scratch_dir, cfg.mol2_file)
        if os.path.isfile(scratch_mol2):
            normalize_residue_roundtrip(scratch_mol2, cfg.residue_name)
        else:
            logging.debug("[prep] scratch mol2 missing for normalization: %s", scratch_mol2)
    except Exception as e:  # pragma: no cover - defensive logging only
        logging.warning("[prep] Residue renaming roundtrip (scratch only) skipped after error: %s", e)
    # Decision logging (before branching)
    try:
        build_initial_decision(
            cfg.protocol,
            cfg.adjust_sym,
            cfg.ignore_sym,
            symmetry_mode=cfg.symmetry_initial_mode,
        ).log('step1')
    except Exception:  # pragma: no cover
        logging.debug('[step1][decisions] logging failed', exc_info=True)
    if cfg.protocol == "opt":
        # Gaussian optimisation + RESP pipeline
        create_gaussian_input(cfg.mol2_file, cfg.molecule_name, cfg.net_charge, scratch_dir, cfg.atom_type)
        gaussian_input = f"{cfg.molecule_name}.gcrt"
        run_gaussian_calculation(gaussian_input, cfg.molecule_name, scratch_dir)
        esp_file = f"{cfg.molecule_name}.esp"
        run_espgen(f"{cfg.molecule_name}.gesp", esp_file, scratch_dir)
        g_out = f"{gaussian_input}.log"
        gaussian_out_to_prepi(g_out, scratch_dir)
        # Validate generated RESP inputs
        respin1 = os.path.join(scratch_dir, "ANTECHAMBER_RESP1.IN")
        respin2 = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
        for f in (respin1, respin2):
            if not os.path.isfile(f):
                raise FileNotFoundError(f"Missing RESP input file: {f}")
        run_python(write_symmetry, respin1, "symmetry_resp.txt", cwd=scratch_dir)
        if cfg.adjust_sym:  # unified helper (shared with conformer path)
            apply_and_optionally_modify(
                run_python,
                find_file_with_extension,
                edit_resp_input,
                write_symmetry,
                RespSymmetryConfig(scratch_dir=scratch_dir, adjust_sym=cfg.adjust_sym, ignore_sym=cfg.ignore_sym),
            )
        # choose possibly modified file
        respin1_mod = os.path.join(scratch_dir, "ANTECHAMBER_RESP1_MOD.IN") if cfg.adjust_sym else respin1
        resp_output1 = f"{cfg.molecule_name}-resp1.out"
        resp_pch1 = f"{cfg.molecule_name}-resp1.pch"
        resp_chg1 = f"{cfg.molecule_name}-resp1.chg"
        run_resp(respin1_mod, esp_file, resp_output1, resp_pch1, resp_chg1, scratch_dir)
        resp_output2 = f"{cfg.molecule_name}-resp2.out"
        resp_pch2 = f"{cfg.molecule_name}-resp2.pch"
        resp_chg2 = f"{cfg.molecule_name}-resp2.chg"
        run_resp(respin2, esp_file, resp_output2, resp_pch2, resp_chg2, scratch_dir, resp_chg1)
        mol2_resp = f"{cfg.molecule_name}_resp.mol2"
        run_python(update_mol2_charges, cfg.mol2_file, resp_chg2, mol2_resp, cwd=scratch_dir)
        run_acpype(mol2_resp, cfg.net_charge, scratch_dir, cfg.atom_type, charges="user")
        acpype_map = _expose_acpype_outputs(mol2_resp, scratch_dir)
        mol2_final = os.path.join(scratch_dir, mol2_resp)
        charges_source = "resp"
    else:  # bcc
        charges_mode = "bcc"
        mol2_filename = cfg.mol2_file
        if cfg.adjust_sym:
            base_name = os.path.splitext(os.path.basename(mol2_filename))[0]
            ac_file = f"{base_name}.ac"
            bcc_ac_file = f"{base_name}.bcc.ac"
            bcc_chg_file = f"{base_name}.bcc.chg"
            eq_chg_file = f"{base_name}.bcc.eq.chg"
            logging.info("[prep][bcc] running Mulliken -> AM1-BCC pipeline with symmetry options (ignore_sym=%s)", cfg.ignore_sym)
            _run_antechamber_mul(mol2_filename, ac_file, cfg.net_charge, scratch_dir)
            _run_am1bcc(ac_file, bcc_ac_file, scratch_dir)
            ac_abs = os.path.join(scratch_dir, bcc_ac_file)
            charges = _extract_charges_from_ac(ac_abs)
            if not charges:
                raise ValueError(f"No charges extracted from {bcc_ac_file}")
            bcc_chg_abs = os.path.join(scratch_dir, bcc_chg_file)
            _write_charges_file(bcc_chg_abs, charges)
            charges_file = bcc_chg_file
            if not cfg.ignore_sym:
                equiv_path = os.path.join(scratch_dir, "equiv_groups.json")
                if os.path.isfile(equiv_path):
                    try:
                        with open(equiv_path, "r") as fh:
                            equivalence_groups = json.load(fh)
                    except Exception as err:
                        logging.warning("[prep][bcc] failed to load symmetry groups %s: %s", equiv_path, err)
                    else:
                        if equivalence_groups:
                            mol2_abs = os.path.join(scratch_dir, mol2_filename)
                            averaged, warnings = apply_symmetry_to_charges(mol2_abs, charges, equivalence_groups)
                            for warning_msg in warnings:
                                logging.warning("[prep][bcc][sym] %s", warning_msg)
                            eq_chg_abs = os.path.join(scratch_dir, eq_chg_file)
                            _write_charges_file(eq_chg_abs, averaged)
                            charges_file = eq_chg_file
                            charges = averaged
                        else:
                            logging.info("[prep][bcc] symmetry groups file is empty; skipping averaging")
                else:
                    logging.warning("[prep][bcc] equiv_groups.json not found in scratch; skipping symmetry averaging")
            run_python(
                update_mol2_charges,
                mol2_filename,
                charges_file,
                mol2_filename,
                cwd=scratch_dir,
            )
            charges_mode = "user"
        run_acpype(mol2_filename, cfg.net_charge, scratch_dir, cfg.atom_type, charges=charges_mode)
        acpype_map = _expose_acpype_outputs(mol2_filename, scratch_dir)
        mol2_final = os.path.join(scratch_dir, os.path.basename(mol2_filename))
        charges_source = "bcc"
    logging.info("Initial structure processing complete (%s)", cfg.molecule_name)
    return InitialResult(
        protocol=cfg.protocol,
        residue_name=cfg.residue_name,
        charges_source=charges_source,
        mol2_final=mol2_final,
        acpype_dir=acpype_map.get('acpype_dir'),
        gmx_gro=acpype_map.get('gro'),
        gmx_itp=acpype_map.get('itp'),
        gmx_top=acpype_map.get('top'),
        posre_itps=acpype_map.get('posre', []),
    )


def _expose_acpype_outputs(mol2_path: str, scratch_dir: str) -> dict[str, object]:
    """Copy a subset of ACPYPE outputs into the scratch root for convenience.

    Returns a mapping with discovered primary GMX artifact paths for structured
    consumption by callers. Missing artifacts simply produce None entries.
    """
    base_name = os.path.splitext(os.path.basename(mol2_path))[0]
    acpype_dir = os.path.join(scratch_dir, f"{base_name}.acpype")
    result: dict[str, object] = {
        'acpype_dir': acpype_dir if os.path.isdir(acpype_dir) else None,
        'gro': None,
        'itp': None,
        'top': None,
        'posre': [],
    }
    if not os.path.isdir(acpype_dir):
        return result
    for fname in os.listdir(acpype_dir):
        copy_needed = False
        key = None
        if fname.endswith("_GMX.gro"):
            key = 'gro'
            copy_needed = True
        elif fname.endswith("_GMX.itp"):
            key = 'itp'
            copy_needed = True
        elif fname.endswith("_GMX.top"):
            key = 'top'
            copy_needed = True
        elif fname.startswith("posre_") and fname.endswith(".itp"):
            key = 'posre'
            copy_needed = True
        if not copy_needed:
            continue
        src = os.path.join(acpype_dir, fname)
        dst = os.path.join(scratch_dir, fname)
        if not os.path.exists(dst):
            try:  # pragma: no cover - defensive copy
                subprocess.run(['cp', src, dst], check=True)
            except Exception:  # pragma: no cover
                logging.debug("[acpype][expose] copy failed for %s", src, exc_info=True)
                continue
        if key == 'posre':
            result['posre'].append(dst)
        else:
            result[key] = dst
    return result

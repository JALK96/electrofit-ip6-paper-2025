"""Conformer charge computation orchestration.

Structured stages (input prep -> Gaussian -> ESP -> RESP prep -> symmetry doc -> RESP fit -> finalize)
with fine‑grained helpers to aid testing and profiling.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
import os
import time

from electrofit.cli.run_commands import (
    run_espgen,
    run_python,
)

from electrofit.io.mol2_ops import pdb_to_mol2
from electrofit.domain.symmetry.write_symmetry import write_symmetry
from electrofit.adapters import gaussian as gaussian_adapter
from electrofit.adapters import resp as resp_adapter
from electrofit.adapters.results import ConformerChargeResult, GaussianResult, RespStageResult
from electrofit.infra.step_logging import log_relevant_config
from electrofit.infra.decisions import build_conformer_decision


@dataclass(slots=True)
class ConformerConfig:
    """Minimal conformer configuration passed into the charge pipeline."""
    molecule_name: str
    pdb_file: str
    net_charge: int
    residue_name: str
    adjust_sym: bool = False
    ignore_sym: bool = False
    symmetry_ensemble_mode: str | None = None
    protocol: str = "bcc"  # or 'opt'


@dataclass(slots=True)
class PreparedInputs:
    """Products of the input preparation stage."""
    mol2_file: str
    gaussian_input: str


@dataclass(slots=True)
class RespPrepResult:
    """RESP preparation outputs (effective inputs + symmetry flag)."""
    resp1: str
    resp2: str
    symmetry_applied: bool


def _time_stage(name: str):  # pragma: no cover - trivial decorator
    def deco(fn):
        def wrapped(*a, **kw):
            t0 = time.perf_counter()
            try:
                return fn(*a, **kw)
            finally:
                dt = time.perf_counter() - t0
                logging.debug(f"[step5][stage] {name} {dt:.3f}s")
        return wrapped
    return deco

__all__ = [
    'ConformerConfig',
    'PreparedInputs',
    'RespPrepResult',
    'prepare_inputs',
    'generate_esp',
    'prepare_resp',
    'write_symmetry_doc',
    'process_conformer',
]


@_time_stage('prepare_inputs')
def prepare_inputs(cfg: ConformerConfig, scratch_dir: str) -> PreparedInputs:
    """Stage 1: Convert PDB to MOL2 and build Gaussian input template.

    Returns a PreparedInputs dataclass; building of final Gaussian run arguments is delegated to
    the Gaussian adapter (which may rebuild inputs if cache hydrated).
    """
    mol2_file = f"{cfg.molecule_name}.mol2"
    pdb_to_mol2(cfg.pdb_file, mol2_file, cfg.residue_name, cwd=scratch_dir)
    gaussian_input = gaussian_adapter.build_input(mol2_file, cfg.molecule_name, cfg.net_charge, scratch_dir)
    return PreparedInputs(mol2_file=mol2_file, gaussian_input=gaussian_input)



@_time_stage('generate_esp')
def generate_esp(cfg: ConformerConfig, gaussian_result: GaussianResult, scratch_dir: str) -> str:
    """Stage 3: Generate ESP file from Gaussian-produced .gesp.

    Preconditions: gaussian_result.gesp_file exists.
    Returns path to generated .esp file.
    """
    gesp_file = gaussian_result.gesp_file
    if not gesp_file or not os.path.isfile(gesp_file):
        raise FileNotFoundError(f"Gaussian stage did not produce .gesp file for {cfg.molecule_name}")
    esp_file = f"{cfg.molecule_name}.esp"
    run_espgen(gesp_file, esp_file, scratch_dir)
    return esp_file


@_time_stage('prepare_resp')
def prepare_resp(cfg: ConformerConfig, scratch_dir: str) -> RespPrepResult:
    """Stage 4: Prepare RESP inputs (BCC protocol) and apply symmetry if requested.

    For non-BCC protocols, we still resolve expected RESP input filenames to keep downstream
    processing uniform (they may belong to previous preparation steps). Symmetry is only applied
    in BCC protocol when adjust_sym is True.
    """
    symmetry_applied = False
    if cfg.protocol == "bcc":
        resp_adapter.prepare_ac(f"{cfg.molecule_name}.mol2", cfg.molecule_name, cfg.net_charge, scratch_dir)
        resp_adapter.generate_inputs(cfg.molecule_name, scratch_dir)
        if cfg.adjust_sym:
            resp_adapter.apply_symmetry(scratch_dir, cfg.adjust_sym, cfg.ignore_sym)
            symmetry_applied = True
    # Resolve effective files (Stage 5 combined logic)
    resp1 = os.path.join(
        scratch_dir,
        "ANTECHAMBER_RESP1_MOD.IN" if (cfg.adjust_sym and cfg.protocol == "bcc") else "ANTECHAMBER_RESP1.IN",
    )
    resp2 = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
    for f in (resp1, resp2):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing RESP input file: {f}")
    return RespPrepResult(resp1=resp1, resp2=resp2, symmetry_applied=symmetry_applied)


@_time_stage('write_symmetry_doc')
def write_symmetry_doc(resp1: str, adjust: bool, scratch_dir: str) -> str | None:
    """Stage 5: Write symmetry documentation file (purely informational).

    Returns the path to the symmetry documentation file (or None if unexpected failure).
    """
    out_name = "symmetry_resp_MOD.txt" if adjust else "symmetry_resp.txt"
    try:
        run_python(write_symmetry, resp1, out_name, cwd=scratch_dir)
        path = os.path.join(scratch_dir, out_name)
        return path if os.path.isfile(path) else None
    except Exception:  # pragma: no cover - logging handled upstream
        logging.debug("[step5][symmetry-doc] generation failed", exc_info=True)
        return None


def process_conformer(cfg: ConformerConfig, scratch_dir: str, original_dir: str, input_files: list[str], defer_finalize: bool) -> ConformerChargeResult:
    """Execute conformer charge pipeline inside prepared scratch directory and return structured result.

    Pipeline (fine-grained stages):
      1. prepare_inputs -> PreparedInputs
      2. run Gaussian (adapter)
      3. generate ESP
      4. prepare RESP (inputs + symmetry) -> RespPrepResult
      5. write symmetry doc file
      6. run RESP two-stage fitting
      7. optional deferred finalize

    Backward compatibility: Signature & behaviour preserved; only internal structuring changed.
    """
    # Stage: selective config + decision logging (side-effect only)
    try:
        log_relevant_config(
            step="step5",
            cfg=cfg,
            fields=[
                "molecule_name",
                "net_charge",
                "protocol",
                "adjust_sym",
                "ignore_sym",
                "symmetry_ensemble_mode",
                "residue_name",
            ],
        )
        try:
            dec = build_conformer_decision(cfg.protocol, cfg.adjust_sym, cfg.ignore_sym)
            if cfg.symmetry_ensemble_mode is not None:
                dec.extra.append(("symmetry.mode.ensemble", cfg.symmetry_ensemble_mode))
            dec.log('step5')
        except Exception:
            logging.debug('[step5][decisions] logging failed', exc_info=True)
    except Exception:
        logging.debug("[step5][log] selective config logging failed", exc_info=True)

    # 1. Input preparation
    prepared = prepare_inputs(cfg, scratch_dir)

    # 2. Gaussian stage (structured)
    gaussian_result: GaussianResult = gaussian_adapter.run_with_result(
        mol2_file=prepared.mol2_file,
        name=cfg.molecule_name,
        net_charge=cfg.net_charge,
        cwd=scratch_dir,
        try_cache=True,
    )

    # 3. ESP generation
    esp_file = generate_esp(cfg, gaussian_result, scratch_dir)

    # 4. RESP preparation (includes resolving inputs + symmetry apply)
    resp_prep = prepare_resp(cfg, scratch_dir)

    # 5. Symmetry documentation
    symmetry_path = write_symmetry_doc(resp_prep.resp1, cfg.adjust_sym, scratch_dir)

    # 6. RESP fitting
    resp_result: RespStageResult = resp_adapter.run_two_stage_with_result(
        name=cfg.molecule_name,
        esp_file=esp_file,
        resp1=resp_prep.resp1,
        resp2=resp_prep.resp2,
        cwd=scratch_dir,
        symmetry_applied=resp_prep.symmetry_applied,
    )

    # 7. Optional deferred finalize
    if defer_finalize:
        from electrofit.infra.scratch_manager import finalize_scratch_directory
        finalize_scratch_directory(
            original_dir,
            scratch_dir,
            input_files,
            output_files=None,
            overwrite=True,
            remove_parent_if_empty=False,
            reason="defer-finalize",
        )
    logging.info("Conformer processing completed for %s", cfg.molecule_name)
    # Log measured wall-clock times collected from adapters so they appear in process logs
    try:
        gtime = gaussian_result.wall_time_s if hasattr(gaussian_result, 'wall_time_s') else None
        rtime = resp_result.wall_time_s if hasattr(resp_result, 'wall_time_s') else None
        if gtime is not None or rtime is not None:
            logging.info(
                "[timing] gaussian=%s s resp=%s s",
                f"{gtime:.3f}" if gtime is not None else "-",
                f"{rtime:.3f}" if rtime is not None else "-",
            )
    except Exception:
        logging.debug("[timing] failed to log wall times", exc_info=True)

    return ConformerChargeResult(
        gaussian=gaussian_result,
        resp=resp_result,
        protocol=cfg.protocol,
        residue_name=cfg.residue_name,
        scratch_dir=scratch_dir,
        symmetry_file=symmetry_path if symmetry_path and os.path.isfile(symmetry_path) else None,
    )

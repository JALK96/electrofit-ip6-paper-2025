"""Extraction helper for step4 isolated from orchestrator for multiprocessing stability.

Provides `_extract_for_molecule` used by both the main step4 orchestrator and
child processes. This avoids circular imports caused by importing the full
`step4` module inside the worker process.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Tuple
import mdtraj as md  # noqa: F401

from electrofit.config.loader import load_config, dump_config, resolve_symmetry_flags
from electrofit.infra.config_snapshot import compose_snapshot
from electrofit.infra.logging import setup_logging, reset_logging
from electrofit.infra.step_logging import log_relevant_config
from electrofit.infra.decisions import build_sampling_decision
from electrofit.domain.sampling import select_frame_indices, prepare_conformer_directory

# (Removed optional tomllib import; unused in current logic.)


def _extract_for_molecule(
    mol_proc_dir: Path,
    project_root: Path,
    sample: int,
    method: str,
    seed: int | None,
    override_cfg: Path | None,
    multi_mol: bool,
    verbose: bool,
) -> Tuple[bool, str]:
    # (Code duplicated from original step4; consider DRY refactor later.)
    sim_dir = mol_proc_dir / "run_gmx_simulation"
    pis_dir = mol_proc_dir / "run_gau_create_gmx_in"
    if not sim_dir.is_dir():
        return False, "no sim dir"
    cfg = load_config(project_root, context_dir=sim_dir, molecule_name=mol_proc_dir.name)
    proj = cfg.project
    molecule_name = proj.molecule_name or mol_proc_dir.name
    residue_name = getattr(proj, 'residue_name', None) or "LIG"
    logging.info(f"[step4][{mol_proc_dir.name}] using residue_name={residue_name}")
    print(f"[step4][debug] residue_name={residue_name}")
    protocol = getattr(proj, "protocol", "bcc")
    adjust_sym, ignore_sym = resolve_symmetry_flags(cfg, "ensemble")

    respin1_file = respin2_file = equiv_groups_file = None
    if protocol == "opt":
        respin1_file = pis_dir / ("ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN")
        respin2_file = pis_dir / "ANTECHAMBER_RESP2.IN"
    elif protocol == "bcc" and adjust_sym:
        json_candidates = sorted(pis_dir.glob("*.json"))
        if json_candidates:
            equiv_groups_file = json_candidates[0]

    input_mol2_file = project_root / "data" / "input" / mol_proc_dir.name / f"{molecule_name}.mol2"
    if not input_mol2_file.is_file():
        logging.warning(f"[step4][{mol_proc_dir.name}] missing input mol2 ({input_mol2_file.name}); continuing (bond insertion skipped)")
    extracted_conforms_dir = mol_proc_dir / "extracted_conforms"
    extracted_conforms_dir.mkdir(exist_ok=True)
    reset_logging()
    setup_logging(str(extracted_conforms_dir / "process.log"), also_console=False)
    existing_snapshot = extracted_conforms_dir / "electrofit.toml"
    if existing_snapshot.is_file():
        try:
            cfg_existing = load_config(project_root, context_dir=extracted_conforms_dir, molecule_name=mol_proc_dir.name)
            logging.info(f"[config] existing extracted_conforms snapshot for {mol_proc_dir.name} detected -> dump below")
            dump_config(cfg_existing, log_fn=logging.debug)
        except Exception:
            logging.debug("[step4] existing snapshot dump failed", exc_info=True)

    try:
        symmetry_json_present = any(pis_dir.glob('*.json'))
        build_sampling_decision(
            protocol=protocol,
            adjust_sym=adjust_sym,
            ignore_sym=ignore_sym,
            sampling_method=method,
            sample_count=sample,
            seed=seed,
            symmetry_json_present=symmetry_json_present,
        ).log('step4')
        # Update project section for logging parity
        proj.adjust_symmetry = adjust_sym  # type: ignore[attr-defined]
        proj.ignore_symmetry = ignore_sym  # type: ignore[attr-defined]
        log_relevant_config('step4', proj, ['molecule_name','residue_name','protocol','adjust_symmetry','ignore_symmetry'])
    except Exception:
        logging.debug('[step4][decisions] logging failed', exc_info=True)

    traj_path = sim_dir / "md_center.xtc"
    gro_path = sim_dir / "md.gro"
    if not traj_path.is_file() or not gro_path.is_file():
        return False, "missing md_center.xtc or md.gro"

    parent_cfg_target = compose_snapshot(
        extracted_conforms_dir,
        project_root,
        mol_proc_dir.name,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        upstream=sim_dir / "electrofit.toml",
        process_cfg=mol_proc_dir / "electrofit.toml",
        molecule_input=project_root / "data" / "input" / mol_proc_dir.name / "electrofit.toml",
        project_defaults=project_root / "electrofit.toml",
        extra_override=override_cfg,
    )

    logging.info(f"Step4: start extraction for molecule_dir={mol_proc_dir.name} method={method} sample={sample}")
    try:
        dump_config(cfg, log_fn=logging.debug)
    except Exception:
        logging.debug("[step4] dump per-molecule config failed", exc_info=True)

    raw_traj = md.load(str(traj_path), top=str(gro_path))
    try:
        # mdtraj's residue.atoms can be an iterator; calling len() on it fails.
        # Count safely without materializing the full list.
        res_counts: dict[str, int] = {}
        for res in raw_traj.topology.residues:
            try:
                count = sum(1 for _ in res.atoms)
            except Exception:
                # Fallback: if the residue provides a direct attribute for atom count
                count = getattr(res, 'n_atoms', 0) or getattr(res, 'atom_count', 0)
            res_counts[res.name] = res_counts.get(res.name, 0) + count
        inv_str = ", ".join(f"{k}:{v}" for k, v in sorted(res_counts.items()))
        logging.info(f"[step4][{mol_proc_dir.name}] residue inventory -> {inv_str}")
    except Exception:
        logging.debug("[step4] residue inventory logging failed", exc_info=True)
    ipl = raw_traj.top.select(f"resname {residue_name}")
    if len(ipl) == 0:
        msg = f"residue '{residue_name}' not in topology"
        logging.error(f"[step4][{mol_proc_dir.name}] {msg}; abort extraction.")
        return False, msg
    traj = raw_traj.atom_slice(ipl)
    if traj.n_atoms == 0:
        logging.warning(f"[step4][{mol_proc_dir.name}] zero atoms after selection; skipping")
        return False, "no atoms after selection"

    total = len(traj)
    n = min(sample, total)
    indices = select_frame_indices(traj, n, method, seed)
    configs = [traj[i] for i in indices]
    logging.info(f"Selected indices (n={len(indices)}) -> {indices}")

    for i, c in enumerate(configs):
        prepare_conformer_directory(
            conform_index=i,
            molecule_name=molecule_name,
            parent_cfg_target=parent_cfg_target,
            override_cfg=override_cfg,
            protocol=protocol,
            respin1_file=respin1_file,
            respin2_file=respin2_file,
            equiv_groups_file=equiv_groups_file,
            pis_dir=pis_dir,
            extracted_conforms_dir=extracted_conforms_dir,
            input_mol2_file=input_mol2_file,
            traj_frame_save_fn=c.save_pdb,
            verbose=verbose,
        )

    logging.info(f"Completed extraction: {len(configs)} conformers (method={method}) for {mol_proc_dir.name}")
    return True, f"Extracted {len(configs)} conformers (method={method}) to {extracted_conforms_dir}"

__all__ = ["_extract_for_molecule"]

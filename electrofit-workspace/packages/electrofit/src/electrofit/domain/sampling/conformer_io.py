"""Utilities for preparing conformer directories during extraction.

Provides a single high-level function that ensures snapshot propagation and
copies protocol-specific auxiliary files (RESP inputs, symmetry groups) before
writing the PDB.
"""
from __future__ import annotations
from pathlib import Path
import shutil
import logging
from electrofit.io.mol2_ops import mol2_to_pdb_with_bonds

__all__ = ["prepare_conformer_directory"]


def prepare_conformer_directory(
    conform_index: int,
    molecule_name: str,
    parent_cfg_target: Path,
    override_cfg: Path | None,
    protocol: str,
    respin1_file: Path | None,
    respin2_file: Path | None,
    equiv_groups_file: Path | None,
    pis_dir: Path,
    extracted_conforms_dir: Path,
    input_mol2_file: Path,
    traj_frame_save_fn,
    verbose: bool,
) -> Path:
    """Create/update a single conformer directory and persist its PDB.

    Parameters are deliberately explicit (rather than passing a large context
    object) to keep test surfaces narrow.
    Returns path to written PDB file.
    """
    conform_dir = extracted_conforms_dir / f"{molecule_name}c{conform_index}"
    conform_dir.mkdir(exist_ok=True)
    logging.info(f"Creating conformer dir {conform_dir.name}")
    snap_local = conform_dir / "electrofit.toml"
    if parent_cfg_target.is_file():
        try:
            refresh = (not snap_local.is_file()) or (
                parent_cfg_target.stat().st_mtime > snap_local.stat().st_mtime
            ) or bool(override_cfg)
            if refresh:
                shutil.copy2(parent_cfg_target, snap_local)
                logging.info(
                    f"[snapshot] {'created' if not snap_local.exists() else 'refreshed'} {snap_local.relative_to(conform_dir)} from parent"
                )
        except Exception:  # pragma: no cover
            logging.debug("[snapshot] copy failure", exc_info=True)

    if protocol == "opt":
        if respin1_file and respin1_file.is_file():
            shutil.copy(str(respin1_file), conform_dir)
        if respin2_file and respin2_file.is_file():
            shutil.copy(str(respin2_file), conform_dir)
    if not equiv_groups_file and (pis_dir / "equiv_groups.json").is_file():
        equiv_groups_file = pis_dir / "equiv_groups.json"  # type: ignore
    if equiv_groups_file and equiv_groups_file.is_file():
        tgt = conform_dir / "equiv_groups.json"
        if not tgt.exists():
            try:
                shutil.copy(str(equiv_groups_file), tgt)
            except Exception:  # pragma: no cover
                logging.debug("[symmetry] copy failed", exc_info=True)

    conform_name = f"{molecule_name}c{conform_index}.pdb"
    conform_path = conform_dir / conform_name
    traj_frame_save_fn(str(conform_path))
    if input_mol2_file.is_file():
        mol2_to_pdb_with_bonds(
            input_file=str(input_mol2_file),
            existing_pdb_file=str(conform_path),
            verbose=verbose,
        )
    logging.info(f"Wrote conformer {conform_name}")
    return conform_path

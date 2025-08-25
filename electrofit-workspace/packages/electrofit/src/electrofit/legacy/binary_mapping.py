"""Legacy helpers for historical binary->IP directory naming & mol2 renaming.

Formerly in ``electrofit.io.files``; retained temporarily for potential
one-off migration scripts. Planned removal once confirmed unused externally.
"""
from __future__ import annotations

import os
import logging

__all__ = [
    'binary_to_ip',
    'rename_mol2_binary',
]

# Static mapping preserved verbatim
binary_to_ip = {
    "010101": "IP21",
    "101010": "IP42",
    "101100": "IP44",
    "111000": "IP56",
    "010111": "IP23",
    "101101": "IP45",
    "111001": "IP57",
    "011111": "IP31",
    "111011": "IP59",
    "101111": "IP47",
    "111101": "IP61",
}

def rename_mol2_binary(base_dir: str, binary: str) -> None:
    """Rename a *.mol2 inside ``run_gau_create_gmx_in`` folder using the binary mapping.

    Legacy one-off helper. Logs warnings instead of printing to stdout.
    """
    ip_name = binary_to_ip.get(binary)
    if ip_name is None:
        logging.warning("[legacy.binary_mapping] binary '%s' not in mapping", binary)
        return
    folder_path = os.path.join(base_dir, "run_gau_create_gmx_in")
    if not os.path.isdir(folder_path):  # pragma: no cover
        logging.warning("[legacy.binary_mapping] folder missing: %s", folder_path)
        return
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.mol2'):
            old_path = os.path.join(folder_path, file_name)
            new_path = os.path.join(folder_path, f"{ip_name}.mol2")
            try:
                os.rename(old_path, new_path)
                logging.info("[legacy.binary_mapping] Renamed '%s' -> '%s'", old_path, new_path)
            except OSError as e:  # pragma: no cover
                logging.warning("[legacy.binary_mapping] rename failed %s -> %s: %s", old_path, new_path, e)
            break

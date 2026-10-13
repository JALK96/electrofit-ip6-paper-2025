from __future__ import annotations

from pathlib import Path

import pytest

from electrofit_analysis.cli.coordination.summerize_nap_dist_count_ip6 import (
    _resolve_coord_bool_file,
)


def _coord_dir(root: Path, pattern: str) -> Path:
    path = root / f"IP_{pattern}" / "analyze_final_sim" / "IonP_coordination"
    path.mkdir(parents=True, exist_ok=True)
    return path


def test_resolve_coord_bool_prefers_canonical_basis_specific_file(tmp_path: Path) -> None:
    coord = _coord_dir(tmp_path, "101101")
    canonical = coord / "IonP_coordination_bool_phosphorus_0p3200nm.npy"
    alias = coord / "IonP_nearest_bool_phosphorus_0p3200nm.npy"
    canonical.write_bytes(b"x")
    alias.write_bytes(b"y")

    resolved = _resolve_coord_bool_file(
        str(tmp_path),
        "101101",
        "analyze_final_sim",
        coord_distance_basis="phosphorus",
    )
    assert resolved == str(canonical)


def test_resolve_coord_bool_ambiguity_requires_override(tmp_path: Path) -> None:
    coord = _coord_dir(tmp_path, "101101")
    (coord / "IonP_coordination_bool_oxygen_0p2200nm.npy").write_bytes(b"x")
    (coord / "IonP_coordination_bool_oxygen_0p3200nm.npy").write_bytes(b"y")

    with pytest.raises(FileNotFoundError, match="Ambiguous coordination bool files"):
        _resolve_coord_bool_file(
            str(tmp_path),
            "101101",
            "analyze_final_sim",
            coord_distance_basis="oxygen",
        )


def test_resolve_coord_bool_override_works_with_multiple_candidates(tmp_path: Path) -> None:
    coord = _coord_dir(tmp_path, "101101")
    target = coord / "IonP_coordination_bool_oxygen_0p3200nm.npy"
    (coord / "IonP_coordination_bool_oxygen_0p2200nm.npy").write_bytes(b"x")
    target.write_bytes(b"y")

    resolved = _resolve_coord_bool_file(
        str(tmp_path),
        "101101",
        "analyze_final_sim",
        coord_distance_basis="oxygen",
        coord_bool_name=target.name,
    )
    assert resolved == str(target)


def test_resolve_coord_bool_falls_back_to_legacy_oxygen_name(tmp_path: Path) -> None:
    coord = _coord_dir(tmp_path, "101101")
    legacy = coord / "IonP_coordination_bool.npy"
    legacy.write_bytes(b"x")

    resolved = _resolve_coord_bool_file(
        str(tmp_path),
        "101101",
        "analyze_final_sim",
        coord_distance_basis="oxygen",
    )
    assert resolved == str(legacy)

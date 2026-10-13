from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electrofit_analysis.cli.coordination.Na_IP6_coordination import (
    _assignment_from_min_dists,
    _with_basis_and_cutoff,
    _write_legacy_unsuffixed_aliases,
    main as coordination_main,
)


def test_basis_suffix_filename_generation() -> None:
    out = _with_basis_and_cutoff(
        path=Path("IonP_coordination_bool.npy"),
        basis="phosphorus",
        cutoff_nm=0.32,
    )
    assert out.name == "IonP_coordination_bool_phosphorus_0p3200nm.npy"


def test_assignment_can_differ_between_basis_geometry() -> None:
    cutoff_a = 2.5

    min_dists_oxygen = np.asarray(
        [
            [2.0, 3.0],  # ion0 assigned to phosphate 0
            [3.0, 2.0],  # ion1 assigned to phosphate 1
        ],
        dtype=float,
    )
    min_dists_phosphorus = np.asarray(
        [
            [3.0, 2.0],  # ion0 assigned to phosphate 1
            [2.0, 3.0],  # ion1 assigned to phosphate 0
        ],
        dtype=float,
    )

    idx_o, mask_o, contacts_o = _assignment_from_min_dists(min_dists_oxygen, cutoff_a)
    idx_p, mask_p, contacts_p = _assignment_from_min_dists(min_dists_phosphorus, cutoff_a)

    assert idx_o.tolist() == [0, 1]
    assert idx_p.tolist() == [1, 0]
    assert mask_o.tolist() == [True, True]
    assert mask_p.tolist() == [True, True]
    assert contacts_o.shape == (2, 2)
    assert contacts_p.shape == (2, 2)


def test_phosphorus_mode_rejects_per_oxygen_rdf_flags(tmp_path) -> None:
    with pytest.raises(ValueError, match="coord_distance_basis='phosphorus'"):
        coordination_main(
            project_dir=str(tmp_path),
            project_mode="project",
            coord_distance_basis="phosphorus",
            rdf_oxygen_mode="both",
            coord_cutoff_nm=0.32,
            skip_rdf=True,
        )


def test_phosphorus_mode_rejects_per_phosphate_scale(tmp_path) -> None:
    with pytest.raises(ValueError, match="coord_distance_basis='phosphorus'"):
        coordination_main(
            project_dir=str(tmp_path),
            project_mode="project",
            coord_distance_basis="phosphorus",
            rdf_scale="per-phosphate",
            coord_cutoff_nm=0.32,
            skip_rdf=True,
        )


def test_legacy_unsuffixed_alias_written_only_for_oxygen(tmp_path) -> None:
    dest = tmp_path / "IonP_coordination"
    dest.mkdir()

    oxygen_src_bool = dest / "IonP_coordination_bool_oxygen_0p3200nm.npy"
    oxygen_src_contacts = dest / "IonP_contacts_bool_oxygen_0p3200nm.npy"
    np.save(oxygen_src_bool, np.zeros((1, 6, 1), dtype=bool))
    np.save(oxygen_src_contacts, np.ones((1, 6, 1), dtype=bool))

    _write_legacy_unsuffixed_aliases(
        dest_dir=dest,
        basis="oxygen",
        cutoff_values=[0.32],
        bool_path_by_cutoff={0.32: oxygen_src_bool},
        contacts_path_by_cutoff={0.32: oxygen_src_contacts},
        overwrite=True,
    )
    assert (dest / "IonP_coordination_bool.npy").exists()
    assert (dest / "IonP_nearest_bool.npy").exists()
    assert (dest / "IonP_contacts_bool.npy").exists()

    for legacy in ("IonP_coordination_bool.npy", "IonP_nearest_bool.npy", "IonP_contacts_bool.npy"):
        (dest / legacy).unlink()

    phosphorus_src_bool = dest / "IonP_coordination_bool_phosphorus_0p3200nm.npy"
    phosphorus_src_contacts = dest / "IonP_contacts_bool_phosphorus_0p3200nm.npy"
    np.save(phosphorus_src_bool, np.zeros((1, 6, 1), dtype=bool))
    np.save(phosphorus_src_contacts, np.ones((1, 6, 1), dtype=bool))

    _write_legacy_unsuffixed_aliases(
        dest_dir=dest,
        basis="phosphorus",
        cutoff_values=[0.32],
        bool_path_by_cutoff={0.32: phosphorus_src_bool},
        contacts_path_by_cutoff={0.32: phosphorus_src_contacts},
        overwrite=True,
    )
    assert not (dest / "IonP_coordination_bool.npy").exists()
    assert not (dest / "IonP_nearest_bool.npy").exists()
    assert not (dest / "IonP_contacts_bool.npy").exists()

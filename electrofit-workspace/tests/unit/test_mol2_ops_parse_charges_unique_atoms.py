from electrofit.io.mol2_ops import parse_charges_from_mol2


def test_parse_charges_from_mol2_does_not_collapse_duplicate_atom_names(tmp_path):
    mol2 = tmp_path / "dup_names.mol2"
    mol2.write_text(
        "\n".join(
            [
                "@<TRIPOS>MOLECULE",
                "TEST",
                " 4 3 0 0 0",
                "SMALL",
                "GASTEIGER",
                "",
                "@<TRIPOS>ATOM",
                "      1 O        0.0000    0.0000    0.0000   O.3   1 MOL1         0.1000",
                "      2 O        0.0000    0.0000    0.0000   O.3   1 MOL1         0.2000",
                "      3 C1       0.0000    0.0000    0.0000   C.3   1 MOL1        -0.3000",
                "      4 H        0.0000    0.0000    0.0000     H   1 MOL1         0.0000",
                "@<TRIPOS>BOND",
                "      1 1 3 1",
                "      2 2 3 1",
                "      3 3 4 1",
                "",
            ]
        )
        + "\n"
    )

    atoms = parse_charges_from_mol2(str(mol2))
    # Should return one entry per atom line, even when MOL2 atom names are reused.
    assert list(atoms.keys()) == ["O1", "O2", "C1", "H1"]
    assert [atoms[k]["charges"] for k in atoms.keys()] == [[0.1], [0.2], [-0.3], [0.0]]
    assert abs(sum(atoms[k]["charges"][0] for k in atoms.keys())) < 1e-12


import json
from pathlib import Path
import pytest

from electrofit.domain.aggregation import average_charges as avg_mod


@pytest.mark.parametrize(
    "remove_outlier,adjust_sym,calc_group_avg",
    [
        (True, False, False),              # cleaned_average_charges path
        (True, True, True),                # cleaned_adjusted_group_average_charges path
        (False, False, False),             # baseline simple average path
    ],
)
def test_step6_acpype_invocation_matrix(tmp_path, monkeypatch, remove_outlier, adjust_sym, calc_group_avg):
    project_root = tmp_path
    (project_root / 'data' / 'input' / 'MOL').mkdir(parents=True)
    # Root project TOML
    (project_root / 'electrofit.toml').write_text('[project]\ncharge = 0\n')

    # Process molecule directories
    mol_dir = project_root / 'process' / 'MOL'
    extracted = mol_dir / 'extracted_conforms'
    pis_dir = mol_dir / 'run_gau_create_gmx_in'
    acp_dir = pis_dir / 'dummy.acpype'
    acp_dir.mkdir(parents=True)
    extracted.mkdir(parents=True)

    # Minimal mol2 file required by directory scan (content unused due to monkeypatch)
    (acp_dir / 'test_gaff2.mol2').write_text('@<TRIPOS>MOLECULE\nTEST\n')

    # Per-molecule config enabling flags
    mol_cfg = ['[project]']
    if adjust_sym:
        mol_cfg.append('adjust_symmetry = true')
    if calc_group_avg:
        mol_cfg.append('calculate_group_average = true')
    (mol_dir / 'electrofit.toml').write_text('\n'.join(mol_cfg) + '\n')

    # Symmetry JSON if needed
    if adjust_sym:
        sym = {"C1": ["O1"]}
        (extracted / 'equiv_groups.json').write_text(json.dumps(sym))

    # Monkeypatch I/O heavy functions inside imported module namespace
    calls = {"acpype": 0, "mol2_updates": []}

    def fake_parse_charges_from_mol2(_path):
        return {"C1": {"charges": [0.1]}, "O1": {"charges": [-0.1]}}

    def fake_extract_charges_from_subdirectories(_base, _res):
        # Provide an obvious outlier in last conformer when remove_outlier=True
        if remove_outlier:
            return {
                "C1": {"charges": [0.10, 0.11, 0.09, 1.50]},
                "O1": {"charges": [-0.10, -0.11, -0.09, -2.00]},
            }
        return {
            "C1": {"charges": [0.10, 0.11, 0.09, 0.10]},
            "O1": {"charges": [-0.10, -0.11, -0.09, -0.10]},
        }

    def fake_update_mol2_charges(inp, chg, out):  # noqa: D401
        calls["mol2_updates"].append((chg, out))
        Path(out).write_text('dummy mol2')

    def fake_run_acpype(mol2_file, charge, outdir, atom_type, charges):  # noqa: D401
        calls["acpype"] += 1
        # simulate acpype directory creation
        Path(outdir).mkdir(exist_ok=True)

    monkeypatch.setattr(avg_mod, 'parse_charges_from_mol2', fake_parse_charges_from_mol2)
    monkeypatch.setattr(avg_mod, 'extract_charges_from_subdirectories', fake_extract_charges_from_subdirectories)
    monkeypatch.setattr(avg_mod, 'update_mol2_charges', fake_update_mol2_charges)
    monkeypatch.setattr(avg_mod, 'run_acpype', fake_run_acpype)
    if adjust_sym:
        monkeypatch.setattr(avg_mod, 'load_symmetry_groups', lambda path: {"C1": ["O1"]})

    ok, msg = avg_mod.process_molecule_average_charges(
        mol_dir=mol_dir,
        project_root=project_root,
        override_cfg=None,
        multi_mol=False,
        remove_outlier=remove_outlier,
        plot_histograms=False,
        hist_combine_groups=False,
        hist_bins=10,
        outlier_iqr_factor=1.5,
    )
    assert ok, msg
    # Exactly one acpype invocation expected in all scenarios
    assert calls["acpype"] == 1, f"acpype calls={calls['acpype']} scenario={(remove_outlier,adjust_sym,calc_group_avg)}"

    results_dir = mol_dir / 'results'
    if remove_outlier:
        # Expect cleaned charge artefacts; group-averaged cleaned may be skipped if config flag not effectively applied
        if adjust_sym and calc_group_avg and (results_dir / 'cleaned_adjusted_group_average_charges.chg').is_file():
            assert (results_dir / 'cleaned_adjusted_group_average_charges.chg').is_file()
        else:
            assert (results_dir / 'cleaned_average_charges.chg').is_file()
    else:
        assert (results_dir / 'average_charges.chg').is_file()

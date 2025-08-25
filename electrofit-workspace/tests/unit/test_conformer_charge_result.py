import os
import pathlib
import pytest

from electrofit.domain.charges.process_conformer import ConformerConfig, process_conformer

@pytest.mark.parametrize("adjust_sym", [False, True])
def test_conformer_charge_result_cache_and_symmetry(monkeypatch, tmp_path, adjust_sym):
    """Variant A: fast test using Gaussian cache + monkeypatched espgen & RESP to avoid heavy external runs.

    Validates structured result invariants without invoking real Gaussian/RESP.
    """
    # Use absolute path anchored at repository root for robustness across cwd changes
    repo_root = pathlib.Path(__file__).resolve().parent.parent  # tests/
    fixtures = repo_root / 'fixtures'
    cache_dir = fixtures / 'gaussian_cache'
    input_pdb = fixtures / 'conformer_min' / 'input.pdb'
    assert input_pdb.is_file(), "Missing minimal input PDB fixture"
    assert (cache_dir / 'TESTCONF.gesp').is_file(), "Cache .gesp missing"

    # Enable Gaussian cache
    monkeypatch.setenv('ELECTROFIT_DEBUG_GAUSSIAN_CACHE', str(cache_dir))

    # Scratch setup
    scratch_dir = tmp_path / 'scratch'
    scratch_dir.mkdir()
    pdb_target = scratch_dir / 'TESTCONF.pdb'
    pdb_target.write_text(input_pdb.read_text())

    # Monkeypatch run_espgen inside process_conformer module (imported symbol) -> create trivial .esp
    def fake_run_espgen(gesp_file, esp_file, cwd):
        assert (pathlib.Path(cwd) / gesp_file).is_file()
        (pathlib.Path(cwd) / esp_file).write_text("ESP DUMMY")
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.run_espgen', fake_run_espgen)

    # Monkeypatch run_python to avoid spawning external interpreter for symmetry file generation
    def fake_run_python(func, *args, cwd=None):
        # Directly call the function with cwd applied via chdir context
        prev = os.getcwd()
        try:
            if cwd:
                os.chdir(cwd)
            func(*args)
        finally:
            os.chdir(prev)
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.run_python', fake_run_python)

    # Monkeypatch RESP two-stage
    class DummyRespResult:
        def __init__(self, name, resp1, resp2, symmetry):
            self.name = name
            self.resp1_in = resp1
            self.resp2_in = resp2
            self.resp1_out = f"{name}-resp1.out"
            self.resp2_out = f"{name}-resp2.out"
            self.final_mol2 = f"{name}_resp.mol2"
            self.charges_file = f"{name}-resp2.chg"
            self.symmetry_applied = symmetry
            self.wall_time_s = 0.0
            self.warnings = []

    def fake_run_two_stage_with_result(name, esp_file, resp1, resp2, cwd, symmetry_applied):
        p = pathlib.Path(cwd)
        for f in [f"{name}-resp1.out", f"{name}-resp2.out", f"{name}-resp2.chg", f"{name}_resp.mol2"]:
            (p / f).write_text("placeholder")
        return DummyRespResult(name, resp1, resp2, symmetry_applied)
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.resp_adapter.run_two_stage_with_result', fake_run_two_stage_with_result)
    # Bypass real symmetry adjustment which expects JSON file
    def fake_apply_symmetry(scratch_dir, adjust_sym, ignore_sym):
        return None
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.resp_adapter.apply_symmetry', fake_apply_symmetry)
    # Provide prepare_ac and generate_inputs replacements
    def fake_prepare_ac(mol2_file, name, net_charge, cwd):
        p = pathlib.Path(cwd) / f"{name}.ac"
        p.write_text("AC DUMMY")
        return p.name
    def fake_generate_inputs(name, cwd):
        p = pathlib.Path(cwd)
        (p / 'ANTECHAMBER_RESP1.IN').write_text('RESP1')
        (p / 'ANTECHAMBER_RESP2.IN').write_text('RESP2')
        # If adjust_sym True, Domain expects MOD variant later; create both to simplify.
        (p / 'ANTECHAMBER_RESP1_MOD.IN').write_text('RESP1MOD')
        return ('ANTECHAMBER_RESP1.IN', 'ANTECHAMBER_RESP2.IN')
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.resp_adapter.prepare_ac', fake_prepare_ac)
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.resp_adapter.generate_inputs', fake_generate_inputs)
    # Dummy symmetry writer
    def fake_write_symmetry(input_filename, output_filename, cwd=None, atomic_number_mapping=None):
        target = pathlib.Path(cwd or '.') / output_filename
        target.write_text('SYMMETRY PLACEHOLDER')
        return str(target)
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.write_symmetry', fake_write_symmetry)

    cfg = ConformerConfig(
        molecule_name='TESTCONF',
        pdb_file=str(pdb_target),
        net_charge=0,
        residue_name='LIG',
        adjust_sym=adjust_sym,
        ignore_sym=False,
        protocol='bcc',
    )
    # Optional explicit cache hydration (diagnostic)
    from electrofit.adapters.gaussian import maybe_use_cache as _maybe
    _maybe('TESTCONF', str(scratch_dir))
    # Strenge Prüfung: Adapter muss .gesp liefern, sonst Test abbrechen
    assert (scratch_dir / 'TESTCONF.gesp').is_file(), 'Gaussian cache hydration failed (.gesp missing)'

    result = process_conformer(cfg, str(scratch_dir), str(tmp_path), ['TESTCONF.pdb'], defer_finalize=False)

    assert result.gaussian.used_cache is True
    assert result.gaussian.warnings == []
    assert result.gaussian.gesp_file and os.path.isfile(result.gaussian.gesp_file)
    assert result.resp.symmetry_applied == adjust_sym
    assert result.symmetry_file is not None and os.path.isfile(result.symmetry_file)
    if adjust_sym:
        assert result.symmetry_file.endswith('symmetry_resp_MOD.txt')
    else:
        assert result.symmetry_file.endswith('symmetry_resp.txt')
    p = pathlib.Path(scratch_dir)
    for f in [f'TESTCONF-resp1.out', f'TESTCONF-resp2.out', f'TESTCONF-resp2.chg', f'TESTCONF_resp.mol2']:
        assert (p / f).is_file(), f"Missing {f}"

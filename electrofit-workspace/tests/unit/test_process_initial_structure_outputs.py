import os, sys, pathlib

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))

from electrofit.domain.prep.process_initial import InitialPrepConfig, process_initial


def test_process_initial_structure_bcc_creates_gmx_files(tmp_path, monkeypatch):
    work = tmp_path / 'work'
    work.mkdir()
    mol2 = work / 'lig.mol2'
    mol2.write_text("@<TRIPOS>MOLECULE\nLIG\n 1 0 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C 1 LIG 0.0\n")
    equiv = work / 'equiv_groups.json'; equiv.write_text('{}')

    # Monkeypatch run_acpype used by domain function to fabricate outputs
    def _fake_run_acpype(mol2_file, net_charge, scratch_dir, atom_type, charges="bcc"):
        base = os.path.splitext(os.path.basename(mol2_file))[0]
        ac_dir = os.path.join(scratch_dir, f"{base}.acpype")
        os.makedirs(ac_dir, exist_ok=True)
        for name, content in {
            f'{base}_GMX.gro': ';gro',
            f'{base}_GMX.itp': ';itp',
            f'{base}_GMX.top': '[ system ]\nX\n\n[ molecules ]\nX 1\n',
            f'posre_{base}.itp': '[ position_restraints ]\n',
        }.items():
            (pathlib.Path(ac_dir)/name).write_text(content)
    # No monkeypatch needed now; we'll pass _fake_run_acpype directly

    old = os.getcwd(); os.chdir(work)
    try:
        # Simulate what legacy shim did: setup scratch manually
        scratch_base = tmp_path/'scratch_base'
        scratch_base.mkdir()
        # Minimal imitation of setup_scratch_directory (copy inputs)
        scratch_dir = scratch_base / 'scratch'
        scratch_dir.mkdir()
        for f in (mol2, equiv):
            (scratch_dir / f.name).write_text(f.read_text())
        cfg = InitialPrepConfig(
            molecule_name='lig',
            mol2_file=str(mol2.name),  # relative inside scratch
            net_charge=0,
            residue_name='LIG',
            atom_type='gaff2',
            adjust_sym=False,
            ignore_sym=False,
            protocol='bcc'
        )
        process_initial(cfg, str(scratch_dir), str(work), [str(mol2), str(equiv)], run_acpype=_fake_run_acpype)
    finally:
        os.chdir(old)

    # Outputs are exposed inside scratch_dir (as per _expose_acpype_outputs logic)
    scratch_outputs = list(os.listdir(scratch_dir))
    assert 'lig_GMX.gro' in scratch_outputs
    assert 'lig_GMX.itp' in scratch_outputs
    assert 'lig_GMX.top' in scratch_outputs
    assert any(n.startswith('posre_lig') and n.endswith('.itp') for n in scratch_outputs)

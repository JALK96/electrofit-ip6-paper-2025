import sys, textwrap, subprocess
from pathlib import Path
import numpy as np
import mdtraj as md

sys.path.append(str(Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))


def _write_minimal_traj(sim_dir: Path, n_frames=4, n_atoms=3):
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue('MOLX', chain)
    for i in range(n_atoms):
        top.add_atom(f'C{i}', element=md.element.carbon, residue=res)
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    for f in range(n_frames):
        xyz[f, :, 0] = f * 0.05
    traj = md.Trajectory(xyz, top)
    traj[0].save_gro(str(sim_dir / 'md.gro'))
    traj.save_xtc(str(sim_dir / 'md_center.xtc'))


def test_step4_reads_sampling_section(tmp_path: Path):
    proj = tmp_path
    sim_dir = proj / 'process' / 'MOLX' / 'run_gmx_simulation'
    sim_dir.mkdir(parents=True)
    _write_minimal_traj(sim_dir)
    # Project-level sampling config
    (proj / 'electrofit.toml').write_text('[sampling]\nmethod = "random"\ncount = 3\nseed = 7\n')
    # Per-molecule input TOML to supply molecule properties
    mol_input_dir = proj / 'data' / 'input' / 'MOLX'
    mol_input_dir.mkdir(parents=True, exist_ok=True)
    (mol_input_dir / 'electrofit.toml').write_text(textwrap.dedent('''\
        [project]
        molecule_name = "MOLX"
        residue_name  = "MOLX"
        protocol      = "bcc"
    '''))

    cmd = [sys.executable, '-m', 'electrofit', 'step4', '--project', str(proj)]
    out = subprocess.run(cmd, text=True, capture_output=True, cwd=proj)
    assert out.returncode == 0, out.stderr + out.stdout
    assert 'method=random' in out.stdout

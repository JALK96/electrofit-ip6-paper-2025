"""CLI und Scratch Management Tests.

Enthält:
- test_cli_help_has_usage_sections: schneller Smoke-Test der CLI-Hilfe.
- test_scratch_finalize_on_exception: stellt sicher, dass ensure_finalized bei Exceptions aufräumt.

Runtime-Flag Tests wurden nach integration/test_runtime_gmx_flags.py verschoben.
"""
import os, sys, subprocess, pathlib

# Paketpfad hinzufügen (lokaler Quellcode)
sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))


def test_cli_help_has_usage_sections():
    out = subprocess.run([sys.executable,'-m','electrofit','--help'], text=True, capture_output=True)
    assert out.returncode == 0
    assert 'usage:' in out.stdout.lower()
    step_out = subprocess.run([sys.executable,'-m','electrofit','step4','--help'], text=True, capture_output=True)
    assert step_out.returncode == 0
    assert 'usage:' in step_out.stdout.lower()


def test_scratch_finalize_on_exception(tmp_path, monkeypatch):
    from electrofit.scratch.manager import setup_scratch_directory
    from electrofit.cli.safe_run import ensure_finalized

    work = tmp_path/'work'; work.mkdir()
    (work/'in.txt').write_text('data')

    old = os.getcwd(); os.chdir(work)
    try:
        scratch, orig = setup_scratch_directory(['in.txt'], base_scratch_dir=str(tmp_path/'scratch_base'))
        try:
            with ensure_finalized(original_dir=orig, scratch_dir=scratch, input_files=['in.txt']):
                raise RuntimeError('boom')
        except RuntimeError:
            pass
    finally:
        os.chdir(old)

    assert not os.path.isdir(scratch)
    assert (work/'in.txt').exists()

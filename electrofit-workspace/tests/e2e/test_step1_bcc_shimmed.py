import os, subprocess, sys
from pathlib import Path

def test_step1_bcc_shimmed(tmp_path, shim_bin, monkeypatch):
    # Put shim on PATH so run_command finds it
    monkeypatch.setenv("PATH", str(shim_bin)+os.pathsep+os.environ["PATH"])

    # minimal input files
    work = tmp_path/"run"
    work.mkdir()
    (work/"IP_000001.mol2").write_text("@<TRIPOS>MOLECULE\nm\n")
    (work/"electrofit.toml").write_text("[compute]\nremote_host='' # local\n[step1]\nname='IP_000001'\nnet_charge=0\nresidue='LIG'\nprotocol='bcc'\n")

    # run CLI (replace with your actual command)
    cmd = [sys.executable, "-m", "electrofit", "step1", "--project", str(tmp_path)]
    out = subprocess.run(cmd, cwd=work, text=True, capture_output=True)
    assert out.returncode == 0, out.stderr + out.stdout

    # check that finalize moved outputs back (your finalize may write to cwd)
    assert any(p.name.endswith((".itp",".top",".gro")) for p in work.iterdir())
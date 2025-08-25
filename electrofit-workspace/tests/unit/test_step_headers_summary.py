import os, subprocess, sys, json, re, shutil
from pathlib import Path

# Minimal integration-ish unit test that each step writes exactly one header line per invocation to step.log
# (Header de-dup active). We invoke lightweight steps only (0 and 6) with controlled dirs.

def _run(cmd, cwd):
    cp = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    assert cp.returncode == 0, cp.stderr + cp.stdout


def test_step0_and_step6_single_headers(tmp_path, monkeypatch):
    monkeypatch.setenv("ELECTROFIT_DEDUP_HEADERS", "1")
    proj = tmp_path / "proj"
    (proj / "data" / "input" / "MOL1").mkdir(parents=True)
    # Pre-create empty process dir (CLI prints creation message else, but we ensure it exists to simplify assertions)
    (proj / "process").mkdir(parents=True, exist_ok=True)
    # Run step0 for single molecule
    cmd_base = [sys.executable, "-m", "electrofit"]
    _run(cmd_base + ["step0", "--project", str(proj), "--molecule", "MOL1"], proj)
    # fabricate minimal structure required for step6 aggregator context (results directories rely on earlier steps; we fake structure)
    (proj/"process"/"MOL1"/"results").mkdir(parents=True, exist_ok=True)
    # Also create extracted_conforms dir + dummy .acpype dir and mol2 so step6 returns ok
    extracted = proj/"process"/"MOL1"/"extracted_conforms"
    extracted.mkdir(parents=True, exist_ok=True)
    (extracted/"electrofit.toml").write_text("[dummy]\n")
    acp = proj/"process"/"MOL1"/"run_gau_create_gmx_in"/"foo.acpype"
    acp.mkdir(parents=True, exist_ok=True)
    (acp/"MOL1_gaff2.mol2").write_text("@<TRIPOS>MOLECULE\nMOL1\n")
    _run(cmd_base + ["step6", "--project", str(proj), "--molecule", "MOL1"], proj)
    step_log = (proj/"step.log")
    assert step_log.is_file()
    lines = step_log.read_text().splitlines()
    assert len([l for l in lines if "step=step0" in l]) == 1, lines
    assert len([l for l in lines if "step=step6" in l]) == 1, lines

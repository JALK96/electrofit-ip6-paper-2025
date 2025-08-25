import subprocess, sys
def test_cli_help():
    out = subprocess.check_output([sys.executable, "-m", "electrofit", "--help"], text=True)
    assert "Usage" in out or "help" in out.lower()
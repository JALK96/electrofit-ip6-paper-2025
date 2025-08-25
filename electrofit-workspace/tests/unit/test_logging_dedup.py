import os, tempfile, logging
from pathlib import Path

def test_header_dedup_basic(monkeypatch):
    # Force dedup enabled via env
    monkeypatch.setenv("ELECTROFIT_DEDUP_HEADERS", "1")
    from electrofit.infra.logging import setup_logging, log_run_header, reset_logging
    with tempfile.TemporaryDirectory() as td:
        log_path = Path(td)/"test.log"
        setup_logging(str(log_path), also_console=False, suppress_initial_message=True)
        log_run_header("stepX")
        log_run_header("stepX")  # should be suppressed
        reset_logging()
        setup_logging(str(log_path), also_console=False, suppress_initial_message=True)
        log_run_header("stepX")  # after reset allowed again
        content = log_path.read_text().splitlines()
        emitted = [l for l in content if "step=stepX" in l]
        # Expect 2 (first before reset, one after reset)
        assert len(emitted) == 2, content


"""Waiting for running commands before scratch finalization on a signal."""
import subprocess
import sys
import time

from electrofit.cli import run_commands
from electrofit.cli.safe_run import _stop_running_commands

CHILD = """
import signal, sys, time
def stop(*_):
    open(sys.argv[1], "w").write("checkpoint")
    sys.exit(0)
signal.signal(signal.SIGTERM, stop)
print("x" * 200000, flush=True)
time.sleep(float(sys.argv[2]))
"""


def _start(tmp_path, sleep):
    flag = tmp_path / "flag"
    proc = subprocess.Popen([sys.executable, "-c", CHILD, str(flag), str(sleep)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with run_commands._ACTIVE_LOCK:
        run_commands._ACTIVE.add(proc)
    return proc, flag


def _forget(proc):
    with run_commands._ACTIVE_LOCK:
        run_commands._ACTIVE.discard(proc)


def test_child_that_ignores_the_signal_gets_one_term(tmp_path):
    proc, flag = _start(tmp_path, sleep=60)
    try:
        start = time.monotonic()
        _stop_running_commands(grace=0.5, term_wait=10)
        assert proc.poll() == 0
        assert flag.read_text() == "checkpoint"
        assert time.monotonic() - start < 10
    finally:
        proc.kill()
        _forget(proc)


def test_child_that_exits_by_itself_is_not_signalled(tmp_path):
    proc, flag = _start(tmp_path, sleep=1)
    try:
        _stop_running_commands(grace=10, term_wait=10)
        assert proc.poll() == 0
        assert not flag.exists()
    finally:
        proc.kill()
        _forget(proc)


SLOW_CHILD = """
import signal, sys, time
def stop(*_):
    time.sleep(2)
    open(sys.argv[1], "w").write("checkpoint")
    sys.exit(0)
signal.signal(signal.SIGTERM, stop)
time.sleep(60)
"""


def test_wait_covers_the_child_of_a_shell(tmp_path):
    """The command runs as sh -c (no exec); the wait must last until the child has finished."""
    flag = tmp_path / "flag"
    script = tmp_path / "child.py"
    script.write_text(SLOW_CHILD)
    proc = subprocess.Popen(["sh", "-c", f"'{sys.executable}' '{script}' '{flag}'; true"],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with run_commands._ACTIVE_LOCK:
        run_commands._ACTIVE.add(proc)
    try:
        time.sleep(0.5)
        _stop_running_commands(grace=0.5, term_wait=10)
        assert flag.read_text() == "checkpoint"
    finally:
        proc.kill()
        _forget(proc)

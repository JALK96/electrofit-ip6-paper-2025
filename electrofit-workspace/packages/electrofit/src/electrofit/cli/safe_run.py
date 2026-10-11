"""
Utilities to ensure scratch directories are finalized even when a run is
interrupted (Ctrl+C) or terminated (SIGTERM).

Supports multiple concurrent/serial molecule runs in the same Python process:
each call to `register_scratch(...)` pushes a new registration onto a stack.
On shutdown/interrupt we finalize them in LIFO order.

Typical usage in a step module after creating the scratch directory:

    from electrofit.cli.safe_run import register_scratch, ensure_finalized

    # Option A: safety net via atexit/signal handlers only
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
    register_scratch(
        original_dir=original_dir,
        scratch_dir=scratch_dir,
        input_files=input_files,
        # output_files=None,      # optional
        # overwrite=True,         # optional
        # remove_parent_if_empty=True,  # optional
    )

    # Option B: scoped guarantee (finalize immediately when the block exits)
    with ensure_finalized(
        original_dir=original_dir,
        scratch_dir=scratch_dir,
        input_files=input_files,
        # output_files=None,
        # overwrite=True,
        # remove_parent_if_empty=True,
    ):
        ... do the work ...

The global cleanup is invoked automatically on normal exit, Ctrl+C, or SIGTERM.
On a signal it first waits for commands started by run_command to exit (GROMACS
writes its checkpoint), so the copied-back files are complete. It is idempotent
and safe to call multiple times.
"""
from __future__ import annotations

import atexit
import signal
import sys
import threading
import traceback
from contextlib import contextmanager
from typing import Iterable, Optional, List

from electrofit.infra.scratch_manager import finalize_scratch_directory


class _Reg:
    __slots__ = (
        "original_dir",
        "scratch_dir",
        "input_files",
        "output_files",
        "overwrite",
        "remove_parent_if_empty",
    )

    def __init__(
        self,
        *,
        original_dir: str,
        scratch_dir: str,
        input_files: Iterable[str],
        output_files: Optional[Iterable[str]],
        overwrite: bool,
        remove_parent_if_empty: bool,
    ) -> None:
        self.original_dir = str(original_dir)
        self.scratch_dir = str(scratch_dir)
        self.input_files = list(input_files)
        self.output_files = list(output_files) if output_files is not None else None
        self.overwrite = bool(overwrite)
        self.remove_parent_if_empty = bool(remove_parent_if_empty)


class _ScratchState:
    """Holds a stack of scratch registrations for the current process."""
    __slots__ = ("_regs", "_lock", "_cleaned")

    def __init__(self) -> None:
        self._regs: List[_Reg] = []
        self._lock = threading.Lock()
        self._cleaned = False

    def push(self, reg: _Reg) -> None:
        with self._lock:
            # if we've already globally cleaned (after an earlier shutdown path),
            # allow further registrations and reset the cleaned flag
            if self._cleaned:
                self._cleaned = False
            self._regs.append(reg)

    def pop_last(self) -> Optional[_Reg]:
        with self._lock:
            if self._regs:
                return self._regs.pop()
            return None

    def snapshot(self) -> List[_Reg]:
        with self._lock:
            return list(self._regs)

    def clear(self) -> None:
        with self._lock:
            self._regs.clear()

    def mark_cleaned(self) -> None:
        with self._lock:
            self._cleaned = True

    def was_cleaned(self) -> bool:
        with self._lock:
            return self._cleaned

    def is_empty(self) -> bool:
        with self._lock:
            return not self._regs


_state = _ScratchState()


def register_scratch(
    *,
    original_dir: str,
    scratch_dir: str,
    input_files: Iterable[str],
    output_files: Optional[Iterable[str]] = None,
    overwrite: bool = True,
    remove_parent_if_empty: bool = True,
) -> None:
    """Register the current run's scratch info for safe finalization.

    You may call this once per run after `setup_scratch_directory(...)`.
    Multiple calls (for multiple molecules/runs) are supported; they will be
    finalized in LIFO order on shutdown/interrupt.
    """
    _state.push(
        _Reg(
            original_dir=original_dir,
            scratch_dir=scratch_dir,
            input_files=input_files,
            output_files=output_files,
            overwrite=overwrite,
            remove_parent_if_empty=remove_parent_if_empty,
        )
    )


def _finalize_one(reg: _Reg, reason: str) -> None:
    try:
        msg = f"[electrofit] finalize scratch (reason={reason}): {reg.scratch_dir} -> {reg.original_dir}\n"
        try:
            sys.stderr.write(msg)
        except Exception:
            pass

        finalize_scratch_directory(
            original_dir=reg.original_dir,
            scratch_dir=reg.scratch_dir,
            input_files=reg.input_files,
            output_files=reg.output_files,
            overwrite=reg.overwrite,
            remove_parent_if_empty=reg.remove_parent_if_empty,
        )
        try:
            sys.stderr.write(f"[electrofit][debug] finalize complete for {reg.scratch_dir}\n")
        except Exception:
            pass
    except Exception:
        # Avoid raising exceptions during interpreter shutdown; print a brief report.
        try:
            sys.stderr.write("[electrofit] WARNING: finalize_scratch_directory raised an exception during cleanup.\n")
            traceback.print_exc()
        except Exception:
            pass


def _cleanup(reason: str = "atexit") -> None:
    # Idempotent global cleanup; safe to run multiple times.
    if _state.was_cleaned():
        return
    regs = list(reversed(_state.snapshot()))
    for reg in regs:
        _finalize_one(reg, reason)
    _state.clear()
    _state.mark_cleaned()


def _descendants(proc) -> list:
    """psutil handles of the processes below a child started by run_command.

    Commands run through ``/bin/sh -c``, so the GROMACS process can be a grandchild.
    """
    try:
        import psutil

        return psutil.Process(proc.pid).children(recursive=True)
    except Exception:
        return []


def _alive(process) -> bool:
    try:
        import psutil

        return process.is_running() and process.status() != psutil.STATUS_ZOMBIE
    except Exception:
        return False


def _signal_command(proc, descendants, signum) -> None:
    """Signal the running descendants of a command, or the command itself if it has none.

    The shell is not signalled: it would exit at once and leave its GROMACS child
    writing the checkpoint unobserved.
    """
    targets = [d for d in descendants if _alive(d)]
    if not targets:
        targets = [proc] if proc.poll() is None else []
    for target in targets:
        try:
            target.send_signal(signum)
        except Exception:
            pass


def _wait_draining(procs, descendants, timeout: float) -> None:
    """Wait until procs and their descendants exit, discarding pipe output meanwhile."""
    import os
    import select
    import time

    def running():
        return any(p.poll() is None for p in procs) or any(_alive(d) for d in descendants)

    fds = set()
    for p in procs:
        for stream in (p.stdout, p.stderr):
            try:
                if stream is not None and not stream.closed:
                    fds.add(stream.fileno())
            except Exception:
                pass
    deadline = time.monotonic() + timeout
    while running() and time.monotonic() < deadline:
        if not fds:
            time.sleep(0.5)
            continue
        try:
            readable = select.select(list(fds), [], [], 0.5)[0]
        except Exception:
            readable = []
        for fd in readable:
            try:
                if not os.read(fd, 65536):
                    fds.discard(fd)
            except Exception:
                fds.discard(fd)


def _stop_running_commands(grace: float = 15.0, term_wait: float = 120.0) -> None:
    """Let commands started by run_command exit before scratch is finalized.

    A terminal interrupt also reaches the children (same process group); GROMACS
    then writes a checkpoint and exits. Commands still running after ``grace``
    seconds, for example when only this process was signalled, receive one
    SIGTERM, and SIGKILL after ``term_wait`` more seconds. Waiting covers the
    descendants recorded at the start, so a shell that exits early does not end
    the wait while its GROMACS child is still writing.
    """
    try:
        from electrofit.cli.run_commands import active_processes

        procs = active_processes()
    except Exception:
        return
    if not procs:
        return
    tree = {proc: _descendants(proc) for proc in procs}
    everything = [d for ds in tree.values() for d in ds]
    _wait_draining(procs, everything, grace)
    for proc, descendants in tree.items():
        _signal_command(proc, descendants, signal.SIGTERM)
    _wait_draining(procs, everything, term_wait)
    for proc, descendants in tree.items():
        for target in [d for d in descendants if _alive(d)] + ([proc] if proc.poll() is None else []):
            try:
                target.send_signal(signal.SIGKILL)
            except Exception:
                pass


_handling_signal = False
# A signal that arrives while ensure_finalized copies scratch back is handled after
# the copy completes, so the run directory never holds half-copied files.
_finalizing = False
_pending_signal = None


def _signal_handler(signum, frame):  # type: ignore[no-untyped-def]
    # Wait for running commands, perform cleanup, then exit with proper signal
    # semantics. Further signals during this bounded wait are ignored.
    global _handling_signal, _pending_signal
    if _finalizing:
        if _pending_signal is None:
            _pending_signal = signum
        return
    if _handling_signal:
        return
    _handling_signal = True
    name = {getattr(signal, k): k for k in ("SIGINT", "SIGTERM")}.get(signum, str(signum))
    _stop_running_commands()
    _cleanup(reason=name)
    # Re-raise default behavior after cleanup for proper exit status
    signal.signal(signum, signal.SIG_DFL)
    try:
        signal.raise_signal(signum)
    except Exception:
        sys.exit(1)


@contextmanager
def ensure_finalized(
    *,
    original_dir: str,
    scratch_dir: str,
    input_files: Iterable[str],
    output_files: Optional[Iterable[str]] = None,
    overwrite: bool = True,
    remove_parent_if_empty: bool = True,
):
    """Context manager to finalize *this* registration on block exit.

    This also registers the scratch with the global handlers as a safety net.
    When the block exits, we pop and finalize this one explicitly, so the
    global cleanup (if it later runs) won't re-finalize it. SIGINT/SIGTERM
    received during that copy-back is handled once the copy has completed.
    """
    global _finalizing
    try:
        sys.stderr.write(f"[ensure-finalized-enter] scratch={scratch_dir} orig={original_dir}\n")
    except Exception:
        pass
    register_scratch(
        original_dir=original_dir,
        scratch_dir=scratch_dir,
        input_files=input_files,
        output_files=output_files,
        overwrite=overwrite,
        remove_parent_if_empty=remove_parent_if_empty,
    )
    try:
        try:
            sys.stderr.write(f"[ensure-finalized-before-yield] scratch={scratch_dir}\n")
        except Exception:
            pass
        yield
    finally:
        try:
            sys.stderr.write(f"[ensure-finalized-after-yield] scratch={scratch_dir}\n")
        except Exception:
            pass
        reg = _state.pop_last()
        if reg is not None:
            try:
                sys.stderr.write(f"[ensure-finalized-finalize-start] scratch={scratch_dir}\n")
            except Exception:
                pass
            _finalizing = True
            try:
                _finalize_one(reg, reason="context-exit")
            finally:
                _finalizing = False
            try:
                sys.stderr.write(f"[ensure-finalized-finalize-done] scratch={scratch_dir}\n")
            except Exception:
                pass
            if _pending_signal is not None:
                _signal_handler(_pending_signal, None)


# Register atexit cleanup once on import
atexit.register(_cleanup)

# Attach signal handlers for interactive/termination signals
for _sig in (signal.SIGINT, signal.SIGTERM):
    try:
        signal.signal(_sig, _signal_handler)
    except Exception:
        # Some environments may not allow setting handlers; ignore
        pass
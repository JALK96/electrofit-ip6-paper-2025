"""Isolated single conformer runner (pipeline version).

It executes exactly one conformer directory using the domain batch logic
and emits a single parsable RESULT:<json> line identical to the legacy
implementation for backwards compatibility.

CLI (intended for internal invocation by step5 orchestrator):
  --conf <path to conformer dir>
  --project <project root>
  [--override <override cfg toml>]
  [--multi-mol]  (flag)
  [--mock]       (flag)
  [--verbose]    (flag)

Exit code: 0 on success, 1 on failure (mirrors ok flag).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from electrofit.domain.charges.conformer_batch import process_conformer_dir

__all__ = ["main"]


def main(argv: list[str] | None = None):  # pragma: no cover - exercised indirectly
    ap = argparse.ArgumentParser(description="Isolated single conformer processor (pipeline).")
    ap.add_argument("--conf", required=True, help="Path to conformer directory")
    ap.add_argument("--project", required=True, help="Project root (contains process/)")
    ap.add_argument("--override", help="Optional override config TOML path")
    ap.add_argument("--multi-mol", action="store_true", help="Indicate multi-molecule project context")
    ap.add_argument("--mock", action="store_true", help="Mock mode (skip heavy Gaussian/RESP)")
    ap.add_argument("--verbose", action="store_true", help="Verbose logging passthrough")
    args = ap.parse_args(argv)

    result = process_conformer_dir(
        Path(args.conf),
        Path(args.project),
        Path(args.override).resolve() if args.override else None,
        args.multi_mol,
        args.mock,
        args.verbose,
    )
    # result expected shape: (relative_path_str, ok_bool, message_str)
    if not (isinstance(result, tuple) and len(result) == 3):  # defensive
        rec = {"rel": str(args.conf), "ok": False, "msg": f"unexpected result shape {result!r}"}
        print("RESULT:" + json.dumps(rec), flush=True)
        sys.exit(1)
    rel, ok, msg = result
    rec = {"rel": rel, "ok": bool(ok), "msg": msg}
    print("RESULT:" + json.dumps(rec), flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":  # pragma: no cover
    main()

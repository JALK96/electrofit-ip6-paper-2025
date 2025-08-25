"""CLI entry point for writing RESP symmetry mapping.

Wraps `electrofit.domain.symmetry.write_symmetry`.
"""
from __future__ import annotations
import argparse
import sys
from electrofit.domain.symmetry import write_symmetry


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate symmetry mapping from RESP input file.")
    p.add_argument("input", help="Path to RESP input file (e.g. ANTECHAMBER_RESP1.IN)")
    p.add_argument("output", help="Output text file path for symmetry mapping")
    p.add_argument("--cwd", help="Execute within this directory (optional)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        out = write_symmetry(args.input, args.output, cwd=args.cwd)
        print(f"Wrote symmetry file to {out}")
        return 0
    except Exception as e:  # pragma: no cover
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

"""CLI entry point for charge plotting.

Reads JSON dictionaries with initial and (one or two) sampled charge sets and
invokes visualization helper functions. Keeps plotting logic in `electrofit.viz`.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any


from electrofit.io.symmetry import load_symmetry_groups
from electrofit.viz.helpers import (
    plot_charges_by_atom,
    plot_charges_by_symmetry,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot charges per atom or symmetry groups (if provided)."
    )
    parser.add_argument(
        "-ic", "--initial-charges", required=True, help="Path to initial_charges_dict.json"
    )
    parser.add_argument("-c", "--charges", required=True, help="Path to charges_dict.json")
    parser.add_argument(
        "-c2", "--charges2", help="Optional path to second charges_dict.json for comparison"
    )
    parser.add_argument(
        "-s", "--symmetry-groups", help="Optional path to symmetry_groups.json"
    )
    parser.add_argument(
        "-d", "--directory", default=".", help="Output directory for plots/files"
    )
    parser.add_argument(
        "--format", choices=["pdf", "png"], default="pdf", help="Output figure format"
    )
    return parser.parse_args(argv)


def _load_json(path: str | Path, desc: str) -> Any:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"ERROR: {desc} file not found: {path}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"ERROR: failed to parse {desc} JSON {path}: {e}", file=sys.stderr)
        sys.exit(1)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    initial_charges_dict = _load_json(args.initial_charges, "initial_charges")
    atoms_dict1 = _load_json(args.charges, "charges")
    atoms_dict2 = _load_json(args.charges2, "charges2") if args.charges2 else None

    base_dir = args.directory
    os.makedirs(base_dir, exist_ok=True)

    if args.symmetry_groups:
        symmetry_groups = load_symmetry_groups(args.symmetry_groups)
        plot_charges_by_symmetry(
            atoms_dict1, initial_charges_dict, base_dir, symmetry_groups, atoms_dict2
        )
    else:
        plot_charges_by_atom(atoms_dict1, initial_charges_dict, base_dir, atoms_dict2)

    # Optionally convert produced PDFs to PNG if requested
    if args.format == "png":
        # For each produced pdf in directory, create png sibling (simple heuristic)
        for name in ("charges.pdf", "charges_comparison.pdf", "charges_by_symmetry.pdf", "charges_by_symmetry_comparison.pdf"):
            pdf_path = Path(base_dir) / name
            if pdf_path.exists():
                try:
                    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: F401 (ensures backend)
                    # Re-load figure is non-trivial; simplest is to advise user to regenerate with format directly.
                    # Placeholder: could integrate saving alternate format in viz helpers.
                except Exception:
                    pass  # Skip silent for now
        # NOTE: real multi-format support can be added by extending viz helpers with a format param.

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

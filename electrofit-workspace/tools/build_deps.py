#!/usr/bin/env python3
# tools/build_deps.py
from __future__ import annotations
import argparse
import ast
import sys
from pathlib import Path

WORKSPACE_PKGS = {"electrofit", "electrofit_fep", "electrofit_analysis"}

# Heuristic mapping: import name -> PyPI package
NAME_MAP = {
    "sklearn": "scikit-learn",
    "PIL": "pillow",
    "bs4": "beautifulsoup4",
    "yaml": "PyYAML",
    "Bio": "biopython",
    "MDAnalysis": "MDAnalysis",
    "MDAnalysisTests": None,  # don't depend on tests
    "ParmEd": "parmed",
    "pymbar": "pymbar",
    "rdkit": "rdkit",  # conda-first; leave as hint if needed
    "openbabel": "openbabel-wheel",  # wheels via PyPI; alt: conda openbabel
    "pybel": "pybel",
    "alchemlyb": "alchemlyb",
    "mdtraj": "mdtraj",
    "pytng": "pytng",
    "pyedr": "pyedr",
    "gsd": "gsd",
    "GridDataFormats": "GridDataFormats",
    "mmtf": "mmtf-python",
    "mmtf_python": "mmtf-python",
    "seaborn": "seaborn",
    "matplotlib": "matplotlib",
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "typer": "typer",
    "rich": "rich",
    "plotly": "plotly",
    "pyscf": "pyscf",
    "jax": "jax",
    "jaxlib": "jaxlib",
    "psycopg2": "psycopg2-binary",
    "SQLAlchemy": "SQLAlchemy",
    "requests": "requests",
}

# Obvious external/CLI ecosystems we don't want in [project.dependencies]
EXTERNAL_HINTS = {
    "antechamber": "AmberTools CLI",
    "parmchk2": "AmberTools CLI",
    "resp": "AmberTools CLI",
    "espgen": "AmberTools CLI",
    "sander": "AmberTools CLI",
    "MMPBSA": "AmberTools CLI",
    "pytraj": "pytraj (often conda; optional)",
    "gmx": "GROMACS CLI",
    "Gaussian": "Gaussian CLI",
}

def top_name(mod: str | None) -> str | None:
    if not mod:
        return None
    return mod.split(".", 1)[0]

def is_stdlib(name: str) -> bool:
    if hasattr(sys, "stdlib_module_names"):
        return name in sys.stdlib_module_names
    # fallback: a few common stdlib modules
    return name in {
        "os", "sys", "re", "math", "json", "pathlib", "subprocess", "typing",
        "argparse", "itertools", "functools", "collections", "dataclasses",
        "logging", "time", "datetime", "shutil", "gzip", "pickle", "hashlib",
        "base64", "statistics", "inspect", "enum", "tempfile", "glob",
    }

def guess_pypi(import_name: str) -> tuple[str | None, str]:
    """Return (pypi_name_or_None, note)"""
    if import_name in WORKSPACE_PKGS:
        return None, "workspace"
    if is_stdlib(import_name):
        return None, "stdlib"
    if import_name in EXTERNAL_HINTS:
        return None, f"external:{EXTERNAL_HINTS[import_name]}"
    if import_name in NAME_MAP:
        mapped = NAME_MAP[import_name]
        if mapped is None:
            return None, "ignore"
        return mapped, "mapped"
    # default: assume same name on PyPI
    return import_name, "direct"

def scan_py_imports(py: Path) -> set[str]:
    mods = set()
    try:
        tree = ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
    except Exception:
        return mods
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for n in node.names:
                t = top_name(n.name)
                if t:
                    mods.add(t)
        elif isinstance(node, ast.ImportFrom):
            # Skip relative imports and empty module (from . import x)
            if node.level and node.module is None:
                continue
            t = top_name(node.module)
            if t:
                mods.add(t)
    return mods

def collect_package_roots(workspace: Path) -> dict[str, Path]:
    roots: dict[str, Path] = {}
    for pkg_dir in (workspace / "packages").glob("*"):
        src = pkg_dir / "src"
        if src.is_dir():
            for top in src.iterdir():
                if top.is_dir():
                    roots[top.name] = top
    return roots

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workspace", default=".", help="Workspace root (has packages/)")
    ap.add_argument("--write", action="store_true", help="(future) write to pyproject")
    args = ap.parse_args()

    ws = Path(args.workspace).resolve()
    roots = collect_package_roots(ws)
    if not roots:
        print("No packages/*/src/<pkg> found.", file=sys.stderr)
        sys.exit(2)

    print(f"# Scanning workspace: {ws}")
    for pkg_name, root in sorted(roots.items()):
        imports: set[str] = set()
        for py in root.rglob("*.py"):
            imports |= scan_py_imports(py)

        # classify
        pypi: set[str] = set()
        external: set[str] = set()
        ignored: set[str] = set()
        for imp in sorted(imports):
            p, note = guess_pypi(imp)
            if note.startswith("external:"):
                external.add(f"{imp}  # {note.split(':',1)[1]}")
            elif note in {"workspace", "stdlib", "ignore"}:
                ignored.add(imp)
            elif p:
                pypi.add(p)

        # minimal tidy
        # if matplotlib used, you probably also want seaborn as optional
        # nothing automatic here; just the raw set.
        print(f"\n## Package: {pkg_name}")
        if pypi:
            print("[project.dependencies]")
            for d in sorted(pypi):
                print(f'#   "{d}"')
            # give a copy-paste friendly TOML block
            deps_array = ",\n  ".join(f'"{d}"' for d in sorted(pypi))
            print("\n# Suggested TOML snippet:")
            print("[project]")
            print(f'name = "{pkg_name.replace("_", "-")}"')
            print('version = "0.1.0"')
            print('requires-python = ">=3.12,<3.13"')
            print("dependencies = [")
            print(f"  {deps_array}")
            print("]")
        else:
            print("No third-party dependencies detected.")

        if external:
            print("\n# Likely external/system tools or conda-first libs (keep out of deps):")
            for e in sorted(external):
                print(f"#   {e}")

        # Always show a few ignored/imports for transparency
        if ignored:
            print("\n# Ignored (stdlib/workspace/tests):")
            print("#  " + ", ".join(sorted(ignored)))

if __name__ == "__main__":
    main()
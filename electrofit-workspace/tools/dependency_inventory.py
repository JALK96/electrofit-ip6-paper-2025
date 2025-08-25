#!/usr/bin/env python3
"""Erzeuge einen Import- und Abhängigkeits-Index für alle Workspace-Pakete.

Ziel: Während des Refactors eine stets aktualisierbare Übersicht haben:
  * Modul -> direkte Imports (intern / extern klassifiziert)
  * Reverse Dependencies
  * Paket-zu-Paket Kanten (electrofit, electrofit_analysis, electrofit_fep)
  * Externe Third‑Party Nutzung (Top-Level Package Namen + Häufigkeit)
  * Aus pyproject.toml gelesene deklarierte Dependencies

Ausgabe:
  tools/inventory/dependency_graph.json (Grunddaten)
  tools/inventory/dependency_report.md (optionaler Markdown Bericht)

Nutzung:
  python tools/dependency_inventory.py --format both

Später bei Refactor-Schritten erneut ausführen und ins Repo committen.

Einschränkungen / Fallstricke:
  - Dynamische / bedingte Imports (inside functions / if TYPE_CHECKING) werden erfasst,
    aber nicht nach Kontext unterschieden.
  - "from X import *" wird nur als Modul X gezählt, Symbole nicht aufgelöst.
  - Relative Imports werden aufgelöst soweit eindeutig.
  - Symbol-Nutzung (ob importierter Name wirklich verwendet wird) wird nicht analysiert.
"""
from __future__ import annotations

import argparse
import ast
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from collections import defaultdict, Counter
import datetime as _dt

try:
    import tomllib  # Python 3.11+
except ImportError:  # pragma: no cover
    tomllib = None  # type: ignore

WORKSPACE_ROOT = Path(__file__).resolve().parent.parent
SRC_PACKAGE_ROOTS = [
    WORKSPACE_ROOT / "packages" / "electrofit" / "src",
    WORKSPACE_ROOT / "packages" / "electrofit-fep" / "src",
    WORKSPACE_ROOT / "packages" / "electrofit-analysis" / "src",
]
INTERNAL_TOP_LEVEL = {"electrofit", "electrofit_fep", "electrofit_analysis"}


@dataclass
class ModuleRecord:
    name: str
    file: str
    imports: list[str]


def iter_python_files() -> list[Path]:
    files: list[Path] = []
    for root in SRC_PACKAGE_ROOTS:
        if not root.exists():
            continue
        for py in root.rglob("*.py"):
            if "__pycache__" in py.parts:
                continue
            files.append(py)
    return files


def module_name_from_path(path: Path) -> str:
    for root in SRC_PACKAGE_ROOTS:
        try:
            rel = path.relative_to(root)
        except ValueError:
            continue
        parts = rel.parts
        if not parts:
            raise RuntimeError(f"Unexpected empty relative path for {path}")
        if parts[-1] == "__init__.py":
            mod = ".".join(parts[:-1])
        else:
            mod = ".".join(parts).removesuffix(".py")
        return mod
    raise RuntimeError(f"File not under known roots: {path}")


def resolve_relative(base_module: str, level: int, name: str | None) -> str | None:
    if level == 0:
        return name
    parts = base_module.split(".")
    if parts[-1] == "__init__":
        parts = parts[:-1]
    if level > len(parts):
        return None
    prefix = parts[:-level]
    if name:
        prefix.append(name)
    return ".".join(prefix)


def collect_imports(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    base_mod = module_name_from_path(path)
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                top = alias.name
                found.append(top)
        elif isinstance(node, ast.ImportFrom):
            mod_name = node.module
            if node.level:
                mod_name = resolve_relative(base_mod, node.level, mod_name)
            if mod_name:
                found.append(mod_name)
    return sorted(set(found))


def build_graph() -> dict[str, ModuleRecord]:
    graph: dict[str, ModuleRecord] = {}
    for py in iter_python_files():
        mod = module_name_from_path(py)
        imps = collect_imports(py)
        graph[mod] = ModuleRecord(name=mod, file=str(py.relative_to(WORKSPACE_ROOT)), imports=imps)
    return graph


def load_pyproject_meta() -> dict[str, dict]:
    data: dict[str, dict] = {}
    for pkg in ["electrofit", "electrofit-fep", "electrofit-analysis"]:
        pyproject = WORKSPACE_ROOT / "packages" / pkg / "pyproject.toml"
        if not pyproject.exists() or tomllib is None:
            continue
        with pyproject.open("rb") as f:
            doc = tomllib.load(f)
        proj = doc.get("project", {})
        data[pkg] = {
            "name": proj.get("name"),
            "version": proj.get("version"),
            "dependencies": proj.get("dependencies", []),
            "scripts": proj.get("scripts", {}),
        }
    return data


def classify_import(mod: str) -> str:
    top = mod.split(".")[0]
    if top in INTERNAL_TOP_LEVEL:
        return "internal"
    return "external"


def build_package_edges(graph: dict[str, ModuleRecord]):
    from collections import Counter
    pkg_edges = Counter()
    for rec in graph.values():
        src_top = rec.name.split(".")[0]
        for imp in rec.imports:
            dst_top = imp.split(".")[0]
            if src_top in INTERNAL_TOP_LEVEL and dst_top in INTERNAL_TOP_LEVEL and src_top != dst_top:
                pkg_edges[(src_top, dst_top)] += 1
    return pkg_edges


def invert_graph(graph: dict[str, ModuleRecord]):
    reverse = defaultdict(list)
    for rec in graph.values():
        for imp in rec.imports:
            reverse[imp].append(rec.name)
    return {k: sorted(v) for k, v in reverse.items()}


def generate_markdown(graph: dict[str, ModuleRecord], meta: dict, out_json_path: Path) -> str:
    pkg_edges = build_package_edges(graph)
    external_counter = Counter()
    for rec in graph.values():
        for imp in rec.imports:
            if classify_import(imp) == "external":
                external_counter[imp.split(".")[0]] += 1
    lines: list[str] = []
    lines.append("# Dependency Inventory\n")
    lines.append(f"Erzeugt: {_dt.datetime.now(_dt.timezone.utc).isoformat()}  ")
    lines.append(f"Quelle JSON: `{out_json_path.relative_to(WORKSPACE_ROOT)}`\n")
    lines.append("## Pakete\n")
    for k, v in meta.items():
        lines.append(f"### {v.get('name')} ({k})\n")
        lines.append(f"Version: {v.get('version')}\n")
        if v.get("dependencies"):
            lines.append("Dependencies (deklariert):\n")
            for d in v["dependencies"]:
                lines.append(f"- {d}")
        if v.get("scripts"):
            lines.append("Scripts:\n")
            for sname, target in v["scripts"].items():
                lines.append(f"- {sname} -> {target}")
        lines.append("")
    lines.append("## Paket-Kanten (interne Cross-Package Imports)\n")
    if pkg_edges:
        for (src, dst), cnt in sorted(pkg_edges.items()):
            lines.append(f"- {src} -> {dst}: {cnt}")
    else:
        lines.append("(Keine)")
    lines.append("")
    lines.append("## Externe Top-Level Imports (Häufigkeit)\n")
    for name, cnt in external_counter.most_common():
        lines.append(f"- {name}: {cnt}")
    lines.append("")
    lines.append("## Modulzählung\n")
    lines.append(f"Anzahl Module: {len(graph)}\n")
    lines.append("## Hinweise\n")
    lines.append("- Dynamische / optionale Imports sind nicht speziell markiert.")
    lines.append("- Nutzung (Symbole) nicht aufgelöst; es handelt sich nur um Modulreferenzen.")
    lines.append("- Bei Refactor erneut Skript ausführen und Änderungen committen.")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--format", choices=["json", "md", "both"], default="both")
    ap.add_argument("--output-dir", default="tools/inventory", help="Ausgabeverzeichnis relativ zum Root")
    args = ap.parse_args()

    graph = build_graph()
    reverse = invert_graph(graph)
    meta = load_pyproject_meta()

    out_dir = WORKSPACE_ROOT / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    json_path = out_dir / "dependency_graph.json"
    if args.format in ("json", "both"):
        payload = {
            "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
            "modules": {k: asdict(v) for k, v in sorted(graph.items())},
            "reverse": reverse,
            "internal_top_level": sorted(INTERNAL_TOP_LEVEL),
            "package_meta": meta,
        }
        json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        print(f"Wrote {json_path}")

    if args.format in ("md", "both"):
        md = generate_markdown(graph, meta, json_path)
        md_path = out_dir / "dependency_report.md"
        md_path.write_text(md, encoding="utf-8")
        print(f"Wrote {md_path}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

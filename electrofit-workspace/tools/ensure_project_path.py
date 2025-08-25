from __future__ import annotations
from pathlib import Path
import io, os, re

root = Path("packages")
files = list(root.rglob("*.py"))

uses_rx   = re.compile(r"\bproject_path\b")
assign_rx = re.compile(r"^\s*project_path\s*=", re.MULTILINE)

def ensure_import_os(lines: list[str]) -> None:
    # If no "import os" present, inject it after the first block of imports.
    has_os = any(l.startswith("import os") or l.startswith("from os ") for l in lines)
    if has_os:
        return
    # find last top-of-file import
    last_imp = -1
    for i, l in enumerate(lines):
        if l.startswith(("import ", "from ")):
            last_imp = i
        elif l.strip() == "" or l.lstrip().startswith("#"):
            continue
        else:
            break
    insert_at = last_imp + 1 if last_imp >= 0 else 0
    lines.insert(insert_at, "import os\n")

def insert_project_path(lines: list[str]) -> None:
    # Insert right after the import block
    insert_at = 0
    for i, l in enumerate(lines):
        if l.startswith(("import ", "from ")):
            insert_at = i + 1
        elif l.strip() == "" or l.lstrip().startswith("#"):
            continue
        else:
            break
    snippet = [
        "PROJECT_PATH = os.environ.get(\"ELECTROFIT_PROJECT_PATH\", os.getcwd())\n",
        "project_path = PROJECT_PATH\n",
    ]
    # Avoid duplicate insertion
    body = "".join(lines)
    if "ELECTROFIT_PROJECT_PATH" not in body and "project_path = PROJECT_PATH" not in body:
        lines[insert_at:insert_at] = snippet

for py in files:
    text = py.read_text(encoding="utf-8")
    if not uses_rx.search(text):
        continue
    if assign_rx.search(text):
        continue  # already defines project_path
    lines = text.splitlines(keepends=True)
    ensure_import_os(lines)
    insert_project_path(lines)
    py.write_text("".join(lines), encoding="utf-8")
    print("inserted project_path shim:", py)
print("Done.")

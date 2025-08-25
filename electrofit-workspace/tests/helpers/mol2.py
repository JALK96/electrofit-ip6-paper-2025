# tests/helpers/mol2.py
from pathlib import Path

_MINIMAL_MOL2 = """@<TRIPOS>MOLECULE
{title}
 1 0 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 C1  0.0000  0.0000  0.0000  C.ar  1  LIG   -0.1000

@<TRIPOS>BOND
"""

def write_minimal_mol2(path: Path, title: str = "MOL1") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_MINIMAL_MOL2.format(title=title))
    return path
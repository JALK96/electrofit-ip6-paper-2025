import os, stat, pytest
from pathlib import Path

@pytest.fixture(scope="session")
def shim_bin(tmp_path_factory):
    """Provide a directory containing shim executables (gmx, acpype) for tests.

    Previously only available to e2e tests; now promoted to top-level so
    integration tests verifying logging isolation can reuse the same fake
    binaries.
    """
    d = tmp_path_factory.mktemp("shims")

    acpype = r"""#!/usr/bin/env bash
set -euo pipefail
in=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) in="$2"; shift 2;;
    *)  shift;;
  esac
done
base="$(basename "${in%.mol2}")"
out="${base}.acpype"
mkdir -p "$out"
touch "$out/${base}_GMX.gro" "$out/${base}_GMX.itp" "$out/${base}.top" "$out/posre_${base}.itp"
echo "acpype shim wrote $out"
"""
    (d/"acpype").write_text(acpype)
    os.chmod(d/"acpype", os.stat(d/"acpype").st_mode | stat.S_IEXEC)

    gmx = r"""#!/usr/bin/env bash
set -euo pipefail
sub="${1:-}"
shift || true
case "$sub" in
  editconf)
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) touch "$2"; shift 2;;
        *) shift;;
      esac
    done
    ;;
  solvate)
    out=""
    top=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) out="$2"; shift 2;;
        -p) top="$2"; shift 2;;
        *) shift;;
      esac
    done
    [[ -n "${out}" ]] && touch "${out}"
    ;;
  grompp)
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) touch "$2"; shift 2;;
        *) shift;;
      esac
    done
    ;;
  mdrun)
    deffnm=""
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -deffnm) deffnm="$2"; shift 2;;
        *) shift;;
      esac
    done
    [[ -n "${deffnm}" ]] && touch "${deffnm}.edr" "${deffnm}.gro" "${deffnm}.cpt" "${deffnm}.xtc"
    ;;
  energy|trjconv|mindist|gyrate)
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) touch "$2"; shift 2;;
        *) shift;;
      esac
    done
    ;;
  *) : ;;
esac
"""
    (d/"gmx").write_text(gmx)
    os.chmod(d/"gmx", os.stat(d/"gmx").st_mode | stat.S_IEXEC)

    return d

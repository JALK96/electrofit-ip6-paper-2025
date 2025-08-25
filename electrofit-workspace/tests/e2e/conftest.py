# tests/e2e/conftest.py
import os, stat, pytest
from pathlib import Path


@pytest.fixture(scope="session")
def shim_bin(tmp_path_factory):
    d = tmp_path_factory.mktemp("shims")

    # acpype: create a directory with expected GMX files
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

    # gmx: create the expected outputs per subcommand used in gromacs.py
    gmx = r"""#!/usr/bin/env bash
set -euo pipefail
sub="${1:-}"
shift || true
case "$sub" in
  editconf)
    # -o <name> present; touch it
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) touch "$2"; shift 2;;
        *) shift;;
      esac
    done
    ;;
  solvate)
    # -o <gro> and update -p <top> by just leaving it as-is
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
    # -o <tpr>
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) touch "$2"; shift 2;;
        *) shift;;
      esac
    done
    ;;
  mdrun)
    # -deffnm <name> -> produce <name>.edr, <name>.gro, <name>.cpt, <name>.xtc
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
    # create the -o output if present
    while [[ $# -gt 0 ]]; do
      case "$1" in
        -o) touch "$2"; shift 2;;
        *) shift;;
      esac
    done
    ;;
  *)
    # no-op for other subcommands in tests
    :
    ;;
esac
"""
    (d/"gmx").write_text(gmx)
    os.chmod(d/"gmx", os.stat(d/"gmx").st_mode | stat.S_IEXEC)

    return d
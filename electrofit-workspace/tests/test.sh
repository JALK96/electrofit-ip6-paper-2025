# Put the same shim on PATH manually
TMP=$(mktemp -d)
cat > "$TMP/acpype" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
in=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) in="$2"; shift 2;;
    *)  shift;;
  esac
done
base="$(basename "${in%.mol2}")"
touch "${base}_GMX.itp" "${base}_GMX.gro" "${base}.top"
SH
chmod +x "$TMP/acpype"
export PATH="$TMP:$PATH"

# Minimal project
RUN=$(mktemp -d)
mkdir -p "$RUN/run"
printf '@<TRIPOS>MOLECULE\nm\n' > "$RUN/run/IP_000001.mol2"
cat > "$RUN/run/electrofit.toml" <<'TOML'
[compute]
remote_host='' # local
[step1]
name='IP_000001'
net_charge=0
residue='LIG'
protocol='bcc'
TOML

python -m electrofit step1 --project "$RUN"
ls -1 "$RUN/run"
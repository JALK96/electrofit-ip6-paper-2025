import os, sys, pathlib, subprocess, re
import pytest

sys.path.append(str(pathlib.Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))
from tests.helpers.project import make_project_tree, install_fake_gmxdata


def _write_config(path: pathlib.Path, extra: str):
    cfg = path/'electrofit.toml'
    txt = cfg.read_text()
    cfg.write_text(txt + '\n' + extra + '\n')


def _install_shims(proj: pathlib.Path, with_runtime=False, threads=4, pin=True):
    if with_runtime:
        _write_config(proj, f'[gmx.runtime]\nthreads={threads}\npin={str(pin).lower()}\n')
    shim_dir = proj/'shim'; shim_dir.mkdir(exist_ok=True)
    (shim_dir/'acpype').write_text("""#!/usr/bin/env bash
base=IP_011101
outdir=${base}.acpype
mkdir -p "$outdir"
cat > "$outdir/${base}_GMX.gro" <<EOF
; gro
EOF
cat > "$outdir/${base}_GMX.itp" <<EOF
; itp
EOF
cat > "$outdir/${base}_GMX.top" <<EOF
[ system ]
X

[ molecules ]
X 1
EOF
echo '[ position_restraints ]' > "$outdir/posre_${base}.itp"
exit 0
""")
    (shim_dir/'acpype').chmod(0o755)
    (shim_dir/'gmx').write_text(r"""#!/usr/bin/env bash
cmd=$1; shift
capture_mdrun(){ echo "CMD:gmx mdrun $@" >> mdrun_calls.txt; touch em_steep.edr; echo '0 0' > em_steep.gro; }
find_out_arg(){ local flag=$1; shift; local prev=""; for a in "$@"; do if [[ $prev == $flag ]]; then echo "$a"; return 0; fi; prev=$a; done; return 1; }
case "$cmd" in
  editconf|solvate) out=$(find_out_arg -o "$@"); [[ -n $out ]] && echo X > "$out" ;;
  grompp) out=$(find_out_arg -o "$@"); [[ -n $out ]] && touch "$out" ;;
  genion) out=$(find_out_arg -o "$@"); [[ -n $out ]] && echo IONS > "$out" ;;
  mdrun) capture_mdrun "$@" ;;
  energy) out=$(find_out_arg -o "$@"); [[ -n $out ]] && echo -e '0 0\n1 0' > "$out" ;;
  mindist) echo -e '0 0\n1 0' > mindist.xvg ;;
  gyrate) echo -e '0 0\n1 0' > gyrate.xvg ;;
  trjconv) out=$(find_out_arg -o "$@"); [[ -n $out ]] && touch "$out" ;;
  *) ;;
 esac
exit 0
""")
    (shim_dir/'gmx').chmod(0o755)
    return shim_dir


def _run_steps(proj: pathlib.Path):
    for step in (1,2,3):
        r = subprocess.run([sys.executable,'-m','electrofit',f'step{step}','--project',str(proj)], cwd=proj, capture_output=True, text=True)
        assert r.returncode == 0, f"step{step} failed: {r.stderr}\n{r.stdout}"


def test_runtime_config_overrides(tmp_path, monkeypatch):
    proj = tmp_path
    make_project_tree(proj)
    shim_dir = _install_shims(proj, with_runtime=True, threads=4, pin=True)
    monkeypatch.setenv('PATH', str(shim_dir)+os.pathsep+os.environ['PATH'])
    # Ensure forcefield validation passes in step3
    gmxdata = install_fake_gmxdata(proj)
    monkeypatch.setenv('GMXDATA', str(gmxdata))
    _run_steps(proj)
    calls_file = proj/'process'/'IP_011101'/'run_gmx_simulation'/'mdrun_calls.txt'
    content = calls_file.read_text()
    assert re.search(r'-nt 4(\s|$)', content)
    assert '-pin on' in content


def test_runtime_defaults_when_unset(tmp_path, monkeypatch):
    proj = tmp_path
    make_project_tree(proj)
    shim_dir = _install_shims(proj, with_runtime=False)
    monkeypatch.setenv('PATH', str(shim_dir)+os.pathsep+os.environ['PATH'])
    # Ensure forcefield validation passes in step3
    gmxdata = install_fake_gmxdata(proj)
    monkeypatch.setenv('GMXDATA', str(gmxdata))
    _run_steps(proj)
    calls_file = proj/'process'/'IP_011101'/'run_gmx_simulation'/'mdrun_calls.txt'
    content = calls_file.read_text()
    nts = re.findall(r'-nt (\d+)', content)
    assert nts, 'expected thread flags present'
    assert len(set(nts)) == 1
    assert '-pin on' in content

import pytest

from electrofit.io.ff import forcefield_digest, is_forcefield_path, stage_forcefield


def _make_ff(root, name="test.ff", na_epsilon="3.65846e-01"):
    ff = root / name
    ff.mkdir(parents=True)
    (ff / "forcefield.itp").write_text('#include "ffnonbonded.itp"\n')
    (ff / "ffnonbonded.itp").write_text(f"[ atomtypes ]\nNa 11 22.99 0.0 A 2.43928e-01 {na_epsilon}\n")
    (ff / "tip3p.itp").write_text("; water\n")
    (ff / "ions.itp").write_text("; ions\n")
    return ff


def test_path_detection():
    assert is_forcefield_path("/abs/amber14sb_corrected-Na.ff")
    assert is_forcefield_path("data/amber14sb_corrected-Na.ff")
    assert is_forcefield_path("./local.ff")
    assert not is_forcefield_path("amber14sb.ff")


def test_stage_copies_folder_next_to_topology(tmp_path):
    src = _make_ff(tmp_path / "src")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    assert stage_forcefield(str(src), str(run_dir)) == "test.ff"
    assert forcefield_digest(str(run_dir / "test.ff")) == forcefield_digest(str(src))
    # Staging again with identical content reuses the copy.
    assert stage_forcefield(str(src), str(run_dir)) == "test.ff"


def test_stage_refuses_a_different_existing_copy(tmp_path):
    src = _make_ff(tmp_path / "src")
    run_dir = tmp_path / "run"
    _make_ff(run_dir, na_epsilon="3.65846e-02")
    with pytest.raises(ValueError):
        stage_forcefield(str(src), str(run_dir))


def test_stage_requires_core_files(tmp_path):
    src = _make_ff(tmp_path / "src")
    (src / "ions.itp").unlink()
    with pytest.raises(FileNotFoundError):
        stage_forcefield(str(src), str(tmp_path))


def test_digest_depends_on_content(tmp_path):
    a = _make_ff(tmp_path / "a")
    b = _make_ff(tmp_path / "b", na_epsilon="3.65846e-02")
    assert forcefield_digest(str(a)) != forcefield_digest(str(b))

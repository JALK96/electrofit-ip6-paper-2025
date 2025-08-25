import json
import os
from pathlib import Path
import subprocess

import pytest
from electrofit.viz.hist_manifest import load_hist_manifest, validate_hist_manifest


@pytest.mark.integration
def test_step6_generates_hist_outputs(tmp_path: Path):
    # Skip if matplotlib truly unavailable at import time (binary deps missing in CI container),
    # since production code now requires it hard; this preserves test suite portability.
    try:  # pragma: no cover
        import matplotlib  # noqa: F401
    except Exception:
        pytest.skip("matplotlib not importable in test environment")
    """Run a minimal synthetic project through step6 with histogram+outlier flags
    and assert expected diagnostic artefacts are created.

    This test purposely does not validate image *content* (expensive / brittle),
    only the presence of the files and a few structural JSON fields.
    """
    # Arrange: create minimal project scaffold matching expected layout
    project = tmp_path / "proj"
    (project / "process" / "MOLX" / "extracted_conforms").mkdir(parents=True)
    (project / "process" / "MOLX" / "run_gau_create_gmx_in").mkdir(parents=True)
    # symmetry json (simple one group with two atoms)
    sym_json = project / "process" / "MOLX" / "extracted_conforms" / "equiv_groups.json"
    sym_json.write_text('{"A1": ["A2"]}')

    # minimal config snapshot: project root electrofit.toml
    (project / "electrofit.toml").write_text("""[project]\nmolecule_name = \"MOLX\"\ncharge = 0\natom_type = \"gaff2\"\nadjust_symmetry = true\ncalculate_group_average = true\n""")

    # Provide dummy mol2 and charges directories expected by step6
    ac_dir = project / "process" / "MOLX" / "run_gau_create_gmx_in" / "dummy.acpype"
    ac_dir.mkdir()
    mol2_file = ac_dir / "MOLX_gaff2.mol2"
    # Create trivial mol2 with two atoms and placeholder charges
    mol2_file.write_text("""@<TRIPOS>MOLECULE\nMOLX\n 2 0 0 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n      1 A1    0.0 0.0 0.0 C.3 1 MOLX  0.0000\n      2 A2    0.0 0.0 0.0 C.3 1 MOLX  0.0000\n""")

    # Create fake conformer subdirectories with simple charge lists
    # step6 expects subdirs each with mol2-like extracted charge JSON; we mimic helper behaviour by
    # simulating directory structure used by extract_charges_from_subdirectories: <extracted_conforms>/<CONF>/charges.json
    for i in range(5):
        cdir = project / "process" / "MOLX" / "extracted_conforms" / f"CONF{i}"
        cdir.mkdir()
        # Two atoms, produce slightly varying charges with an outlier
        if i == 4:
            charges = {"A1": 0.8, "A2": -0.8}  # outlier extreme
        else:
            charges = {"A1": 0.1 + 0.01 * i, "A2": -0.1 - 0.01 * i}
        (cdir / "charges.json").write_text(json.dumps(charges))

    # Act: invoke step6 via module (ensuring project env var)
    env = os.environ.copy()
    env["ELECTROFIT_PROJECT_PATH"] = str(project)
    # We call through CLI to exercise argparse & flag wiring
    proc = subprocess.run([
        "python", "-m", "electrofit.cli.app", "step6", "--project", str(project),
        "--plot-histograms", "--remove-outlier", "--hist-combine-groups", "--hist-bins", "10"
    ], env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        stderr = proc.stderr or ""
        if "matplotlib" in stderr or "CXXABI" in stderr:
            pytest.skip("matplotlib runtime import failed in subprocess: " + stderr.splitlines()[-1])
        else:
            raise AssertionError(f"step6 subprocess failed: rc={proc.returncode}\nSTDERR:\n{stderr}\nSTDOUT:\n{proc.stdout}")

    results_dir = project / "process" / "MOLX" / "results"
    # Assert: base histograms (initial + after removal) must exist (hard dependency on matplotlib).
    hist_pdf = results_dir / "hist.pdf"
    hist_no_outlier_pdf = results_dir / "hist_no_outlier.pdf"
    assert hist_pdf.is_file(), "hist.pdf missing"
    assert hist_no_outlier_pdf.is_file(), "hist_no_outlier.pdf missing"
    # Manifest records which optional plots were expected/created
    manifest_path = results_dir / "hist_manifest.json"
    assert manifest_path.is_file(), "hist_manifest.json missing"
    manifest = load_hist_manifest(manifest_path)
    validate_hist_manifest(manifest)
    # Adjusted/group hist optional: if expected=True in manifest then created must be True
    for key in ("adjusted", "group_combined"):
        if key in manifest and manifest[key]["expected"]:
            assert manifest[key]["created"], f"Expected histogram '{key}' not created: reason={manifest[key].get('reason')}"
    # Outlier summary JSON
    summary_file = results_dir / "hist_summary.json"
    assert summary_file.is_file(), "hist_summary.json missing"
    summary = json.loads(summary_file.read_text())
    assert "total_conformers" in summary and summary["total_conformers"] == 5
    assert summary.get("removed_conformers", 0) >= 1, "Expected at least one outlier conformer removed"
    assert "per_atom_outliers" in summary and set(summary["per_atom_outliers"].keys()) == {"A1", "A2"}

from __future__ import annotations

import os
from pathlib import Path
import types
import pytest

from electrofit.domain.charges.process_conformer import (
    ConformerConfig,
    prepare_inputs,
    generate_esp,
    prepare_resp,
    write_symmetry_doc,
)
from electrofit.adapters.results import GaussianResult


@pytest.fixture()
def scratch(tmp_path: Path):
    d = tmp_path / "scratch"
    d.mkdir()
    return d


def _make_cfg(tmp_path: Path) -> ConformerConfig:
    pdb = tmp_path / "mol.pdb"
    pdb.write_text("ATOM\nEND\n")
    return ConformerConfig(
        molecule_name="mol",
        pdb_file=str(pdb),
        net_charge=0,
        residue_name="LIG",
        adjust_sym=False,
        ignore_sym=False,
        protocol="bcc",
    )


def test_prepare_inputs_basic(scratch: Path, monkeypatch):
    cfg = _make_cfg(scratch)
    # monkeypatch gaussian build_input to avoid external complexity
    import electrofit.adapters.gaussian as ga
    def fake_build_input(mol2, name, net_charge, cwd):  # noqa: D401
        # create dummy input marker file
        p = Path(cwd)/f"{name}.gjf"; p.write_text("%chk\n# HF\n")
        return str(p)
    monkeypatch.setattr(ga, "build_input", fake_build_input)
    res = prepare_inputs(cfg, str(scratch))
    assert res.mol2_file.endswith('.mol2')
    assert Path(str(scratch))/res.mol2_file.split('/')[-1]


def test_generate_esp_requires_gesp(scratch: Path):
    cfg = _make_cfg(scratch)
    # GaussianResult with missing gesp triggers FileNotFoundError
    gr = GaussianResult(name="mol", input_file="in.gjf", gesp_file=None, log_file=None, used_cache=False, wall_time_s=0.0)
    with pytest.raises(FileNotFoundError):
        generate_esp(cfg, gr, str(scratch))


def test_prepare_resp_missing_inputs_errors(scratch: Path, monkeypatch):
    cfg = _make_cfg(scratch)
    # Force absence of generated RESP inputs by not creating them
    import electrofit.domain.charges.process_conformer as pc
    def fake_generate_inputs(name, cwd):  # noqa: D401
        return ("RESP1", "RESP2")
    monkeypatch.setattr('electrofit.domain.charges.process_conformer.resp_adapter', types.SimpleNamespace(
        prepare_ac=lambda *a, **k: None,
        generate_inputs=fake_generate_inputs,
        apply_symmetry=lambda *a, **k: None,
    ))
    # Intentionally do not create files -> should raise
    with pytest.raises(FileNotFoundError):
        prepare_resp(cfg, str(scratch))


def test_write_symmetry_doc_creates_file(scratch: Path, monkeypatch):
    # create dummy resp1
    resp1 = scratch / 'ANTECHAMBER_RESP1.IN'
    resp1.write_text('X')
    # monkeypatch run_python to just write output
    import electrofit.domain.charges.process_conformer as pc
    def fake_run_python(func, *args, **kwargs):  # noqa: D401
        # args[1] is output file name per write_symmetry_doc convention
        out = Path(kwargs.get('cwd', '.'))/args[1]
        out.write_text('sym')
    monkeypatch.setattr(pc, 'run_python', fake_run_python)
    path = write_symmetry_doc(str(resp1), adjust=False, scratch_dir=str(scratch))
    assert path is not None and os.path.isfile(path)

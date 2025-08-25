from pathlib import Path
from electrofit.io.files import replace_charge_in_ac_file

def test_replace_charge_in_ac_file(tmp_path: Path):
    f = tmp_path/"x.ac"
    f.write_text("FOO\nCHARGE     0.00 ( 0 )\nBAR\n")
    replace_charge_in_ac_file(str(f), new_charge_float=-1.0, cwd=None)
    assert "CHARGE     -1.00 ( -1 )" in f.read_text()
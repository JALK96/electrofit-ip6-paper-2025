from pathlib import Path
from electrofit.io.files import modify_gaussian_input

def test_modify_gaussian_input(tmp_path: Path):
    f = tmp_path/"in.gcrt"
    f.write_text("%mem=8GB\n#p b3lyp/6-31g* opt scf=tight\n\nTitle\n\n0 1\n")
    modify_gaussian_input(str(f))
    txt = f.read_text()
    assert " opt" not in txt
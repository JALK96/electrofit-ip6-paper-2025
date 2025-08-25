import pathlib, shutil
from electrofit.io.ff import remove_defaults_section_lines

FIX = pathlib.Path(__file__).parent / 'fixtures'

def _read(p):
    return (FIX/p).read_text().splitlines()

def test_remove_defaults_section_present(tmp_path):
    src = FIX/'min_top_with_defaults.top'
    dst = tmp_path/'copy.top'
    shutil.copy2(src, dst)
    remove_defaults_section_lines(str(dst))
    txt = dst.read_text()
    assert '[ defaults ]' not in txt
    # system section preserved
    assert '[ system ]' in txt


def test_remove_defaults_section_absent_idempotent(tmp_path):
    src = FIX/'min_top_without_defaults.top'
    dst = tmp_path/'copy.top'
    shutil.copy2(src, dst)
    before = dst.read_text()
    remove_defaults_section_lines(str(dst))
    after = dst.read_text()
    assert before == after

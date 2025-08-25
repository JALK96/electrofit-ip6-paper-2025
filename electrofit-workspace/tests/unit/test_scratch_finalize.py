import os
from pathlib import Path
import sys

# add package path
sys.path.append(str(Path(__file__).resolve().parents[2] / 'packages' / 'electrofit' / 'src'))

from electrofit.scratch.manager import finalize_scratch_directory  # type: ignore
from electrofit.cli.safe_run import ensure_finalized  # type: ignore


def test_finalize_scratch_file_overwrite(tmp_path: Path):
    orig = tmp_path / 'orig'
    orig.mkdir()
    (orig / 'input.txt').write_text('ORIGINAL')

    scratch = tmp_path / 'scratch'
    scratch.mkdir()
    # copy original file (simulate setup)
    (scratch / 'input.txt').write_text('MODIFIED')

    finalize_scratch_directory(
        original_dir=str(orig),
        scratch_dir=str(scratch),
        input_files=['input.txt'],
        output_files=None,
        overwrite=True,
        remove_parent_if_empty=True,
    )

    # original dir now should have input.txt replaced and renamed backup
    replaced = (orig / 'input.txt').read_text()
    assert replaced == 'MODIFIED'
    backups = list(orig.glob('input.input_file.txt'))
    assert backups, 'Expected renamed backup of original input file'
    assert not scratch.exists(), 'Scratch directory should be removed'


def test_ensure_finalized_context(tmp_path: Path):
    orig = tmp_path / 'orig2'
    orig.mkdir()
    (orig / 'a.dat').write_text('A')
    scratch = tmp_path / 'scratch2'
    scratch.mkdir()
    (scratch / 'a.dat').write_text('A')

    with ensure_finalized(
        original_dir=str(orig),
        scratch_dir=str(scratch),
        input_files=['a.dat'],
        output_files=None,
        overwrite=True,
        remove_parent_if_empty=True,
    ):
        # modify scratch copy
        (scratch / 'a.dat').write_text('CHANGED')

    assert (orig / 'a.dat').read_text() == 'CHANGED'
    assert list(orig.glob('a.input_file.dat')), 'Backup of original should exist'
    assert not scratch.exists()

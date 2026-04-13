from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electrofit_analysis.cli.h_bonds.make_pp_matrix_ip6 import build_matrix, parse_hbond_log
from electrofit_analysis.structure.util.hbond_io import parse_xpm


def _write_minimal_xpm(path: Path, rows: list[str], *, one_symbol: str = "o") -> None:
    width = len(rows[0])
    height = len(rows)
    for row in rows:
        if len(row) != width:
            raise ValueError("All synthetic XPM rows must have equal width")
    path.write_text(
        "\n".join(
            [
                "/* title: synthetic */",
                "static char *gromacs_xpm[] = {",
                f"\"{width} {height} 2 1\",",
                "\"   c #FFFFFF\",",
                f"\"{one_symbol}  c #FF0000\",",
                *(f"\"{row}\"," for row in rows),
                "};",
            ]
        )
        + "\n"
    )


def _write_minimal_log(path: Path, pairs: list[tuple[str, str]]) -> None:
    path.write_text("\n".join(f"{d} - {a}" for d, a in pairs) + "\n")


def test_parse_xpm_align_rows_to_log_reverses_image_row_order(tmp_path: Path) -> None:
    xpm = tmp_path / "hb.xpm"
    _write_minimal_xpm(
        xpm,
        rows=[
            "o  ",  # image row 0
            " o ",  # image row 1
            "  o",  # image row 2
        ],
    )

    raw, _ = parse_xpm(xpm, align_rows_to_log=False)
    aligned, _ = parse_xpm(xpm, align_rows_to_log=True)

    assert np.array_equal(aligned, raw[::-1, :])
    assert np.array_equal(raw[0], np.array([1, 0, 0], dtype=np.uint8))
    assert np.array_equal(aligned[0], np.array([0, 0, 1], dtype=np.uint8))


def test_ip101101_regression_row_alignment_lowers_implausible_p4_to_p2() -> None:
    test_file = Path(__file__).resolve()
    candidate_roots = [test_file.parents[4], test_file.parents[3]]
    repo_root = next((root for root in candidate_roots if (root / "projects").is_dir()), None)
    if repo_root is None:
        pytest.skip("Cannot locate repository root with projects/ directory")

    hb_dir = (
        repo_root
        / "projects/ellas_results/ellas_meta-project/K/100mM/process/IP_101101/analyze_final_sim/h_bonds"
    )
    xpm_path = hb_dir / "intra_hb_matrix.xpm"
    log_path = hb_dir / "intra_hb.log"

    if not (xpm_path.is_file() and log_path.is_file()):
        pytest.skip("IP_101101 regression data not present in this checkout")

    pairs = parse_hbond_log(log_path)
    no_align, _ = parse_xpm(xpm_path, align_rows_to_log=False)
    aligned, _ = parse_xpm(xpm_path, align_rows_to_log=True)

    m_old = build_matrix(no_align, pairs, mode="union", oxygen_mode="all")
    m_new = build_matrix(aligned, pairs, mode="union", oxygen_mode="all")

    p4_to_p2_old = float(m_old[3, 1])
    p4_to_p2_new = float(m_new[3, 1])

    assert p4_to_p2_old > 0.2
    assert p4_to_p2_new < 0.01

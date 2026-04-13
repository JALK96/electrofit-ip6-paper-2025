from __future__ import annotations

import numpy as np

from electrofit_analysis.cli.h_bonds.make_pp_matrix_ip6 import build_matrix


def test_terminal_mode_excludes_bridging_oxygen_rows() -> None:
    # rows: (donor O token, acceptor O token)
    pairs = [
        ("O", "O2"),      # bridging donor (excluded in terminal mode)
        ("O6", "O12"),    # terminal donor+acceptor (kept)
        ("O10", "O15"),   # another terminal pair (kept)
    ]
    # 3 bonds x 5 frames
    # row0/row1 both map to P1->P3; row0 is bridging, row1 is terminal
    xpm = np.array(
        [
            [1, 1, 0, 0, 0],  # mean 0.4 (bridging)
            [0, 1, 0, 1, 0],  # mean 0.4 (terminal)
            [1, 0, 1, 0, 1],  # mean 0.6 (terminal)
        ],
        dtype=np.uint8,
    )

    m_all = build_matrix(xpm, pairs, mode="union", oxygen_mode="all")
    m_terminal = build_matrix(xpm, pairs, mode="union", oxygen_mode="terminal")

    # P1->P3 (index [0,2]): union over row0+row1 => [1,1,0,1,0] => 3/5 = 0.6
    assert np.isclose(m_all[0, 2], 0.6)
    # terminal excludes row0 => row1 only => 2/5 = 0.4
    assert np.isclose(m_terminal[0, 2], 0.4)

    # P2->P4 from row2 stays unchanged between modes
    assert np.isclose(m_all[1, 3], 0.6)
    assert np.isclose(m_terminal[1, 3], 0.6)

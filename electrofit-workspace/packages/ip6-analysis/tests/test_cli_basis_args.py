from __future__ import annotations

import sys
import types

from electrofit_analysis.cli.app import build_parser


def test_coordination_parser_default_basis_is_oxygen() -> None:
    parser = build_parser()
    args = parser.parse_args(["coordination", "-p", "/tmp/project"])
    assert args.coord_distance_basis == "oxygen"


def test_summarize_parser_accepts_basis_and_bool_override() -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "summarize-nap",
            "-p",
            "/tmp/project",
            "-o",
            "/tmp/out",
            "--coord-distance-basis",
            "phosphorus",
            "--coord-bool-name",
            "IonP_coordination_bool_phosphorus_0p3200nm.npy",
        ]
    )
    assert args.coord_distance_basis == "phosphorus"
    assert args.coord_bool_name == "IonP_coordination_bool_phosphorus_0p3200nm.npy"


def test_coordination_cli_forwards_basis(monkeypatch) -> None:
    parser = build_parser()
    captured: dict[str, object] = {}

    fake_mod = types.ModuleType("electrofit_analysis.cli.coordination.Na_IP6_coordination")

    def _fake_main(**kwargs):
        captured.update(kwargs)

    fake_mod.main = _fake_main  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "electrofit_analysis.cli.coordination.Na_IP6_coordination", fake_mod)

    args = parser.parse_args(
        [
            "coordination",
            "-p",
            "/tmp/project",
            "--coord-distance-basis",
            "phosphorus",
            "--coord-cutoff-nm",
            "0.32",
            "--skip-rdf",
        ]
    )
    args.func(args)
    assert captured["coord_distance_basis"] == "phosphorus"


def test_summarize_cli_forwards_basis_and_bool_override(monkeypatch) -> None:
    parser = build_parser()
    captured: dict[str, object] = {}

    fake_mod = types.ModuleType(
        "electrofit_analysis.cli.coordination.summerize_nap_dist_count_ip6"
    )

    def _fake_main(**kwargs):
        captured.update(kwargs)

    fake_mod.main = _fake_main  # type: ignore[attr-defined]
    monkeypatch.setitem(
        sys.modules,
        "electrofit_analysis.cli.coordination.summerize_nap_dist_count_ip6",
        fake_mod,
    )

    args = parser.parse_args(
        [
            "summarize-nap",
            "-p",
            "/tmp/project",
            "-o",
            "/tmp/out",
            "--coordination",
            "--coord-distance-basis",
            "phosphorus",
            "--coord-bool-name",
            "IonP_coordination_bool_phosphorus_0p3200nm.npy",
        ]
    )
    args.func(args)
    assert captured["coord_distance_basis"] == "phosphorus"
    assert captured["coord_bool_name"] == "IonP_coordination_bool_phosphorus_0p3200nm.npy"

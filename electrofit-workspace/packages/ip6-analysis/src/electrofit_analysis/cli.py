"""Legacy CLI shim.

This module remains for backward compatibility with callers importing
`electrofit_analysis.cli:main`. The authoritative CLI implementation is
`electrofit_analysis.cli.app:main`.
"""

from __future__ import annotations

from .cli.app import main as _app_main


def main(argv: list[str] | None = None) -> None:
    _app_main(argv)


if __name__ == "__main__":
    main()

"""Pipeline Step 0: Initialize process workspace.

Copies input molecule subdirectories from ``data/input/*`` into ``process/*`` as
working copies. Pure orchestration: delegates filesystem copying to
``electrofit.io.files.copy_and_rename_folders`` and logging to infra layer.

CLI equivalent (legacy): ``electrofit step0 --project <path>``
"""
from __future__ import annotations

from pathlib import Path
import os
import argparse
import logging

from electrofit.io.files import copy_and_rename_folders
from electrofit.infra.logging import setup_logging, log_run_header

__all__ = ["main"]

def main():  # pragma: no cover (CLI orchestration)
    parser = argparse.ArgumentParser(description="Step0: initialize process/ workspace by copying data/input")
    parser.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    parser.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    parser.add_argument("--molecule", help="If set, only copy this molecule subfolder (data/input/<name>)", default=None)
    args = parser.parse_args()

    project_path = Path(args.project).resolve()
    src = project_path / "data" / "input"
    dst = project_path / "process"
    dst.mkdir(exist_ok=True)
    setup_logging(str(project_path / "step.log"), also_console=args.log_console)
    log_run_header("step0")
    logging.info("[step0] Source: %s", src)
    logging.info("[step0] Destination: %s", dst)
    if not src.exists():
        raise FileNotFoundError(f"Source directory '{src}' does not exist.")
    if args.molecule:
        # Konkurrenzsicher nur das eine Molekül kopieren (kein späteres Pruning -> vermeidet Race Conditions bei parallelen Läufen).
        mol = args.molecule
        sub = src / mol
        if not sub.is_dir():
            raise SystemExit(f"[step0][abort] requested molecule '{mol}' not found under {src}")
        nested = dst / mol / "run_gau_create_gmx_in"
        lock_path = dst / f".{mol}.step0.lock"
        # Einfacher Lock um gleichzeitiges Schreiben in dasselbe Molekül zu verhindern.
        import errno
        import time
        fd = None
        try:
            try:
                fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, b"lock\n")
            except OSError as e:
                if e.errno == errno.EEXIST:
                    # Kurzes Warten falls anderer Prozess gerade kopiert.
                    for _ in range(50):  # bis ~5s
                        if not lock_path.exists():
                            break
                        time.sleep(0.1)
                    if lock_path.exists():
                        raise SystemExit(f"[step0][abort] lock for molecule '{mol}' still present: {lock_path}")
                else:
                    raise
            # Falls schon vorhanden und nicht leer: wir gehen davon aus, dass ein anderer Lauf abgeschlossen hat.
            if nested.is_dir() and any(nested.iterdir()):
                logging.info("[step0] Molecule '%s' already prepared (reusing).", mol)
            else:
                from shutil import copy2, copytree
                nested.mkdir(parents=True, exist_ok=True)
                for item in sub.iterdir():
                    dest_path = nested / item.name
                    if item.is_file():
                        copy2(item, dest_path)
                    elif item.is_dir():
                        # dirs_exist_ok ab Python 3.8
                        copytree(item, dest_path, dirs_exist_ok=True)
                logging.info("[step0] Copied molecule '%s' -> %s", mol, nested)
        finally:
            if fd is not None:
                try:
                    os.close(fd)
                except Exception:
                    pass
            # Lock entfernen (ignorieren, falls schon weg)
            try:
                if lock_path.exists():
                    os.unlink(lock_path)
            except Exception:
                pass
    else:
        copy_and_rename_folders(source=str(src), destination=str(dst))
    logging.info("[step0] Copy complete.")

if __name__ == "__main__":  # pragma: no cover
    main()

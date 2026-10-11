import hashlib
import os
import shutil
import logging
from typing import Iterable, Optional, Tuple





def replace_posres_in_file(file_path):
    """Replace legacy token POSRES_LIG with POSRES in topology if present.

    Idempotent: on subsequent calls when already replaced it logs at INFO and leaves file unchanged.
    """
    try:
        with open(file_path, "r") as fh:
            content = fh.read()
        if "POSRES_LIG" not in content:
            logging.info("POSRES replacement already applied (no 'POSRES_LIG' in %s)", file_path)
            return
        updated_content = content.replace("POSRES_LIG", "POSRES")
        with open(file_path, "w") as fh:
            fh.write(updated_content)
        logging.info("Replaced 'POSRES_LIG' with 'POSRES' in %s", file_path)
    except FileNotFoundError:
        logging.error("The file %s does not exist.", file_path)
    except Exception as e:
        logging.error("An error occurred: %s", e)

def include_tip3p(file_path, forcefield="amber14sb.ff"):
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if '[ system ]' is present
        if "[ system ]" in content:
            # Avoid duplicate insertion
            include_line = f'#include "{forcefield}/tip3p.itp"'
            if include_line in content:
                logging.info(f"tip3p.itp already included in {file_path}")
                return

            # Insert the tip3p topology block before [ system ]
            updated_content = content.replace(
                "[ system ]",
                f"""
; Include water topology
#include "{forcefield}/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

[ system ]
""",
            )
        else:
            logging.error(
                "No '[ system ]' section found in the file to include tip3p topology (tip3p.itp)."
            )
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(
            f"Successfully included tip3p.itp and position restraints in {file_path}"
        )

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def include_ions(file_path, forcefield="amber14sb.ff"):
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if '[ system ]' is present
        if "[ system ]" in content:
            # Avoid duplicate insertion
            include_line = f'#include "{forcefield}/ions.itp"'
            if include_line in content:
                logging.info(f"ions.itp already included in {file_path}")
                return

            # Insert the ions topology block before [ system ]
            updated_content = content.replace(
                "[ system ]",
                f"""
; Include ion topology
#include "{forcefield}/ions.itp"

[ system ]
""",
            )
        else:
            logging.error(
                "No '[ system ]' section found in the file to include ions.itp."
            )
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(f"Successfully included ions.itp in {file_path}")

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def remove_defaults_section_lines(file_path):
    """Remove a legacy '[ defaults ]' block (and next two lines) if present.

    Idempotent: if no such block exists nothing is written.
    """
    try:
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")
        with open(file_path, "r") as fh:
            lines = fh.readlines()
        new_lines = []
        skip_next = 0
        found = False
        for index, line in enumerate(lines):
            if line.strip() == "[ defaults ]" and skip_next == 0:
                logging.info("Removing obsolete '[ defaults ]' section in %s", file_path)
                skip_next = 3  # remove this + next two lines
                found = True
                continue
            if skip_next > 0:
                skip_next -= 1
                continue
            new_lines.append(line)
        if not found:
            logging.debug("No '[ defaults ]' section present in %s (already clean)", file_path)
            return
        with open(file_path, "w") as fh:
            fh.writelines(new_lines)
        logging.info("Removed '[ defaults ]' section from %s", file_path)
    except FileNotFoundError as e:
        logging.error(e)
        raise
    except Exception as e:
        logging.error("Error while processing '[ defaults ]' section in %s: %s", file_path, e)
        raise

def include_ff(file_path, forcefield="amber14sb.ff"):
    try:
        # Open the file and read all lines
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Define the lines to insert
        include_comment = "; Include forcefield\n"
        include_line = f'#include "{forcefield}/forcefield.itp"\n'
        include_target = f"{forcefield}/forcefield.itp"

        def _extract_target(line: str) -> str | None:
            if "forcefield.itp" not in line:
                return None
            parts = line.replace("#include", "").strip().split()
            if not parts:
                return None
            # last token should be the path, possibly quoted
            return parts[-1].strip('"')

        # If a forcefield include exists but points elsewhere, replace it
        existing_idx = next((i for i, line in enumerate(lines) if "forcefield.itp" in line), None)
        if existing_idx is not None:
            current = _extract_target(lines[existing_idx])
            if current == include_target:
                logging.info(
                    f'"{forcefield}/forcefield.itp" is already included in {file_path}.'
                )
                return
            lines[existing_idx] = include_line
            if existing_idx == 0 or "Include forcefield" not in lines[existing_idx - 1]:
                lines.insert(existing_idx, include_comment)
            logging.info(
                f"Replaced existing forcefield include with '{forcefield}/forcefield.itp' in {file_path}."
            )
        else:
            # Insert the include lines as the second and third lines (or at start if empty)
            if lines:
                lines.insert(1, include_comment)
                lines.insert(2, include_line)
            else:
                lines.append(include_comment)
                lines.append(include_line)

        with open(file_path, "w") as file:
            file.writelines(lines)

        logging.info(
            f'Successfully included "{forcefield}/forcefield.itp" in {file_path}.'
        )

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


# ---- Forcefield discovery / validation utilities ---------------------------------

def _candidate_gmxdata_paths() -> Iterable[str]:
    """Yield plausible GMXDATA roots in descending priority.

    We prefer an explicit GMXDATA env, then common conda layout, then /usr/local.
    """
    env = os.environ.get("GMXDATA")
    if env:
        yield env
    # Conda layout: $CONDA_PREFIX/share/gromacs
    cprefix = os.environ.get("CONDA_PREFIX")
    if cprefix:
        share_candidate = os.path.join(cprefix, "share", "gromacs")
        yield share_candidate
    # Fallback typical system prefix
    yield "/usr/local/gromacs/share/gromacs"


def locate_gmx_topdir() -> Optional[str]:
    """Return the resolved path to the GROMACS 'top' directory or None if not found."""
    for root in _candidate_gmxdata_paths():
        top = os.path.join(root, "top")
        if os.path.isdir(top):
            return top
    return None


def forcefield_exists(forcefield: str) -> Tuple[bool, Optional[str]]:
    """Check whether a forcefield directory (<name>.ff) exists under an accessible top/.

    Returns (exists, top_dir_used). If top_dir could not be resolved returns (False, None).
    """
    top = locate_gmx_topdir()
    if not top:
        return False, None
    path = os.path.join(top, forcefield)
    return os.path.isdir(path), top


def validate_forcefield(forcefield: str, explicit: bool) -> None:
    """Validate presence of forcefield directory, emit clear diagnostics.

    Parameters
    ----------
    forcefield : str
        Name of the forcefield directory (e.g. 'amber14sb.ff').
    explicit : bool
        Whether user explicitly specified it (controls warning wording).
    """
    exists, top = forcefield_exists(forcefield)
    if not exists:
        if top is None:
            logging.warning(
                "[ff] Could not resolve GMXDATA/top path (GMXDATA env unset & no standard location). "
                "Cannot validate forcefield '%s'. If grompp fails, ensure GROMACS environment sourced.",
                forcefield,
            )
            return
        msg = (
            f"[ff] Forcefield directory '{forcefield}' not found under {top}. "
            "Available entries: "
            + ", ".join(sorted(d for d in os.listdir(top) if d.endswith('.ff')))  # type: ignore[arg-type]
        )
        # escalate to error only if user explicitly requested a non-existent ff
        if explicit:
            raise FileNotFoundError(msg)
        logging.warning(msg + " (using default fallback anyway)")
    else:
        logging.info("[ff] Validated forcefield '%s' in %s", forcefield, top)


# ---- Run-local staging of path-provided forcefields ------------------------------

_REQUIRED_FF_FILES = ("forcefield.itp", "tip3p.itp", "ions.itp")


def is_forcefield_path(forcefield: str) -> bool:
    """True if the configured forcefield is a path rather than a bare '<name>.ff'."""
    return isinstance(forcefield, str) and (
        os.path.isabs(forcefield) or os.sep in forcefield or forcefield.startswith(".")
    )


def forcefield_digest(ff_dir: str) -> str:
    """SHA-256 over the relative path and content hash of every file in a '.ff' directory."""
    digest = hashlib.sha256()
    for root, dirs, files in os.walk(ff_dir):
        dirs.sort()
        for name in sorted(files):
            path = os.path.join(root, name)
            digest.update(os.path.relpath(path, ff_dir).encode() + b"\0")
            with open(path, "rb") as handle:
                digest.update(hashlib.sha256(handle.read()).digest())
    return digest.hexdigest()


def stage_forcefield(ff_path: str, dest_dir: str) -> str:
    """Copy a '.ff' directory next to the topology in dest_dir and return its name.

    grompp resolves ``#include "<name>.ff/..."`` in the topology's directory before
    the GROMACS installation, so the staged copy defines the parameters even if a
    folder of the same name is installed with GROMACS. The copy travels with the
    run inputs into scratch. An existing copy in dest_dir is reused only if its
    content is identical to the source.
    """
    src = os.path.expanduser(os.path.expandvars(ff_path))
    if not os.path.isdir(src):
        raise FileNotFoundError(f"[ff] Source forcefield path not found: {src}")
    ff_name = os.path.basename(src.rstrip(os.sep))
    if not ff_name.endswith(".ff"):
        raise ValueError(f"[ff] Expected a '.ff' directory, got: {ff_name}")
    missing = [p for p in _REQUIRED_FF_FILES if not os.path.isfile(os.path.join(src, p))]
    if missing:
        raise FileNotFoundError(f"[ff] Missing files in {src}: {', '.join(missing)}")

    digest = forcefield_digest(src)
    dest = os.path.join(dest_dir, ff_name)
    if os.path.isdir(dest):
        if forcefield_digest(dest) != digest:
            raise ValueError(
                f"[ff] {dest} differs from {src}; remove it or use a new run directory "
                "so one run does not mix force-field versions"
            )
    else:
        shutil.copytree(src, dest)
    logging.info("[ff] Staged '%s' from %s into %s (sha256 %s)", ff_name, src, dest_dir, digest)
    return ff_name

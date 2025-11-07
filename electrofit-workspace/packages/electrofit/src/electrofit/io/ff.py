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

        # Check if the include line already exists to prevent duplication
        if any(forcefield in line and "forcefield.itp" in line for line in lines):
            logging.info(
                f'"{forcefield}/forcefield.itp" is already included in {file_path}.'
            )
            return

        # Insert the include lines as the second and third lines
        if lines:
            lines.insert(1, include_comment)
            lines.insert(2, include_line)
        else:
            # If the file is empty, add the include lines at the beginning
            lines.append(include_comment)
            lines.append(include_line)

        # Write the updated content back to the file
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


# ---- Installation of path-provided forcefields into GMX top ---------------------

def _conda_top_dir() -> Optional[str]:
    """Return $CONDA_PREFIX/share/gromacs/top if it exists; else fall back to locate_gmx_topdir()."""
    cprefix = os.environ.get("CONDA_PREFIX")
    if cprefix:
        top = os.path.join(cprefix, "share", "gromacs", "top")
        if os.path.isdir(top):
            return top
    return locate_gmx_topdir()


def ensure_forcefield_installed(ff_path: str) -> str:
    """Ensure a '.ff' directory from a path is available under the active GMX 'top' directory.

    - ff_path may be absolute or relative ("~/...", "./...", or project-local path).
    - Verifies required files exist in the source folder.
    - Copies into <top>/<ff_name> if missing. Returns the short name (e.g. 'amber14sb.ff').

    Raises FileNotFoundError, ValueError, or PermissionError with helpful messages.
    """
    src = os.path.expanduser(os.path.expandvars(ff_path))
    if not os.path.isdir(src):
        raise FileNotFoundError(f"[ff] Source forcefield path not found: {src}")
    ff_name = os.path.basename(src.rstrip(os.sep))
    if not ff_name.endswith(".ff"):
        raise ValueError(f"[ff] Expected a '.ff' directory, got: {ff_name}")

    # Basic sanity: ensure core files exist
    required = ("forcefield.itp", "tip3p.itp", "ions.itp")
    missing = [p for p in required if not os.path.isfile(os.path.join(src, p))]
    if missing:
        raise FileNotFoundError(f"[ff] Missing files in {src}: {', '.join(missing)}")

    top = _conda_top_dir()
    if not top:
        raise FileNotFoundError(
            "[ff] Could not resolve GROMACS 'top' directory. Ensure CONDA_PREFIX is set or GMXDATA configured."
        )
    dest = os.path.join(top, ff_name)
    if os.path.isdir(dest):
        logging.info("[ff] '%s' already available in %s", ff_name, top)
        return ff_name
    try:
        shutil.copytree(src, dest)
        logging.info("[ff] Copied '%s' -> %s", src, dest)
    except Exception as e:
        raise PermissionError(f"[ff] Failed to install '{ff_name}' into '{top}': {e}") from e
    return ff_name

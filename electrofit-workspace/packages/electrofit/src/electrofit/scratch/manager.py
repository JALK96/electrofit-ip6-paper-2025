import filecmp
import logging
import os
import shutil
import hashlib

from electrofit.infra.logging import setup_logging


def setup_scratch_directory(input_files, base_scratch_dir="/scratch/johannal96/tmp"):
    """
    Sets up a scratch directory and copies essential input files.

    Parameters:
    - input_files (list): List of input file and directory names to copy into the scratch directory.
    - base_scratch_dir (str): Base directory for scratch space.

    Returns:
    - scratch_dir (str): Path to the created scratch directory.
    - original_dir (str): Path to the original working directory.
    """
    base_scratch_dir = base_scratch_dir or "/tmp/electrofit_scratch"
    # Expand env vars and user (~); allow ${USER} etc.
    base_scratch_dir = os.path.expanduser(os.path.expandvars(base_scratch_dir))
    fullpath = os.getcwd()
    calcdir = os.path.basename(fullpath)
    parent_dir = os.path.dirname(fullpath)
    parent_folder_name = os.path.basename(parent_dir)
    scratch_dir = os.path.join(base_scratch_dir, parent_folder_name, calcdir)

    # Initialize logging with log file in original_dir
    log_file_path = os.path.join(fullpath, "process.log")
    suppress = any(h for h in logging.getLogger().handlers if isinstance(h, logging.FileHandler))
    setup_logging(log_file_path, suppress_initial_message=suppress)

    # Create scratch directory
    os.makedirs(scratch_dir, exist_ok=True)
    logging.info(f"Created scratch directory: {scratch_dir}")

    # Copy input files
    for file in input_files:
        src = os.path.join(fullpath, file)
        dst = os.path.join(scratch_dir, file)
        if os.path.exists(src):
            if os.path.isfile(src):
                # Avoid SameFileError if source already in scratch (e.g., test uses scratch under original dir)
                if os.path.abspath(src) == os.path.abspath(dst):
                    logging.debug(f"Skip copying '{file}' – source and destination identical.")
                else:
                    shutil.copy2(src, dst)
                    logging.info(f"Copied file '{file}' to scratch directory.")
            elif os.path.isdir(src):
                # If destination dir already exists from a previous aborted run, refresh it
                if os.path.isdir(dst):
                    try:
                        shutil.rmtree(dst)
                        logging.info(
                            f"Removed pre-existing scratch subdirectory '{file}' (stale)."
                        )
                    except Exception as e:
                        logging.warning(
                            f"Could not remove existing scratch subdirectory '{file}': {e}. Attempting to proceed."
                        )
                try:
                    shutil.copytree(src, dst)
                    logging.info(
                        f"Copied directory '{file}' to scratch directory."
                    )
                except FileExistsError:
                    # Fallback: Python <3.8 without dirs_exist_ok or race condition
                    logging.warning(
                        f"Directory '{file}' already existed after removal attempt; leaving as-is."
                    )
            else:
                logging.warning(
                    f"Input item '{file}' is neither a file nor a directory. Skipping copy."
                )
        else:
            logging.warning(
                f"Input file or directory '{file}' not found. Skipping copy."
            )

    # Log the contents of the scratch directory
    scratch_contents = os.listdir(scratch_dir)
    logging.info(f"Scratch directory contents after setup: {scratch_contents}")

    return scratch_dir, fullpath



def directories_differ(src, dst):
    """
    Recursively check if two directories differ.

    Returns True if:
      - There are files or subdirectories present in one but not the other,
      - There are files with different contents,
      - Or any comparison issues are detected.
    Otherwise, returns False.
    """
    dcmp = filecmp.dircmp(src, dst)
    if dcmp.left_only or dcmp.right_only or dcmp.diff_files or dcmp.funny_files:
        return True
    # Recurse into each subdirectory present in both directories
    for subdir in dcmp.subdirs:
        new_src = os.path.join(src, subdir)
        new_dst = os.path.join(dst, subdir)
        if directories_differ(new_src, new_dst):
            return True
    return False


# Helper: Remove directory if empty, for scratch cleanup
def _rmdir_if_empty(path: str) -> bool:
    """Try to remove a directory if it is empty.

    Returns True if removed, False otherwise.
    """
    try:
        # os.rmdir only removes empty dirs; raises OSError otherwise
        os.rmdir(path)
        logging.info(f"Removed empty parent scratch directory: {path}")
        return True
    except OSError:
        # Not empty or not removable; that's fine
        logging.debug(f"Parent scratch directory not empty or not removable: {path}")
        return False


def finalize_scratch_directory(
    original_dir,
    scratch_dir,
    input_files,
    output_files=None,
    overwrite=True,
    remove_parent_if_empty=True,
    reason: str | None = None,
):
    """Finalize a scratch directory by copying back changed inputs & desired outputs then cleaning up.
    Recursively checks input files and directories for changes and copies updated content
    back to the original directory. Also processes output files/directories from scratch_dir.
    If an item appears in both input_files and output_files, it will only be processed as an input.
    Finally, it cleans up the scratch directory.

    Parameters:
    - original_dir (str): Path to the original directory.
    - scratch_dir (str): Path to the scratch directory.
    - input_files (list): List of input file and directory names to check for changes.
                          For changed items, the original is renamed and the modified copy is copied over.
    - output_files (list, optional): Specific list of output files/directories to copy back.
                                     If None, all items in scratch_dir excluding input_files are processed.
    - overwrite (bool, optional):
          If True, updated input items will be copied back (after renaming the original).
          If False, even if differences are detected the original content will be preserved.
    - reason (str | None, optional): Diagnostic tag for caller context (e.g. "defer-finalize").
    """
    print(f"[scratch-debug] enter finalize original={original_dir} scratch={scratch_dir}", flush=True)
    if reason:
        logging.debug(f"Finalize scratch directory reason={reason}")
    # Short-circuit if scratch dir vanished earlier (already finalized)
    if not os.path.isdir(scratch_dir):
        logging.debug(f"Scratch directory already removed: {scratch_dir}")
        print(f"[scratch-debug] skip finalize (missing) {scratch_dir}", flush=True)
        return

    # --- Process Input Files/Directories ---
    for item in input_files:
        src = os.path.join(scratch_dir, item)
        dst = os.path.join(original_dir, item)

        if not os.path.exists(src):
            logging.warning(f"Input item '{item}' does not exist in scratch directory.")
            continue

        if os.path.isfile(src):
            if os.path.exists(dst):
                if not filecmp.cmp(src, dst, shallow=False):
                    if overwrite:
                        base, ext = os.path.splitext(dst)
                        renamed_dst = f"{base}.input_file{ext}"
                        os.rename(dst, renamed_dst)
                        logging.info(
                            f"Original input file '{item}' renamed to '{os.path.basename(renamed_dst)}'."
                        )
                        shutil.copy2(src, dst)
                        logging.info(
                            f"Modified input file '{item}' copied back to original directory."
                        )
                    else:
                        logging.info(
                            f"Input file '{item}' has changes but overwrite is disabled. No action taken."
                        )
                else:
                    logging.info(f"Input file '{item}' unchanged. No action taken.")
            else:
                shutil.copy2(src, dst)
                logging.info(
                    f"Input file '{item}' does not exist in original directory. Copied from scratch."
                )
        elif os.path.isdir(src):
            if os.path.exists(dst):
                if directories_differ(src, dst):
                    if overwrite:
                        renamed_dst = f"{dst}.input_file"
                        os.rename(dst, renamed_dst)
                        logging.info(
                            f"Original input directory '{item}' renamed to '{os.path.basename(renamed_dst)}'."
                        )
                        shutil.copytree(src, dst)
                        logging.info(
                            f"Modified input directory '{item}' copied back to original directory."
                        )
                    else:
                        logging.info(
                            f"Input directory '{item}' has changes but overwrite is disabled. No action taken."
                        )
                else:
                    logging.info(
                        f"Input directory '{item}' unchanged. No action taken."
                    )
            else:
                shutil.copytree(src, dst)
                logging.info(
                    f"Input directory '{item}' does not exist in original directory. Copied from scratch."
                )
        else:
            logging.warning(
                f"Input item '{item}' is neither a file nor a directory. Skipping."
            )

    # --- Process Output Files/Directories ---
    # If output_files is None, compute all items excluding input_files.
    if output_files is None:
        output_files = [f for f in os.listdir(scratch_dir) if f not in input_files]
    else:
        # Filter out any items that were already processed as input_files.
        output_files = [f for f in output_files if f not in input_files]

    logging.info(f"Output files to be copied back: {output_files}")

    # Helper: streaming SHA256 (avoids loading large files fully into RAM)
    def _sha256(path: str, _buf_size: int = 1024 * 1024) -> str:
        h = hashlib.sha256()
        try:
            with open(path, "rb") as f:
                while True:
                    chunk = f.read(_buf_size)
                    if not chunk:
                        break
                    h.update(chunk)
        except FileNotFoundError:
            return ""  # signal missing
        return h.hexdigest()

    for item in output_files:
        src = os.path.join(scratch_dir, item)
        dst = os.path.join(original_dir, item)

        if not os.path.exists(src):
            logging.warning(
                f"Output item '{item}' does not exist in scratch directory."
            )
            continue

        if os.path.isfile(src):
            if os.path.exists(dst):
                # Hash-Vergleich: nur kopieren wenn Inhalt verschieden
                try:
                    src_hash = _sha256(src)
                    dst_hash = _sha256(dst)
                except Exception as e:
                    logging.debug(f"Hash compare failed for '{item}': {e}; falling back to duplicate copy policy.")
                    src_hash = dst_hash = None  # force copy branch below
                if src_hash and dst_hash and src_hash == dst_hash:
                    logging.info(f"File '{item}' unchanged (SHA256 match). No copy generated.")
                else:
                    base, ext = os.path.splitext(dst)
                    counter = 1
                    new_dst = f"{base}_copy{counter}{ext}"
                    while os.path.exists(new_dst):
                        counter += 1
                        new_dst = f"{base}_copy{counter}{ext}"
                    shutil.copy2(src, new_dst)
                    logging.info(
                        f"File '{item}' already exists and differs (hash/content). Copied as '{os.path.basename(new_dst)}'."
                    )
            else:
                shutil.copy2(src, dst)
                logging.info(f"Copied file '{item}' back to original directory.")
        elif os.path.isdir(src):
            if os.path.exists(dst):
                base = dst
                counter = 1
                new_dst = f"{base}_copy{counter}"
                while os.path.exists(new_dst):
                    counter += 1
                    new_dst = f"{base}_copy{counter}"
                shutil.copytree(src, new_dst)
                logging.info(
                    f"Directory '{item}' already exists. Renamed to '{os.path.basename(new_dst)}'."
                )
            else:
                shutil.copytree(src, dst)
                logging.info(f"Copied directory '{item}' back to original directory.")
        else:
            logging.warning(
                f"Output item '{item}' is neither a file nor a directory. Skipping."
            )

    # --- Cleanup: Remove the scratch directory ---
    try:
        shutil.rmtree(scratch_dir)
        logging.info(f"Removed scratch directory: {scratch_dir}")
        if remove_parent_if_empty:
            parent_dir = os.path.dirname(scratch_dir)
            _rmdir_if_empty(parent_dir)
    except Exception as e:
        logging.error(f"Failed to remove scratch directory '{scratch_dir}': {e}")
    finally:
        print(f"[scratch-debug] leave finalize original={original_dir} scratch={scratch_dir}", flush=True)

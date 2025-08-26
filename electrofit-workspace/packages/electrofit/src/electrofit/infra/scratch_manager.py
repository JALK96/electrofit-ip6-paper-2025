"""Scratch Directory Management (verschoben aus electrofit.scratch.manager)."""
from __future__ import annotations
import filecmp
import logging
import os
import shutil
import hashlib

from electrofit.infra.logging import setup_logging

__all__ = [
    "setup_scratch_directory",
    "directories_differ",
    "finalize_scratch_directory",
]


def setup_scratch_directory(input_files, base_scratch_dir="/scratch/johannal96/tmp"):
    base_scratch_dir = base_scratch_dir or "/tmp/electrofit_scratch"
    base_scratch_dir = os.path.expanduser(os.path.expandvars(base_scratch_dir))
    fullpath = os.getcwd()
    calcdir = os.path.basename(fullpath)
    parent_dir = os.path.dirname(fullpath)
    parent_folder_name = os.path.basename(parent_dir)
    scratch_dir = os.path.join(base_scratch_dir, parent_folder_name, calcdir)
    log_file_path = os.path.join(fullpath, "process.log")
    suppress = any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers)
    setup_logging(log_file_path, suppress_initial_message=suppress)
    os.makedirs(scratch_dir, exist_ok=True)
    logging.info(f"Created scratch directory: {scratch_dir}")
    for file in input_files:
        src = os.path.join(fullpath, file)
        dst = os.path.join(scratch_dir, file)
        if os.path.exists(src):
            if os.path.isfile(src):
                if os.path.abspath(src) != os.path.abspath(dst):
                    shutil.copy2(src, dst)
                    logging.info(f"Copied file '{file}' to scratch directory.")
            elif os.path.isdir(src):
                if os.path.isdir(dst):
                    try:
                        shutil.rmtree(dst)
                        logging.info(f"Removed pre-existing scratch subdirectory '{file}' (stale).")
                    except Exception as e:
                        logging.warning(f"Could not remove existing scratch subdirectory '{file}': {e}. Proceeding.")
                try:
                    shutil.copytree(src, dst)
                    logging.info(f"Copied directory '{file}' to scratch directory.")
                except FileExistsError:
                    logging.warning(f"Directory '{file}' already existed; leaving as-is.")
            else:
                logging.warning(f"Input item '{file}' is neither file nor directory. Skipping copy.")
        else:
            logging.warning(f"Input file or directory '{file}' not found. Skipping copy.")
    logging.info(f"Scratch directory contents after setup: {os.listdir(scratch_dir)}")
    return scratch_dir, fullpath


def directories_differ(src, dst):
    dcmp = filecmp.dircmp(src, dst)
    if dcmp.left_only or dcmp.right_only or dcmp.diff_files or dcmp.funny_files:
        return True
    for subdir in dcmp.subdirs:
        if directories_differ(os.path.join(src, subdir), os.path.join(dst, subdir)):
            return True
    return False


def _rmdir_if_empty(path: str) -> bool:
    try:
        os.rmdir(path)
        logging.info(f"Removed empty parent scratch directory: {path}")
        return True
    except OSError:
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
    print(f"[scratch-debug] enter finalize original={original_dir} scratch={scratch_dir}", flush=True)
    if reason:
        logging.debug(f"Finalize scratch directory reason={reason}")
    if not os.path.isdir(scratch_dir):
        logging.debug(f"Scratch directory already removed: {scratch_dir}")
        print(f"[scratch-debug] skip finalize (missing) {scratch_dir}", flush=True)
        return
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
                        logging.info(f"Original input file '{item}' renamed to '{os.path.basename(renamed_dst)}'.")
                        shutil.copy2(src, dst)
                        logging.info(f"Modified input file '{item}' copied back to original directory.")
                    else:
                        logging.info(f"Input file '{item}' has changes but overwrite disabled.")
                else:
                    logging.info(f"Input file '{item}' unchanged. No action taken.")
            else:
                shutil.copy2(src, dst)
                logging.info(f"Input file '{item}' missing in original. Copied from scratch.")
        elif os.path.isdir(src):
            if os.path.exists(dst):
                if directories_differ(src, dst):
                    if overwrite:
                        renamed_dst = f"{dst}.input_file"
                        os.rename(dst, renamed_dst)
                        logging.info(f"Original input directory '{item}' renamed to '{os.path.basename(renamed_dst)}'.")
                        shutil.copytree(src, dst)
                        logging.info(f"Modified input directory '{item}' copied back.")
                    else:
                        logging.info(f"Input directory '{item}' changed but overwrite disabled.")
                else:
                    logging.info(f"Input directory '{item}' unchanged.")
            else:
                shutil.copytree(src, dst)
                logging.info(f"Input directory '{item}' missing in original. Copied from scratch.")
        else:
            logging.warning(f"Input item '{item}' is neither a file nor a directory. Skipping.")
    if output_files is None:
        output_files = [f for f in os.listdir(scratch_dir) if f not in input_files]
    else:
        output_files = [f for f in output_files if f not in input_files]
    logging.info(f"Output files to be copied back: {output_files}")
    def _sha256(path: str, _buf_size: int = 1024 * 1024) -> str:
        h = hashlib.sha256()
        try:
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(_buf_size), b""):
                    h.update(chunk)
        except FileNotFoundError:
            return ""
        return h.hexdigest()
    for item in output_files:
        src = os.path.join(scratch_dir, item)
        dst = os.path.join(original_dir, item)
        if not os.path.exists(src):
            logging.warning(f"Output item '{item}' does not exist in scratch directory.")
            continue
        if os.path.isfile(src):
            if os.path.exists(dst):
                try:
                    src_hash = _sha256(src)
                    dst_hash = _sha256(dst)
                except Exception as e:
                    logging.debug(f"Hash compare failed for '{item}': {e}; forcing copy")
                    src_hash = dst_hash = None
                if src_hash and dst_hash and src_hash == dst_hash:
                    logging.info(f"File '{item}' unchanged (SHA256 match). No copy.")
                else:
                    base, ext = os.path.splitext(dst)
                    counter = 1
                    new_dst = f"{base}_copy{counter}{ext}"
                    while os.path.exists(new_dst):
                        counter += 1
                        new_dst = f"{base}_copy{counter}{ext}"
                    shutil.copy2(src, new_dst)
                    logging.info(f"File '{item}' differs. Copied as '{os.path.basename(new_dst)}'.")
            else:
                shutil.copy2(src, dst)
                logging.debug(f"Copied file '{item}' back to original directory.") # set the level to debug to reduce log length. 
        elif os.path.isdir(src):
            if os.path.exists(dst):
                base = dst
                counter = 1
                new_dst = f"{base}_copy{counter}"
                while os.path.exists(new_dst):
                    counter += 1
                    new_dst = f"{base}_copy{counter}"
                shutil.copytree(src, new_dst)
                logging.info(f"Directory '{item}' already exists. Copied as '{os.path.basename(new_dst)}'.")
            else:
                shutil.copytree(src, dst)
                logging.info(f"Copied directory '{item}' back to original directory.")
        else:
            logging.warning(f"Output item '{item}' is neither a file nor a directory. Skipping.")
    try:
        shutil.rmtree(scratch_dir)
        logging.info(f"Removed scratch directory: {scratch_dir}")
        if remove_parent_if_empty:
            _rmdir_if_empty(os.path.dirname(scratch_dir))
    except Exception as e:
        logging.error(f"Failed to remove scratch directory '{scratch_dir}': {e}")
    finally:
        print(f"[scratch-debug] leave finalize original={original_dir} scratch={scratch_dir}", flush=True)

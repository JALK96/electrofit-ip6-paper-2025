from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


def load_hb_num_xvg(filename: str) -> np.ndarray:
    """Load data from an XVG file, skipping comment/header lines."""
    data = []
    with open(filename, "r") as file:
        for line in file:
            if line.startswith(("@", "#")):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    floats = [float(part) for part in parts]
                    data.append(floats)
                except ValueError:
                    continue
    return np.array(data)


def refine_atom_name(atom: str) -> str:
    """Convert GROMACS O labels to one-based convention used by plotting tables."""
    if re.fullmatch(r"[A-Za-z]+", atom):
        return atom + "1"
    match = re.fullmatch(r"([A-Za-z]+)(\d+)", atom)
    if match:
        name = match.group(1)
        number = int(match.group(2)) + 1
        return f"{name}{number}"
    return atom


def _extract_o_token(token: str) -> str | None:
    """Return canonical oxygen label token from a donor/acceptor string."""
    if not token:
        return None
    match = re.search(r"O(\d{1,2})(?!\d)", token)
    if match:
        return f"O{match.group(1)}"
    if re.search(r"O(?!\d)", token):
        return "O"
    return None


def parse_hbond_log_pairs(file_path: str | Path) -> list[tuple[str | None, str | None]]:
    """Parse hb.log into ordered donor/acceptor oxygen tokens (O, O1, ...)."""
    pairs: list[tuple[str | None, str | None]] = []
    line_pattern = re.compile(r"^\s*(\S+)\s*-\s*(\S+)\s*$")

    with open(file_path, "r") as file:
        for line in file:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "*", '"')):
                continue
            match = line_pattern.match(stripped)
            if not match:
                continue
            donor_full, acceptor_full = match.groups()
            pairs.append((_extract_o_token(donor_full), _extract_o_token(acceptor_full)))

    return pairs


def parse_hbond_log_to_dataframe(
    file_path: str | Path,
    *,
    one_based_labels: bool = True,
) -> pd.DataFrame:
    """Parse hb.log into a DataFrame with [idx, donor, acceptor]."""
    hbond_pairs = []
    line_pattern = re.compile(r"^\s*(\S+)\s+-\s+(\S+)\s*$")
    atom_pattern = re.compile(r"^[A-Za-z]+\d+([A-Za-z]+\d*)$")

    with open(file_path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            stripped = line.strip()
            if (
                not stripped
                or stripped.startswith("#")
                or stripped.startswith('"""')
                or stripped.startswith("*")
            ):
                continue

            match = line_pattern.match(stripped)
            if not match:
                continue

            donor_full = match.group(1)
            acceptor_full = match.group(2)
            donor_match = atom_pattern.match(donor_full)
            acceptor_match = atom_pattern.match(acceptor_full)
            if not (donor_match and acceptor_match):
                continue

            donor_atom = donor_match.group(1)
            acceptor_atom = acceptor_match.group(1)
            if one_based_labels:
                donor_atom = refine_atom_name(donor_atom)
                acceptor_atom = refine_atom_name(acceptor_atom)

            hbond_pairs.append({"donor": donor_atom, "acceptor": acceptor_atom})

    df = pd.DataFrame(hbond_pairs)
    df.reset_index(inplace=True)
    df.rename(columns={"index": "idx"}, inplace=True)
    df["idx"] = df.index
    return df


def parse_xpm_binary(
    file_path: str | Path,
    *,
    align_rows_to_log: bool = True,
) -> np.ndarray:
    """Parse an XPM file into a binary matrix (rows=bonds, cols=frames)."""
    data_matrix, _ = parse_xpm(file_path, align_rows_to_log=align_rows_to_log)
    return data_matrix


def parse_xpm(
    file_path: str | Path,
    *,
    align_rows_to_log: bool = True,
) -> tuple[np.ndarray, dict]:
    """Parse XPM into a binary matrix and metadata.

    Notes
    -----
    GROMACS hbond XPM files are rendered with image-row order opposite to the
    hb.log/hb_idx logical bond index order. By default, this parser returns rows
    aligned to hb.log order (`align_rows_to_log=True`).
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    metadata: dict[str, str] = {}
    for raw in lines:
        stripped = raw.strip()
        if not stripped.startswith("/*"):
            continue
        comment = stripped.strip("/* ").strip(" */")
        if ":" in comment:
            key, value = comment.split(":", 1)
            metadata[key.strip().lower()] = value.strip().strip('"')

    header_index = None
    for i, line in enumerate(lines):
        match = re.match(r'^\s*"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
        if match:
            width, height, num_colors, chars_per_pixel = map(int, match.groups())
            header_index = i
            break
    if header_index is None:
        raise ValueError("XPM header not found")

    color_map: dict[str, str] = {}
    i = header_index + 1
    for _ in range(num_colors):
        line = lines[i].rstrip("\n")
        i += 1
        q1 = line.find('"')
        q2 = line.rfind('"')
        if q1 < 0 or q2 <= q1:
            raise ValueError(f"Bad XPM color line: {line}")
        body = line[q1 + 1 : q2]
        symbol = body[:chars_per_pixel]
        rest = body[chars_per_pixel:]
        match = re.search(r'\bc\s+([^ \t"]+)', rest, flags=re.I)
        color = (match.group(1) if match else "").lower()
        color_map[symbol] = color

    present_symbols = {s for s, c in color_map.items() if c in ("#ff0000", "red", "#f00")}
    if not present_symbols:
        present_symbols = {s for s, c in color_map.items() if c not in ("#ffffff", "white")}
        if not present_symbols:
            present_symbols = {next(iter(color_map))}

    bg_symbols = [s for s, c in color_map.items() if c in ("#ffffff", "white")]
    bg_symbol = bg_symbols[0] if bg_symbols else next(iter(color_map))

    data = []
    collected_rows = 0
    needed_chars = width * chars_per_pixel
    while collected_rows < height and i < len(lines):
        line = lines[i].rstrip("\n")
        i += 1
        if "/*" in line or not line.strip().startswith('"'):
            continue
        q1 = line.find('"')
        q2 = line.rfind('"')
        if q1 < 0 or q2 <= q1:
            continue
        body = line[q1 + 1 : q2]
        if len(body) < needed_chars:
            pad_n = (needed_chars - len(body)) // max(1, chars_per_pixel)
            body = body + (bg_symbol * pad_n)
        elif len(body) > needed_chars:
            body = body[:needed_chars]

        if chars_per_pixel == 1:
            tokens = list(body)
        else:
            tokens = [body[k : k + chars_per_pixel] for k in range(0, needed_chars, chars_per_pixel)]
        if len(tokens) != width:
            raise ValueError(f"XPM row length mismatch: got {len(tokens)}, expected {width}")
        data.append([1 if tok in present_symbols else 0 for tok in tokens])
        collected_rows += 1

    if collected_rows != height:
        raise ValueError(f"Collected {collected_rows} XPM rows, expected {height}")

    matrix = np.array(data, dtype=np.uint8)
    if align_rows_to_log:
        matrix = matrix[::-1, :]

    return matrix, metadata


def analyze_hydrogen_bonds(data_matrix: np.ndarray, metadata: dict | None = None) -> dict:
    """Compute common summaries from a binary H-bond matrix."""
    _ = metadata  # kept for API compatibility

    hbonds_over_time = np.sum(data_matrix, axis=0)
    hbonds_per_index = np.sum(data_matrix, axis=1)

    lifetimes = []
    for bond in data_matrix:
        current_lifetime = 0
        bond_lifetimes = []
        for state in bond:
            if state == 1:
                current_lifetime += 1
            else:
                if current_lifetime > 0:
                    bond_lifetimes.append(current_lifetime)
                    current_lifetime = 0
        if current_lifetime > 0:
            bond_lifetimes.append(current_lifetime)
        lifetimes.append(bond_lifetimes)

    return {
        "hbonds_over_time": hbonds_over_time,
        "hbonds_per_index": hbonds_per_index,
        "lifetimes": lifetimes,
    }


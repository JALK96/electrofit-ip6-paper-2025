import re
import numpy as np
import pandas as pd
import logging
from pathlib import Path
from typing import Iterable, List
from electrofit_analysis.structure.util.hbond_io import (
    analyze_hydrogen_bonds,
    load_hb_num_xvg,
    parse_hbond_log_to_dataframe,
    parse_xpm,
    refine_atom_name,
)
from electrofit_analysis.cli.h_bonds.make_pp_matrix_ip6 import (
    draw_phosphorus_diagram as _draw_phosphorus_diagram_shared,
)


logger = logging.getLogger()

def count_donor_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Sums frames for 'donor' => avoid double-count. Returns Series: donor->frames
    """
    logger.info("Counting DONOR-based frames only.")
    df=hbond_df.copy()
    if 'count' in df.columns:
        df.drop(columns=['count'],inplace=True)
    df['count']=hbonds_per_index[df['idx'].values]
    melted=df.melt(
        id_vars=['idx','count'],
        value_vars=['donor'],
        value_name='oxygen'
    ).dropna(subset=['oxygen'])
    grouped=melted.groupby('oxygen')['count'].sum()
    logger.info(f"Donor-based occurrences => {len(grouped)} donors.")
    return grouped

def count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Summation for both donor & acceptor => used for final existence heatmap
    """
    logger.info("Counting oxygen occurrences (both donor & acceptor).")
    df=hbond_df.copy()
    if 'count' in df.columns:
        df.drop(columns=['count'],inplace=True)
    df['count']=hbonds_per_index[df['idx'].values]
    melted=df.melt(
        id_vars=['idx','count'],
        value_vars=['donor','acceptor'],
        value_name='oxygen'
    ).dropna(subset=['oxygen'])
    grouped=melted.groupby('oxygen')['count'].sum()
    logger.info(f"Found {len(grouped)} total oxygen atoms.")
    return grouped

def expand_phosphate_data_for_violins(df_phosphate):
    """
    Expand each row of df_phosphate (which has columns e.g.:
       [species, phosphate_group, engaged_ns]
    into multiple rows, effectively simulating a distribution
    for each (species, phosphate_group).

    Returns a new DataFrame with columns [species, pgroup_index].
    'pgroup_index' is repeated ~ engaged_ns times for each row.
    """
    # 1) Map each phosphate group to an integer index
    #    e.g. P1->0, P2->1, ... P6->5
    pgroup_index_map = {
        'P1': 0,
        'P2': 1,
        'P3': 2,
        'P4': 3,
        'P5': 4,
        'P6': 5
    }

    expanded_rows = []

    for row in df_phosphate.itertuples(index=False):
        species_id = row.species
        pgroup     = row.phosphate_group
        engaged_ns = row.engaged_ns * 1000

        # choose how many times to repeat. E.g. round it:
        repeat_count = int(round(engaged_ns))  
        # If engaged_ns is large, repeat_count could be big => large DataFrame.

        idx_val = pgroup_index_map.get(pgroup, None)
        if idx_val is None:
            continue  # skip unknown group

        # expand
        for _ in range(repeat_count):
            expanded_rows.append({
                'species': species_id,
                'pgroup_index': idx_val
            })

    expanded_df = pd.DataFrame(expanded_rows)
    return expanded_df

def draw_phosphorus_diagram(species_id, top_three_map, edge_weights=None):
    """Backward-compatible wrapper using the shared phosphorus graph renderer."""
    _draw_phosphorus_diagram_shared(
        species_bits=species_id,
        species_id=species_id,
        top_three_map=top_three_map,
        edge_weights=edge_weights,
        width_mode="continuous",
        width_ref="absolute",
    )




# -----------------------
# Utilities (dihedrals_ip6)
# -----------------------
def shift_atom_number(name: str) -> str:
    """
    If 'name' ends with a number, shift that number by 1. Example: C1 -> C2.
    If no number is present, append '1'.
    """
    match = re.match(r"([A-Za-z]+)(\d+)?$", name)
    if not match:
        return name
    base, num = match.group(1), match.group(2)
    return f"{base}{int(num)+1}" if num is not None else f"{base}1"


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def chunks(lst: List, n: int) -> Iterable[List]:
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import networkx as nx
from pathlib import Path
from typing import Iterable, List
from electrofit_analysis.structure.util.hbond_io import (
    analyze_hydrogen_bonds as _shared_analyze_hydrogen_bonds,
    load_hb_num_xvg as _shared_load_hb_num_xvg,
    parse_hbond_log_to_dataframe as _shared_parse_hbond_log_to_dataframe,
    parse_xpm as _shared_parse_xpm,
    refine_atom_name as _shared_refine_atom_name,
)


logger = logging.getLogger()

# ------------------------------------------------------------------------------
# Common Utility Functions (hb_comparison.py)
# ------------------------------------------------------------------------------
def refine_atom_name(atom):
    """Thin wrapper over the shared H-bond atom-name canonicalizer."""
    return _shared_refine_atom_name(atom)

def load_hb_num_xvg(filename):
    """
    (#H-bonds vs time) from .xvg => Nx2 array.
    """
    logger.info(f"Loading H-bond number data from: {filename}")
    arr = _shared_load_hb_num_xvg(filename)
    logger.info(f"XVG shape: {arr.shape}")
    return arr

def parse_xpm(file_path, align_rows_to_log=True):
    """Parse H-bond XPM and align rows to hb.log order by default."""
    logger.info(f"Parsing XPM: {file_path}")
    arr, metadata = _shared_parse_xpm(
        file_path,
        align_rows_to_log=align_rows_to_log,
    )
    logger.info(f"XPM matrix => shape={arr.shape}")
    return arr, metadata

def analyze_hydrogen_bonds(data_matrix, metadata):
    """Summaries => over_time, per_index, lifetimes."""
    logger.info("Analyzing hydrogen bond matrix.")
    return _shared_analyze_hydrogen_bonds(data_matrix, metadata)

def parse_hbond_log_to_dataframe(file_path):
    """
    GROMACS .log => DataFrame [idx, donor, acceptor], refining numeric bits
    """
    logger.info(f"Parsing H-bond log: {file_path}")
    df = _shared_parse_hbond_log_to_dataframe(file_path)
    logger.info(f"Found {len(df)} H-bond pairs.")
    return df

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
    """
    Draws a directed graph (P1..P6) showing donor->target relationships,
    with the species label centered and pure‐acceptor nodes in white.
    """
    if edge_weights is None:
        edge_weights = {}

    # 1) Build graph
    G = nx.DiGraph()
    for i in range(1, 7):
        G.add_node(i)
    for donor, ranks in top_three_map.items():
        for rank, target in ranks.items():
            G.add_edge(donor, target, rank=rank)

    # 2) Custom circular layout (rotated so P1 sits where P3 normally is)
    oldpos = nx.circular_layout(G)
    pos = {1: oldpos[3], 2: oldpos[4], 3: oldpos[5],
           4: oldpos[6], 5: oldpos[1], 6: oldpos[2]}

    # 3) Determine node colors: 
    #    white for nodes with in_degree>0 and out_degree==0 ("pure acceptors"), else lightgray
    node_colors = []
    for n in G.nodes():
        if G.in_degree(n) > 0 and G.out_degree(n) == 0:
            node_colors.append("white")
        else:
            node_colors.append("lightgray")

    # 4) Prepare figure
    fig, ax = plt.subplots(figsize=(8,6))

    # 5) Draw nodes
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors,
        edgecolors="black", node_size=1300, ax=ax
    )

    # 6) Split edges by rank
    edges_first  = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="first"]
    edges_second = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="second"]
    edges_third  = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="third"]

    #  – edge-width scaling
    min_w, max_w = 1.0, 6.0
    if edge_weights:
        max_val = max(edge_weights.values()) or 1.0
        def scale(t):
            return min_w + (max_w - min_w) * (t / max_val)
    else:
        def scale(t):
            return min_w

    # 7) Draw first‐tier (red), handling reciprocals
    mutual = [(u,v) for (u,v) in edges_first if (v,u) in edges_first and u<v]
    mutual_set = set(mutual) | set((v,u) for (u,v) in mutual)
    nonmutual = [e for e in edges_first if e not in mutual_set]

    nx.draw_networkx_edges(
        G, pos, edgelist=nonmutual,
        width=[scale(edge_weights.get(e,0.0)) for e in nonmutual],
        edge_color="red", arrows=True, arrowstyle='-|>', arrowsize=20,
        connectionstyle="arc3,rad=0.0",
        min_source_margin=45, min_target_margin=45, ax=ax
    )
    for u,v in mutual:
        w_uv = scale(edge_weights.get((u,v),0.0))
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u,v)],
            width=[w_uv], edge_color="red", arrows=True,
            arrowstyle='-|>', arrowsize=20,
            connectionstyle="arc3,rad=+0.2",
            min_source_margin=45, min_target_margin=45, ax=ax
        )
        w_vu = scale(edge_weights.get((v,u),0.0))
        nx.draw_networkx_edges(
            G, pos, edgelist=[(v,u)],
            width=[w_vu], edge_color="red", arrows=True,
            arrowstyle='-|>', arrowsize=20,
            connectionstyle="arc3,rad=+0.2",
            min_source_margin=45, min_target_margin=45, ax=ax
        )

    # 8) Second‐tier (black)
    nx.draw_networkx_edges(
        G, pos, edgelist=edges_second,
        width=[scale(edge_weights.get(e,0.0)) for e in edges_second],
        edge_color="black", arrows=True, arrowstyle='-|>', arrowsize=20,
        connectionstyle="arc3,rad=+0.2",
        min_source_margin=45, min_target_margin=45, ax=ax
    )

    # 9) Third‐tier (blue)
    nx.draw_networkx_edges(
        G, pos, edgelist=edges_third,
        width=[scale(edge_weights.get(e,0.0)) for e in edges_third],
        edge_color="blue", arrows=True, arrowstyle='-|>', arrowsize=20,
        connectionstyle="arc3,rad=-0.2",
        min_source_margin=45, min_target_margin=45, ax=ax
    )

    # 10) Draw node labels
    nx.draw_networkx_labels(
        G, pos, labels={i: f"P{i}" for i in G.nodes()}, font_size=14, ax=ax
    )

    # 11) Add centered species label in the middle of the plot
    #    Use axis‐coordinates (0.5,0.5) so it's always at the visual center
    ax.text(
        0.5, 0.5, f"({species_id})",
        fontsize=18, fontweight="bold",
        ha="center", va="center",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", alpha=0)
    )

    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(f"species_{species_id}.pdf")
    plt.close()




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



import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import networkx as nx
from pathlib import Path
from typing import Iterable, List


logger = logging.getLogger()

# ------------------------------------------------------------------------------
# Common Utility Functions (hb_comparison.py)
# ------------------------------------------------------------------------------
def refine_atom_name(atom):
    """
    If we see e.g. 'O' => 'O1', or 'O10' => 'O11'.
    """
    if re.fullmatch(r'[A-Za-z]+', atom):
        return atom + '1'
    match = re.fullmatch(r'([A-Za-z]+)(\d+)', atom)
    if match:
        name = match.group(1)
        number = int(match.group(2)) + 1
        return f"{name}{number}"
    return atom

def load_hb_num_xvg(filename):
    """
    (#H-bonds vs time) from .xvg => Nx2 array.
    """
    logger.info(f"Loading H-bond number data from: {filename}")
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('@', '#')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    floats = list(map(float, parts[:2]))
                    data.append(floats)
                except ValueError:
                    pass
    arr = np.array(data)
    logger.info(f"XVG shape: {arr.shape}")
    return arr

def parse_xpm(file_path):
    """Parses XPM => binary array (#FF0000 =>1). Returns (arr, metadata)."""
    logger.info(f"Parsing XPM: {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    metadata={}
    data_lines=[]
    color_map={}
    header_found=False
    data_started=False
    width=height=num_colors=chars_per_pixel=None

    for idx, raw_line in enumerate(lines, start=1):
        line=raw_line.strip()
        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ':' in comment:
                key, value = comment.split(":",1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        if line.startswith('static char'):
            continue

        if (not header_found) and line.startswith('"'):
            line_str=line.strip(',').lstrip('"').rstrip('"')
            tokens=line_str.split()
            if len(tokens)>=4:
                try:
                    width,height,num_colors,chars_per_pixel=map(int,tokens[:4])
                    header_found=True
                    logger.info(f"XPM header => w={width},h={height},colors={num_colors},cpp={chars_per_pixel}")
                    continue
                except ValueError:
                    logger.error(f"Invalid header line {idx}: {raw_line}")
                    raise
            else:
                logger.error(f"Header line {idx} missing tokens: {raw_line}")
                raise ValueError(f"Bad XPM header line: {raw_line}")
        if header_found and not data_started:
            cdef=line.strip(',').lstrip('"').rstrip('"')
            pattern=rf'(.{{{chars_per_pixel}}})\s+c\s+(\S+)'
            match=re.match(pattern, cdef)
            if match:
                symbol,color_val=match.groups()
                color_map[symbol]=color_val
            if len(color_map)==num_colors:
                data_started=True
            continue
        if data_started and line.startswith('"'):
            row_data=line.strip(',').lstrip('"').rstrip('"')
            data_lines.append(row_data)

    if width is None or height is None:
        logger.error("No valid XPM header found.")
        raise ValueError("XPM header missing.")
    if len(data_lines)!=height:
        logger.warning(f"XPM mismatch: expected {height} lines,found {len(data_lines)}")

    arr=np.zeros((height,width), dtype=int)
    for y,row_str in enumerate(data_lines):
        for x,ch in enumerate(row_str):
            color=color_map.get(ch,None)
            if color=='#FF0000':
                arr[y,x]=1
    logger.info(f"XPM matrix => shape={arr.shape}")
    return arr,metadata

def analyze_hydrogen_bonds(data_matrix, metadata):
    """Summaries => over_time, per_index, lifetimes."""
    logger.info("Analyzing hydrogen bond matrix.")
    hbonds_over_time=np.sum(data_matrix,axis=0)
    hbonds_per_index=np.sum(data_matrix,axis=1)
    lifetimes=[]
    for row in data_matrix:
        cur=0
        run_list=[]
        for val in row:
            if val==1:
                cur+=1
            else:
                if cur>0:
                    run_list.append(cur)
                    cur=0
        if cur>0:
            run_list.append(cur)
        lifetimes.append(run_list)
    return {
        'hbonds_over_time':hbonds_over_time,
        'hbonds_per_index':hbonds_per_index,
        'lifetimes':lifetimes
    }

def parse_hbond_log_to_dataframe(file_path):
    """
    GROMACS .log => DataFrame [idx, donor, acceptor], refining numeric bits
    """
    logger.info(f"Parsing H-bond log: {file_path}")
    pairs=[]
    line_pattern=re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')
    atom_pattern=re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')
    with open(file_path,'r') as f:
        for line_no,raw_line in enumerate(f,start=1):
            line=raw_line.strip()
            if not line or line.startswith('#') or line.startswith('"""') or line.startswith('*'):
                continue
            match=line_pattern.match(line)
            if match:
                d_full,a_full=match.groups()
                d_match=atom_pattern.match(d_full)
                a_match=atom_pattern.match(a_full)
                if d_match and a_match:
                    donor_atom=refine_atom_name(d_match.group(1))
                    acceptor_atom=refine_atom_name(a_match.group(1))
                    pairs.append({'donor':donor_atom,'acceptor':acceptor_atom})
                else:
                    logger.warning(f"Line {line_no}: parse fail => {line}")
            else:
                logger.warning(f"Line {line_no}: no match => {line}")
    df=pd.DataFrame(pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index':'idx'},inplace=True)
    df['idx']=df.index
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




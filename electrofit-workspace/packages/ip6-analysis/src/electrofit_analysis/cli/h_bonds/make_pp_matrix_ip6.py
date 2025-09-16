#!/usr/bin/env python3
import os
import re
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
from matplotlib.patches import FancyArrowPatch


# MODE notes:
#  - 'union' (default): fraction of time with at least one O–O bond for each P_i->P_j (per-frame OR)
#  - 'sum'           : sum of per-O–O fractions (can be >1.0)
#  - 'both'          : writes both outputs; union keeps base name, sum uses suffix '_sum'

#
# ---------- O-atom -> phosphate mapping ----------

# GROMACS mapping (as found in hb.log): oxygen token -> GROMACS phosphate label (P, P1..P5)
GRO_O2P = {
    'O':'P',  'O6':'P',  'O7':'P',  'O8':'P',
    'O1':'P1','O9':'P1','O10':'P1','O11':'P1',
    'O2':'P2','O12':'P2','O13':'P2','O14':'P2',
    'O3':'P3','O15':'P3','O16':'P3','O17':'P3',
    'O4':'P4','O18':'P4','O19':'P4','O20':'P4',
    'O5':'P5','O21':'P5','O22':'P5','O23':'P5'
}
# Relabel GROMACS phosphate labels to paper convention P1..P6
RELABEL_P = {'P':'P1','P1':'P2','P2':'P3','P3':'P4','P4':'P5','P5':'P6'}

Lauren = False
if Lauren: 
    # Lauren has already correct labels - her map here: 
    GRO_O2P = {
        'O1': 'P1',  'O7': 'P1',  'O8': 'P1',  'O9': 'P1',
        'O2': 'P2', 'O10': 'P2', 'O11': 'P2', 'O12': 'P2',
        'O3': 'P3', 'O13': 'P3', 'O14': 'P3', 'O15': 'P3',
        'O4': 'P4', 'O16': 'P4', 'O17': 'P4', 'O18': 'P4',
        'O5': 'P5', 'O19': 'P5', 'O20': 'P5', 'O21': 'P5',
        'O6': 'P6', 'O22': 'P6', 'O23': 'P6', 'O24': 'P6'
    }
    RELABEL_P = {'P1':'P1','P2':'P2','P3':'P3','P4':'P4','P5':'P5', 'P6':'P6'}

P2IDX = {f'P{i}': i-1 for i in range(1,7)}

# Map an oxygen token (GROMACS) to a paper phosphate label (P1..P6)
def map_oxygen_to_paper_pg(otoken: str):
    pg = GRO_O2P.get(otoken)
    if not pg:
        return None
    return RELABEL_P.get(pg, pg)


def parse_xpm_binary(xpm_path):
    import numpy as np
    import re

    with open(xpm_path, "r") as f:
        lines = f.readlines()

    # --- header ---
    hdr_i = None
    for i, L in enumerate(lines):
        m = re.match(r'^\s*"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', L)
        if m:
            w, h, ncolors, cpp = map(int, m.groups())
            hdr_i = i
            break
    if hdr_i is None:
        raise ValueError("XPM header not found")

    # --- colors ---
    cmap = {}  # symbol -> color (lowercase)
    i = hdr_i + 1
    for _ in range(ncolors):
        L = lines[i].rstrip("\n")
        i += 1
        # keep everything inside the first quoted string verbatim
        q1 = L.find('"')
        q2 = L.rfind('"')
        if q1 < 0 or q2 <= q1:
            raise ValueError(f"Bad color line: {L}")
        body = L[q1+1:q2]              # do NOT strip; preserves spaces
        symbol = body[:cpp]
        rest = body[cpp:]
        mcol = re.search(r'\bc\s+([^ \t"]+)', rest, flags=re.I)
        color = (mcol.group(1) if mcol else "").lower()
        cmap[symbol] = color

    present_syms = {s for s,c in cmap.items() if c in ("#ff0000", "red", "#f00")}
    if not present_syms:
        # fallback: anything not white
        present_syms = {s for s,c in cmap.items() if c not in ("#ffffff", "white")}
        if not present_syms:
            present_syms = {next(iter(cmap))}

    bg_syms = [s for s,c in cmap.items() if c in ("#ffffff", "white")]
    bg = bg_syms[0] if bg_syms else next(s for s in cmap if s not in present_syms)

    # --- pixel rows ---
    data = []
    rows = 0
    need = w * cpp
    while rows < h and i < len(lines):
        L = lines[i].rstrip("\n")
        i += 1
        if "/*" in L or not L.strip().startswith('"'):
            continue  # skip comments and non-row lines
        q1 = L.find('"')
        q2 = L.rfind('"')
        if q1 < 0 or q2 <= q1:
            continue
        body = L[q1+1:q2]  # keep spaces!
        if len(body) < need:
            # pad with background symbol (cpp-safe)
            body = body + bg * ((need - len(body)) // max(1, cpp))
        elif len(body) > need:
            body = body[:need]

        if cpp == 1:
            tokens = list(body)  # single chars, includes spaces
        else:
            tokens = [body[k:k+cpp] for k in range(0, need, cpp)]
        if len(tokens) != w:
            raise ValueError(f"Row {rows}: got {len(tokens)} tokens, expected {w}")
        data.append([1 if t in present_syms else 0 for t in tokens])
        rows += 1

    if rows != h:
        raise ValueError(f"Collected {rows} rows, expected {h}")

    return np.array(data, dtype=np.uint8)  # shape (h=bonds, w=frames)

# --- helper: always treat GROMACS O labels as 1-based (O, O1, O2, ...)
def _extract_O_token(token: str) -> str | None:
    """Return canonical oxygen label used for mapping.
    Rules:
    - Prefer an explicit numeric oxygen label anywhere in the token: O<1-2 digits> (e.g., "IP61O15" -> "O15").
    - Otherwise, accept a bare 'O' that is **not** followed by a digit, even if it is glued to residue text
      (e.g., "IP61O" -> "O").
    - If no oxygen-like token is present, return None (we still keep row alignment; mapping may skip it).
    """
    if not token:
        return None
    # 1) Numeric oxygen like O7, O14, up to two digits; allow trailing non-digits
    m = re.search(r"O(\d{1,2})(?!\d)", token)
    if m:
        return f"O{m.group(1)}"
    # 2) Bare 'O' not followed by a digit (e.g., 'IP61O')
    if re.search(r"O(?!\d)", token):
        return "O"
    return None

# ---------- Parse hb.log into ordered donor/acceptor pairs (keep ALL rows) ----------
def parse_hbond_log(log_path):
    pairs = []
    # capture two non-empty tokens around a single dash
    pat = re.compile(r'^\s*(\S+)\s*-\s*(\S+)\s*$')
    with open(log_path,'r') as f:
        for ln in f:
            if not ln.strip() or ln.lstrip().startswith(('#','*','"')):
                continue
            m = pat.match(ln)
            if not m:
                continue
            d_raw, a_raw = m.groups()
            d_tok = _extract_O_token(d_raw)
            a_tok = _extract_O_token(a_raw)
            # keep alignment even if one side isn't an oxygen; mapping will skip it later
            pairs.append((d_tok, a_tok))
    return pairs  # length should match XPM height

# ---------- Build 6x6 matrix of per-time fractions
# mode:
#   - 'union' (default): per-frame OR across all O–O rows for a P_i->P_j pair
#   - 'sum'            : sum of per-row fractions (can exceed 1.0 if multiple O–O bonds coexist)
#                        This mirrors the comparison script that sums lifetimes per O–O bond.
def build_matrix(xpm_mat, pairs, mode: str = 'union'):
    nbonds, nframes = xpm_mat.shape
    assert nbonds == len(pairs), "XPM rows must match log pairs count."
    # collect row indices per phosphate pair
    rows_by_pair = {(i,j): [] for i in range(6) for j in range(6)}
    for row_idx, (dO, aO) in enumerate(pairs):
        dP = map_oxygen_to_paper_pg(dO)
        aP = map_oxygen_to_paper_pg(aO)
        if dP and aP:
            rows_by_pair[(P2IDX[dP], P2IDX[aP])].append(row_idx)
    mapped_rows = sum(len(v) for v in rows_by_pair.values())
    if mapped_rows == 0:
        print("[WARN] No O->O rows mapped to phosphate pairs. Check O-label extraction and GRO_O2P mapping for this system.")

    M = np.zeros((6,6), dtype=float)
    for (i,j), rows in rows_by_pair.items():
        if not rows:
            continue
        sub = xpm_mat[rows, :]  # shape (n_rows_for_pair, nframes)
        if mode == 'union':
            # fraction of frames with at least one O–O bond present
            present_any = sub.max(axis=0)  # 0/1 per frame
            M[i, j] = present_any.mean()
        elif mode == 'sum':
            # sum of per-row fractions (average number of simultaneous O–O bonds)
            # equals sum over rows of (frames_with_bond / nframes)
            M[i, j] = sub.mean(axis=1).sum()
        else:
            raise ValueError("mode must be 'union' or 'sum'")
    return M

# ---------- Save helpers ----------
def save_matrix_csv_npy(M, out_stem):
    np.save(out_stem + ".npy", M)
    header = ",".join([f"P{k}" for k in range(1,7)])
    with open(out_stem + ".csv","w") as f:
        f.write("," + header + "\n")
        for i in range(6):
            row = ",".join(f"{M[i,j]:.6f}" for j in range(6))
            f.write(f"P{i+1}," + row + "\n")

# --------- Diagram drawing and helper functions ----------
def draw_phosphorus_diagram(
        species_bits,
        species_id, top_three_map, edge_weights=None,
        width_mode="continuous",        # "continuous" or "coarse"
        width_ref="auto",               # "auto" (relative to max in this graph) or "absolute" (fraction 0..1)
        min_w=1.0, max_w=20.0,          # continuous scaling bounds (up to 20 pt)
        coarse_step=0.01,               # 1% bins by default
        coarse_min=1.0, coarse_max=20.0 # 1..20 pt for coarse
):
    """
    Draws a directed graph (P1..P6) showing donor->target relationships,
    with the species label centered and pure‐acceptor nodes in white.
    """
    # --- layout helpers ---
    def draw_self_loop(ax, center, outward_vec, width,
                   color="orangered", arrowsize=40,
                   r_off=0.2, t_off=0.15, rad=1.5):
        """
        Draw a self-loop outside the hexagon at a node.
        center: (x,y) of node
        outward_vec: radial outward direction (will be normalized)
        width: line width (points)
        r_off: radial offset (data units) from node center
        t_off: tangential half-span (data units)
        rad: curvature for arc3 (positive -> bulge to the left from start->end)
        """
        v = np.array(outward_vec, dtype=float)
        nv = v / (np.linalg.norm(v) if np.linalg.norm(v) else 1.0)
        # tangential direction = +90° rotate
        t = np.array([-nv[1], nv[0]])

        start = center + nv * r_off - t * t_off
        end   = center + nv * r_off + t * t_off

        patch = FancyArrowPatch(
            start, end,
            arrowstyle='-|>',
            mutation_scale=arrowsize,
            lw=width,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            zorder=5,
            clip_on=False,
            joinstyle="miter",
            capstyle="butt",
        )
        ax.add_patch(patch)

    def freeze_equal_hex_bounds(ax, pos, pad=0.2, extra=0.0):
        # radius of your layout
        R = max(np.linalg.norm(v) for v in pos.values())
        # pad is a fractional margin; extra is an absolute extra radius
        Rpad = R * (1 + pad) + extra
        ax.set_xlim(-Rpad, Rpad)
        ax.set_ylim(-Rpad, Rpad)
        ax.set_aspect('equal', adjustable='box')
        ax.set_anchor('C')  # center the plot

    def hexagon_pos_clockwise(order=(1,2,3,4,5,6), start_angle_deg=90, radius=1.0):
        """Return positions at the vertices of a REGULAR HEXAGON, placed CLOCKWISE
        in the given order. start_angle_deg is where the first element in `order`
        sits (0°=+x axis, 90°=+y axis)."""
        theta0 = np.deg2rad(start_angle_deg)
        step = 2*np.pi/6
        pos = {}
        for k, node in enumerate(order):
            theta = theta0 - k*step  # minus => clockwise
            pos[node] = np.array([radius*np.cos(theta), radius*np.sin(theta)], dtype=float)
        return pos

    def rotate_layout(pos, deg):
        """Rotate all points around origin by +deg degrees."""
        t = np.deg2rad(deg)
        R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
        return {k: (R @ v) for k, v in pos.items()}

    if edge_weights is None:
        edge_weights = {}

    # 1) Build graph
    G = nx.DiGraph()
    for i in range(1, 7):
        G.add_node(i)
    for donor, ranks in top_three_map.items():
        for rank, target in ranks.items():
            G.add_edge(donor, target, rank=rank)

    # 2) Layout: clockwise P1..P6 on a REGULAR HEXAGON, then flip across P2–P5 axis and rotate +90°
    pos = hexagon_pos_clockwise(order=[2,3,4,5,6,1], start_angle_deg=0)
    pos = rotate_layout(pos, 90)          # rotate whole network by +90°

    # 3) Determine node colors:


    def _parse_species_bits(species_id: str):
        """
        Extract the last 6 binary digits from the species label.
        Accepts formats like '010101' or 'IP_010101'.
        Returns a dict {1..6 -> 0/1} or None if not found.
        """
        bits = re.findall(r'[01]', species_id)
        if len(bits) >= 6:
            last6 = bits[-6:]
            return {i+1: (1 if last6[i] == '1' else 0) for i in range(6)}
        return None

    # --- inside draw_phosphorus_diagram(...) just before drawing nodes ---
    bitmap = _parse_species_bits(species_bits)

    if bitmap is not None:
        # color by microstate: 1=protonated (blue), 0=deprotonated (white)
        blue_rgb = (0.7, 0.7, 1.0)  # your current non-transparent blue
        node_colors = [blue_rgb if bitmap[n] == 1 else "white" for n in G.nodes()]
    else:
        # fallback (old behavior): white if pure acceptor, else blue
        blue_rgb = (0.7, 0.7, 1.0)
        node_colors = []
        for n in G.nodes():
            if G.in_degree(n) > 0 and G.out_degree(n) == 0:
                node_colors.append("white")
            else:
                node_colors.append(blue_rgb)

    # 4) Set up the plot
    fig, ax = plt.subplots(figsize=(8,8))  # square figure helps too
    # Make figure and axes fully transparent
    fig.patch.set_facecolor('none')
    ax.set_facecolor('none')
    freeze_equal_hex_bounds(ax, pos, pad=0.30, extra=0.25)
    # --- thin dotted scaffold: connect adjacent P nodes on the hexagon ---
    scaffold_edges = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1)]
    nx.draw_networkx_edges(
        G, pos, edgelist=scaffold_edges,
        edge_color="black", style="dotted", width=1.2,
        arrows=False, ax=ax
    )

    # 5) Draw nodes
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors,
        edgecolors="black", node_size=3600, ax=ax
    )

    # 6) Split edges by rank
    edges_first  = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="first"]
    edges_second = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="second"]
    edges_third  = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="third"]

    # --- carve out self-loops (u==v) ---
    def _split_loops(edgelist):
        loops = [e for e in edgelist if e[0] == e[1]]
        non   = [e for e in edgelist if e[0] != e[1]]
        return non, loops

    edges_first,  loops1 = _split_loops(edges_first)
    edges_second, loops2 = _split_loops(edges_second)
    edges_third,  loops3 = _split_loops(edges_third)
    self_loops = loops1 + loops2 + loops3

    #  – edge-width scaling (continuous vs coarse-grained)
    def _scale_continuous(f):
        """
        Continuous scaling.
        - width_ref == "auto": scale relative to the maximum edge weight in this graph.
        - width_ref == "absolute": scale by the absolute fraction f (0..1), so 1.0 -> max_w.
        """
        if width_ref == "absolute":
            ff = max(0.0, min(1.0, float(f)))
            return min_w + (max_w - min_w) * ff
        # default: auto/relative to the maximum in this graph
        mv = max(edge_weights.values()) if edge_weights else 0.0
        mv = mv if mv > 0 else 1.0
        return min_w + (max_w - min_w) * (float(f) / mv)

    def _scale_coarse(f):
        """
        Coarse-grained widths with binning.

        - width_ref == "auto": normalize by the maximum edge weight in this graph,
          then bin the normalized value in [0,1].
        - width_ref == "absolute": bin the absolute fraction f in [0,1].
        """
        # choose base in [0,1]
        if width_ref == "auto":
            mv = max(edge_weights.values()) if edge_weights else 0.0
            base = (float(f) / mv) if mv > 0 else 0.0
        else:
            base = float(f)
        base = 0.0 if base < 0 else (1.0 if base > 1.0 else base)

        step = max(coarse_step, 1e-9)
        nbins = max(1, int(round(1.0 / step)))
        # bin index 1..nbins
        width_idx = int(np.ceil(base / step))
        width_idx = 1 if width_idx < 1 else (nbins if width_idx > nbins else width_idx)

        if nbins == 1:
            return coarse_max
        return coarse_min + (coarse_max - coarse_min) * (width_idx - 1) / (nbins - 1)

    scale = _scale_continuous if width_mode == "continuous" else _scale_coarse

    # --- arrowhead sizing + connection style helpers ---
    def _scaled_arrowsize(w, amin=20.0, amax=30.0):
        """Scale arrowhead size with the visual line width."""
        denom = (max_w - min_w) if (max_w - min_w) > 1e-9 else 1.0
        frac = (w - min_w) / denom
        frac = 0.0 if frac < 0 else (1.0 if frac > 1.0 else frac)
        return amin + frac * (amax - amin)

    def _connstyle(u, v):
        """Exception: P2<->P5 edges should be straight lines (rad=0)."""
        if (u, v) in ((2,5), (5,2)):
            return "arc3,rad=0.0"
        return "arc3,rad=+0.4"

    def _draw_edge_list(edgelist):
        """Draw edges one-by-one with FancyArrowPatch so we can set
        join/cap styles and custom curvature. Also shrink arrow at the
        node boundaries to avoid overlapping node circles.
        """
        for (u, v) in edgelist:
            w  = scale(edge_weights.get((u, v), 0.0))
            az = _scaled_arrowsize(w)
            # straight line for P2<->P5
            rad = 0.0 if (u, v) in ((2, 5), (5, 2)) else 0.4

            p1 = np.asarray(pos[u], dtype=float)
            p2 = np.asarray(pos[v], dtype=float)

            patch = FancyArrowPatch(
                p1, p2,
                arrowstyle='-|>',
                mutation_scale=az,
                lw=w,
                color='orangered',
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=45, shrinkB=45,  # roughly match previous margins
                joinstyle='miter',
                capstyle='butt',
                zorder=5,
                clip_on=False,
            )
            ax.add_patch(patch)

    # 7) Draw first‐tier (red), handling reciprocals
    mutual = [(u,v) for (u,v) in edges_first if (v,u) in edges_first and u<v]
    mutual_set = set(mutual) | set((v,u) for (u,v) in mutual)
    nonmutual = [e for e in edges_first if e not in mutual_set]

    _draw_edge_list(nonmutual)
    for u, v in mutual:
        _draw_edge_list([(u, v)])
        _draw_edge_list([(v, u)])

    # 8) Second‐tier (red — per your current code)
    _draw_edge_list(edges_second)

    # 9) Third‐tier (red — per your current code)
    _draw_edge_list(edges_third)
    # --- self-loops outside each node, aligned with outward vector ---
    for (n, _n) in self_loops:
        center = np.asarray(pos[n], dtype=float)
        outward = center.copy()  # hexagon is centered at (0,0), so pos is radial
        w = scale(edge_weights.get((n, n), 0.0))
        az = _scaled_arrowsize(w)
        draw_self_loop(ax, center, outward, width=w, color="orangered", arrowsize=az)
        
    # 10) Draw node labels
    nx.draw_networkx_labels(
        G, pos, labels={i: f"P{i}" for i in G.nodes()}, font_size=25, font_family="Nimbus Roman", ax=ax
    )

    # 11) Add centered species label
    ax.text(
        0.5, 0.0, f"{species_bits}",
        fontsize=26, #fontweight="bold",
        fontdict={"family": "Nimbus Roman"},
        ha="center", va="center",
        transform=ax.transAxes
    )

    ax.set_axis_off()
    #plt.tight_layout()
    #plt.show()
    plt.savefig(f"species_{species_id}.pdf", transparent=True, bbox_inches='tight', pad_inches=0.15)
    plt.close()

# --- layout helpers ---

def matrix_to_top_three_and_weights(M):
    """
    Translate a 6x6 P->P matrix of fractions into:
      - top_three_map: {donor_int: {'first': acc_int, 'second': acc_int, 'third': acc_int}}
      - edge_weights : {(donor_int, acc_int): weight_float}
    Only non-zero entries are considered.
    """
    top_three_map = {}
    edge_weights = {}

    for i in range(6):
        row = M[i, :]
        pairs = [(j+1, float(row[j])) for j in range(6) if row[j] > 0.0]
        # record weights (for scaling)
        for j, val in pairs:
            edge_weights[(i+1, j)] = val
        if not pairs:
            continue
        # choose up to top 3 acceptors
        pairs_sorted = sorted(pairs, key=lambda x: x[1], reverse=True)
        ranks = {}
        if len(pairs_sorted) >= 1:
            ranks['first'] = pairs_sorted[0][0]
        if len(pairs_sorted) >= 2:
            ranks['second'] = pairs_sorted[1][0]
        if len(pairs_sorted) >= 3:
            ranks['third'] = pairs_sorted[2][0]
        if ranks:
            top_three_map[i+1] = ranks

    return top_three_map, edge_weights

# ---------- Main ----------
def main(root="process", hbond_kind="intra", mode: str = "union", width_mode: str = "continuous", width_ref: str = "auto"):
    root = Path(root)
    for ip_dir in sorted(root.glob("IP_*")):
        hbdir = ip_dir / "analyze_final_sim" / "h_bonds"
        xpm = hbdir / f"{hbond_kind}_hb_matrix.xpm"
        log = hbdir / f"{hbond_kind}_hb.log"
        if not (xpm.is_file() and log.is_file()):
            continue
        xmat = parse_xpm_binary(xpm)
        pairs = parse_hbond_log(log)
        out_stem = str(hbdir / f"{ip_dir.name}_PtoP_matrix")

        if mode == 'both':
            # union (default naming)
            M_union = build_matrix(xmat, pairs, mode='union')
            save_matrix_csv_npy(M_union, out_stem)
            print(f"[OK] {ip_dir.name}: saved {out_stem}.npy / .csv (union)")
            # Draw diagram (union) using the original function; save into hbdir
            species_bits = ip_dir.name.replace("IP_", "")
            top3_u, ew_u = matrix_to_top_three_and_weights(M_union)
            old_cwd = os.getcwd()
            try:
                os.chdir(hbdir)
                draw_phosphorus_diagram(species_bits, f"{species_bits}_union", top3_u, edge_weights=ew_u, width_mode=width_mode, width_ref=width_ref)
            finally:
                os.chdir(old_cwd)
            # sum-of-bonds (suffix)
            M_sum = build_matrix(xmat, pairs, mode='sum')
            save_matrix_csv_npy(M_sum, out_stem + "_sum")
            print(f"[OK] {ip_dir.name}: saved {out_stem}_sum.npy / .csv (sum-of-bonds)")
            # Draw diagram (sum) using the original function; save into hbdir
            top3_s, ew_s = matrix_to_top_three_and_weights(M_sum)
            old_cwd = os.getcwd()
            try:
                os.chdir(hbdir)
                draw_phosphorus_diagram(species_bits, f"{species_bits}_sum", top3_s, edge_weights=ew_s, width_mode=width_mode, width_ref=width_ref)
            finally:
                os.chdir(old_cwd)
        else:
            M = build_matrix(xmat, pairs, mode=mode)
            suffix = "" if mode == 'union' else "_sum"
            save_matrix_csv_npy(M, out_stem + suffix)
            tag = "union" if mode == 'union' else "sum-of-bonds" # remove tag for "union"
            print(f"[OK] {ip_dir.name}: saved {out_stem}{suffix}.npy / .csv ({tag})")
            # Draw diagram for the selected mode
            species_bits = ip_dir.name.replace("IP_", "")
            top3, ew = matrix_to_top_three_and_weights(M)
            label = f"{species_bits}_{tag.replace(' ', '_')}"
            old_cwd = os.getcwd()
            try:
                os.chdir(hbdir)
                draw_phosphorus_diagram(species_bits, label, top3, edge_weights=ew, width_mode=width_mode, width_ref=width_ref)
            finally:
                os.chdir(old_cwd)

if __name__ == "__main__":
    # Usage: python make_pp_matrix.py [process_root] [intra|inter] [union|sum|both] [continuous|coarse] [auto|absolute]
    root = sys.argv[1] if len(sys.argv) > 1 else "process"
    kind = sys.argv[2] if len(sys.argv) > 2 else "intra"
    mode = sys.argv[3].lower() if len(sys.argv) > 3 else "union"
    if mode not in ("union", "sum", "both"):
        print(f"[WARN] Unknown mode '{mode}', defaulting to 'union'.")
        mode = "union"

    width_mode = sys.argv[4].lower() if len(sys.argv) > 4 else "continuous"
    if width_mode not in ("continuous", "coarse"):
        print(f"[WARN] Unknown width mode '{width_mode}', defaulting to 'continuous'.")
        width_mode = "continuous"

    width_ref = sys.argv[5].lower() if len(sys.argv) > 5 else "auto"
    if width_ref not in ("auto", "absolute"):
        print(f"[WARN] Unknown width_ref '{width_ref}', defaulting to 'auto'.")
        width_ref = "auto"

    main(root, kind, mode, width_mode, width_ref)
#!/usr/bin/env python3
"""
Combined Script for H-Bond Analysis with:
 1) A 2×2 Figure (Lifetime, Distance, Time/Atom, H-Bond Count)
 2) A Summary Plot (Occurrence bars + Existence Heatmap)
 3) A Phosphate Group Violin Plot (Donor occurrences by P1..P6)
 4) A textual table summarizing donor→acceptor phosphate usage.
    - additionally, a directed graph (network) showing donor→target relationships.

We assume IP6 microprotonation states, with donor→phosphate mapping:
  - P1: O1, O7, O8, O9
  - P2: O2, O10, O11, O12
  - P3: O3, O13, O14, O15
  - P4: O4, O16, O17, O18
  - P5: O5, O19, O20, O21
  - P6: O6, O22, O23, O24

Author: Arthur Laux
Date:   2025-01-27
"""

import os
import re
import sys
import logging
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import argparse
from electrofit.infra.logging import setup_logging

# ------------------------------------------------------------------------------
# Common Utility Functions (PLACEHOLDER)
# ------------------------------------------------------------------------------

from electrofit_analysis.structure.util.common_util import (
    load_hb_num_xvg,
    parse_xpm,
    analyze_hydrogen_bonds,
    parse_hbond_log_to_dataframe,
    count_donor_occurrences_from_matrix,
    count_oxygen_occurrences_from_matrix, 
    draw_phosphorus_diagram
)

from electrofit_analysis.viz.hbond_helpers import (
    three_by_two_plots_box_violin,
    generate_summary_plot
)

# ------------------------------------------------------------------------------
# Logging & Font Setup (configured in main)
# ------------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# Use Times New Roman + "talk" context from Seaborn
sns.set_context("talk")
sns.set_style({'font.family':'serif', 'font.serif':'Times New Roman'})

# ------------------------------------------------------------------------------
# Donor & Acceptor → Phosphate Mapping
# ------------------------------------------------------------------------------
donor_to_phosphate = {
    'O1': 'P1','O7': 'P1','O8': 'P1','O9': 'P1',
    'O2': 'P2','O10': 'P2','O11': 'P2','O12': 'P2',
    'O3': 'P3','O13': 'P3','O14': 'P3','O15': 'P3',
    'O4': 'P4','O16': 'P4','O17': 'P4','O18': 'P4',
    'O5': 'P5','O19': 'P5','O20': 'P5','O21': 'P5',
    'O6': 'P6','O22': 'P6','O23': 'P6','O24': 'P6'
}
acceptor_to_phosphate = donor_to_phosphate  # same mapping if acceptors are in [O1..O24]

# ------------------------------------------------------------------------------
# Textual Table: Donor→Acceptor usage
# ------------------------------------------------------------------------------

def build_donor_acceptor_summary(species_id, bond_df, hbonds_per_index, time_per_frame_ns, donor_map, acceptor_map):
    """
    For each bond => map donor->pgroup, acceptor->pgroup => accumulate times
    """
    result = {}
    hbonds_length = len(hbonds_per_index)  # Cache length for efficiency
    for i, row in bond_df.iterrows():
        b_idx = row['idx']
        d_atom = row['donor']
        a_atom = row['acceptor']
        
        # Safe access to hbonds_per_index
        if 0 <= b_idx < hbonds_length:
            frames = hbonds_per_index[b_idx]
        else:
            frames = 0  # Default to 0 if index is out of bounds
        
        engaged_ns = frames * time_per_frame_ns

        donor_pgroup = donor_map.get(d_atom, None)
        acceptor_pgroup = acceptor_map.get(a_atom, None)
        if donor_pgroup is None or acceptor_pgroup is None:
            continue

        if donor_pgroup not in result:
            result[donor_pgroup] = {"sum_donor_time": 0.0, "targets": {}}
        result[donor_pgroup]["sum_donor_time"] += engaged_ns

        key_t = (acceptor_pgroup, a_atom)
        if key_t not in result[donor_pgroup]["targets"]:
            result[donor_pgroup]["targets"][key_t] = 0.0
        result[donor_pgroup]["targets"][key_t] += engaged_ns
    return result

def print_donor_acceptor_table(species_id, data, logger, summary_format='arrow', draw_figure=True):
    """
    Generate a summary table for each donor pgroup, log it, and (optionally) draw a figure
    for the top 3 phosphorus-target relationships (1st=red, 2nd=black, 3rd=blue).
    Now we also produce separate textual summaries for first, second, and third P-targets.

    Parameters:
    - species_id (str): Identifier for the species (e.g. "101111").
    - data (dict): Summary data from build_donor_acceptor_summary.
    - logger (logging.Logger): Configured logger instance.
    - summary_format (str): Format of the summary line ('arrow', 'graphical', etc.).
    - draw_figure (bool): If True, generate a color-coded figure for each species.
    """
    table_lines = []
    table_lines.append(f"Species: {species_id}\n")

    # We'll store up to 3 P-targets for each donor in top_three_map => used by the figure
    # Also store them in rank_p_targets => used by textual summary
    top_three_map   = {}  # e.g. {1: {"first":2, "second":4, "third":5}, ...}
    rank_p_targets  = { "first": {}, "second": {}, "third": {} }  
    # e.g. rank_p_targets["first"] = { "P2":"P3", "P4":"P5", ... }

    donor_sorted = sorted(data.keys())  # e.g. ['P1','P2','P3','P4','P5','P6'] if they exist
    for dpgroup in donor_sorted:
        sum_time = data[dpgroup]["sum_donor_time"]
        table_lines.append(f"  Donor {dpgroup} => total time: {sum_time:.2f} ns")
        
        # Sort all acceptor targets by descending time
        t_map = data[dpgroup]["targets"]  # e.g. {(acc_pg, acc_atom): time_ns, ...}
        t_list = sorted(t_map.items(), key=lambda x: x[1], reverse=True)

        if not t_list:
            table_lines.append("    No acceptor targets found.")
            table_lines.append("    No Phosphorus acceptor targets found.\n")
            continue
        
        # Print each target with rank (1st, 2nd, 3rd, etc.) for logging
        for i, ((acc_pg, acc_atom), val_ns) in enumerate(t_list, start=1):
            # ordinal suffix
            if 10 <= i % 100 <= 20:
                suffix = 'th'
            else:
                suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(i % 10, 'th')
            target_label = f"{i}{suffix} target"
            percentage = (val_ns / sum_time * 100) if sum_time > 0 else 0.0
            table_lines.append(f"    {target_label:15s} => {acc_pg}({acc_atom}) = {percentage:.1f} %")

        # Aggregate phosphorus targets only
        p_group_summary = {}
        for ((acc_pg, acc_atom), val_ns) in t_list:
            if acc_pg.startswith('P'):
                p_group_summary[acc_pg] = p_group_summary.get(acc_pg, 0.0) + val_ns

        if p_group_summary:
            # Sort P-targets by time, descending
            p_targets_sorted = sorted(p_group_summary.items(), key=lambda x: x[1], reverse=True)
            table_lines.append("    Phosphorus Targets:")

            # Donor numeric index (P3 => 3)
            donor_num = int(dpgroup.replace('P',''))
            top_three_map[donor_num] = {}

            # We’ll store up to 3 in top_three_map: 'first','second','third'
            rank_to_key = {1: 'first', 2: 'second', 3: 'third'}

            for i, (acc_pg, total_ns) in enumerate(p_targets_sorted, start=1):
                if 10 <= i % 100 <= 20:
                    suffix = 'th'
                else:
                    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(i % 10, 'th')
                p_target_label = f"{i}{suffix} P target"

                p_percentage = (total_ns / sum_time * 100) if sum_time > 0 else 0.0
                table_lines.append(f"      {p_target_label:15s} => {acc_pg} = {p_percentage:.1f} %")

                # If i <= 3 => store in top_three_map for the figure
                if i <= 3:
                    top_three_map[donor_num][rank_to_key[i]] = int(acc_pg.replace('P',''))

                # Also store in rank_p_targets => for the textual summary
                if i <= 3:
                    # e.g. rank_p_targets["first"][dpgroup] = acc_pg
                    rank_str = rank_to_key[i]
                    rank_p_targets[rank_str][dpgroup] = acc_pg
        else:
            table_lines.append("    No Phosphorus acceptor targets found.")
        
        table_lines.append("")  # spacing

    # Now we produce three separate summaries (first, second, third)
    # e.g. "Summary of First Phosphorus Targets:\n  P2 -> P3 | P4 -> P5 | ..."
    # If any donors have a first, second, third P-target.

    # Helper to produce lines like "P2 -> P3 | P4 -> P5 | ..."
    def build_summary_line(mapping):
        # mapping e.g. rank_p_targets["first"] => { "P2":"P3", "P4":"P5", ... }
        pairs = []
        for donor, first_p in sorted(mapping.items()):
            pairs.append(f"{donor} -> {first_p}")
        return " | ".join(pairs)

    # 1) First
    if rank_p_targets["first"]:
        line_str = build_summary_line(rank_p_targets["first"])
        table_lines.append("Summary of First Phosphorus Targets:")
        table_lines.append(f"  {line_str}\n")

    # 2) Second
    if rank_p_targets["second"]:
        line_str = build_summary_line(rank_p_targets["second"])
        table_lines.append("Summary of Second Phosphorus Targets:")
        table_lines.append(f"  {line_str}\n")

    # 3) Third
    if rank_p_targets["third"]:
        line_str = build_summary_line(rank_p_targets["third"])
        table_lines.append("Summary of Third Phosphorus Targets:")
        table_lines.append(f"  {line_str}\n")

    table_lines.append("-" * 70)
    table_lines.append("")

    # Combine lines + log
    table_str = "\n".join(table_lines)
    logger.info(table_str)

    # Finally, if desired, draw the figure
    if draw_figure:
        # Build a simple phosphate→phosphate weight map in ns
        edge_weights = {}
        for donor_pg, info in data.items():
            u = int(donor_pg.replace('P',''))
            for (acceptor_pg, _atom), t_ns in info['targets'].items():
                if acceptor_pg.startswith('P'):
                    v = int(acceptor_pg.replace('P',''))
                    edge_weights[(u,v)] = edge_weights.get((u,v), 0.0) + t_ns
        # Now call the figure, passing this map:
        draw_phosphorus_diagram(species_id, top_three_map, edge_weights)
        #draw_phosphorus_diagram(species_id, top_three_map)

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main(hbond_type: str, project_root: str) -> None:
    """Run H-bond comparison analysis across all species in <project_root>/process.

    Parameters
    ----------
    hbond_type : str
        The H-bond analysis type prefix used for file discovery (e.g., 'intra').
    project_root : str
        Path to the project root containing a 'process' directory.
    """

   
    # 1) Find process directory
    process_dir = os.path.join(project_root, "process")
    log_path = os.path.join(process_dir, "hbond_comparison.log")
    # Configure logging for this run
    setup_logging(log_path=log_path, also_console=True)
    logger.info("Starting combined analysis with 3x2 figure, summary plot, phosphate violin, textual table...")

    if not os.path.isdir(process_dir):
        logger.error(f"Invalid process directory => {process_dir}")
        sys.exit(1)
    logger.info(f"Using process directory => {process_dir}")

    process_path = Path(process_dir)
    time_per_frame_ns = 0.01

    # 2) Lists/dicts for data accumulation
    folder_list = []
    all_species_ids = []
    hbonds_data = []      # for 2x2 hbond count
    lifetimes_data = []   # for 2x2 lifetime
    distance_data = []    # for 2x2 distance
    engaged_data = []     # for 2x2 time/atom
    df_phosphate_rows = []  # for phosphate group violin
    donor_atoms_by_species = {}
    occurrence_data = {}
    existence_data = {}
    donor_acceptor_summaries = {}

    # 3) Identify species subfolders
    for fname in sorted(os.listdir(process_path)):
        if not fname.startswith("IP_"):
            continue
        folder_path = process_path / fname
        if not folder_path.is_dir():
            continue
        species_id = fname.replace("IP_", "")
        n_ones = species_id.count('1')
        if n_ones not in [3,4,5]:
            continue

        folder_list.append((fname, species_id, n_ones))
        all_species_ids.append(species_id)

    # sort 5->4->3
    folder_list.sort(key=lambda x: (-x[2], x[1]))
    sorted_species_list = [x[1] for x in folder_list]

    # 4) Loop over species => parse data
    for folder_name, species_id, n_ones in folder_list:
        folder_path = process_path / folder_name

        # Required files
        hb_num_file = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_num.xvg"
        xpm_file    = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_matrix.xpm"
        dist_file   = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_dist.xvg"
        log_file    = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb.log"

        # Check existence
        if not (hb_num_file.is_file() and xpm_file.is_file() and dist_file.is_file() and log_file.is_file()):
            logger.warning(f"Skipping species {species_id}: missing hbond files.")
            continue

        # (a) H-bond count vs time
        data_num = load_hb_num_xvg(str(hb_num_file))
        if data_num.size<2:
            logger.warning(f"No hbond count data => {species_id}")
            continue
        hbond_counts = data_num[:,1]
        for val in hbond_counts:
            hbonds_data.append((species_id, n_ones, val))

        # (b) XPM => analyze => lifetimes
        xpm_matrix, meta = parse_xpm(str(xpm_file))
        analysis_res = analyze_hydrogen_bonds(xpm_matrix, meta)
        lf_frames = [lf for bond_lf in analysis_res['lifetimes'] for lf in bond_lf]
        lf_ns = [f*time_per_frame_ns for f in lf_frames]
        for val_ns in lf_ns:
            lifetimes_data.append((species_id, n_ones, val_ns))

        # (c) Distances
        d_list = []
        with open(dist_file, 'r') as df_:
            for line in df_:
                if not line.strip() or line.startswith(('#','@')):
                    continue
                parts = line.split()
                if len(parts)==2:
                    try:
                        dx=float(parts[0])
                        dy=float(parts[1])
                        d_list.append((dx,dy))
                    except ValueError:
                        pass
        arr_d = np.array(d_list)
        if arr_d.size<2:
            logger.warning(f"No distance data => {species_id}")
            continue
        arr_d = arr_d[arr_d[:,0]>=0.2]
        for row_val in arr_d:
            distance_data.append((species_id,n_ones,row_val[0], row_val[1]))

        # (d) bond_df => parse log
        bond_df = parse_hbond_log_to_dataframe(str(log_file))
        if bond_df.empty:
            logger.warning(f"No hbond pairs => {species_id}")
            continue
        # store donor sets
        donor_atoms_by_species[species_id] = set(bond_df['donor'].unique())

        # from analysis => frames/bond
        hpi = analysis_res['hbonds_per_index']
        # donor-based sum => time/atom
        donor_counts = count_donor_occurrences_from_matrix(bond_df, hpi)
        total_donor_frames = donor_counts.sum()
        engaged_ns = (total_donor_frames * time_per_frame_ns)/n_ones
        engaged_data.append((species_id, n_ones, engaged_ns))

        # Build phosphate group data => for violin
        for donor_atom, frames in donor_counts.items():
            engaged_time_ns = frames*time_per_frame_ns
            pgroup = donor_to_phosphate.get(donor_atom, None)
            if pgroup:
                df_phosphate_rows.append({
                    'species': species_id,
                    'num_ones': n_ones,
                    'donor_atom': donor_atom,
                    'phosphate_group': pgroup,
                    'engaged_ns': engaged_time_ns
                })

        # Build donor→acceptor summary => textual table
        summ_data = build_donor_acceptor_summary(
            species_id, bond_df, hpi, time_per_frame_ns,
            donor_to_phosphate, acceptor_to_phosphate
        )
        donor_acceptor_summaries[species_id] = summ_data

    # 5) Build DataFrames for the 2×2 figure
    df_hbonds    = pd.DataFrame(hbonds_data, columns=["species","num_ones","hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species","num_ones","lifetime_ns"])
    df_distance  = pd.DataFrame(distance_data, columns=["species","num_ones","distance_nm","frequency"])
    df_engaged   = pd.DataFrame(engaged_data,  columns=["species","num_ones","engaged_time_ns"])

    # Build color_dict
    all_species_ids = sorted(set(all_species_ids))
    palette = sns.color_palette("tab10", len(all_species_ids))
    color_dict = {}
    for i, sp in enumerate(all_species_ids):
        color_dict[sp] = palette[i % len(palette)]

    # Build a DataFrame for the phosphate group violin
    df_phosphate = pd.DataFrame(df_phosphate_rows)
    # 6) Optionally create the 3×2 figure
    three_by_two_plots_box_violin(
        df_lifetimes=df_lifetimes,
        df_distance=df_distance,
        df_engaged=df_engaged,
        df_hbonds=df_hbonds,
        df_phosphate=df_phosphate,
        color_dict=color_dict,
        output_prefix=f"{hbond_type}_hbonds",
        figure_size=(16,18)
    )

    # 7) Create data for summary existence plot
    occurrence_data.clear()
    existence_data.clear()

    # We'll also collect all oxygens across these species
    all_oxys = set()
    for fname, species_id, n_ones in folder_list:
        log_path = process_path / fname / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb.log"
        if not log_path.is_file():
            continue
        bdf = parse_hbond_log_to_dataframe(str(log_path))
        if bdf.empty:
            continue
        all_oxys.update(bdf['donor'].unique())
        all_oxys.update(bdf['acceptor'].unique())

    for fname, species_id, n_ones in folder_list:
        fpath = process_path / fname
        xpm_file = fpath / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_matrix.xpm"
        log_file = fpath / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb.log"
        if not (xpm_file.is_file() and log_file.is_file()):
            continue
        bond_df = parse_hbond_log_to_dataframe(str(log_file))
        if bond_df.empty:
            continue
        mat, meta = parse_xpm(str(xpm_file))
        an = analyze_hydrogen_bonds(mat, meta)
        hpi = an['hbonds_per_index']

        # Sum donor+acceptor → occurrence_data
        from_both = count_oxygen_occurrences_from_matrix(bond_df, hpi)
        occurrence_data[species_id] = from_both

        # Build existence map
        bin_size_ns = 0.2
        frames_per_bin = int(bin_size_ns / time_per_frame_ns)
        if frames_per_bin < 1:
            continue
        n_frames = mat.shape[1]
        nbins = n_frames // frames_per_bin
        if n_frames % frames_per_bin != 0:
            nbins += 1

        binned = np.zeros((mat.shape[0], nbins))
        for i in range(nbins):
            st = i * frames_per_bin
            en = min((i + 1) * frames_per_bin, n_frames)
            chunk = mat[:, st:en]
            binned[:, i] = np.mean(chunk, axis=1)

        # Build a mapping from oxygen to bond indices
        oxy_map = {}
        for idx, row in bond_df.iterrows():
            d = row['donor']
            a = row['acceptor']
            oxy_map.setdefault(d, []).append(idx)
            oxy_map.setdefault(a, []).append(idx)

        # Filter out "OW1" from the global oxygen list
        sorted_oxy = sorted([o for o in all_oxys if o != "OW1"],
                            key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)
        
        agg = np.zeros((len(sorted_oxy), nbins))
        for o_i, oxy in enumerate(sorted_oxy):
            b_list = oxy_map.get(oxy, [])
            if b_list:
                for bn in range(nbins):
                    agg[o_i, bn] = np.max(binned[b_list, bn])
        existence_data[species_id] = agg

    
    # 8) Once we have occurrence_data + existence_data:
    # This is time consuming...
    generate_summary_plot(
        occurrence_data=occurrence_data,
        existence_data=existence_data,
        donor_atoms_by_species=donor_atoms_by_species,
        color_dict=color_dict,
        time_per_frame_ns=time_per_frame_ns,
        output_file=f"{hbond_type}_oxygen_occurrences_summary.pdf",
        folder_order=sorted_species_list
    )

    # 10) Print textual table
    for sp, data in donor_acceptor_summaries.items():
        print_donor_acceptor_table(sp, data, logger, "arrow")

    logger.info("All analysis steps completed successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="H-bond comparison analysis: generate 2x2 figure, summary plots, and donor→acceptor tables."
    )
    parser.add_argument(
        "--project",
        "-p",
        required=True,
        help="Path to the project root directory (must contain 'process').",
    )
    parser.add_argument(
        "--type",
        "-t",
        default="intra",
        help="H-bond analysis type prefix used in filenames (default: 'intra').",
    )
    args = parser.parse_args()
    main(args.type, os.path.abspath(args.project))
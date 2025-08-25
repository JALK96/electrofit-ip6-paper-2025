#!/usr/bin/env python3
import re
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from electrofit_analysis.structure.util.common_util import (
    expand_phosphate_data_for_violins,
)

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
# 1) 2×2 Figure (Lifetime, Distance, Time/Atom, H-Bond Count)
# ------------------------------------------------------------------------------
def three_by_two_plots_box_violin(
    df_lifetimes,   
    df_distance,    
    df_engaged,     
    df_hbonds,      
    df_phosphate,   # [species, phosphate_group, engaged_ns, ...]
    color_dict=None,
    output_prefix="intra_hbonds",
    figure_size=(16, 22),   
    spacer_prefix="~space_"
):
    """
    Creates a 3×2 layout:

      Row 0: (0,0) => (a) Lifetime (Box)  
              (0,1) => (b) Distance (Box)

      Row 1: (1,0) => (c) Time/Atom (Stacked)
              (1,1) => (d) H-bond Count (Violin)

      Row 2: (2,0) => (e) Horizontal Phosphate Violin
              (2,1) => Empty (white)

    We do not alter the 4 subplots from your 2×2 figure, 
    but add a 5th plot (phosphate violin) in row=2, col=0.
    Spacers are injected so all y-axis labels match across the subplots.
    """

    logger.info("Generating a 3×2 figure: the original 4 subplots + a 5th horizontal violin below them.")

    # ------------------------------------------------------
    # 0) Check empties for the top-4 data
    # ------------------------------------------------------
    all_empty = True
    for df_check in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df_check is not None and not df_check.empty:
            all_empty = False
            break
    if all_empty:
        logger.warning("All top-4 DataFrames empty => skipping top-4 plots.")
        # We still allow the 5th plot if df_phosphate is non-empty

    # ------------------------------------------------------
    # 1) Gather species from the top-4 DataFrames
    # ------------------------------------------------------
    species_in_any = set()
    for df_check in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df_check is not None and not df_check.empty:
            species_in_any.update(df_check["species"].unique())

    # Helper: figure out how many '1' bits for a species
    def get_num_ones(sp, df_list):
        for d in df_list:
            if d is not None and not d.empty:
                row = d.loc[d["species"] == sp]
                if not row.empty:
                    return row.iloc[0]["num_ones"]
        return None

    # Keep only species with 3,4,5
    all_with_num = []
    for sp in species_in_any:
        val = get_num_ones(sp, [df_lifetimes, df_distance, df_engaged, df_hbonds])
        if val in [3, 4, 5]:
            all_with_num.append((sp, val))

    if not all_with_num:
        logger.warning("No species with num_ones in [3,4,5]. Skipping all top-4 plots + 5th.")
        return

    # Sort => 5->top, 4->middle, 3->bottom
    all_with_num.sort(key=lambda x: (-x[1], x[0]))

    # Insert spacers for transitions (5->4->3)
    grouped_species_order = []
    prev_n_ones = None
    for i, (sp, n_ones) in enumerate(all_with_num):
        if i > 0 and n_ones != prev_n_ones:
            grouped_species_order.append((f"{spacer_prefix}{prev_n_ones}to{n_ones}", None))
        grouped_species_order.append((sp, n_ones))
        prev_n_ones = n_ones
    # final list => species + placeholders
    full_order = [t[0] for t in grouped_species_order]

    # A helper to inject spacers into a DataFrame
    def inject_spacers(df, measure_col):
        """
        For each spacer in grouped_species_order, add a row with measure_col=NaN.
        """
        if df is None or df.empty:
            return pd.DataFrame(columns=["species", measure_col, "num_ones"])
        new_df = df.copy()
        row_spacers = []
        for sp, val in grouped_species_order:
            if val is None:  # a dummy spacer
                row_spacers.append({
                    "species": sp,
                    measure_col: np.nan,
                    "num_ones": -1
                })
        if row_spacers:
            dummy_df = pd.DataFrame(row_spacers)
            combined = pd.concat([new_df, dummy_df], ignore_index=True)
            return combined
        return new_df

    # ------------------------------------------------------
    # 2) For df_distance, expand by frequency if needed
    # ------------------------------------------------------
    if df_distance is not None and not df_distance.empty:
        df_distance = df_distance.copy()
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        df_dist_expanded = df_distance.loc[df_distance.index.repeat(df_distance['frequency'])].reset_index(drop=True)
    else:
        df_dist_expanded = pd.DataFrame()

    # Now inject spacers for each measure
    df_life2  = inject_spacers(df_lifetimes, "lifetime_ns")
    df_dist2  = inject_spacers(df_dist_expanded, "distance_nm")
    df_eng2   = inject_spacers(df_engaged, "engaged_time_ns")
    df_hbond2 = inject_spacers(df_hbonds, "hbond_count")

    # ------------------------------------------------------
    # 3) Summarize df_phosphate => row=species, col=P1..P6, so we can do time/atom stacked
    #    (unchanged from your code)
    # ------------------------------------------------------
    grouped_phos = df_phosphate.groupby(['species','phosphate_group'])['engaged_ns'].sum().reset_index()
    pivot_phos = grouped_phos.pivot(index='species', columns='phosphate_group', values='engaged_ns').fillna(0)

    phosphate_cols = ['P1','P2','P3','P4','P5','P6']
    for c in phosphate_cols:
        if c not in pivot_phos.columns:
            pivot_phos[c] = 0.0
    pivot_phos = pivot_phos[phosphate_cols]

    # If no color_dict => create an empty one
    if color_dict is None:
        color_dict = {}
    else:
        color_dict = dict(color_dict)  # copy to avoid mutating the original

    # Ensure spacers have a color
    spacer_color = "#D3D3D3"
    for sp, val in grouped_species_order:
        if val is None and sp not in color_dict:
            color_dict[sp] = spacer_color

    # ------------------------------------------------------
    # 4) Create a 3×2 figure 
    # ------------------------------------------------------
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=figure_size)

    # Ax references
    ax_life   = axes[0][0]  # (a)
    ax_dist   = axes[0][1]  # (b)
    ax_stack  = axes[1][0]  # (c)
    ax_hbond  = axes[1][1]  # (d)
    ax_violin = axes[2][0]  # (e)
    ax_blank  = axes[2][1]  
    ax_blank.axis('off')   # White cell

    # ~~~~~~~~~~~~~ (a) Lifetime (Box) ~~~~~~~~~~~~~
    if df_life2 is not None and not df_life2.empty:
        sns.boxplot(
            data=df_life2,
            x="lifetime_ns",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            showmeans=True,
            meanprops={'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8},
            showfliers=True,
            ax=ax_life
        )
        ax_life.set_title("(a) Lifetime (Box)")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("Species")
        ax_life.set_xscale("log")

        # Hide spacer row labels
        new_lbls = []
        for lbl in ax_life.get_yticklabels():
            txt = lbl.get_text()
            if txt.startswith(spacer_prefix):
                new_lbls.append("")
            else:
                new_lbls.append(txt)
        ax_life.set_yticklabels(new_lbls)
    else:
        ax_life.set_title("No Lifetime Data")
        ax_life.set_xlabel("")
        ax_life.set_ylabel("")

    # ~~~~~~~~~~~~~ (b) Distance (Box) ~~~~~~~~~~~~~
    if df_dist2 is not None and not df_dist2.empty:
        sns.boxplot(
            data=df_dist2,
            x="distance_nm",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            showmeans=True,
            meanprops={'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8},
            showfliers=False,
            ax=ax_dist
        )
        ax_dist.set_title("(b) Distance (Box)")
        ax_dist.set_xlabel("Distance (nm)")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])
    else:
        ax_dist.set_title("No Distance Data")
        ax_dist.set_xlabel("")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])

    # Reference df_eng2 to avoid unused-variable lint when top-4 are empty
    _ = df_eng2
    # ~~~~~~~~~~~~~ (d) H-bond Count (Violin) ~~~~~~~~~~~~~
    if df_hbond2 is not None and not df_hbond2.empty:
        sns.violinplot(
            data=df_hbond2,
            x="hbond_count",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            cut=0,
            inner=None,
            ax=ax_hbond
        )
        ax_hbond.set_title("(d) H-bond Count (Violin)")
        ax_hbond.set_xlabel("H-bond Count")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

        # Overlay means
        df_means = df_hbond2.groupby('species')['hbond_count'].mean().reset_index()
        df_means = df_means[~df_means['species'].str.startswith(spacer_prefix)]
        ax_hbond.scatter(
            df_means['hbond_count'],
            df_means['species'],
            color='black',
            marker='o',
            s=80,
            label='Mean',
            zorder=5
        )
    else:
        ax_hbond.set_title("No H-bond Data")
        ax_hbond.set_xlabel("")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

    # ~~~~~~~~~~~~~ (c) Time/Atom (Stacked) ~~~~~~~~~~~~~
    ax_stack.set_title("(c) Time / H-Atom")
    ax_stack.set_xlabel("Engaged Time (ns)")
    ax_stack.set_ylabel("Species")

    # We'll map each species/spacer -> a y-position
    y_positions = {}
    real_index = 0
    for sp, val in grouped_species_order:
        y_positions[sp] = real_index
        real_index += 1

    # For each species, sum of engaged time in pivot_phos => stacked bars
    for sp in y_positions:
        if sp not in pivot_phos.index:
            continue
        species_color = color_dict.get(sp, "gray")  
        row_vals = pivot_phos.loc[sp]
        left_val = 0.0
        for pgroup in ['P1','P2','P3','P4','P5','P6']:
            seg_val = row_vals[pgroup]
            if seg_val>0:
                ax_stack.barh(
                    y=y_positions[sp],
                    width=seg_val,
                    left=left_val,
                    height=0.8,
                    color=species_color, 
                    edgecolor='k',
                    alpha=0.9
                )
                if seg_val>12.5:
                    ax_stack.annotate(
                        f"{seg_val:.2f}",
                        (left_val + seg_val/2, y_positions[sp]),
                        ha='center', va='center',
                        color='white', fontsize=10
                    )
                left_val += seg_val

    ax_stack.set_ylim(-0.5, real_index-0.5)
    ax_stack.set_yticks(np.arange(real_index))
    # Build final y labels
    y_lbls_stack = []
    for sp, val in grouped_species_order:
        if val is None:
            y_lbls_stack.append("")
        else:
            y_lbls_stack.append(sp)
    ax_stack.set_yticklabels(y_lbls_stack)
    ax_stack.invert_yaxis()

    # ~~~~~~~~~~~~~ (e) Horizontal Phosphate Violin ~~~~~~~~~~~~~
    # We want to inject spacers for the expanded data, so the y-axis matches top plots
    expanded_df = expand_phosphate_data_for_violins(df_phosphate)
    if expanded_df.empty:
        logger.warning("No phosphate-group donor data => skipping subplot (e).")
        ax_violin.set_title("(e) No Phosphate Data")
    else:
        # We'll treat 'pgroup_index' as measure column so we can inject spacers
        # though they won't have pgroup_index in the spacer row. It's still okay.
        # We'll create a new column measure for injection
        expanded_for_injection = expanded_df.copy()
        expanded_for_injection['pgroup_index'] = expanded_for_injection['pgroup_index'].astype(float)

        # Now we do the same injection approach
        # We'll rename measure_col to 'pgroup_index'
        def inject_violin_spacers(df_in):
            if df_in is None or df_in.empty:
                return pd.DataFrame(columns=["species", "pgroup_index"])
            new_df = df_in.copy()
            row_spacers = []
            for sp, val in grouped_species_order:
                if val is None:  
                    row_spacers.append({
                        "species": sp,
                        "pgroup_index": np.nan
                    })
            if row_spacers:
                dummy_df = pd.DataFrame(row_spacers)
                return pd.concat([new_df, dummy_df], ignore_index=True)
            return new_df

        expanded_with_spacers = inject_violin_spacers(expanded_for_injection)

        # Now we can do a horizontal violin
        sns.violinplot(
            data=expanded_with_spacers,
            x='pgroup_index',
            y='species',
            hue='species',
            order=full_order,
            palette=color_dict,
            orient='h',
            cut=0,
            inner=None,
            ax=ax_violin
        )
        ax_violin.set_title("(e) Phosphate Group Horizontal Violin")

        # x-axis ticks => 0..5 => P1..P6
        phosphate_labels = ['P1','P2','P3','P4','P5','P6']
        ax_violin.set_xticks(range(len(phosphate_labels)))
        ax_violin.set_xticklabels(phosphate_labels)
        ax_violin.set_xlabel("Phosphate Group")
        ax_violin.set_ylabel("Species")

        # We'll create the same y ticks as the other subplots => full_order
        ax_violin.set_yticks(np.arange(real_index))
        y_lbls_violin = []
        for sp, val in grouped_species_order:
            if val is None:
                y_lbls_violin.append("")
            else:
                y_lbls_violin.append(sp)
        ax_violin.set_yticklabels(y_lbls_violin)

        # If the legend is large, you may remove or place it differently
        # ax_violin.legend(loc='best')

    # ------------------------------------------------------
    # 5) Final
    # ------------------------------------------------------
    plt.tight_layout()
    out_pdf = f"{output_prefix}_3x2_with_violin.pdf"
    plt.savefig(out_pdf, dpi=300)
    plt.close()
    logger.info(f"Saved 3×2 figure with the 5th plot => '{out_pdf}'!")

# ------------------------------------------------------------------------------
# 2) Summary Plot (Occurrence bars + Existence Heatmap)
# ------------------------------------------------------------------------------
def generate_summary_plot(
    occurrence_data,
    existence_data,
    donor_atoms_by_species=None,
    color_dict=None,
    time_per_frame_ns=0.01,
    output_file='oxygen_occurrences_summary.pdf',
    folder_order=None
):
    """
    Snippet-like summary plot:
      Left: barplot of occurrence times (reversed x-axis),
      Right: existence heatmap.

    If donor_atoms_by_species is given, any oxygen in that set is labeled red on y-axis.
    """
    logger.info("Generating summary plot with existence maps & occurrence bars.")
    num_species = len(occurrence_data)
    if num_species == 0:
        logger.warning("No data in occurrence_data; skipping summary plot.")
        return

    if donor_atoms_by_species is None:
        donor_atoms_by_species = {}

    # ----------------------------------------------------------------------
    # 1) Remove 'OW1' from each species' occurrence_data
    # ----------------------------------------------------------------------
    for sp in list(occurrence_data.keys()):
        ser = occurrence_data[sp]
        if "OW1" in ser.index:
            print(ser)
            print(ser.index)
            occurrence_data[sp] = ser.drop("OW1")

    # Gather all possible oxygens
    all_oxygens = set()
    for ser in occurrence_data.values():
        all_oxygens.update(ser.index)
    all_oxygens = sorted(
        all_oxygens,
        key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0
    )

    if folder_order:
        sorted_species_ids = [sp for sp in folder_order if sp in occurrence_data]
    else:
        sorted_species_ids = sorted(occurrence_data.keys())

    max_occ = max([s.max() for s in occurrence_data.values()]) if occurrence_data else 0
    x_limit = max_occ * time_per_frame_ns * 1.1

    fig_height = max(6, 5.5*num_species)
    fig, axes = plt.subplots(
        nrows=num_species, ncols=2,
        figsize=(16, fig_height),
        constrained_layout=True,
        gridspec_kw={'width_ratios':[1,1]}
    )
    if num_species == 1:
        axes = [axes]

    for i, sp in enumerate(sorted_species_ids):
        ax_occ, ax_exist = axes[i]
        occ_series = occurrence_data[sp].reindex(all_oxygens, fill_value=0)
        time_vals = occ_series * time_per_frame_ns

        bar_color = color_dict.get(sp, "gray") if color_dict else "gray"

        # Left subplot: bars
        sns.barplot(
            x=time_vals.values,
            y=all_oxygens,
            color=bar_color,
            ax=ax_occ
        )
        ax_occ.set_title(f"H-Bond Occurrence Time - {sp}", fontweight='bold')
        ax_occ.set_xlabel("Time (ns)")
        ax_occ.set_xlim(x_limit, 0)

        for p in ax_occ.patches:
            width = p.get_width()
            if width>0:
                ax_occ.annotate(
                    f"{width:.2f}",
                    (width, p.get_y() + p.get_height()/2),
                    ha='right', va='center',
                    xytext=(-5,0),
                    textcoords='offset points',
                    fontsize=13
                )
        ax_occ.grid(True, linestyle='--', alpha=0.5, axis='x')
        ax_occ.set_yticklabels([])
        ax_occ.yaxis.set_ticks_position('right')

        # Right subplot: existence heatmap
        mat = existence_data[sp]
        if mat.shape[0] != len(all_oxygens):
            logger.warning(f"Mismatch for {sp}: existence map {mat.shape[0]} rows vs {len(all_oxygens)} oxygens.")
            ax_exist.set_visible(False)
            continue

        sns.heatmap(
            mat, cmap="Reds", cbar=True, ax=ax_exist,
            linewidths=0, linecolor='white',
            vmin=0, vmax=1,
            cbar_kws={"label":"Fraction of Bond Presence"}
        )
        ax_exist.set_title(f"H-Bond Existence Map - {sp}", fontweight='bold')
        ax_exist.set_xlabel("Time (ns)")

        nbins = mat.shape[1]
        bin_size_ns = 0.2
        tvals = np.arange(nbins)*bin_size_ns
        ticks_count = min(6, nbins)
        tick_positions = np.linspace(0, nbins-1, ticks_count, dtype=int)
        tick_labels = [f"{int(tvals[pos])}" for pos in tick_positions]

        ax_exist.set_xticks(tick_positions + 0.5)
        ax_exist.set_xticklabels(tick_labels, rotation=0, ha='right')
        ax_exist.set_yticks(np.arange(len(all_oxygens)) + 0.5)
        ax_exist.set_yticklabels(all_oxygens, rotation=0)

        # Color donor y-labels red
        donorset = donor_atoms_by_species.get(sp, set())
        for txt_label in ax_exist.get_yticklabels():
            oxy_txt = txt_label.get_text()
            if oxy_txt in donorset:
                txt_label.set_color("red")

        for spine in ax_exist.spines.values():
            spine.set_visible(True)

    plt.savefig(output_file, dpi=300)
    plt.close()
    logger.info(f"Summary plot saved as '{output_file}'")
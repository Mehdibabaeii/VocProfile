#!/usr/bin/env python3
#
# VocProfile - 12_plot_filled_distributions.py
#
# Violin and box plot of every compound across accessions after the accession means
# without detection are filled.
#
# in   MEANS_FILE, DATASETS, SUBCLASS_ORDER_TABLE, COLOUR_TABLE
# out  figures/12_plot_filled_distributions/distributions_filled_<dataset>/violin_box_filled_<dataset>.png
#      results/12_plot_filled_distributions/12_qc.txt
#
# Usage:
#   python scripts/12_plot_filled_distributions.py <config.txt>
#
# Notes: docs/NOTES.md, step 12
#
import gc
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, save_figure

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')

cfg, step, args = load_config("12")
start_step(cfg, 'means')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'accession means without detection: set to NONDETECT_FRACTION of the smallest positive accession mean of the compound',
)
ABUNDANCE_SHEET = cfg["abundance_sheet"]
DPI = step["dpi"]


def build_palette(df_cols):
    """Subclass and main-class colours, same order as in step 2."""
    all_sg = list(cfg["subclass_order"]["distributions"])
    present = list(dict.fromkeys(df_cols.get_level_values('Subclass')))
    unique_sg = [sg for sg in all_sg if sg in present]
    for sg in present:
        if sg not in unique_sg:
            unique_sg.append(sg)
    colors = sns.color_palette("husl", len(unique_sg))
    sg_pal = dict(zip(unique_sg, colors))
    mc_pal = {}
    for mc in list(dict.fromkeys(df_cols.get_level_values('Main_Class'))):
        first_sg = next((sg for mc2, sg, _ in df_cols if mc2 == mc), None)
        mc_pal[mc] = sg_pal.get(first_sg, (0.5, 0.5, 0.5))
    return sg_pal, mc_pal, unique_sg


file_path = input_path(cfg, step["input"])
sheets_to_run = step["sheets"]
require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets_to_run])
output_folder_base = os.path.join(output_dir(cfg, "12_plot_filled_distributions"), step["output_folder"])

for sheet_name in sheets_to_run:
    print(f"Dataset: {sheet_name}")

    output_folder = f"{output_folder_base}_{sheet_name}"
    os.makedirs(output_folder, exist_ok=True)
    unit_label = "Abundance of VOCs" if sheet_name == ABUNDANCE_SHEET else "Estimated Concentration (µg/g)"

    try:
        df = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], header=cfg["header_rows"], index_col=0, engine='openpyxl')
        df.columns.names = ['Main_Class', 'Subclass', 'Compound']
        print(f"Loaded: {df.shape}")
    except Exception as e:
        print(f"Error: {e}")
        continue
    df = df.dropna(axis=0, how='all')

    numeric_df = df.select_dtypes(include=[np.number])

    numeric_df_filled = numeric_df.copy()
    for col in numeric_df_filled.columns:
        col_min = numeric_df_filled[col].min(skipna=True)
        if pd.notna(col_min) and col_min > 0:
            numeric_df_filled[col] = numeric_df_filled[col].fillna(col_min * cfg["settings"].value("NONDETECT_FRACTION"))
        else:
            numeric_df_filled[col] = numeric_df_filled[col].fillna(0)

    na_before = numeric_df.isna().sum().sum()
    print(f"Filled {na_before} NAs with LOD/2")

    subclass_colors_palette, main_class_color_map, unique_subclasss = build_palette(numeric_df_filled.columns)
    main_class_color_map = dict(cfg["main_class_colours"])
    main_classes = list(dict.fromkeys(numeric_df_filled.columns.get_level_values('Main_Class')))
    print(f"Palette: {len(unique_subclasss)} subclasses, {len(main_classes)} main classes")

    compound_columns = sorted(numeric_df_filled.columns.tolist(), key=lambda x: (x[0], x[1]))
    n_compounds = len(compound_columns)
    n_cols = 8
    n_rows = (n_compounds + n_cols - 1) // n_cols

    fig_width = n_cols * 3.5
    fig_height = n_rows * 4.0

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height),
                              sharex=False, sharey=False)
    axes = axes.flatten()

    for i, compound in enumerate(compound_columns):
        ax = axes[i]
        main_class = compound[0]
        color = main_class_color_map.get(main_class, (0.5, 0.5, 0.5))
        data = numeric_df_filled[compound]

        if not data.empty:
            sns.violinplot(y=data, ax=ax, color=color, inner=None, cut=0, linewidth=0)
            sns.boxplot(y=data, ax=ax, width=0.2, color='white',
                        boxprops=dict(alpha=0.5, edgecolor='black', linewidth=1.5),
                        medianprops=dict(color='black', linewidth=2),
                        capprops=dict(color='black', linewidth=1.5),
                        whiskerprops=dict(color='black', linewidth=1.5))
            ax.set_title(compound[2], fontsize=20, fontweight='bold', pad=4)
            ax.set_xlabel("")
            ax.tick_params(axis='x', bottom=False, labelbottom=False)
            ax.tick_params(axis='y', labelsize=25)
            if i % n_cols == 0:
                ax.set_ylabel(unit_label, fontsize=18, fontweight='bold')
            else:
                ax.set_ylabel("")
            ax.set_ylim(bottom=0)
            ax.grid(False)
            for spine in ax.spines.values():
                spine.set_edgecolor('black')
                spine.set_linewidth(1.2)
        else:
            ax.set_visible(False)

    empty_axes = []
    for i in range(n_compounds, len(axes)):
        empty_axes.append(axes[i])
        axes[i].set_visible(False)

    handles = [mpatches.Patch(color=main_class_color_map[mc], label=mc) for mc in main_classes]
    if empty_axes:
        leg_ax = empty_axes[-1]
        leg_ax.set_visible(True)
        leg_ax.axis('off')
        leg_ax.legend(handles=handles, title='Main Class', loc='center',
                      fontsize=11, title_fontsize=12, frameon=False)
    else:
        fig.legend(handles=handles, title='Main Class',
                   bbox_to_anchor=(1.01, 0.5), loc='center left', frameon=False, fontsize=10)

    plt.tight_layout(rect=[0, 0, 1, 1])
    save_path = os.path.join(output_folder, f"violin_box_filled_{sheet_name}.png")
    save_figure(save_path, dpi=DPI, bbox_inches='tight')
    plt.close('all')
    print(f"Saved: {save_path}")

    gc.collect()
    print(f"\n{sheet_name} complete.")

finish_step(cfg, TREATMENT)

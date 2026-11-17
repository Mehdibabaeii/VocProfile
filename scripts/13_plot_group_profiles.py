#!/usr/bin/env python3
#
# VocProfile - 13_plot_group_profiles.py
#
# Stacked profiles of the groups by main class, subclass and compound, once per input
# set: K groups on the accession means of every dataset, clusters on the BLUP
# estimates for all compounds, and for the 12 markers. The three layouts differ in
# group order, palette and figure sizes, as in the figures of the paper.
#
# in   MEANS_FILE, BLUP_CLUSTER_FILE, INPUT_SETS_TABLE (step 13), MARKER_TABLE, SUBCLASS_ORDER_TABLE
# out  figures/13_plot_group_profiles/<input set>/group_profiles_<dataset>/
#      results/13_plot_group_profiles/13_qc.txt
#
# Usage:
#   python scripts/13_plot_group_profiles.py <config.txt>
#
# Notes: docs/NOTES.md, step 13
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

from lib_config import load_config, input_path, require_sheets, output_dir, selected_sets
from lib_run import start_step, finish_step, save_figure, set_scope

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')

cfg, step, args = load_config("13")
start_step(cfg, 'means', 'blup_clusters')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'group means use the accessions with a value; a compound without any value in a group counts as zero in the stack',
)
DPI = step["dpi"]
ABUNDANCE_SHEET = cfg["abundance_sheet"]


def draw_k_groups(set_name, input_set):
    file_path = input_path(cfg, input_set["input"])
    sheets_to_run = input_set["sheets"]
    require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets_to_run])
    output_folder_base = os.path.join(output_dir(cfg, "13_plot_group_profiles", set_name), input_set["output_folder"])

    for sheet_name in sheets_to_run:
        print(f"Dataset: {sheet_name}")
        output_folder = f"{output_folder_base}_{sheet_name}"
        os.makedirs(output_folder, exist_ok=True)
        unit_label = "Abundance of VOCs" if sheet_name == ABUNDANCE_SHEET else "Estimated Concentration of VOCs (µg/g)"

        try:
            df_full = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name],
                                    header=cfg["header_rows"], index_col=[0, 1], engine="openpyxl")
            df_full.index.names = ["Acc", "K"]
            df_full.columns.names = ["Main_Class", "Subclass", "Compound"]
            print(f"Loaded: {df_full.shape}")
        except Exception as e:
            print(f"Error: {e}")
            continue
        df_full = df_full.dropna(axis=0, how="all")

        k_series = df_full.index.get_level_values("K").astype(str)
        df = df_full.reset_index(level="K", drop=True)
        df = df.apply(pd.to_numeric, errors="coerce")
        print(f"   K groups: {sorted(set(k_series))}")

        df_grouped = df.groupby(k_series).mean()
        k_group_names = sorted(df_grouped.index.tolist())
        df_grouped = df_grouped.loc[k_group_names]
        print(f"K groups: {k_group_names} | samples per group: {pd.Series(k_series).value_counts().to_dict()}")

        mc_order = list(cfg["main_class_colours"])
        mc_order_idx = {mc: i for i, mc in enumerate(mc_order)}

        sg_by_mc = {}
        for mc, sg, cpd in df.columns:
            if mc not in sg_by_mc: sg_by_mc[mc] = []
            if sg not in sg_by_mc[mc]: sg_by_mc[mc].append(sg)

        unique_subclasss = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                if sg not in unique_subclasss:
                    unique_subclasss.append(sg)
        sg_order_idx = {sg: i for i, sg in enumerate(unique_subclasss)}

        cpd_order = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                for cpd in sorted([c for m2,s2,c in df.columns if m2==mc and s2==sg]):
                    cpd_order.append((mc, sg, cpd))

        sg_colors = sns.color_palette("husl", len(unique_subclasss))
        sg_palette = dict(zip(unique_subclasss, sg_colors))
        mc_palette = dict(cfg["main_class_colours"])

        width_bar = 0.85
        FS_TITLE, FS_LABEL, FS_TICK = 20, 16, 14

        def style_ax(ax, title="", ylabel="", xlabel="K Group"):
            ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", color="black", pad=15)
            ax.set_ylabel(ylabel, fontsize=FS_LABEL, fontweight="bold", color="black")
            ax.set_xlabel(xlabel, fontsize=FS_LABEL, fontweight="bold", color="black")
            ax.tick_params(axis="x", colors="black", labelsize=FS_TICK)
            ax.tick_params(axis="y", colors="black", labelsize=FS_TICK)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("black")
            ax.spines["bottom"].set_color("black")
            ax.grid(False)

        def get_vals(cols, grouped):
            vals = np.zeros(len(k_group_names))
            for col in cols:
                if col in grouped.columns:
                    vals += grouped[col].fillna(0).values
            return vals

        fig4, ax4 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for mc in mc_order:
            mc_cols = [(m,s,c) for m,s,c in cpd_order if m==mc]
            vals = get_vals(mc_cols, df_grouped)
            ax4.bar(k_group_names, vals, bottom=bottom, color=mc_palette[mc],
                    alpha=0.9, width=width_bar, label=mc)
            bottom += vals
        style_ax(ax4, "K-Group Profiles: Main Class Abundance", unit_label)
        ax4.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
                   fontsize=12, title="Main Class", title_fontsize=13)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"4_{sheet_name}_K_main_class_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 4 saved.")

        fig5, ax5 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        totals = get_vals(cpd_order, df_grouped)
        for mc in mc_order:
            mc_cols = [(m,s,c) for m,s,c in cpd_order if m==mc]
            vals = get_vals(mc_cols, df_grouped)
            pct = np.where(totals > 0, vals / totals * 100, 0)
            ax5.bar(k_group_names, pct, bottom=bottom, color=mc_palette[mc],
                    alpha=0.9, width=width_bar, label=mc)
            bottom += pct
        ax5.set_ylim(0, 100)
        style_ax(ax5, "K-Group Profiles: Main Class Distribution", "Percentage (%)")
        ax5.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
                   fontsize=12, title="Main Class", title_fontsize=13)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"5_{sheet_name}_K_main_class_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 5 saved.")

        fig6, ax6 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for sg in unique_subclasss:
            sg_cols = [(m,s,c) for m,s,c in cpd_order if s==sg]
            vals = get_vals(sg_cols, df_grouped)
            ax6.bar(k_group_names, vals, bottom=bottom, color=sg_palette[sg],
                    alpha=0.9, width=width_bar, label=sg)
            bottom += vals
        style_ax(ax6, "K-Group Profiles: Subclass Abundance", unit_label)
        legend_patches6 = [mpatches.Patch(color=sg_palette[sg], label=sg) for sg in unique_subclasss]
        ax6.legend(handles=legend_patches6, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"6_{sheet_name}_K_subclass_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 6 saved.")

        fig7, ax7 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for sg in unique_subclasss:
            sg_cols = [(m,s,c) for m,s,c in cpd_order if s==sg]
            vals = get_vals(sg_cols, df_grouped)
            pct = np.where(totals > 0, vals / totals * 100, 0)
            ax7.bar(k_group_names, pct, bottom=bottom, color=sg_palette[sg],
                    alpha=0.9, width=width_bar, label=sg)
            bottom += pct
        ax7.set_ylim(0, 100)
        style_ax(ax7, "K-Group Profiles: Subclass Distribution", "Percentage (%)")
        legend_patches7 = [mpatches.Patch(color=sg_palette[sg], label=sg) for sg in unique_subclasss]
        ax7.legend(handles=legend_patches7, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"7_{sheet_name}_K_subclass_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 7 saved.")

        cpd_shade_palette = {}
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                cols_in_sg = [(mc, sg, c) for mc2, sg2, c in cpd_order if mc2==mc and sg2==sg]
                n = len(cols_in_sg)
                base = np.array(sg_colors[unique_subclasss.index(sg)])
                for i, col in enumerate(cols_in_sg):
                    scale = 1 - (i / max(1, n-1)) * 0.5
                    shade = base * scale + (1 - scale) * np.array([1,1,1])
                    cpd_shade_palette[col] = shade

        fig8, ax8 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in df_grouped.columns: continue
            vals = df_grouped[col].fillna(0).values
            color = cpd_shade_palette.get(col, sg_palette.get(sg,(0.5,0.5,0.5)))
            ax8.bar(k_group_names, vals, bottom=bottom, color=color, alpha=0.9, width=width_bar)
            bottom += vals
        style_ax(ax8, "K-Group Profiles: Compound Abundance", unit_label)
        legend_patches8 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                legend_patches8.append(mpatches.Patch(color=sg_palette[sg], label=f"── {sg}"))
                for mc2, sg2, cpd in cpd_order:
                    if mc2==mc and sg2==sg:
                        col = (mc2,sg2,cpd)
                        legend_patches8.append(mpatches.Patch(
                            color=cpd_shade_palette.get(col,(0.5,0.5,0.5)), label=f"  {cpd}", alpha=0.85))
        ax8.legend(handles=legend_patches8, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=8, title="Subclass → Compound", title_fontsize=10, ncol=2)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"8_{sheet_name}_K_compound_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 8 saved.")

        fig9, ax9 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in df_grouped.columns: continue
            vals = df_grouped[col].fillna(0).values
            pct = np.where(totals > 0, vals / totals * 100, 0)
            color = cpd_shade_palette.get(col, sg_palette.get(sg,(0.5,0.5,0.5)))
            ax9.bar(k_group_names, pct, bottom=bottom, color=color, alpha=0.9, width=width_bar)
            bottom += pct
        ax9.set_ylim(0, 100)
        style_ax(ax9, "K-Group Profiles: Compound Distribution", "Percentage (%)")
        legend_patches9 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                legend_patches9.append(mpatches.Patch(color=sg_palette[sg], label=f"── {sg}"))
                for mc2, sg2, cpd in cpd_order:
                    if mc2==mc and sg2==sg:
                        col = (mc2,sg2,cpd)
                        legend_patches9.append(mpatches.Patch(
                            color=cpd_shade_palette.get(col,(0.5,0.5,0.5)), label=f"  {cpd}", alpha=0.85))
        ax9.legend(handles=legend_patches9, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=8, title="Subclass → Compound", title_fontsize=10, ncol=2)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"9_{sheet_name}_K_compound_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 9 saved.")

        gc.collect()
        print(f"\n{sheet_name} complete.")


def draw_clusters_all(set_name, input_set):
    file_path = input_path(cfg, input_set["input"])
    sheets_to_run = [input_set["output_label"]]
    output_folder_base = os.path.join(output_dir(cfg, "13_plot_group_profiles", set_name), input_set["output_folder"])

    for sheet_name in sheets_to_run:
        print(f"Dataset: {sheet_name}")
        output_folder = f"{output_folder_base}_{sheet_name}"
        os.makedirs(output_folder, exist_ok=True)
        unit_label = "Abundance of VOCs" if sheet_name == ABUNDANCE_SHEET else "Estimated Concentration (µg/g)"

        try:
            df_full = pd.read_csv(file_path, header=cfg["header_rows"], index_col=[0, 1])
            df_full.index.names = ["Acc", "K"]
            df_full.columns.names = ["Main_Class", "Subclass", "Compound"]
            print(f"Loaded: {df_full.shape}")
        except Exception as e:
            print(f"Error: {e}")
            continue
        df_full = df_full.dropna(axis=0, how="all")

        k_series = df_full.index.get_level_values("K").astype(str)
        df = df_full.reset_index(level="K", drop=True)
        df = df.apply(pd.to_numeric, errors="coerce")
        print(f"   K groups: {sorted(set(k_series))}")

        df_grouped = df.groupby(k_series).mean()
        _K_ORDER = cfg["clusters"]
        k_group_names = [k for k in _K_ORDER if k in df_grouped.index]
        df_grouped = df_grouped.loc[k_group_names]
        print(f"K groups: {k_group_names} | samples per group: {pd.Series(k_series).value_counts().to_dict()}")

        mc_order = list(cfg["main_class_colours"])
        mc_order_idx = {mc: i for i, mc in enumerate(mc_order)}

        sg_by_mc = {}
        for mc, sg, cpd in df.columns:
            if mc not in sg_by_mc: sg_by_mc[mc] = []
            if sg not in sg_by_mc[mc]: sg_by_mc[mc].append(sg)

        unique_subclasss = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                if sg not in unique_subclasss:
                    unique_subclasss.append(sg)
        sg_order_idx = {sg: i for i, sg in enumerate(unique_subclasss)}

        cpd_order = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                for cpd in sorted([c for m2,s2,c in df.columns if m2==mc and s2==sg]):
                    cpd_order.append((mc, sg, cpd))

        _all_sg = list(cfg["subclass_order"]["profiles"])
        _all_colors = sns.color_palette("husl", len(_all_sg))
        sg_colors = [_all_colors[_all_sg.index(sg)] if sg in _all_sg else (0.5,0.5,0.5) for sg in unique_subclasss]
        sg_palette = dict(zip(unique_subclasss, sg_colors))
        mc_palette = dict(cfg["main_class_colours"])

        width_bar = 0.85
        FS_TITLE, FS_LABEL, FS_TICK = 20, 16, 14

        def style_ax(ax, title="", ylabel="", xlabel=input_set["group_axis_label"]):
            ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", color="black", pad=15)
            ax.set_ylabel(ylabel, fontsize=FS_LABEL, fontweight="bold", color="black")
            ax.set_xlabel(xlabel, fontsize=FS_LABEL, fontweight="bold", color="black")
            ax.tick_params(axis="x", colors="black", labelsize=FS_TICK)
            ax.tick_params(axis="y", colors="black", labelsize=FS_TICK)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("black")
            ax.spines["bottom"].set_color("black")
            ax.grid(False)

        def get_vals(cols, grouped):
            vals = np.zeros(len(k_group_names))
            for col in cols:
                if col in grouped.columns:
                    vals += grouped[col].fillna(0).values
            return vals

        fig4, ax4 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for mc in mc_order:
            mc_cols = [(m,s,c) for m,s,c in cpd_order if m==mc]
            vals = get_vals(mc_cols, df_grouped)
            ax4.bar(k_group_names, vals, bottom=bottom, color=mc_palette[mc],
                    alpha=0.9, width=width_bar, label=mc)
            bottom += vals
        style_ax(ax4, "K-Group Profiles: Main Class Abundance", unit_label)
        ax4.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
                   fontsize=12, title="Main Class", title_fontsize=13)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"4_{sheet_name}_K_main_class_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 4 saved.")

        fig5, ax5 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        totals = get_vals(cpd_order, df_grouped)
        for mc in mc_order:
            mc_cols = [(m,s,c) for m,s,c in cpd_order if m==mc]
            vals = get_vals(mc_cols, df_grouped)
            pct = np.where(totals > 0, vals / totals * 100, 0)
            ax5.bar(k_group_names, pct, bottom=bottom, color=mc_palette[mc],
                    alpha=0.9, width=width_bar, label=mc)
            bottom += pct
        ax5.set_ylim(0, 100)
        style_ax(ax5, "K-Group Profiles: Main Class Distribution", "Percentage (%)")
        ax5.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
                   fontsize=12, title="Main Class", title_fontsize=13)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"5_{sheet_name}_K_main_class_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 5 saved.")

        fig6, ax6 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for sg in unique_subclasss:
            sg_cols = [(m,s,c) for m,s,c in cpd_order if s==sg]
            vals = get_vals(sg_cols, df_grouped)
            ax6.bar(k_group_names, vals, bottom=bottom, color=sg_palette[sg],
                    alpha=0.9, width=width_bar, label=sg)
            bottom += vals
        style_ax(ax6, "K-Group Profiles: Subclass Abundance", unit_label)
        legend_patches6 = [mpatches.Patch(color=sg_palette[sg], label=sg) for sg in unique_subclasss]
        ax6.legend(handles=legend_patches6, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"6_{sheet_name}_K_subclass_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 6 saved.")

        fig7, ax7 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for sg in unique_subclasss:
            sg_cols = [(m,s,c) for m,s,c in cpd_order if s==sg]
            vals = get_vals(sg_cols, df_grouped)
            pct = np.where(totals > 0, vals / totals * 100, 0)
            ax7.bar(k_group_names, pct, bottom=bottom, color=sg_palette[sg],
                    alpha=0.9, width=width_bar, label=sg)
            bottom += pct
        ax7.set_ylim(0, 100)
        style_ax(ax7, "K-Group Profiles: Subclass Distribution", "Percentage (%)")
        legend_patches7 = [mpatches.Patch(color=sg_palette[sg], label=sg) for sg in unique_subclasss]
        ax7.legend(handles=legend_patches7, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"7_{sheet_name}_K_subclass_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 7 saved.")

        cpd_shade_palette = {}
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                cols_in_sg = [(mc, sg, c) for mc2, sg2, c in cpd_order if mc2==mc and sg2==sg]
                n = len(cols_in_sg)
                base = np.array(sg_colors[unique_subclasss.index(sg)])
                for i, col in enumerate(cols_in_sg):
                    scale = 1 - (i / max(1, n-1)) * 0.5
                    shade = base * scale + (1 - scale) * np.array([1,1,1])
                    cpd_shade_palette[col] = shade

        fig8, ax8 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in df_grouped.columns: continue
            vals = df_grouped[col].fillna(0).values
            color = cpd_shade_palette.get(col, sg_palette.get(sg,(0.5,0.5,0.5)))
            ax8.bar(k_group_names, vals, bottom=bottom, color=color, alpha=0.9, width=width_bar)
            bottom += vals
        style_ax(ax8, "K-Group Profiles: Compound Abundance", unit_label)
        legend_patches8 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                legend_patches8.append(mpatches.Patch(color=sg_palette[sg], label=f"── {sg}"))
                for mc2, sg2, cpd in cpd_order:
                    if mc2==mc and sg2==sg:
                        col = (mc2,sg2,cpd)
                        legend_patches8.append(mpatches.Patch(
                            color=cpd_shade_palette.get(col,(0.5,0.5,0.5)), label=f"  {cpd}", alpha=0.85))
        ax8.legend(handles=legend_patches8, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=8, title="Subclass → Compound", title_fontsize=10, ncol=2)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"8_{sheet_name}_K_compound_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 8 saved.")

        fig9, ax9 = plt.subplots(figsize=(12, 7))
        bottom = np.zeros(len(k_group_names))
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in df_grouped.columns: continue
            vals = df_grouped[col].fillna(0).values
            pct = np.where(totals > 0, vals / totals * 100, 0)
            color = cpd_shade_palette.get(col, sg_palette.get(sg,(0.5,0.5,0.5)))
            ax9.bar(k_group_names, pct, bottom=bottom, color=color, alpha=0.9, width=width_bar)
            bottom += pct
        ax9.set_ylim(0, 100)
        style_ax(ax9, "K-Group Profiles: Compound Distribution", "Percentage (%)")
        legend_patches9 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                legend_patches9.append(mpatches.Patch(color=sg_palette[sg], label=f"── {sg}"))
                for mc2, sg2, cpd in cpd_order:
                    if mc2==mc and sg2==sg:
                        col = (mc2,sg2,cpd)
                        legend_patches9.append(mpatches.Patch(
                            color=cpd_shade_palette.get(col,(0.5,0.5,0.5)), label=f"  {cpd}", alpha=0.85))
        ax9.legend(handles=legend_patches9, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=8, title="Subclass → Compound", title_fontsize=10, ncol=2)
        plt.xticks(rotation=0)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"9_{sheet_name}_K_compound_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 9 saved.")

        gc.collect()
        print(f"\n{sheet_name} complete.")


def draw_clusters_markers(set_name, input_set):
    file_path = input_path(cfg, input_set["input"])
    sheets_to_run = [input_set["output_label"]]
    _feat12 = cfg["markers"]
    _all_sg = list(cfg["subclass_order"]["profiles"])
    K_ORDER   = cfg["clusters"]
    K_COLORS  = dict(cfg["cluster_colours"])
    output_folder_base = os.path.join(output_dir(cfg, "13_plot_group_profiles", set_name), input_set["output_folder"])

    for sheet_name in sheets_to_run:
        print(f"Dataset: {sheet_name}")
        output_folder = f"{output_folder_base}_{sheet_name}"
        os.makedirs(output_folder, exist_ok=True)
        unit_label = "Estimated Concentration of VOCs (µg/g)"

        try:
            df_full = pd.read_csv(file_path, header=cfg["header_rows"], index_col=[0, 1])
            df_full.index.names = ["Acc", "K"]
            df_full.columns.names = ["Main_Class", "Subclass", "Compound"]
            print(f"Loaded: {df_full.shape}")
        except Exception as e:
            print(f"Error: {e}")
            continue
        df_full = df_full.dropna(axis=0, how="all")

        k_series = df_full.index.get_level_values("K").astype(str)
        df = df_full.reset_index(level="K", drop=True)
        df = df.apply(pd.to_numeric, errors="coerce")
        print(f"   K groups found: {sorted(set(k_series))}")

        mc_order = list(cfg["main_class_colours"])
        mc_order_idx = {mc: i for i, mc in enumerate(mc_order)}

        sg_by_mc = {}
        for mc, sg, cpd in df.columns:
            if mc not in sg_by_mc: sg_by_mc[mc] = []
            if sg not in sg_by_mc[mc]: sg_by_mc[mc].append(sg)

        unique_subclasss = [sg for sg in _all_sg if any(sg == s for _, s, _ in df.columns)]

        cpd_order = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                for cpd in sorted([c for m2,s2,c in df.columns if m2==mc and s2==sg]):
                    cpd_order.append((mc, sg, cpd))

        _all_colors = sns.color_palette("husl", len(_all_sg))
        sg_palette  = {sg: _all_colors[_all_sg.index(sg)] for sg in unique_subclasss}
        mc_palette  = dict(cfg["main_class_colours"])

        cpd_shade_palette = {}
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                cols_in_sg = [(mc,sg,c) for mc2,sg2,c in cpd_order if mc2==mc and sg2==sg]
                n = len(cols_in_sg)
                base = np.array(sg_palette[sg])
                for i, col in enumerate(cols_in_sg):
                    scale = 1 - (i / max(1, n-1)) * 0.5
                    cpd_shade_palette[col] = base * scale + (1-scale)*np.array([1,1,1])

        cpd_order     = [(mc,sg,cpd) for mc,sg,cpd in cpd_order if cpd in _feat12]
        _feat12_sgs   = {sg for _,sg,_ in cpd_order}
        unique_subclasss = [sg for sg in unique_subclasss if sg in _feat12_sgs]
        sg_palette    = {sg: v for sg,v in sg_palette.items() if sg in _feat12_sgs}
        sg_by_mc      = {mc: [sg for sg in sgs if sg in _feat12_sgs] for mc,sgs in sg_by_mc.items()}
        print(f"Filtered to {len(cpd_order)} compounds in {len(unique_subclasss)} subclasses")

        df_grouped = df.groupby(k_series).mean()
        k_group_names = [k for k in K_ORDER if k in df_grouped.index]
        df_grouped = df_grouped.loc[k_group_names]
        k_bar_colors = [K_COLORS[k] for k in k_group_names]
        print(f"K groups: {k_group_names}")

        width_bar = 0.85
        FS_TITLE, FS_LABEL, FS_TICK = 22, 18, 16

        def style_ax(ax, title="", ylabel="", xlabel="Cluster"):
            ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", color="black", pad=15)
            ax.set_ylabel(ylabel, fontsize=FS_LABEL, fontweight="bold", color="black")
            ax.set_xlabel(xlabel, fontsize=FS_LABEL, fontweight="bold", color="black")
            ax.tick_params(axis="x", colors="black", labelsize=FS_TICK, width=1.5)
            ax.tick_params(axis="y", colors="black", labelsize=FS_TICK, width=1.5)
            for lbl in ax.get_xticklabels(): lbl.set_fontweight("bold")
            for lbl in ax.get_yticklabels(): lbl.set_fontweight("bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("black")
            ax.spines["bottom"].set_color("black")
            ax.grid(False)

        def get_vals(cols, grouped):
            vals = np.zeros(len(k_group_names))
            for col in cols:
                if col in grouped.columns:
                    vals += grouped[col].fillna(0).values
            return vals

        totals = get_vals(cpd_order, df_grouped)

        fig4, ax4 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for mc in mc_order:
            mc_cols = [(m,s,c) for m,s,c in cpd_order if m==mc]
            vals = get_vals(mc_cols, df_grouped)
            ax4.bar(k_group_names, vals, bottom=bottom, color=mc_palette[mc], alpha=0.9, width=width_bar, label=mc)
            bottom += vals
        style_ax(ax4, "Cluster Profiles: Main Class Abundance", unit_label)
        ax4.legend(loc="upper left", bbox_to_anchor=(1.01,1), frameon=False, fontsize=12, title="Main Class", title_fontsize=13)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"4_K_main_class_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all"); print("Plot 4 saved.")

        fig5, ax5 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for mc in mc_order:
            mc_cols = [(m,s,c) for m,s,c in cpd_order if m==mc]
            vals = get_vals(mc_cols, df_grouped)
            pct = np.where(totals>0, vals/totals*100, 0)
            ax5.bar(k_group_names, pct, bottom=bottom, color=mc_palette[mc], alpha=0.9, width=width_bar, label=mc)
            bottom += pct
        ax5.set_ylim(0, 100)
        style_ax(ax5, "Cluster Profiles: Main Class Distribution", "Percentage (%)")
        ax5.legend(loc="upper left", bbox_to_anchor=(1.01,1), frameon=False, fontsize=12, title="Main Class", title_fontsize=13)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"5_K_main_class_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all"); print("Plot 5 saved.")

        fig6, ax6 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for sg in unique_subclasss:
            sg_cols = [(m,s,c) for m,s,c in cpd_order if s==sg]
            vals = get_vals(sg_cols, df_grouped)
            ax6.bar(k_group_names, vals, bottom=bottom, color=sg_palette[sg], alpha=0.9, width=width_bar, label=sg)
            bottom += vals
        style_ax(ax6, "Cluster Profiles: Subclass Abundance", unit_label)
        ax6.legend(handles=[mpatches.Patch(color=sg_palette[sg], label=sg) for sg in unique_subclasss],
                   loc="upper left", bbox_to_anchor=(1.01,1), frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"6_K_subclass_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all"); print("Plot 6 saved.")

        fig7, ax7 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for sg in unique_subclasss:
            sg_cols = [(m,s,c) for m,s,c in cpd_order if s==sg]
            vals = get_vals(sg_cols, df_grouped)
            pct = np.where(totals>0, vals/totals*100, 0)
            ax7.bar(k_group_names, pct, bottom=bottom, color=sg_palette[sg], alpha=0.9, width=width_bar, label=sg)
            bottom += pct
        ax7.set_ylim(0, 100)
        style_ax(ax7, "Cluster Profiles: Subclass Distribution", "Percentage (%)")
        ax7.legend(handles=[mpatches.Patch(color=sg_palette[sg], label=sg) for sg in unique_subclasss],
                   loc="upper left", bbox_to_anchor=(1.01,1), frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"7_K_subclass_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all"); print("Plot 7 saved.")

        fig8, ax8 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in df_grouped.columns: continue
            vals = df_grouped[col].fillna(0).values
            color = cpd_shade_palette.get(col, sg_palette.get(sg,(0.5,0.5,0.5)))
            ax8.bar(k_group_names, vals, bottom=bottom, color=color, alpha=0.9, width=width_bar)
            bottom += vals
        style_ax(ax8, "Cluster Profiles: Compound Abundance", unit_label)
        leg8 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc,[]):
                leg8.append(mpatches.Patch(color=sg_palette[sg], label=f"── {sg}"))
                for mc2,sg2,cpd in cpd_order:
                    if mc2==mc and sg2==sg:
                        col=(mc2,sg2,cpd)
                        leg8.append(mpatches.Patch(color=cpd_shade_palette.get(col,(0.5,0.5,0.5)), label=f"  {cpd}", alpha=0.85))
        ax8.legend(handles=leg8, loc="upper left", bbox_to_anchor=(1.01,1), frameon=False,
                   fontsize=9, title="Subclass → Compound", title_fontsize=11, ncol=1)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"8_K_compound_abundance.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all"); print("Plot 8 saved.")

        fig9, ax9 = plt.subplots(figsize=(10, 7))
        bottom = np.zeros(len(k_group_names))
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in df_grouped.columns: continue
            vals = df_grouped[col].fillna(0).values
            pct = np.where(totals>0, vals/totals*100, 0)
            color = cpd_shade_palette.get(col, sg_palette.get(sg,(0.5,0.5,0.5)))
            ax9.bar(k_group_names, pct, bottom=bottom, color=color, alpha=0.9, width=width_bar)
            bottom += pct
        ax9.set_ylim(0, 100)
        style_ax(ax9, "Cluster Profiles: Compound Distribution", "Percentage (%)")
        leg9 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc,[]):
                leg9.append(mpatches.Patch(color=sg_palette[sg], label=f"── {sg}"))
                for mc2,sg2,cpd in cpd_order:
                    if mc2==mc and sg2==sg:
                        col=(mc2,sg2,cpd)
                        leg9.append(mpatches.Patch(color=cpd_shade_palette.get(col,(0.5,0.5,0.5)), label=f"  {cpd}", alpha=0.85))
        ax9.legend(handles=leg9, loc="upper left", bbox_to_anchor=(1.01,1), frameon=False,
                   fontsize=9, title="Subclass → Compound", title_fontsize=11, ncol=1)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"9_K_compound_percentage.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all"); print("Plot 9 saved.")

        gc.collect()
        print(f"\n{sheet_name} complete.")


LAYOUTS = {"k_groups": draw_k_groups, "clusters_all": draw_clusters_all, "clusters_markers": draw_clusters_markers}
for set_name, input_set in selected_sets(step, args).items():
    set_scope(input_set["grouping"])
    LAYOUTS[input_set["layout"]](set_name, input_set)

finish_step(cfg, TREATMENT)

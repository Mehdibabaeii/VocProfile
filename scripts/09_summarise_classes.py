#!/usr/bin/env python3
#
# VocProfile - 09_summarise_classes.py
#
# Detection by main class and subclass and stacked profiles per accession, once per
# input set: accession means of every dataset, BLUP estimates with clusters for all
# compounds, and for the 12 markers (Fig. 2d, Fig. 3c, Figs. S3-S4).
#
# in   MEANS_FILE, BLUP_CLUSTER_FILE, INPUT_SETS_TABLE (step 09), MARKER_TABLE
# out  figures/09_summarise_classes/<input set>/class_summaries_<dataset>/
#      results/09_summarise_classes/09_qc.txt
#
# Usage:
#   python scripts/09_summarise_classes.py <config.txt>
#
# Notes: docs/NOTES.md, step 09
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
from matplotlib.legend import Legend

from lib_config import load_config, input_path, require_sheets, output_dir, selected_sets, term_function
from lib_run import start_step, finish_step, save_figure, set_scope, revised, save_legend, deferred

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')

cfg, step, args = load_config("09")
start_step(cfg, 'means', 'blup_clusters')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'detection counts use the values as they are',
    'in the stacked profiles a value without detection counts as zero',
)
ABUNDANCE_SHEET = cfg["abundance_sheet"]
DPI = step["dpi"]
_feat12 = cfg["markers"]
TT = deferred(term_function(cfg, "revised_terms_in_titles"))

for set_name, input_set in selected_sets(step, args).items():
    file_path = input_path(cfg, input_set["input"])
    set_scope(input_set["grouping"])
    FROM_CSV = file_path.lower().endswith(".csv")
    CLUSTER_LEGEND = input_set["cluster_legend"]
    MARKERS_ONLY = input_set["compounds"] == "markers"
    if FROM_CSV:
        sheets = [input_set["output_label"]]
    else:
        sheets = input_set["sheets"]
        require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets])
    output_root = output_dir(cfg, "09_summarise_classes", set_name)

    for sheet_name in sheets:
        print(f"\nInput set {set_name}, dataset {sheet_name}")
        if sheet_name == ABUNDANCE_SHEET:
            unit_label = "Abundance of VOCs"
            unit_short = "Normalized Abundance"
        else:
            unit_label = "Estimated Concentration of VOCs (µg/g)"
            unit_short = "Concentration (µg/g)"

        output_folder = os.path.join(output_root, f"{input_set['output_folder']}_{sheet_name}")
        os.makedirs(output_folder, exist_ok=True)
        _kcolors = dict(cfg["cluster_colours"])
        if CLUSTER_LEGEND:
            try:
                _raw2 = pd.read_csv(file_path, header=cfg["header_rows"], index_col=[0,1])
                class_mapping = {str(acc): str(k) for acc, k in _raw2.index}
                print(f"Cluster mapping: {len(class_mapping)} accessions, clusters={set(class_mapping.values())}")
            except Exception as _ce:
                class_mapping = {}
                print(f"No K column found: {_ce}")

        def _color_labels(ax):
            for lbl in ax.get_xticklabels():
                c = class_mapping.get(lbl.get_text())
                if c in _kcolors:
                    lbl.set_color(_kcolors[c])
                    lbl.set_fontweight('bold')

        def _cluster_legend(ax):
            if revised():
                # the legend would cover bars; step 06 writes it as an image of its own
                return
            handles = [mpatches.Patch(color=col, label=f'Cluster {k}')
                       for k, col in _kcolors.items()]
            leg = Legend(ax, handles, [h.get_label() for h in handles],
                         title='Clusters', loc='upper right',
                         frameon=True, fontsize=12, title_fontsize=12)
            ax.add_artist(leg)

        try:
            if FROM_CSV:
                _df2 = pd.read_csv(file_path, header=cfg["header_rows"], index_col=[0, 1])
                _df2.columns.names = ["Main_Class", "Subclass", "Compound"]
                df = _df2.copy()
                df.index = _df2.index.get_level_values(0)
                df.index.name = None
                df = df.apply(pd.to_numeric, errors='coerce')
                print(f"CSV loaded: {df.shape[0]} accessions × {df.shape[1]} compounds")
            else:
                df = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], header=cfg["header_rows"], index_col=0)
                df.columns.names = ["Main_Class", "Subclass", "Compound"]
                print(f"Sheet '{sheet_name}' loaded.")
        except Exception as e:
            print(f"Error loading '{file_path}': {e}")
            continue

        df = df.dropna(axis=0, how="all")

        numeric_df = df.select_dtypes(include=[np.number])
        total_samples = len(numeric_df)
        accessions = df.index.to_list()

        class_missing_counts = numeric_df.isnull().sum(axis=0).groupby(level=0).sum()
        class_present_counts = (total_samples * numeric_df.columns.get_level_values(0).value_counts()) - class_missing_counts
        class_plot_df = pd.DataFrame({
            "Missing_Count": class_missing_counts,
            "Present_Count": class_present_counts
        }).reset_index()
        class_plot_df.columns = ["Main_Class", "Missing_Count", "Present_Count"]

        missing_counts = numeric_df.isnull().sum(axis=0)
        missing_df_plot = missing_counts.to_frame(name="Missing_Count").reset_index()
        level_cols = [c for c in missing_df_plot.columns if c != "Missing_Count"]
        if len(level_cols) == 3:
            missing_df_plot.columns = ["Main_Class", "Subclass", "Compound", "Missing_Count"]
        elif len(level_cols) == 2:
            missing_df_plot.columns = ["Main_Class", "Subclass", "Missing_Count"]
        missing_df_plot["Present_Count"] = total_samples - missing_df_plot["Missing_Count"]

        mc_order = list(cfg["main_class_colours"])
        mc_order_idx = {mc: i for i, mc in enumerate(mc_order)}

        sg_by_mc = {}
        for mc, sg, cpd in numeric_df.columns:
            if mc not in sg_by_mc:
                sg_by_mc[mc] = []
            if sg not in sg_by_mc[mc]:
                sg_by_mc[mc].append(sg)
        unique_subclasss = []
        for mc in mc_order:
            if mc in sg_by_mc:
                for sg in sg_by_mc[mc]:
                    if sg not in unique_subclasss:
                        unique_subclasss.append(sg)
        sg_order_idx = {sg: i for i, sg in enumerate(unique_subclasss)}

        cpd_order = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                cpds = sorted([cpd for mc2, sg2, cpd in numeric_df.columns if mc2 == mc and sg2 == sg])
                for cpd in cpds:
                    cpd_order.append((mc, sg, cpd))
        cpd_order_idx = {cpd: i for i, (mc, sg, cpd) in enumerate(cpd_order)}

        missing_df_plot["mc_order"] = missing_df_plot["Main_Class"].map(mc_order_idx)
        missing_df_plot["sg_order"] = missing_df_plot["Subclass"].map(sg_order_idx)
        if MARKERS_ONLY:
            missing_df_plot = missing_df_plot[missing_df_plot["Compound"].isin(_feat12)].copy()
            missing_df_plot = missing_df_plot.sort_values(["mc_order","sg_order"]).reset_index(drop=True)
        if "Compound" in missing_df_plot.columns:
            missing_df_plot["cpd_order"] = missing_df_plot["Compound"].map(cpd_order_idx)
            missing_df_plot = missing_df_plot.sort_values(["mc_order", "sg_order", "cpd_order"]).reset_index(drop=True)
        else:
            missing_df_plot = missing_df_plot.sort_values(["mc_order", "sg_order"]).reset_index(drop=True)

        if MARKERS_ONLY:
            cpd_order_plot = [(mc,sg,cpd) for mc,sg,cpd in cpd_order if cpd in _feat12]
        else:
            cpd_order_plot = cpd_order

        subclass_base_colors = sns.color_palette("husl", len(unique_subclasss))
        subclass_colors_palette = dict(zip(unique_subclasss, subclass_base_colors))

        main_class_color_map = dict(cfg["main_class_colours"])

        if MARKERS_ONLY:
            cpd_order = [(mc, sg, cpd) for mc, sg, cpd in cpd_order if cpd in _feat12]
            _feat12_cols = [col for col in numeric_df.columns if col[2] in _feat12]
            numeric_df = numeric_df[_feat12_cols]
            _feat12_sgs = {sg for _, sg, _ in cpd_order}
            unique_subclasss = [sg for sg in unique_subclasss if sg in _feat12_sgs]
            subclass_colors_palette = {sg: v for sg, v in subclass_colors_palette.items() if sg in _feat12_sgs}
            sg_by_mc = {mc: [sg for sg in sgs if sg in _feat12_sgs] for mc, sgs in sg_by_mc.items()}
            sg_order_idx = {sg: i for i, sg in enumerate(unique_subclasss)}
            print(f"Filtered to {len(cpd_order)} compounds in {len(unique_subclasss)} subclasses")

        class_plot_df["mc_order"] = class_plot_df["Main_Class"].map(mc_order_idx)
        class_plot_df = class_plot_df.sort_values("mc_order").reset_index(drop=True)
        class_plot_df["Total_Count"] = class_plot_df["Missing_Count"] + class_plot_df["Present_Count"]
        class_plot_df["Present_Percentage"] = class_plot_df["Present_Count"] / class_plot_df["Total_Count"] * 100
        class_plot_df["Missing_Percentage"] = class_plot_df["Missing_Count"] / class_plot_df["Total_Count"] * 100

        fig1, ax1 = plt.subplots(figsize=(8, 6))
        mc_colors = [main_class_color_map[mc] for mc in class_plot_df["Main_Class"]]
        ax1.bar(class_plot_df["Main_Class"], class_plot_df["Present_Percentage"],
                color=mc_colors, alpha=0.9, width=0.9)
        ax1.bar(class_plot_df["Main_Class"], class_plot_df["Missing_Percentage"],
                bottom=class_plot_df["Present_Percentage"], color="lightgray", width=0.9)
        for i, (p, m) in enumerate(zip(class_plot_df["Present_Percentage"], class_plot_df["Missing_Percentage"])):
            if p > 0:
                ax1.text(i, p/2, f"{p:.1f}%", ha="center", va="center", fontweight="bold", fontsize=14, color="white")
            if m > 0:
                ax1.text(i, p+m/2, f"{m:.1f}%", ha="center", va="center", fontweight="bold", fontsize=14, color="white")
        ax1.set_ylabel("Percentage (%)", fontsize=16, fontweight="bold", color="black")
        ax1.set_xlabel("Main Chemical Classes", fontsize=16, fontweight="bold", color="black")
        ax1.tick_params(axis="x", colors="black", labelsize=16)
        ax1.tick_params(axis="y", colors="black", labelsize=16)
        plt.xticks(rotation=45, ha="right")
        ax1.grid(False)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.spines["left"].set_color("black")
        ax1.spines["bottom"].set_color("black")
        legend_patches1 = [mpatches.Patch(color=main_class_color_map[mc], label=mc) for mc in mc_order if mc in class_plot_df["Main_Class"].values]
        legend_patches1.append(mpatches.Patch(color="lightgray", label="Missing"))
        ax1.legend(handles=legend_patches1, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False, fontsize=12)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"1_{sheet_name}_main_class_distribution.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 1 saved.")

        subclass_stats = missing_df_plot.groupby(["Main_Class", "Subclass"], sort=False).agg(
            {"Present_Count": "sum", "Missing_Count": "sum"}
        ).reset_index()
        subclass_stats["sg_order"] = subclass_stats["Subclass"].map(sg_order_idx)
        subclass_stats["mc_order"] = subclass_stats["Main_Class"].map(mc_order_idx)
        subclass_stats = subclass_stats.sort_values(["mc_order", "sg_order"]).reset_index(drop=True)
        subclass_stats["Total_Count"] = subclass_stats["Present_Count"] + subclass_stats["Missing_Count"]
        subclass_stats["Present_Percentage"] = subclass_stats["Present_Count"] / subclass_stats["Total_Count"] * 100
        subclass_stats["Missing_Percentage"] = subclass_stats["Missing_Count"] / subclass_stats["Total_Count"] * 100

        fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(16, 6))
        bw = 1
        sg_colors_2a = [subclass_colors_palette[sg] for sg in subclass_stats["Subclass"]]
        ax2a.bar(range(len(subclass_stats)), subclass_stats["Present_Percentage"],
                 color=sg_colors_2a, alpha=0.9, width=bw)
        ax2a.bar(range(len(subclass_stats)), subclass_stats["Missing_Percentage"],
                 bottom=subclass_stats["Present_Percentage"], color="lightgray", width=bw)
        for i, (p, m) in enumerate(zip(subclass_stats["Present_Percentage"], subclass_stats["Missing_Percentage"])):
            if p > 0:
                ax2a.text(i, p/2, f"{p:.1f}%", ha="center", va="center", fontweight="bold", fontsize=10, color="white")
            if m > 0:
                ax2a.text(i, p+m/2, f"{m:.1f}%", ha="center", va="center", fontweight="bold", fontsize=10, color="white")
        ax2a.set_ylabel("Presence Percentage (%)", fontsize=16, fontweight="bold", color="black")
        ax2a.set_xlabel("Subclass", fontsize=16, fontweight="bold", color="black")
        ax2a.set_xticks(range(len(subclass_stats)))
        ax2a.set_xticklabels(subclass_stats["Subclass"], rotation=45, ha="right", fontsize=8)
        ax2a.tick_params(axis="x", colors="black", labelsize=16)
        ax2a.tick_params(axis="y", colors="black", labelsize=16)
        ax2a.grid(False)
        ax2a.spines["top"].set_visible(False)
        ax2a.spines["right"].set_visible(False)
        ax2a.spines["left"].set_color("black")
        ax2a.spines["bottom"].set_color("black")

        sg_counts_2b = [(sg, len(missing_df_plot[missing_df_plot["Subclass"] == sg]))
                        for sg in unique_subclasss if len(missing_df_plot[missing_df_plot["Subclass"] == sg]) > 0]
        sg_labels_2b = [x[0] for x in sg_counts_2b]
        sg_vals_2b   = [x[1] for x in sg_counts_2b]
        sg_colors_2b = [subclass_colors_palette[sg] for sg in sg_labels_2b]
        ax2b.bar(range(len(sg_labels_2b)), sg_vals_2b, color=sg_colors_2b, alpha=0.9, width=bw)
        for i, count in enumerate(sg_vals_2b):
            ax2b.text(i, count/2, str(count), ha="center", va="center", fontweight="bold", fontsize=14, color="white")
        ax2b.set_ylabel("Number of Compounds", fontsize=16, fontweight="bold", color="black")
        ax2b.set_xlabel("Subclass", fontsize=16, fontweight="bold", color="black")
        ax2b.set_xticks(range(len(sg_labels_2b)))
        ax2b.set_xticklabels(sg_labels_2b, rotation=45, ha="right", fontsize=8)
        ax2b.tick_params(axis="x", colors="black", labelsize=16)
        ax2b.tick_params(axis="y", colors="black", labelsize=16)
        ax2b.grid(False)
        ax2b.spines["top"].set_visible(False)
        ax2b.spines["right"].set_visible(False)
        ax2b.spines["left"].set_color("black")
        ax2b.spines["bottom"].set_color("black")

        legend_patches2 = [mpatches.Patch(color=subclass_colors_palette[sg], label=sg) for sg in unique_subclasss if sg in subclass_stats["Subclass"].values]
        legend_patches2.append(mpatches.Patch(color="lightgray", label="Missing"))
        fig2.legend(handles=legend_patches2, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"2_{sheet_name}_subclass_analysis.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 2 saved.")

        main_class_counts = missing_df_plot.groupby("Main_Class", sort=False)["Compound"].count()
        main_class_counts = main_class_counts.reindex([mc for mc in mc_order if mc in main_class_counts.index])
        subclass_counts = missing_df_plot.groupby(["Main_Class", "Subclass"], sort=False)["Compound"].count()
        ordered_idx = [(mc, sg) for mc in mc_order for sg in unique_subclasss
                       if (mc, sg) in subclass_counts.index]
        subclass_counts = subclass_counts.reindex(ordered_idx)

        inner_colors = [main_class_color_map[mc] for mc in main_class_counts.index]
        outer_colors  = [subclass_colors_palette[sg] for mc, sg in subclass_counts.index]

        fig_d, ax_d = plt.subplots(figsize=(12, 10))
        wedges1, texts1, autotexts1 = ax_d.pie(
            main_class_counts, labels=None, radius=0.6, colors=inner_colors,
            autopct="%1.1f%%", startangle=90, pctdistance=0.6,
            textprops={"fontsize": 16, "fontweight": "bold", "color": "white"},
            wedgeprops={"edgecolor": "white", "linewidth": 0.8}
        )
        wedges2, texts2, autotexts2 = ax_d.pie(
            subclass_counts, labels=None, radius=1.0, colors=outer_colors,
            autopct="%1.1f%%", startangle=90, pctdistance=0.85,
            textprops={"fontsize": 16, "fontweight": "bold", "color": "white"},
            wedgeprops={"width": 0.35, "edgecolor": "white", "linewidth": 1}
        )
        for at in autotexts1 + autotexts2:
            at.set_fontsize(16)
            at.set_fontweight("bold")
            at.set_color("white")
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"donut_{sheet_name}.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Donut plot saved.")

        missing_df_sorted = missing_df_plot.copy()
        compound_colors3 = [subclass_colors_palette.get(sg, (0.5,0.5,0.5)) for sg in missing_df_sorted["Subclass"]]

        fig3, ax3 = plt.subplots(figsize=(28, 10))
        bar_w3 = 0.92
        ax3.bar(range(len(missing_df_sorted)), missing_df_sorted["Present_Count"],
                color=compound_colors3, alpha=0.9, width=bar_w3)
        ax3.bar(range(len(missing_df_sorted)), missing_df_sorted["Missing_Count"],
                bottom=missing_df_sorted["Present_Count"], color="lightgray", alpha=0.6, width=bar_w3)
        ax3.set_ylabel("Number of Cultivars", fontsize=20, fontweight="bold", color="black")
        ax3.set_xlabel("Compounds", fontsize=20, fontweight="bold", color="black")
        ax3.set_xticks(range(len(missing_df_sorted)))
        ax3.set_xticklabels(missing_df_sorted["Compound"], rotation=90, ha="center", fontsize=20, fontweight="bold")
        ax3.tick_params(axis="y", colors="black", labelsize=20)
        ax3.tick_params(axis="x", colors="black", pad=2)
        ax3.set_xlim(-0.5, len(missing_df_sorted) - 0.5)
        ax3.grid(False)
        ax3.spines["top"].set_visible(False)
        ax3.spines["right"].set_visible(False)
        ax3.spines["left"].set_color("black")
        ax3.spines["left"].set_linewidth(1.5)
        ax3.spines["bottom"].set_color("black")
        ax3.spines["bottom"].set_linewidth(1.5)

        current_class = missing_df_sorted["Main_Class"].iloc[0]
        class_positions = [0]
        for i in range(1, len(missing_df_sorted)):
            if missing_df_sorted["Main_Class"].iloc[i] != current_class:
                ax3.axvline(x=i-0.5, color="black", linestyle="--", linewidth=2, alpha=0.7)
                class_positions.append(i)
                current_class = missing_df_sorted["Main_Class"].iloc[i]
        class_positions.append(len(missing_df_sorted))

        legend_patches3 = [mpatches.Patch(color=subclass_colors_palette[sg], label=sg) for sg in unique_subclasss]
        legend_patches3.append(mpatches.Patch(color="lightgray", label="Missing"))
        ax3.legend(handles=legend_patches3, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        plt.tight_layout(rect=[0, 0, 0.88, 0.95])
        save_figure(os.path.join(output_folder, f"3_{sheet_name}_detailed_compound_view.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 3 saved.")

        n_accessions = len(df)
        sample_accessions = df.iloc[:n_accessions]
        accession_names = sample_accessions.index.to_list()
        width_bar = 0.9

        compounds_data_simple = {}
        for mc, sg, cpd in cpd_order:
            col = (mc, sg, cpd)
            if col not in numeric_df.columns:
                continue
            compounds_data_simple[cpd] = {
                "values": [sample_accessions.loc[idx, col] if not pd.isna(sample_accessions.loc[idx, col]) else 0
                           for idx in sample_accessions.index],
                "color": subclass_colors_palette.get(sg, (0.5,0.5,0.5)),
                "subgroup": sg,
                "main_class": mc,
                "column": col,
            }

        compounds_data_detailed = {}
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                cols_in_sg = [(mc2, sg2, cpd2) for mc2, sg2, cpd2 in cpd_order if mc2 == mc and sg2 == sg]
                n_cpd = len(cols_in_sg)
                base_color = np.array(subclass_colors_palette[sg])
                for i, col in enumerate(cols_in_sg):
                    mc2, sg2, cpd2 = col
                    if col not in numeric_df.columns:
                        continue
                    scale = 1 - (i / max(1, n_cpd-1)) * 0.5
                    shade = base_color * scale + (1 - scale) * np.array([1,1,1])
                    compounds_data_detailed[cpd2] = {
                        "values": [sample_accessions.loc[idx, col] if not pd.isna(sample_accessions.loc[idx, col]) else 0
                                   for idx in sample_accessions.index],
                        "color": shade,
                        "subgroup": sg2,
                        "main_class": mc2,
                        "column": col,
                    }

        fig4, ax4 = plt.subplots(figsize=(35, 10))
        bottom_values = np.zeros(n_accessions)
        prev_sg = None
        for mc, sg, cpd in cpd_order_plot:
            if cpd not in compounds_data_simple:
                continue
            cinfo = compounds_data_simple[cpd]
            if prev_sg is not None and sg != prev_sg:
                pass
            ax4.bar(accession_names, cinfo["values"], bottom=bottom_values,
                    color=cinfo["color"], alpha=0.8, width=width_bar, label=sg)
            bottom_values += np.array(cinfo["values"])
            prev_sg = sg
        ax4.set_title(TT("Genotype Profiles: Compound Abundance"), fontsize=20, fontweight="bold", color="black", pad=20)
        ax4.set_ylabel(unit_label, fontsize=25, fontweight="bold", color="black")
        ax4.set_xlabel("Cultivars", fontsize=25, fontweight="bold", color="black")
        ax4.set_xlim(-0.5, n_accessions-0.5)
        ax4.spines["top"].set_visible(False)
        ax4.spines["right"].set_visible(False)
        ax4.spines["left"].set_color("black")
        ax4.spines["bottom"].set_color("black")
        ax4.tick_params(axis="x", colors="black", labelsize=16)
        ax4.tick_params(axis="y", colors="black", labelsize=20)
        plt.xticks(rotation=90)
        if CLUSTER_LEGEND:
            _color_labels(ax4)
            _cluster_legend(ax4)
        ax4.grid(False)
        handles4, labels4 = ax4.get_legend_handles_labels()
        by_label4 = dict(zip(labels4, handles4))
        ax4.legend(by_label4.values(), by_label4.keys(), loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9, title="Subclass", title_fontsize=10)
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        save_figure(os.path.join(output_folder, f"4_{sheet_name}_accession_compound_profiles_stacked.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 4 saved.")

        fig5, ax5 = plt.subplots(figsize=(35, 10))
        totals = np.zeros(n_accessions)
        for mc, sg, cpd in cpd_order_plot:
            if cpd in compounds_data_simple:
                totals += np.array(compounds_data_simple[cpd]["values"])
        bottom_values = np.zeros(n_accessions)
        for mc, sg, cpd in cpd_order_plot:
            if cpd not in compounds_data_simple:
                continue
            cinfo = compounds_data_simple[cpd]
            pct = np.where(totals > 0, np.array(cinfo["values"]) / totals * 100, 0)
            ax5.bar(accession_names, pct, bottom=bottom_values,
                    color=cinfo["color"], alpha=0.8, width=width_bar)
            bottom_values += pct
        ax5.set_title(TT("Genotype Profiles: Compound Distribution"), fontsize=20, fontweight="bold", color="black", pad=20)
        ax5.set_ylabel("Percentage (%)", fontsize=25, fontweight="bold", color="black")
        ax5.set_xlabel("Cultivars", fontsize=25, fontweight="bold", color="black")
        ax5.set_ylim(0, 100)
        ax5.set_xlim(-0.5, n_accessions-0.5)
        ax5.spines["top"].set_visible(False)
        ax5.spines["right"].set_visible(False)
        ax5.spines["left"].set_color("black")
        ax5.spines["bottom"].set_color("black")
        ax5.tick_params(axis="x", colors="black", labelsize=16)
        ax5.tick_params(axis="y", colors="black", labelsize=20)
        plt.xticks(rotation=90)
        if CLUSTER_LEGEND:
            _color_labels(ax5)
            _cluster_legend(ax5)
        ax5.grid(False)
        handles5, labels5 = ax5.get_legend_handles_labels()
        by_label5 = dict(zip(labels5, handles5))
        ax5.legend(by_label5.values(), by_label5.keys(), loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9, title="Subclass", title_fontsize=10)
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        save_figure(os.path.join(output_folder, f"5_{sheet_name}_accession_compound_percentages_stacked.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 5 saved.")

        main_class_data = {mc: {"values": np.zeros(n_accessions), "color": main_class_color_map[mc]}
                           for mc in mc_order}
        for i, (idx, row) in enumerate(df.iterrows()):
            for mc, sg, cpd in cpd_order_plot:
                col = (mc, sg, cpd)
                if col in numeric_df.columns:
                    val = row[col] if not pd.isna(row[col]) else 0
                    main_class_data[mc]["values"][i] += val

        fig6, ax6 = plt.subplots(figsize=(35, 10))
        bottom_values = np.zeros(n_accessions)
        for mc in mc_order:
            vals = main_class_data[mc]["values"]
            ax6.bar(accession_names, vals, bottom=bottom_values,
                    color=main_class_data[mc]["color"], alpha=0.8, width=width_bar, label=mc)
            bottom_values += vals
        ax6.set_title(TT("Genotype Profiles: Main Class Abundance"), fontsize=20, fontweight="bold", color="black", pad=20)
        ax6.set_ylabel(unit_label, fontsize=25, fontweight="bold", color="black")
        ax6.set_xlabel("Cultivars", fontsize=25, fontweight="bold", color="black")
        ax6.set_xlim(-0.5, n_accessions-0.5)
        ax6.spines["top"].set_visible(False)
        ax6.spines["right"].set_visible(False)
        ax6.spines["left"].set_color("black")
        ax6.spines["bottom"].set_color("black")
        ax6.tick_params(axis="x", colors="black", labelsize=16)
        ax6.tick_params(axis="y", colors="black", labelsize=20)
        plt.xticks(rotation=90)
        if CLUSTER_LEGEND:
            _color_labels(ax6)
            _cluster_legend(ax6)
        ax6.grid(False)
        ax6.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=11, title="Main Class", title_fontsize=12)
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        save_figure(os.path.join(output_folder, f"6_{sheet_name}_accession_main_class_profiles.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 6 saved.")

        fig7, ax7 = plt.subplots(figsize=(35, 10))
        totals_mc = sum(main_class_data[mc]["values"] for mc in mc_order)
        bottom_values = np.zeros(n_accessions)
        for mc in mc_order:
            pct = np.where(totals_mc > 0, main_class_data[mc]["values"] / totals_mc * 100, 0)
            ax7.bar(accession_names, pct, bottom=bottom_values,
                    color=main_class_data[mc]["color"], alpha=0.8, width=width_bar, label=mc)
            bottom_values += pct
        ax7.set_title(TT("Genotype Profiles: Main Class Distribution"), fontsize=20, fontweight="bold", color="black", pad=20)
        ax7.set_ylabel("Percentage (%)", fontsize=25, fontweight="bold", color="black")
        ax7.set_xlabel("Cultivars", fontsize=25, fontweight="bold", color="black")
        ax7.set_ylim(0, 100)
        ax7.set_xlim(-0.5, n_accessions-0.5)
        ax7.spines["top"].set_visible(False)
        ax7.spines["right"].set_visible(False)
        ax7.spines["left"].set_color("black")
        ax7.spines["bottom"].set_color("black")
        ax7.tick_params(axis="x", colors="black", labelsize=16)
        ax7.tick_params(axis="y", colors="black", labelsize=20)
        plt.xticks(rotation=90)
        if CLUSTER_LEGEND:
            _color_labels(ax7)
            _cluster_legend(ax7)
        ax7.grid(False)
        ax7.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=11, title="Main Class", title_fontsize=12)
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        save_figure(os.path.join(output_folder, f"7_{sheet_name}_accession_main_class_percentages.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 7 saved.")

        fig8, ax8 = plt.subplots(figsize=(35, 10))
        bottom_values = np.zeros(n_accessions)
        for mc, sg, cpd in cpd_order_plot:
            if cpd not in compounds_data_detailed:
                continue
            cinfo = compounds_data_detailed[cpd]
            ax8.bar(accession_names, cinfo["values"], bottom=bottom_values,
                    color=cinfo["color"], alpha=0.8, width=width_bar)
            bottom_values += np.array(cinfo["values"])
        ax8.set_title(TT("Genotype Profiles: Compound Abundance"), fontsize=20, fontweight="bold", color="black", pad=20)
        ax8.set_ylabel(unit_label, fontsize=16, fontweight="bold", color="black")
        ax8.set_xlabel("Cultivars", fontsize=16, fontweight="bold", color="black")
        ax8.set_xlim(-0.5, n_accessions-0.5)
        ax8.spines["top"].set_visible(False)
        ax8.spines["right"].set_visible(False)
        ax8.spines["left"].set_color("black")
        ax8.spines["bottom"].set_color("black")
        ax8.tick_params(axis="x", colors="black", labelsize=16)
        ax8.tick_params(axis="y", colors="black", labelsize=16)
        plt.xticks(rotation=90)
        if CLUSTER_LEGEND:
            _color_labels(ax8)
            _cluster_legend(ax8)
        ax8.grid(False)
        legend_patches8 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                legend_patches8.append(mpatches.Patch(color=subclass_colors_palette[sg], label=f"── {sg}", alpha=1.0))
                for mc2, sg2, cpd2 in cpd_order:
                    if mc2 == mc and sg2 == sg and cpd2 in compounds_data_detailed:
                        legend_patches8.append(mpatches.Patch(
                            color=compounds_data_detailed[cpd2]["color"], label=f"  {cpd2}", alpha=0.85))
        ax8.legend(handles=legend_patches8, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=8, title=input_set["compound_legend_title"], title_fontsize=10,
                   ncol=2)
        plt.tight_layout(rect=[0, 0, 0.75, 1])
        save_figure(os.path.join(output_folder, f"8_{sheet_name}_accession_all_compounds_profiles.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Plot 8 saved.")

        fig9, ax9 = plt.subplots(figsize=(35, 10))
        totals_cpd = np.zeros(n_accessions)
        for mc, sg, cpd in cpd_order_plot:
            if cpd in compounds_data_detailed:
                totals_cpd += np.array(compounds_data_detailed[cpd]["values"])
        bottom_values = np.zeros(n_accessions)
        for mc, sg, cpd in cpd_order_plot:
            if cpd not in compounds_data_detailed:
                continue
            cinfo = compounds_data_detailed[cpd]
            pct = np.where(totals_cpd > 0, np.array(cinfo["values"]) / totals_cpd * 100, 0)
            ax9.bar(accession_names, pct, bottom=bottom_values,
                    color=cinfo["color"], alpha=0.8, width=width_bar)
            bottom_values += pct
        ax9.set_title(TT("Genotype Profiles: Compound Distribution"), fontsize=20, fontweight="bold", color="black", pad=20)
        ax9.set_ylabel("Percentage (%)", fontsize=16, fontweight="bold", color="black")
        ax9.set_xlabel("Cultivars", fontsize=16, fontweight="bold", color="black")
        ax9.set_ylim(0, 100)
        ax9.set_xlim(-0.5, n_accessions-0.5)
        ax9.spines["top"].set_visible(False)
        ax9.spines["right"].set_visible(False)
        ax9.spines["left"].set_color("black")
        ax9.spines["bottom"].set_color("black")
        ax9.tick_params(axis="x", colors="black", labelsize=16)
        ax9.tick_params(axis="y", colors="black", labelsize=16)
        plt.xticks(rotation=90)
        if CLUSTER_LEGEND:
            _color_labels(ax9)
            _cluster_legend(ax9)
        ax9.grid(False)
        legend_patches9 = []
        for mc in mc_order:
            for sg in sg_by_mc.get(mc, []):
                legend_patches9.append(mpatches.Patch(color=subclass_colors_palette[sg], label=f"── {sg}", alpha=1.0))
                for mc2, sg2, cpd2 in cpd_order:
                    if mc2 == mc and sg2 == sg and cpd2 in compounds_data_detailed:
                        legend_patches9.append(mpatches.Patch(
                            color=compounds_data_detailed[cpd2]["color"], label=f"  {cpd2}", alpha=0.85))
        ax9.legend(handles=legend_patches9, loc="upper left", bbox_to_anchor=(1.01, 1),
                   frameon=False, fontsize=8, title=input_set["compound_legend_title"], title_fontsize=10,
                   ncol=2)
        plt.tight_layout(rect=[0, 0, 0.75, 1])
        save_figure(os.path.join(output_folder, f"9_{sheet_name}_accession_all_compounds_percentages.png"), dpi=DPI, bbox_inches="tight")
        save_legend(os.path.join(output_folder, f"legend_compounds_{sheet_name}.png"), legend_patches9,
                    [patch.get_label() for patch in legend_patches9], input_set["compound_legend_title"], DPI,
                    fontsize=8, title_fontsize=10, ncol=2)
        plt.close("all")
        print("Plot 9 saved.")

        main_classes_order = [mc for mc in mc_order if mc in missing_df_plot["Main_Class"].unique()]
        FS_TITLE, FS_LABEL, FS_TICK, FS_BAR = 22, 18, 16, 14
        bar_width = 0.85

        fig_c, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(20, 8))

        bottom_vals = [0] * len(main_classes_order)
        for sg in unique_subclasss:
            counts = [len(missing_df_plot[(missing_df_plot["Main_Class"] == mc) &
                                          (missing_df_plot["Subclass"] == sg)])
                      for mc in main_classes_order]
            ax_left.bar(range(len(main_classes_order)), counts, bottom=bottom_vals,
                        color=subclass_colors_palette.get(sg, (0.5,0.5,0.5)),
                        edgecolor="white", linewidth=0.8, width=bar_width, alpha=0.9)
            for j, val in enumerate(counts):
                if val > 0:
                    ax_left.text(j, bottom_vals[j]+val/2, str(val), ha="center", va="center",
                                 fontsize=FS_BAR, color="white", fontweight="bold")
            bottom_vals = [sum(x) for x in zip(bottom_vals, counts)]
        ax_left.set_ylabel("Number of Compounds", fontsize=FS_LABEL, fontweight="bold", color="black")
        ax_left.set_xlabel("Main Chemical Classes", fontsize=FS_LABEL, fontweight="bold", color="black")
        ax_left.set_title("Compounds per Main Class and Subclass", fontsize=FS_TITLE, fontweight="bold", color="black")
        ax_left.set_xticks(range(len(main_classes_order)))
        ax_left.set_xticklabels(main_classes_order, rotation=45, ha="right", fontsize=FS_TICK, fontweight="bold")
        ax_left.tick_params(axis="y", colors="black", labelsize=FS_TICK)
        ax_left.spines["top"].set_visible(False)
        ax_left.spines["right"].set_visible(False)
        ax_left.spines["left"].set_color("black")
        ax_left.spines["bottom"].set_color("black")
        ax_left.grid(False)

        sg_counts_r = [(sg, len(missing_df_plot[missing_df_plot["Subclass"] == sg]))
                       for sg in unique_subclasss if len(missing_df_plot[missing_df_plot["Subclass"] == sg]) > 0]
        sg_labels_r = [x[0] for x in sg_counts_r]
        sg_vals_r   = [x[1] for x in sg_counts_r]
        sg_colors_r = [subclass_colors_palette[sg] for sg in sg_labels_r]
        ax_right.bar(range(len(sg_labels_r)), sg_vals_r, color=sg_colors_r, alpha=0.9, width=bar_width)
        for i, count in enumerate(sg_vals_r):
            ax_right.text(i, count/2, str(count), ha="center", va="center",
                          fontweight="bold", fontsize=FS_BAR, color="white")
        ax_right.set_title("Number of Compounds per Subclass", fontsize=FS_TITLE, fontweight="bold", color="black")
        ax_right.set_ylabel("Number of Compounds", fontsize=FS_LABEL, fontweight="bold", color="black")
        ax_right.set_xlabel("Subclass", fontsize=FS_LABEL, fontweight="bold", color="black")
        ax_right.set_xticks(range(len(sg_labels_r)))
        ax_right.set_xticklabels(sg_labels_r, rotation=45, ha="right", fontsize=FS_TICK, fontweight="bold")
        ax_right.tick_params(axis="y", colors="black", labelsize=FS_TICK)
        ax_right.spines["top"].set_visible(False)
        ax_right.spines["right"].set_visible(False)
        ax_right.spines["left"].set_color("black")
        ax_right.spines["bottom"].set_color("black")
        ax_right.grid(False)

        legend_patches_c = [mpatches.Patch(color=subclass_colors_palette[sg], label=sg) for sg in unique_subclasss]
        fig_c.legend(handles=legend_patches_c, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=11, title="Subclass", title_fontsize=12)
        fig_c.align_xlabels([ax_left, ax_right])
        plt.tight_layout()
        save_figure(os.path.join(output_folder, f"combined_subclass_mainclass_{sheet_name}.png"), dpi=DPI, bbox_inches="tight")
        plt.close("all")
        print("Combined figure saved.")

        gc.collect()
        print(f"\nDataset {sheet_name} complete.")

print("\nAll input sets processed.")

finish_step(cfg, TREATMENT)

#!/usr/bin/env python3
#
# VocProfile - 10_plot_detection.py
#
# Number of accessions in which each compound was detected, drawn in the fixed colour
# of each compound (Fig. 2c). A compound counts as detected in an accession when its
# accession mean is not missing.
#
# in   MEANS_FILE, DATASETS, COMPOUND_COLOUR_TABLE
# out  figures/10_plot_detection/detection/3_<dataset>_detailed_compound_view.png
#      results/10_plot_detection/10_qc.txt
#
# Usage:
#   python scripts/10_plot_detection.py <config.txt>
#
# Notes: docs/NOTES.md, step 10
#
import ast
import gc
import os
import re
import unicodedata
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

from lib_config import load_config, input_path, table_path, require_sheets, output_dir, die
from lib_run import start_step, finish_step, save_figure, later

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')

cfg, step, args = load_config("10")
start_step(cfg, 'means')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'values without detection are left out of every statistic; nothing is filled',
    'the step stops if a compound has no colour in COMPOUND_COLOUR_TABLE',
)
file_path = input_path(cfg, step["input"])
color_file = table_path(cfg, "COMPOUND_COLOUR_TABLE")
sheets = step["sheets"]
require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets])
output_folder = output_dir(cfg, "10_plot_detection", step["output_folder"])

mc_order = list(cfg["main_class_colours"])
mc_order_idx = {mc: i for i, mc in enumerate(mc_order)}
GRAY = (0.5, 0.5, 0.5, 1.0)

ALIASES = {
    "Germacrene β": "Germacrene B",
}

DISPLAY_NAMES = {
    "Germacrene β": "Germacrene B",
}


def norm(s):
    """Compound name reduced to a key that ignores Greek/Latin spelling, case and punctuation."""
    s = unicodedata.normalize("NFKC", str(s)).replace("ß", "β").casefold()
    for g, l in {"α": "a", "β": "b", "γ": "g", "δ": "d", "ε": "e"}.items():
        s = s.replace(g, l)
    return re.sub(r"[\s\-_()\[\],;:'’]", "", s)


ALIASES_NORM = {norm(k): norm(v) for k, v in ALIASES.items()}


def parse_color(v):
    """hex string, 'r,g,b' string, or tuple/list -> RGBA tuple"""
    if isinstance(v, str):
        v = v.strip()
        if v.startswith("#") or v.isalpha():
            return to_rgba(v)
        try:
            v = ast.literal_eval(v)
        except Exception:
            v = tuple(float(x) for x in v.strip("()[] ").split(","))
    v = tuple(float(x) for x in v)
    if max(v[:3]) > 1:
        v = tuple(x / 255 for x in v[:3]) + tuple(v[3:])
    return to_rgba(v)


def load_color_map(path):
    d = pd.read_csv(path, sep="\t")
    d.columns = [str(c).strip() for c in d.columns]
    low = {c.lower(): c for c in d.columns}
    col_color = next(low[k] for k in low if any(t in k for t in ("color", "colour", "hex", "rgb")))
    col_cpd = next(low[k] for k in low if "compound" in k)
    cmap = {}
    for _, r in d.dropna(subset=[col_color, col_cpd]).iterrows():
        cmap[norm(r[col_cpd])] = parse_color(r[col_color])
    return cmap


def get_color(cpd, sg, missing_colors):
    k = norm(cpd)
    k = ALIASES_NORM.get(k, k)
    if k in cpd_color_map:
        return cpd_color_map[k]
    missing_colors.add((str(sg).strip(), str(cpd).strip()))
    return GRAY


print("Loading color map...")
cpd_color_map = load_color_map(color_file)


for sheet_name in sheets:
    print(f"\nDataset: {sheet_name}")
    try:
        df = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], header=cfg["header_rows"], index_col=0)
        df.columns.names = ["Main_Class", "Subclass", "Compound"]
        print(f"Sheet '{sheet_name}' loaded successfully.")
    except Exception as e:
        print(f"Error loading sheet '{sheet_name}': {e}")
        continue

    df = df.dropna(axis=0, how="all")
    numeric_df = df.select_dtypes(include=[np.number])
    total_samples = len(numeric_df)

    missing_counts = numeric_df.isnull().sum(axis=0)
    missing_df_plot = missing_counts.to_frame(name="Missing_Count").reset_index()
    missing_df_plot.columns = ["Main_Class", "Subclass", "Compound", "Missing_Count"]
    missing_df_plot["Present_Count"] = total_samples - missing_df_plot["Missing_Count"]

    sg_by_mc = {}
    for mc, sg, cpd in numeric_df.columns:
        sg_by_mc.setdefault(mc, [])
        if sg not in sg_by_mc[mc]:
            sg_by_mc[mc].append(sg)
    unique_subclasss = []
    for mc in mc_order:
        for sg in sg_by_mc.get(mc, []):
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
    missing_df_plot["cpd_order"] = missing_df_plot["Compound"].map(cpd_order_idx)
    missing_df_sorted = missing_df_plot.sort_values(
        ["mc_order", "sg_order", "cpd_order"]).reset_index(drop=True)

    missing_colors = set()
    compound_colors3 = [get_color(c, s, missing_colors) for c, s in
                        zip(missing_df_sorted["Compound"], missing_df_sorted["Subclass"])]
    if missing_colors:
        die(f"{len(missing_colors)} compound(s) have no colour in {color_file}: {sorted(missing_colors)}")

    fig3, ax3 = plt.subplots(figsize=(28, 10))
    bar_w3 = 0.92
    ax3.bar(range(len(missing_df_sorted)), missing_df_sorted["Present_Count"],
            color=compound_colors3, alpha=0.9, width=bar_w3)
    ax3.bar(range(len(missing_df_sorted)), missing_df_sorted["Missing_Count"],
            bottom=missing_df_sorted["Present_Count"], color="lightgray", alpha=0.6, width=bar_w3)
    ax3.set_ylabel(later(cfg["settings"].tables_label("DETECTION_Y_AXIS", "submitted"), cfg["labels"]["fig2c_y_axis"]), fontsize=20, fontweight="bold", color="black")
    ax3.set_xlabel("Compounds", fontsize=20, fontweight="bold", color="black")
    ax3.set_xticks(range(len(missing_df_sorted)))
    ax3.set_xticklabels([DISPLAY_NAMES.get(c, c) for c in missing_df_sorted["Compound"]],
                        rotation=90, ha="center", fontsize=20, fontweight="bold")
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
    for i in range(1, len(missing_df_sorted)):
        if missing_df_sorted["Main_Class"].iloc[i] != current_class:
            ax3.axvline(x=i - 0.5, color="black", linestyle="--", linewidth=2, alpha=0.7)
            current_class = missing_df_sorted["Main_Class"].iloc[i]

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(output_folder, f"3_{sheet_name}_detailed_compound_view.png")
    save_figure(out_path, dpi=step["dpi"], bbox_inches="tight")
    plt.close("all")
    print(f"Plot 3 saved: {out_path}")

    gc.collect()

print("\nAll datasets done.")

finish_step(cfg, TREATMENT)

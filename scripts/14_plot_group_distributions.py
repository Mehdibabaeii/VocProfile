#!/usr/bin/env python3
#
# VocProfile - 14_plot_group_distributions.py
#
# Per input set and dataset: violin and box plots by group with Tukey letters,
# clustered heat map of clipped z-scores (accessions by compounds), PCA, and the
# clustered Pearson correlation matrix between compounds. The input set decides the
# grouping (K groups or clusters) and the transformation: filled and ln(1 + x),
# filled, or none. The correlation matrix of clusters_log1p_fig3b is Fig. 3b.
#
# in   MEANS_FILE, MEANS_CLUSTER_FILE, INPUT_SETS_TABLE (step 14), FIGURE_STYLES_TABLE
#      CLIP_LIMIT, TUKEY_ALPHA, JITTER_SEED, NONDETECT_FRACTION, COLOUR_TABLE, CLUSTER_TABLE
# out  figures/14_plot_group_distributions/<input set>/plots_<dataset>/
#      results/14_plot_group_distributions/14_qc.txt
#
# Usage:
#   python scripts/14_plot_group_distributions.py <config.txt>
#
# Notes: docs/NOTES.md, step 14
#
import gc
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
from matplotlib.transforms import Bbox
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage as scipy_linkage
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from lib_config import load_config, input_path, require_sheets, output_dir, selected_sets, term_function, display_names
from lib_run import start_step, finish_step, save_figure, set_scope, save_legend, note

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')

cfg, step, args = load_config("14")
start_step(cfg, 'means', 'means_clusters')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'log1p and filled input sets: accession means without detection: set to NONDETECT_FRACTION of the smallest positive accession mean of the compound',
    'log1p input sets: values then transformed as ln(1 + x)',
    'input sets without filling: box plots and correlations use the values present; the heat map sets values without detection to zero before the z-scores; the PCA uses accessions without missing values',
    'heat map: z-score per compound, clipped at CLIP_LIMIT',
    'the jitter of the points is seeded with JITTER_SEED before each box-plot figure',
)
ABUNDANCE_SHEET = cfg["abundance_sheet"]
DPI = step["dpi"]
T = term_function(cfg)
DISPLAY_NAMES = display_names(cfg)
JITTER_SEED = cfg["random"]["jitter_seed"]   # None: point jitter differs from run to run


def build_palette(df_cols):
    all_sg = list(cfg["subclass_order"]["distributions"])
    present = list(dict.fromkeys(df_cols.get_level_values('Subclass')))
    unique_sg = [sg for sg in all_sg if sg in present]
    for sg in present:
        if sg not in unique_sg:
            unique_sg.append(sg)
    colors = sns.color_palette("husl", len(unique_sg))
    sg_pal = dict(zip(unique_sg, colors))
    mc_pal = dict(cfg["main_class_colours"])
    return sg_pal, mc_pal, unique_sg


def detect_binary_compounds(df, compound_cols):
    binaries = []
    for col in compound_cols:
        vals = pd.Series(df[col].dropna().unique())
        non_nan = set([v for v in vals if pd.notna(v)])
        if non_nan and non_nan.issubset({0, 1, 0.0, 1.0}):
            binaries.append(col)
    return binaries


def get_tukey_groups(df, col, group_col='K', alpha=step["tukey_alpha"]):
    """Tukey HSD letters per group; groups that do not differ share a letter."""
    temp_df = df[[col, group_col]].dropna()
    if len(temp_df[group_col].unique()) < 2 or len(temp_df) < 3:
        return None
    try:
        tukey = pairwise_tukeyhsd(endog=temp_df[col], groups=temp_df[group_col], alpha=alpha)
    except ValueError:
        return None
    if not np.any(tukey.pvalues < alpha):
        return None
    tukey_summary = pd.DataFrame(tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    non_sig_pairs = set()
    for _, row in tukey_summary.iterrows():
        if row['p-adj'] > alpha:
            non_sig_pairs.add(tuple(sorted([row['group1'], row['group2']])))
    group_means = temp_df.groupby(group_col)[col].mean()
    groups = group_means.sort_values(ascending=False).index.tolist()
    letters = 'abcdefghijklmnopqrstuvwxyz'
    assigned_letters = {g: [] for g in groups}
    letter_idx = 0
    for current_group in groups:
        if not assigned_letters[current_group]:
            current_letter = letters[letter_idx]
            groups_for_this_letter = [current_group]
            for other_group in groups:
                if other_group != current_group:
                    pair = tuple(sorted([current_group, other_group]))
                    if pair in non_sig_pairs:
                        groups_for_this_letter.append(other_group)
            for g in groups_for_this_letter:
                if current_letter not in assigned_letters[g]:
                    assigned_letters[g].append(current_letter)
            letter_idx += 1
            if letter_idx >= len(letters):
                break
    final_assignments = {}
    for group, letter_list in assigned_letters.items():
        final_assignments[group] = ''.join(sorted(set(letter_list)))
    return final_assignments


def perform_all_tukey_tests(df, compound_cols, group_col='K'):
    results = {}
    for col in compound_cols:
        results[col] = get_tukey_groups(df, col, group_col)
    return results


def plot_traits_by_K(df_data, group_traits, group_name, k_colors_map, df_original_multi_index,
                      tukey_results, binary_traits, unit_label, output_folder, group_axis_label):
    group_traits = [t for t in group_traits if t not in binary_traits]
    if not group_traits:
        print("No non-binary traits to plot for violins.")
        return
    n_traits = len(group_traits)
    n_cols = 8
    n_rows = (n_traits + n_cols - 1) // n_cols
    fig_width = n_cols * 5.0
    fig_height = max(5.0, n_rows * 5.0)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()

    df_data = df_data.copy()
    df_data['K'] = pd.Categorical(df_data['K'], categories=sorted(df_data['K'].unique()), ordered=True)
    sorted_k_categories = df_data['K'].cat.categories.tolist()
    palette = {k: k_colors_map.get(k, '#333333') for k in sorted_k_categories}

    for i, trait in enumerate(group_traits):
        if i >= len(axes):
            break
        ax = axes[i]
        if trait in df_data.columns and pd.api.types.is_numeric_dtype(df_data[trait]):
            temp_df = df_data[['K', trait]].dropna(subset=[trait])
            if temp_df.empty:
                fig.delaxes(ax)
                continue
            violin_parts = sns.violinplot(
                x='K', y=trait, data=temp_df, ax=ax,
                inner=None, linewidth=2, order=sorted_k_categories,
                color='white'
            )
            for j, violin in enumerate(violin_parts.collections):
                color = palette.get(sorted_k_categories[j % len(sorted_k_categories)], '#333333')
                violin.set_edgecolor(color)
                violin.set_facecolor('none')
                violin.set_alpha(1)
            sns.boxplot(
                x='K', y=trait, data=temp_df, ax=ax,
                width=0.2,
                boxprops={'facecolor':'none', 'edgecolor':'gray', 'linewidth':1.5},
                medianprops={'color':'black', 'linestyle':'-', 'linewidth':2},
                whiskerprops={'color':'gray', 'linewidth':1.2},
                capprops={'color':'gray', 'linewidth':1.2},
                showcaps=True, showfliers=False,
                order=sorted_k_categories
            )
            sns.stripplot(
                x='K', y=trait, data=temp_df, ax=ax,
                palette=palette, jitter=0.25, alpha=0.7, size=4,
                order=sorted_k_categories
            )
            if trait in tukey_results and tukey_results[trait] is not None:
                y_max = temp_df[trait].max()
                y_pos = y_max * 1.15
                for k_idx, group in enumerate(sorted_k_categories):
                    letters = tukey_results[trait].get(group, '')
                    ax.text(k_idx, y_pos, letters, ha='center', va='bottom', fontsize=14, fontweight='bold')
            ax.grid(False)
            for spine in ax.spines.values():
                spine.set_edgecolor("black")
                spine.set_linewidth(1)
            ax.set_ylim(0, temp_df[trait].max() * 1.25)
            ax.set_ylabel('')
            ax.tick_params(axis='y', labelsize=14)
            if i % n_cols == 0:
                ax.set_ylabel(unit_label, fontsize=16, fontweight='bold')
            ax.set_xlabel('')
            ax.tick_params(axis='x', labelsize=14)
            if i + n_cols >= n_traits:
                ax.set_xlabel(group_axis_label, fontsize=16, fontweight='bold')
            ax.set_title(trait, fontsize=16, fontweight='bold', pad=6)
        else:
            fig.delaxes(ax)

    for i in range(n_traits, len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()
    save_path = os.path.join(output_folder, f'{group_name.replace(" ","_")}_distributions.png')
    save_figure(save_path, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    gc.collect()
    print(f"Violin plots saved: {save_path}")


def plot_compound_heatmap(
    df_data,
    compound_traits,
    k_colors_map,
    df_original_multi_index,
    main_class_colors_palette,
    subclass_colors_palette,
    output_folder,
    style
):
    compound_traits = [c for c in compound_traits if c in df_data.columns]
    temp_df = df_data[compound_traits + ['K']].copy()
    if temp_df.empty:
        print("Skipping heatmap - no data.")
        return

    X = StandardScaler().fit_transform(temp_df[compound_traits].fillna(0))
    scaled_df = pd.DataFrame(X, columns=compound_traits, index=temp_df.index)
    scaled_df = scaled_df.clip(lower=-step["clip_limit"], upper=step["clip_limit"])

    row_colors = temp_df['K'].map(k_colors_map)

    cols_to_use = [col for col in df_original_multi_index.columns if col[2] in compound_traits]
    if not cols_to_use:
        print("Warning: No matching columns found for color mapping.")
        return
    meta_df = pd.DataFrame(cols_to_use, columns=['Main_Class', 'Subclass', 'Compound'])
    col_colors_main = meta_df['Main_Class'].map(main_class_colors_palette)
    col_colors_sub = meta_df['Subclass'].map(subclass_colors_palette)
    col_colors_main.fillna('lightgray', inplace=True)
    col_colors_sub.fillna('lightgray', inplace=True)

    fig_w = max(50, len(compound_traits) * 1.5)
    fig_h = max(70, scaled_df.shape[0] * 0.7)

    row_linkage = scipy_linkage(pdist(scaled_df.values, metric='euclidean'), method='ward')
    col_linkage = scipy_linkage(pdist(scaled_df.values.T, metric='euclidean'), method='ward')

    g = sns.clustermap(
        scaled_df,
        row_cluster=True,
        col_cluster=True,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_colors=row_colors,
        col_colors=[col_colors_main.values, col_colors_sub.values],
        cmap='coolwarm',
        linewidths=0.5,
        cbar_kws={'label': 'Z-score'},
        figsize=(fig_w, fig_h),
        dendrogram_ratio=(0.1, 0.1),
        colors_ratio=(0.02, 0.02)
    )

    for ax in [g.ax_row_dendrogram, g.ax_col_dendrogram]:
        for line in ax.collections:
            line.set_color('black')
            line.set_linewidth(2.5)

    row_colors_pos = g.ax_row_colors.get_position()
    label_x = row_colors_pos.x0 + row_colors_pos.width / 2
    label_y = row_colors_pos.y0 - 0
    g.fig.text(label_x, label_y, "K", ha='center', va='top', fontsize=25, fontweight='bold')

    col_colors_pos = g.ax_col_colors.get_position()
    label_x = col_colors_pos.x1 + 0.01
    main_class_y = col_colors_pos.y0 + col_colors_pos.height * 0.25
    subclass_y = col_colors_pos.y0 + col_colors_pos.height * 0.75
    g.fig.text(label_x, main_class_y, style["subclass_label"], rotation=0, va='center', ha='left',
               fontsize=style["strip_label_size"], transform=g.fig.transFigure)
    g.fig.text(label_x, subclass_y, style["main_class_label"], rotation=0, va='center', ha='left',
               fontsize=style["strip_label_size"], transform=g.fig.transFigure)

    g.ax_cbar.set_position([1.05, 0.3, 0.015, 0.4])
    g.ax_cbar.set_ylabel('Z-score', fontsize=style["colourbar_label_size"], fontweight='bold')
    g.ax_cbar.tick_params(labelsize=style["colourbar_tick_size"])

    compound_labels = g.ax_heatmap.get_xticklabels()
    if style["display_names"]:
        compound_labels = [DISPLAY_NAMES.get(t.get_text(), t.get_text()) for t in compound_labels]
    g.ax_heatmap.set_xticklabels(compound_labels, rotation=45, ha='right', fontsize=style["compound_label_size"])
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=style["accession_label_size"])
    g.ax_heatmap.set_xlabel('Compounds', fontsize=style["axis_label_size"], fontweight='bold')
    g.ax_heatmap.set_ylabel(T(style["row_axis_label"]), fontsize=style["axis_label_size"], fontweight='bold')

    if style["legend"]:
        draw_heatmap_legend(g, k_colors_map, main_class_colors_palette, subclass_colors_palette)

    save_path = os.path.join(output_folder, step["heatmap_file"])
    save_figure(save_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close('all')
    gc.collect()
    print(f"Main heatmap saved: {save_path}")


def draw_heatmap_legend(g, k_colors_map, main_class_colors_palette, subclass_colors_palette):
    leg_k = [mpatches.Patch(facecolor=k_colors_map[k], edgecolor='black', linewidth=1, label=f"Clade {k}")
             for k in sorted(k_colors_map.keys())]
    leg_main = [mpatches.Patch(facecolor=main_class_colors_palette[m], edgecolor='black', linewidth=1, label=f"Main: {m}")
                for m in main_class_colors_palette.keys() if m not in ['Main_Class', None]]
    leg_sub = [mpatches.Patch(facecolor=subclass_colors_palette[s], edgecolor='black', linewidth=1, label=f"Sub: {s}")
               for s in subclass_colors_palette.keys() if s not in ['Subclass', None]]
    spacer = mpatches.Patch(facecolor='none', edgecolor='none')
    all_handles = leg_k + [spacer] + leg_main + [spacer] + leg_sub

    g.fig.legend(
        handles=all_handles, labels=[h.get_label() for h in all_handles],
        title="Classification", bbox_to_anchor=(0, 0.99), loc='upper left',
        fontsize=30, title_fontsize=40, frameon=False, fancybox=False, shadow=False, ncol=1
    )


def two_legends(figure, main_handles, subclass_handles, style):
    """Main classes and subclasses as two legends with their own titles, one under the other."""
    options = dict(fontsize=19, title_fontsize=20, frameon=False, ncol=1)
    first = figure.legend(handles=main_handles, labels=[h.get_label().split(": ", 1)[-1] for h in main_handles],
                          title=style["legend_main_title"], bbox_to_anchor=(0, 1), loc='upper left', **options)
    figure.canvas.draw()
    lower_edge = figure.transFigure.inverted().transform(first.get_window_extent().p0)[1]
    figure.legend(handles=subclass_handles, labels=[h.get_label().split(": ", 1)[-1] for h in subclass_handles],
                  title=style["legend_subclass_title"], bbox_to_anchor=(0, lower_edge), loc='upper left', **options)


def enlarge_compound_names(g, box, style):
    """Largest size up to compound_label_fit_max at which the compound names stay on the canvas of the
    submitted figure and clear of each other, of the axis title and of the colour bar."""
    largest, base = style["compound_label_fit_max"], style["compound_label_size"]
    if not largest or largest <= base or box is None:
        return
    figure, axes = g.fig, g.ax_heatmap
    renderer = figure.canvas.get_renderer()
    title = axes.xaxis.label
    placed = title.get_window_extent(renderer)
    axes.xaxis.set_label_coords(*axes.transAxes.inverted().transform(((placed.x0 + placed.x1) / 2, placed.y1)))
    canvas = Bbox(box.get_points() * figure.dpi)
    pitch = axes.get_window_extent(renderer).height / max(len(axes.get_yticklabels()), 1)
    columns, rows = axes.get_xticklabels(), axes.get_yticklabels()
    for size in range(int(largest), int(base), -1):
        for label in columns + rows:
            label.set_fontsize(size)
        figure.canvas.draw()
        extents = [label.get_window_extent(renderer) for label in columns + rows]
        text_height = size * figure.dpi / 72
        inside = all(canvas.x0 <= e.x0 and e.x1 <= canvas.x1 and canvas.y0 <= e.y0 and e.y1 <= canvas.y1 for e in extents)
        clear_of_title = not any(e.overlaps(title.get_window_extent(renderer)) for e in extents[:len(columns)])
        clear_of_bar = not any(e.overlaps(g.ax_cbar.get_tightbbox(renderer)) for e in extents[len(columns):])
        apart = text_height <= pitch * 0.7071
        if inside and clear_of_title and clear_of_bar and apart:
            note(f"correlation matrix: compound names at {size} pt (submitted: {base} pt)")
            return
    for label in columns + rows:
        label.set_fontsize(base)
    note(f"correlation matrix: compound names stay at {base} pt, a larger size does not fit")


def plot_correlation_heatmap(
    df_data,
    df_original_multi_index,
    main_class_colors_palette,
    subclass_colors_palette,
    output_folder,
    style
):
    if 'K' in df_data.columns:
        df_corr_calc = df_data.drop(columns='K')
    else:
        df_corr_calc = df_data.copy()

    corr_matrix = df_corr_calc.corr(method='pearson')

    cols_to_use = [col for col in df_original_multi_index.columns if col[2] in corr_matrix.columns]
    if not cols_to_use:
        print("Warning: No matching columns found for color mapping.")
        return
    meta_df = pd.DataFrame(cols_to_use, columns=['Main_Class', 'Subclass', 'Compound'])
    col_colors_main = meta_df['Main_Class'].map(main_class_colors_palette)
    col_colors_sub = meta_df['Subclass'].map(subclass_colors_palette)
    col_colors_main.fillna('lightgray', inplace=True)
    col_colors_sub.fillna('lightgray', inplace=True)
    row_colors_main = col_colors_main.copy()
    row_colors_sub = col_colors_sub.copy()

    fig_w = max(50, len(corr_matrix.columns) * 0.8)
    fig_h = max(50, len(corr_matrix.index) * 0.5)

    row_linkage = scipy_linkage(pdist(corr_matrix.values, metric='euclidean'), method='ward')
    col_linkage = scipy_linkage(pdist(corr_matrix.values.T, metric='euclidean'), method='ward')

    g = sns.clustermap(
        corr_matrix,
        row_cluster=True,
        col_cluster=True,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_colors=[row_colors_main.values, row_colors_sub.values],
        col_colors=[col_colors_main.values, col_colors_sub.values],
        cmap='coolwarm',
        vmin=-1, vmax=1,
        linewidths=0.5,
        cbar_kws={'label': 'Pearson Correlation'},
        figsize=(fig_w, fig_h),
        dendrogram_ratio=(0.1, 0.1),
        colors_ratio=(0.02, 0.02)
    )

    for ax in [g.ax_row_dendrogram, g.ax_col_dendrogram]:
        for line in ax.collections:
            line.set_color('black')
            line.set_linewidth(2.5)

    col_colors_pos = g.ax_col_colors.get_position()
    label_x = col_colors_pos.x1 + 0.01
    main_class_y = col_colors_pos.y0 + col_colors_pos.height * 0.25
    subclass_y = col_colors_pos.y0 + col_colors_pos.height * 0.75
    g.fig.text(label_x, main_class_y, style["subclass_label"], rotation=0, va='center', ha='left',
               fontsize=style["strip_label_size"], transform=g.fig.transFigure)
    g.fig.text(label_x, subclass_y, style["main_class_label"], rotation=0, va='center', ha='left',
               fontsize=style["strip_label_size"], transform=g.fig.transFigure)

    g.ax_cbar.set_position([style["colourbar_left"], 0.3, 0.015, 0.4])
    g.ax_cbar.set_ylabel('Pearson Correlation', fontsize=35, fontweight='bold')
    g.ax_cbar.tick_params(labelsize=30)

    x_labels = g.ax_heatmap.get_xticklabels()
    y_labels = g.ax_heatmap.get_yticklabels()
    if style["display_names"]:
        x_labels = [DISPLAY_NAMES.get(t.get_text(), t.get_text()) for t in x_labels]
        y_labels = [DISPLAY_NAMES.get(t.get_text(), t.get_text()) for t in y_labels]
    g.ax_heatmap.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=style["compound_label_size"])
    g.ax_heatmap.set_yticklabels(y_labels, fontsize=style["compound_label_size"])
    g.ax_heatmap.set_xlabel('Compounds', fontsize=style["axis_label_size"], fontweight='bold')

    leg_main = [mpatches.Patch(facecolor=main_class_colors_palette[m], edgecolor='black', linewidth=1, label=f"Main: {m}")
                for m in main_class_colors_palette.keys() if m not in ['Main_Class', None]]
    leg_sub = [mpatches.Patch(facecolor=subclass_colors_palette[s], edgecolor='black', linewidth=1, label=f"Sub: {s}")
               for s in subclass_colors_palette.keys() if s not in ['Subclass', None]]
    spacer = mpatches.Patch(facecolor='none', edgecolor='none')
    all_handles = leg_main + [spacer] + leg_sub

    if style["legend_title"]:
        g.fig.legend(
            handles=all_handles, labels=[h.get_label() for h in all_handles],
            title=style["legend_title"], bbox_to_anchor=(0, 1), loc='upper left',
            fontsize=19, title_fontsize=20, frameon=False, fancybox=False, shadow=False, ncol=1
        )
    else:
        two_legends(g.fig, leg_main, leg_sub, style)
    if style["legend_files"]:
        for name, handles, title in (("main_class", leg_main, style["legend_main_title"] or style["main_class_label"]),
                                     ("subclass", leg_sub, style["legend_subclass_title"] or style["subclass_label"])):
            save_legend(os.path.join(output_folder, f"legend_{name}.png"), handles,
                        [h.get_label().split(": ", 1)[-1] for h in handles], title, DPI, fontsize=19, title_fontsize=20)

    save_path = os.path.join(output_folder, 'compound_correlation_heatmap.png')
    save_figure(save_path, panel=bool(style["panel"]), dpi=DPI, bbox_inches='tight', facecolor='white',
                after_layout=lambda figure, box: enlarge_compound_names(g, box, style))
    plt.close('all')
    gc.collect()
    print(f"Pearson correlation heatmap saved: {save_path}")


def plot_pca_by_K(df_data, compound_traits, k_colors_map, output_folder):
    temp_df = df_data[compound_traits + ['K']].dropna()
    if temp_df.empty:
        print("Skipping PCA - no data.")
        return
    X = StandardScaler().fit_transform(temp_df[compound_traits])
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    pca_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'], index=temp_df.index)
    pca_df['K'] = temp_df['K']
    plt.figure(figsize=(12, 10))
    sns.scatterplot(x='PC1', y='PC2', hue='K', style='K', data=pca_df, palette=k_colors_map, s=130, alpha=0.9)
    plt.title('PCA of Compounds by K', fontsize=20, fontweight='bold')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=16)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=12, title_fontsize=14)
    plt.tight_layout()
    save_path = os.path.join(output_folder, 'compound_pca.png')
    save_figure(save_path, dpi=DPI, bbox_inches='tight')
    plt.close('all')
    gc.collect()
    print(f"PCA plot saved: {save_path}")

for set_name, input_set in selected_sets(step, args).items():
    file_path = input_path(cfg, input_set["input"])
    sheets_to_run = input_set["sheets"]
    require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets_to_run])
    output_folder_base = os.path.join(output_dir(cfg, "14_plot_group_distributions", set_name), input_set["output_folder"])
    k_colors = input_set["group_colours"]
    heatmap_style = input_set["heatmap"]
    transform = input_set["transform"]
    set_scope(input_set["grouping"])

    for sheet_name in sheets_to_run:
        print(f"\nInput set {set_name}, dataset {sheet_name}")
        output_folder = f"{output_folder_base}_{sheet_name}"
        os.makedirs(output_folder, exist_ok=True)
        unit_label = step["unit_labels"][transform]["abundance"] if sheet_name == ABUNDANCE_SHEET else step["unit_labels"][transform]["concentration"]

        try:
            df_raw = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], header=cfg["header_rows"],
                                   index_col=[0, 1], engine='openpyxl')
            df_raw.index.names = ['Acc', 'K']
            df_raw.columns.names = ['Main_Class', 'Subclass', 'Compound']
            print(f"Data loaded: {df_raw.shape}")
        except Exception as e:
            print(f"Error: {e}")
            continue
        df_raw = df_raw.dropna(axis=0, how='all')

        k_groups_series = df_raw.index.get_level_values('K').astype(str)
        df_numeric_only = df_raw.copy()
        df_numeric_only = df_numeric_only.reset_index(level='K', drop=True)

        if transform in ("log1p", "filled"):
            col_mins = df_numeric_only.min(skipna=True)
            fill_vals = (col_mins * cfg["settings"].value("NONDETECT_FRACTION")).where(col_mins > 0, 1e-9)
            df_numeric_only = df_numeric_only.fillna(fill_vals)
            print(f"Means without detection filled. Remaining: {df_numeric_only.isna().sum().sum()}")
        if transform == "log1p":
            df_numeric_only = np.log1p(df_numeric_only)

        subclass_colors_palette, main_class_colors_palette, unique_subclasss = build_palette(df_numeric_only.columns)

        df_plot = df_numeric_only.copy()
        df_plot.columns = df_plot.columns.get_level_values('Compound')
        compound_columns = df_plot.columns.to_list()

        df_combined_for_plotting = df_plot.copy()
        df_combined_for_plotting['K'] = k_groups_series.values

        binary_compounds = detect_binary_compounds(df_combined_for_plotting, compound_columns)
        nonbinary_compounds = [c for c in compound_columns if c not in binary_compounds]
        print(f"Compounds: {len(compound_columns)} ({len(binary_compounds)} binary, {len(nonbinary_compounds)} continuous)")

        if compound_columns:
            tukey_results = perform_all_tukey_tests(df_combined_for_plotting, nonbinary_compounds)

            if JITTER_SEED is not None:
                np.random.seed(JITTER_SEED)
            plot_traits_by_K(df_combined_for_plotting, nonbinary_compounds, "Compounds",
                             k_colors, df_raw, tukey_results, binary_compounds, unit_label, output_folder,
                             input_set["group_axis_label"])

            plot_compound_heatmap(
                df_combined_for_plotting, compound_columns, k_colors, df_raw,
                main_class_colors_palette, subclass_colors_palette,
                output_folder=output_folder, style=heatmap_style
            )

            plot_pca_by_K(df_combined_for_plotting, compound_columns, k_colors, output_folder)

            plot_correlation_heatmap(
                df_combined_for_plotting, df_raw,
                main_class_colors_palette, subclass_colors_palette,
                output_folder=output_folder, style=input_set["correlation"]
            )

        plt.close('all')
        gc.collect()
        print(f"\n{sheet_name} complete.")

print("\nAll datasets processed.")

finish_step(cfg, TREATMENT)

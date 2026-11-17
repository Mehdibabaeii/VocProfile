#!/usr/bin/env python3
#
# VocProfile - 16_compute_odor_index.py
#
# Relative odor index of the compounds with a published threshold: accession mean of
# the runs with a detection divided by the odor threshold, cluster means, the panel-
# wide value weighted by cluster size, tests across clusters, scalings and their
# figures (Figs. S5-S6, Tables S7, S18).
#
# in   REPLICATE_FILE, MAIN_DATASET, BLUP_CLUSTER_FILE (cluster of each accession)
#      ODOR_THRESHOLD_TABLE, CLUSTER_TABLE, SEED
# out  results/16_compute_odor_index/odor_index_*.xlsx and other tables
#      figures/16_compute_odor_index/plot*.png
#      results/16_compute_odor_index/16_qc.txt
#
# Usage:
#   python scripts/16_compute_odor_index.py <config.txt>
#
# Notes: docs/NOTES.md, step 16
#
import os
import warnings
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch as MPatch
from scipy import stats
from scipy.cluster.hierarchy import linkage as scipy_linkage
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

from lib_config import load_config, input_path, require_sheets, output_dir, term_function
from lib_run import start_step, finish_step, save_figure, write_table, check_same_accessions, use_font, read_own_table, revised, deferred

warnings.filterwarnings('ignore')

cfg, step, args = load_config("16")
start_step(cfg, 'replicates', 'blup_clusters')
LEVELS = cfg["significance_levels"]
check_same_accessions(cfg, 'replicates', 'blup_clusters')
TREATMENT = (
    'accession value: mean of the runs with a detection; nothing is filled',
    'cluster mean: accessions with a detection; panel-wide value: cluster means weighted by cluster size',
    'Pareto-scaled values without detection are set to zero, the mean of the detected accessions',
)
SEED = cfg["seed"]
SHEET = step["sheet"]
LABEL_FILE = input_path(cfg, step["labels"])
DPI = step["dpi"]
DPI_HIGH = step["dpi_high"]
OUT_DIR = output_dir(cfg, "16_compute_odor_index")


T = term_function(cfg)                                  # table headings
TT = deferred(term_function(cfg, "revised_terms_in_titles"))      # titles
LABELS = cfg["labels"]
ODOR_WORDING = LABELS["odor_index_wording"] if LABELS["revised_odor_wording"] else []


def out_path(name):
    return os.path.join(OUT_DIR, name)


def revised_odor_wording(text):
    for old, new in ODOR_WORDING:
        text = text.replace(old, new)
    return text


OI = deferred(revised_odor_wording)     # axis labels, colour-bar labels and titles


def bh_significance(stats_results):
    """Significance marks from Benjamini-Hochberg adjusted p-values across the compounds."""
    names = list(stats_results)
    p = np.array([stats_results[c]['p_value'] for c in names])
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    adjusted = {}
    for i, value in zip(order, np.minimum(q, 1)):
        mark = '***' if value < LEVELS[0] else '**' if value < LEVELS[1] else '*' if value < LEVELS[2] else 'ns'
        adjusted[names[i]] = {'p_value': stats_results[names[i]]['p_value'], 'q_value': float(value),
                              'significance': mark}
    return adjusted


CONFIG = {
    'random_seed': SEED,
    'top_n_for_profiles': 15,
    'colors_class': dict(cfg["cluster_colours"])
}

np.random.seed(CONFIG['random_seed'])

plt.rcParams.update({
    'font.size': 15,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'legend.fontsize': 13,
    'font.family': use_font(step["font_family"]),
    'axes.linewidth': 1.5,
})

ODT_DATABASE = cfg["odor_thresholds"]


def place_cluster_legend(g, handles):
    """Cluster legend of a clustered heat map. Submitted: in the upper left corner, over the colour bar.
    Revised: in one row above the column dendrogram, and the strip named above its upper end."""
    options = dict(handles=handles, title='Chemotype', fontsize=13, title_fontsize=14, frameon=False)
    if not revised():
        g.fig.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98), **options)
        return
    heatmap, top = g.ax_heatmap.get_position(), g.ax_col_dendrogram.get_position()
    g.fig.legend(loc='lower left', bbox_to_anchor=(heatmap.x0, top.y1 + 0.004), ncol=len(handles), **options)
    strip = g.ax_row_colors
    strip.xaxis.tick_top()
    strip.set_xticks(strip.get_xticks())
    strip.set_xticklabels(['Cluster'] * len(strip.get_xticks()), rotation=90)
    strip.tick_params(axis='x', length=0)


def clean_compound_name(name):
    name = str(name).strip()
    name = name.replace(', (', ' (')
    name = ' '.join(name.split())
    return name


def load_data(file_path, label_path=LABEL_FILE):
    print("Loading replicate data (non-detects left missing)...")

    raw = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][SHEET],
                       na_values=cfg["nondetect_code"], engine="openpyxl")
    compound_cols = [c for c in raw.columns if c not in (cfg["accession_column"], cfg["replicate_column"])]
    raw[compound_cols] = raw[compound_cols].apply(pd.to_numeric, errors="coerce")

    # accession value: mean of the replicates with a detection; non-detects are not filled
    df = raw.groupby(cfg["accession_column"])[compound_cols].mean()
    df.columns = [clean_compound_name(c) for c in df.columns]

    lab = pd.read_csv(label_path, header=cfg["header_rows"], index_col=[0, 1])
    classes_full = pd.Series(lab.index.get_level_values(1),
                             index=lab.index.get_level_values(0))

    common = df.index.intersection(classes_full.index)
    df = df.loc[common]
    classes = classes_full.loc[common]

    main_classes = list(cfg["clusters"])
    mask = classes.isin(main_classes)

    sample_counts = {cls: (classes[mask] == cls).sum() for cls in main_classes}
    n_missing_cells = int(df[mask].isna().sum().sum())

    print(f"   {mask.sum()} samples, {len(df.columns)} compounds")
    print(f"   {n_missing_cells} accession-compound cells with NO detected "
          f"replicate at all (left as missing, NOT imputed)")
    for cls in main_classes:
        print(f"      {cls}: {sample_counts[cls]} samples")

    return df[mask], classes[mask], sample_counts


def match_compounds_exact(data_df, odt_database):
    print("Exact compound matching...")
    matched = {}

    data_compounds = {clean_compound_name(col): col for col in data_df.columns}
    odt_compounds = {clean_compound_name(k): k for k in odt_database.keys()}

    for data_clean, data_original in data_compounds.items():
        if data_clean in odt_compounds:
            matched[data_original] = odt_compounds[data_clean]

    print(f"   Matched: {len(matched)}/{len(data_df.columns)} compounds")

    if len(matched) < len(data_df.columns):
        unmatched = set(data_df.columns) - set(matched.keys())
        print(f"   Unmatched: {len(unmatched)} compounds")
        if len(unmatched) <= 10:
            for comp in list(unmatched)[:10]:
                print(f"      {comp}")

    return matched


def calculate_odor_index(data_df, matched_compounds, odt_database, classes, sample_counts):
    print("Calculating relative odor index (missing stays missing: NOT filled with 0)...")

    compound_list = list(matched_compounds.keys())
    odor_index_df = pd.DataFrame(index=data_df.index, columns=compound_list, dtype=float)

    for compound, db_name in matched_compounds.items():
        odt_value = odt_database[db_name]['ODT']
        conc = pd.to_numeric(data_df[compound], errors='coerce').values
        odor_index_df[compound] = conc / odt_value

    odor_index_with_class = odor_index_df.copy()
    odor_index_with_class['Class'] = classes
    class_means = odor_index_with_class.groupby('Class').mean()
    class_std = odor_index_with_class.groupby('Class').std()

    detect_pct = odor_index_with_class.groupby('Class').apply(
        lambda g: g[compound_list].notna().mean() * 100)
    print("   Detection frequency by class (%):")
    for cls in list(cfg["clusters"]):
        low = detect_pct.loc[cls][detect_pct.loc[cls] < 80]
        if len(low):
            print(f"      {cls}: <80% detected for "
                  f"{', '.join(f'{c} ({v:.0f}%)' for c, v in low.items())}")

    total_samples = sum(sample_counts.values())
    overall_mean_weighted = pd.Series(0.0, index=compound_list)

    # overall index: cluster means weighted by cluster size; a cluster without detection adds 0
    for cls in list(cfg["clusters"]):
        weight = sample_counts[cls] / total_samples
        overall_mean_weighted += class_means.loc[cls].fillna(0) * weight

    print(f"   relative odor index calculated (class means computed from detected accessions only)")
    return odor_index_df, data_df[compound_list], class_means, class_std, overall_mean_weighted, detect_pct


def calculate_cv_methods(odor_index_df, classes, class_means):
    print("Calculating CV (nan-aware: missing values excluded, not zero)...")

    cv_results = {}

    for compound in odor_index_df.columns:
        class_values = [float(class_means.loc[cls, compound]) for cls in list(cfg["clusters"])]
        cv_between_class = np.nanstd(class_values) / (np.nanmean(class_values) + 1e-10)

        all_sample_values = odor_index_df[compound].values
        cv_within_sample = np.nanstd(all_sample_values) / (np.nanmean(all_sample_values) + 1e-10)

        cv_results[compound] = {
            'cv_between_class': cv_between_class,
            'cv_within_sample': cv_within_sample
        }

    print(f"   CV calculated")
    return cv_results


def perform_statistical_tests(odor_index_df, classes):
    print("Statistical tests...")

    stats_results = {}
    for compound in odor_index_df.columns:
        groups = [odor_index_df.loc[classes == cls, compound].astype(float).dropna()
                 for cls in list(cfg["clusters"])]

        if all(len(g) > 0 for g in groups):
            try:
                stat, p_value = stats.kruskal(*groups)
                if p_value < LEVELS[0]:
                    sig = '***'
                elif p_value < LEVELS[1]:
                    sig = '**'
                elif p_value < LEVELS[2]:
                    sig = '*'
                else:
                    sig = 'ns'
                stats_results[compound] = {'p_value': float(p_value), 'significance': sig}
            except:
                stats_results[compound] = {'p_value': 1.0, 'significance': 'ns'}
        else:
            stats_results[compound] = {'p_value': 1.0, 'significance': 'ns'}

    n_sig = sum(1 for v in stats_results.values() if v['significance'] != 'ns')
    print(f"   Significant: {n_sig}/{len(stats_results)}")
    return stats_results


def calculate_specificity(class_means):
    print("Calculating specificity...")

    specificity = {}

    for compound in class_means.columns:
        sum_all = class_means[compound].sum()
        if sum_all > 0:
            specificity[compound] = {
                cls: class_means.loc[cls, compound] / sum_all
                for cls in list(cfg["clusters"])
            }
        else:
            specificity[compound] = {cls: 0.33 for cls in list(cfg["clusters"])}

    print(f"   Specificity calculated")
    return specificity


def analyze_all_odorants(odor_index_df, conc_df, classes, odt_database, matched_compounds,
                        stats_results, cv_results, overall_mean_weighted, specificity):
    print("Analyzing odorants...")

    odor_index_with_class = odor_index_df.copy()
    odor_index_with_class['Class'] = classes
    class_means = odor_index_with_class.groupby('Class').mean()

    conc_with_class = conc_df.copy()
    conc_with_class['Class'] = classes
    conc_class_means = conc_with_class.groupby('Class').mean()

    all_odorants = {}

    for compound in odor_index_df.columns:
        db_name = matched_compounds.get(compound, compound)
        odor_info = odt_database.get(db_name, {'odor': 'Unknown', 'category': 'Unknown'})

        all_odorants[compound] = {
            'mean_oav': float(overall_mean_weighted[compound]),
            'max_oav': float(class_means[compound].max()),
            'mean_conc': float(conc_class_means[compound].mean()),
            'odor': odor_info['odor'],
            'category': odor_info['category'],
            'odt': odt_database[db_name]['ODT'] if db_name in odt_database else None,
            'p_value': stats_results[compound]['p_value'],
            'significance': stats_results[compound]['significance'],
            'cv_between_class': cv_results[compound]['cv_between_class'],
            'cv_within_sample': cv_results[compound]['cv_within_sample'],
            'class_means': {
                'ZERO': float(class_means.loc['ZERO', compound]),
                'ONE': float(class_means.loc['ONE', compound]),
                'TWO': float(class_means.loc['TWO', compound])
            },
            'specificity': specificity.get(compound, {})
        }

    all_odorants_sorted = dict(sorted(all_odorants.items(),
                                     key=lambda x: x[1]['mean_oav'],
                                     reverse=True))

    print(f"   {len(all_odorants)} odorants analyzed")
    return all_odorants_sorted, class_means


def generate_enhanced_profiles(all_odorants, class_means, specificity):
    print("Generating profiles...")

    profiles = {}

    for cls in list(cfg["clusters"]):
        class_ranking = sorted(all_odorants.items(),
                             key=lambda x: x[1]['class_means'][cls],
                             reverse=True)[:CONFIG['top_n_for_profiles']]

        weighted_categories = {}
        specific_compounds = []

        for rank, (compound, info) in enumerate(class_ranking, 1):
            cat = info['category']
            odor_index_cls = info['class_means'][cls]
            spec = info['specificity'].get(cls, 0.33)

            if spec > 0.5:
                specific_compounds.append((compound, cat, spec))

            weight = odor_index_cls * spec * (1.0 / np.sqrt(rank))
            weighted_categories[cat] = weighted_categories.get(cat, 0) + weight

        sorted_cats = sorted(weighted_categories.items(),
                           key=lambda x: x[1], reverse=True)

        top_cats = [cat for cat, _ in sorted_cats[:3]]

        if len(sorted_cats) >= 2:
            dominance_ratio = sorted_cats[0][1] / (sorted_cats[1][1] + 1e-10)
        else:
            dominance_ratio = 999

        if dominance_ratio > 3.0 and specific_compounds:
            primary_cat = top_cats[0]
            label_short = f"{primary_cat}-Dominant"
            label_detailed = f"{primary_cat}-Dominant ({specific_compounds[0][0][:15]})"
        elif dominance_ratio > 2.0:
            label_short = f"{top_cats[0]}-Rich"
            label_detailed = f"{top_cats[0]}-Rich / {top_cats[1]}"
        elif len(top_cats) >= 3:
            label_short = f"{top_cats[0]}/{top_cats[1]}/{top_cats[2]}"
            label_detailed = f"{top_cats[0]}-{top_cats[1]}-{top_cats[2]} Blend"
        else:
            label_short = f"{top_cats[0]}-{top_cats[1]}"
            label_detailed = f"{top_cats[0]}-{top_cats[1]} Complex"

        profiles[cls] = {
            'label_short': label_short,
            'label_detailed': label_detailed,
            'top_compounds': [comp for comp, _ in class_ranking],
            'weighted_categories': dict(sorted_cats),
            'dominance_ratio': dominance_ratio,
            'specific_compounds': specific_compounds[:3]
        }

        print(f"   {cls}: {label_short}")

    return profiles


def count_category_distribution(all_odorants):
    print("Counting categories...")

    categories = [info['category'] for info in all_odorants.values()]
    category_counts = Counter(categories)

    print(f"   Found {len(category_counts)} categories")
    return category_counts


def run_complete_pipeline(file_path):

    data_df, classes, sample_counts = load_data(file_path)
    matched = match_compounds_exact(data_df, ODT_DATABASE)

    if len(matched) == 0:
        print("\nNo compounds matched!")
        return None

    odor_index_df, conc_df, class_means, class_std, overall_mean_weighted, detect_pct = calculate_odor_index(
        data_df, matched, ODT_DATABASE, classes, sample_counts)

    cv_results = calculate_cv_methods(odor_index_df, classes, class_means)
    stats_results = perform_statistical_tests(odor_index_df, classes)
    specificity = calculate_specificity(class_means)

    all_odorants, class_means = analyze_all_odorants(
        odor_index_df, conc_df, classes, ODT_DATABASE, matched,
        stats_results, cv_results, overall_mean_weighted, specificity)

    profiles = generate_enhanced_profiles(all_odorants, class_means, specificity)

    profiles = {name: {"label_short": text} for name, text in cfg["cluster_descriptors"].items()}
    category_counts = count_category_distribution(all_odorants)

    write_table(odor_index_df, out_path('odor_index_accessions.xlsx'))
    write_table(class_means, out_path('odor_index_cluster_means.xlsx'))

    odorants_df = pd.DataFrame([
        {
            'Rank': i,
            'Compound': clean_compound_name(comp),
            'Mean_OAV': info['mean_oav'],
            'Category': info['category'],
            'Odor': info['odor'],
            'CV_between_class': info['cv_between_class'],
            'CV_within_sample': info['cv_within_sample'],
            'Significance': info['significance'],
            'ZERO_OAV': info['class_means']['ZERO'],
            'ONE_OAV': info['class_means']['ONE'],
            'TWO_OAV': info['class_means']['TWO']
        }
        for i, (comp, info) in enumerate(all_odorants.items(), 1)
    ])
    write_table(odorants_df, out_path('odor_index_compounds.xlsx'), index=False)

    print(f"\nTotal: {len(all_odorants)} compounds with a threshold")

    return {
        'oav': odor_index_df,
        'all_odorants': all_odorants,
        'profiles': profiles,
        'top_15': list(all_odorants.keys())[:15]
    }


FILE_PATH = input_path(cfg, step["input"])
require_sheets(FILE_PATH, [cfg["sheet_of"][SHEET]])
results = run_complete_pipeline(FILE_PATH)

_odor_index_df         = results['oav']
_all_odorants   = results['all_odorants']
_profiles       = results['profiles']
# the pipeline returns no test results; this stays empty and plot 2 computes its own tests
_stats_results  = results.get('stats_results', results.get('stats', {}))
_category_counts = Counter([info['category'] for info in results['all_odorants'].values()])
_class_means    = read_own_table(out_path('odor_index_cluster_means.xlsx'), index_col=0)


all_odorants    = _all_odorants
profiles        = _profiles
stats_results   = _stats_results if _stats_results else {}
category_counts = _category_counts
class_means_odor_index = _class_means
top_15_compounds = list(all_odorants.keys())[:15]
all_compounds    = list(all_odorants.keys())
colors_class    = CONFIG['colors_class']
categories      = sorted(set(info['category'] for info in all_odorants.values()))
cat_colors      = {cat: plt.cm.tab20(i/max(1,len(categories))) for i, cat in enumerate(categories)}
cat_marker_map  = {cat: 'o' for cat in categories}

data_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    for cls in list(cfg["clusters"]):
        try:
            data_matrix.loc[compound, cls] = class_means_odor_index.loc[cls, compound]
        except:
            data_matrix.loc[compound, cls] = 0.0
zscore_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    row = data_matrix.loc[compound].values.astype(float)
    mu, sigma = row.mean(), row.std()
    zscore_matrix.loc[compound, :] = (row - mu) / sigma if sigma > 1e-10 else 0.0
z_min, z_max = zscore_matrix.min().min(), zscore_matrix.max().max()
scaled_matrix = ((zscore_matrix - z_min) / (z_max - z_min + 1e-10)) * 100

fig1, ax1 = plt.subplots(figsize=(22, 6.7))

plotted_cats = set()
texts = []

for compound in all_compounds:
    info     = all_odorants[compound]
    cat      = info['category']
    marker   = cat_marker_map[cat]
    color    = cat_colors[cat]
    mean_pac = info['mean_conc']
    mean_odor_index = info['mean_oav']
    if mean_pac <= 0:
        continue
    log_pac  = np.log10(mean_pac)
    odor_index_plot = mean_odor_index + 1

    label = cat if cat not in plotted_cats else ""
    ax1.scatter(log_pac, odor_index_plot, s=120, c=[color], marker=marker,
                alpha=0.9, edgecolors='black', linewidth=2, label=label, zorder=3)
    plotted_cats.add(cat)

    t = ax1.text(log_pac, odor_index_plot, clean_compound_name(compound),
                 fontsize=12, fontweight='bold', zorder=5)
    texts.append(t)

ax1.set_yscale('log')
ax1.set_ylim(bottom=0.5)

try:
    from adjustText import adjust_text
    adjust_text(texts, ax=ax1,
                arrowprops=dict(arrowstyle='-', color='#888888', lw=0.6),
                expand_text=(1.15, 1.3), expand_points=(1.3, 1.5),
                force_text=(0.4, 0.6), force_points=(0.3, 0.5),
                only_move={'points':'y', 'text':'xy'})
except ImportError:
    print("   adjustText not installed: run: pip install adjustText")

ax1.set_xlabel('Log₁₀ Estimated Concentration of VOCs (µg/g)', fontsize=16, weight='bold')
ax1.set_ylabel(OI('OAV'), fontsize=16, weight='bold')
ax1.set_title('Odor Impact vs Estimated Concentration of VOCs',
              fontsize=17, weight='bold', pad=12)
ax1.tick_params(labelsize=14)
ax1.legend(loc='center left', bbox_to_anchor=(1.01, 0.5),
           fontsize=12, title='Odor Category', title_fontsize=13,
           frameon=True, framealpha=0.9)
ax1.grid(False)

plt.tight_layout()
save_figure(out_path('plot1_concentration_vs_odor_index.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("   Plot 1 saved.")


all_odorants    = _all_odorants
profiles        = _profiles
stats_results   = _stats_results if _stats_results else {}
category_counts = _category_counts
class_means_odor_index = _class_means
top_15_compounds = list(all_odorants.keys())[:15]
all_compounds    = list(all_odorants.keys())
colors_class    = CONFIG['colors_class']
categories      = sorted(set(info['category'] for info in all_odorants.values()))
cat_colors      = {cat: plt.cm.tab20(i/max(1,len(categories))) for i, cat in enumerate(categories)}
cat_markers     = ['o','s','^','D','v','<','>','p','*','h']
cat_marker_map  = {cat: cat_markers[i % len(cat_markers)] for i, cat in enumerate(categories)}

data_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    for cls in list(cfg["clusters"]):
        try:
            data_matrix.loc[compound, cls] = class_means_odor_index.loc[cls, compound]
        except:
            data_matrix.loc[compound, cls] = 0.0
zscore_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    row = data_matrix.loc[compound].values.astype(float)
    mu, sigma = row.mean(), row.std()
    zscore_matrix.loc[compound, :] = (row - mu) / sigma if sigma > 1e-10 else 0.0
z_min, z_max = zscore_matrix.min().min(), zscore_matrix.max().max()
scaled_matrix = ((zscore_matrix - z_min) / (z_max - z_min + 1e-10)) * 100

_lab = pd.read_csv(LABEL_FILE,
                   header=cfg["header_rows"], index_col=[0,1])
_classes_full = pd.Series(_lab.index.get_level_values(1),
                          index=_lab.index.get_level_values(0))
_common = _odor_index_df.index.intersection(_classes_full.index)
_classes_arr = _classes_full.loc[_common]
_classes_arr = _classes_arr[_classes_arr.isin(list(cfg["clusters"]))]
stats_results = perform_statistical_tests(_odor_index_df.loc[_classes_arr.index], _classes_arr.values)
if LABELS["fig_s5_significance"] == "bh":
    stats_results = bh_significance(stats_results)

top_all    = all_compounds[:33]
odor_indexs       = [all_odorants[c]['mean_oav'] for c in top_all]
sigs       = [stats_results.get(c, {}).get('significance', 'ns') for c in top_all]
sig_colors = {'***':'#d62728','**':'#ff7f0e','*':'#2ca02c','ns':'#aaaaaa'}
bar_colors = [sig_colors[s] for s in sigs]

n_sig = sum(1 for s in sigs if s != 'ns')
print(f"   Significant: {n_sig}/{len(sigs)}")

fig2, ax2 = plt.subplots(figsize=(11, 12))
ax2.barh(range(len(odor_indexs)), odor_indexs, color=bar_colors,
         edgecolor='black', linewidth=0.8, height=0.68)
ax2.set_yticks(range(len(odor_indexs)))
ax2.set_yticklabels([clean_compound_name(c) for c in top_all],
                    fontsize=11, fontweight='bold')
ax2.set_xlabel(OI('Mean OAV'), fontsize=13, weight='bold')
ax2.invert_yaxis()
ax2.tick_params(axis='x', labelsize=11)
ax2.set_xscale('log')
ax2.grid(True, alpha=0.25, axis='x')
for i, (val, sig) in enumerate(zip(odor_indexs, sigs)):
    ax2.text(val * 1.1, i, sig, va='center', fontsize=11, weight='bold',
             color=sig_colors[sig])
leg = [MPatch(color=sig_colors[s], label=s) for s in ['***','**','*','ns']]
ax2.legend(handles=leg, title='Significance', fontsize=10,
           title_fontsize=10, loc='lower right')
plt.tight_layout()
save_figure(out_path('plot2_odor_index_ranking.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("   Plot 2 saved.")


all_odorants    = _all_odorants
profiles        = _profiles
stats_results   = _stats_results if _stats_results else {}
category_counts = _category_counts
class_means_odor_index = _class_means
top_15_compounds = list(all_odorants.keys())[:15]
all_compounds    = list(all_odorants.keys())
colors_class    = CONFIG['colors_class']
categories      = sorted(set(info['category'] for info in all_odorants.values()))
cat_colors      = {cat: plt.cm.tab20(i/max(1,len(categories))) for i, cat in enumerate(categories)}
cat_markers     = ['o','s','^','D','v','<','>','p','*','h']
cat_marker_map  = {cat: cat_markers[i % len(cat_markers)] for i, cat in enumerate(categories)}

data_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    for cls in list(cfg["clusters"]):
        try:
            data_matrix.loc[compound, cls] = class_means_odor_index.loc[cls, compound]
        except:
            data_matrix.loc[compound, cls] = 0.0
zscore_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    row = data_matrix.loc[compound].values.astype(float)
    mu, sigma = row.mean(), row.std()
    zscore_matrix.loc[compound, :] = (row - mu) / sigma if sigma > 1e-10 else 0.0
z_min, z_max = zscore_matrix.min().min(), zscore_matrix.max().max()
scaled_matrix = ((zscore_matrix - z_min) / (z_max - z_min + 1e-10)) * 100


cats_sorted = sorted(category_counts.items(), key=lambda x: x[1], reverse=True)
cat_names = [c[0] for c in cats_sorted]
cat_vals  = [c[1] for c in cats_sorted]
colors_cat = [cat_colors[c] for c in cat_names]

fig3, ax3 = plt.subplots(figsize=(11, 8))
ax3.barh(range(len(cat_names)), cat_vals, color=colors_cat,
         edgecolor='black', linewidth=1.0, height=0.7)
ax3.set_yticks(range(len(cat_names)))
ax3.set_yticklabels(cat_names, fontsize=13, fontweight='bold')
ax3.set_xlabel('Number of Compounds', fontsize=15, weight='bold')
ax3.set_title('Odor Category Distribution', fontsize=18, weight='bold')
ax3.tick_params(axis='x', labelsize=13)
ax3.invert_yaxis()
ax3.grid(True, alpha=0.3, axis='x')
for i, val in enumerate(cat_vals):
    ax3.text(val + 0.08, i, str(val), va='center', fontsize=13, weight='bold')
plt.tight_layout()
save_figure(out_path('plot3_category_distribution.png'), panel=True, dpi=DPI_HIGH, bbox_inches='tight', facecolor='white')
plt.close(); print("   Plot 3 saved.")


all_odorants    = _all_odorants
profiles        = _profiles
stats_results   = _stats_results if _stats_results else {}
category_counts = _category_counts
class_means_odor_index = _class_means
top_15_compounds = list(all_odorants.keys())[:15]
all_compounds    = list(all_odorants.keys())
colors_class    = CONFIG['colors_class']
categories      = sorted(set(info['category'] for info in all_odorants.values()))
cat_colors      = {cat: plt.cm.tab20(i/max(1,len(categories))) for i, cat in enumerate(categories)}
cat_markers     = ['o','s','^','D','v','<','>','p','*','h']
cat_marker_map  = {cat: cat_markers[i % len(cat_markers)] for i, cat in enumerate(categories)}

data_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    for cls in list(cfg["clusters"]):
        try:
            data_matrix.loc[compound, cls] = class_means_odor_index.loc[cls, compound]
        except:
            data_matrix.loc[compound, cls] = 0.0
zscore_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    row = data_matrix.loc[compound].values.astype(float)
    mu, sigma = row.mean(), row.std()
    zscore_matrix.loc[compound, :] = (row - mu) / sigma if sigma > 1e-10 else 0.0
z_min, z_max = zscore_matrix.min().min(), zscore_matrix.max().max()
scaled_matrix = ((zscore_matrix - z_min) / (z_max - z_min + 1e-10)) * 100


n_vars  = len(top_15_compounds)
angles  = np.linspace(0, 2*np.pi, n_vars, endpoint=False).tolist()
angles += angles[:1]

fig4, ax4 = plt.subplots(figsize=(13, 13), subplot_kw=dict(projection='polar'))

for cls in list(cfg["clusters"]):
    values  = [float(scaled_matrix.loc[c, cls]) for c in top_15_compounds]
    values += values[:1]
    lbl = f"{cls}: {profiles[cls]['label_short']}"
    ax4.plot(angles, values, 'o-', linewidth=3,
             label=lbl, color=colors_class[cls], markersize=9)
    ax4.fill(angles, values, alpha=0.2, color=colors_class[cls])

ax4.set_xticks(angles[:-1])
label_names = [clean_compound_name(c) for i, c in enumerate(top_15_compounds)]
ax4.set_xticklabels(label_names, fontsize=14, fontweight='bold')
ax4.tick_params(axis='x', pad=35)
ax4.set_ylim(0, 100)
ax4.set_yticks([0, 25, 50, 75, 100])
ax4.set_yticklabels(['0','25','50','75','100'], fontsize=12)
ax4.set_title('Chemotype Aroma Profiles\n(Scaled Relative Intensity 0–100)',
              fontsize=18, weight='bold', pad=40)
ax4.legend(loc='upper right', bbox_to_anchor=(1.35, 1.15),
           fontsize=13, title='Chemotype', title_fontsize=13)
ax4.grid(True, alpha=0.4)
plt.tight_layout()
save_figure(out_path('plot4_spider.png'), panel=True, dpi=DPI_HIGH, bbox_inches='tight', facecolor='white')
plt.close(); print("   Plot 4 saved.")


all_odorants    = _all_odorants
profiles        = _profiles
stats_results   = _stats_results if _stats_results else {}
category_counts = _category_counts
class_means_odor_index = _class_means
top_15_compounds = list(all_odorants.keys())[:15]
all_compounds    = list(all_odorants.keys())
colors_class    = CONFIG['colors_class']
categories      = sorted(set(info['category'] for info in all_odorants.values()))
cat_colors      = {cat: plt.cm.tab20(i/max(1,len(categories))) for i, cat in enumerate(categories)}
cat_markers     = ['o','s','^','D','v','<','>','p','*','h']
cat_marker_map  = {cat: cat_markers[i % len(cat_markers)] for i, cat in enumerate(categories)}

data_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    for cls in list(cfg["clusters"]):
        try:
            data_matrix.loc[compound, cls] = class_means_odor_index.loc[cls, compound]
        except:
            data_matrix.loc[compound, cls] = 0.0
zscore_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    row = data_matrix.loc[compound].values.astype(float)
    mu, sigma = row.mean(), row.std()
    zscore_matrix.loc[compound, :] = (row - mu) / sigma if sigma > 1e-10 else 0.0
z_min, z_max = zscore_matrix.min().min(), zscore_matrix.max().max()
scaled_matrix = ((zscore_matrix - z_min) / (z_max - z_min + 1e-10)) * 100


z_abs_max = max(abs(z_min), abs(z_max))
y_labels  = [clean_compound_name(c) for i, c in enumerate(top_15_compounds)]

fig5, ax5 = plt.subplots(figsize=(8, 11))
sns.heatmap(zscore_matrix.astype(float),
            annot=True, fmt='.2f', cmap='RdYlGn',
            ax=ax5, cbar_kws={'label': 'Z-Score', 'shrink': 0.8},
            center=0, vmin=-z_abs_max, vmax=z_abs_max,
            yticklabels=y_labels,
            xticklabels=list(cfg["clusters"]),
            annot_kws={'size': 13, 'weight': 'bold'})
ax5.set_title('Heatmap – Relative Aroma Intensity (Z-Score)',
              fontsize=16, weight='bold', pad=15)
ax5.set_xlabel('Chemotype', fontsize=14, weight='bold')
ax5.set_ylabel('Rank. Compound', fontsize=14, weight='bold')
ax5.tick_params(axis='y', labelsize=12)
ax5.tick_params(axis='x', labelsize=13)
ax5.figure.axes[-1].tick_params(labelsize=11)
ax5.figure.axes[-1].set_ylabel('Z-Score', fontsize=12)
plt.tight_layout()
save_figure(out_path('plot5_heatmap.png'), panel=True, dpi=DPI_HIGH, bbox_inches='tight', facecolor='white')
plt.close(); print("   Plot 5 saved.")


all_odorants     = _all_odorants
class_means_odor_index  = _class_means
top_15_compounds = list(all_odorants.keys())[:15]
colors_class     = CONFIG['colors_class']

raw_matrix = pd.DataFrame(index=top_15_compounds, columns=list(cfg["clusters"]), dtype=float)
for compound in top_15_compounds:
    for cls in list(cfg["clusters"]):
        try:
            raw_matrix.loc[compound, cls] = class_means_odor_index.loc[cls, compound]
        except:
            raw_matrix.loc[compound, cls] = 0.0

print("Relative odor index per cluster (cluster means)")
print(f"{'Compound':<28} {'Category':<12} {'ZERO':>12} {'ONE':>12} {'TWO':>12}")
for compound in top_15_compounds:
    cat = all_odorants[compound]['category']
    z, o, t = raw_matrix.loc[compound, list(cfg["clusters"])]
    print(f"{clean_compound_name(compound)[:27]:<28} {cat:<12} {z:>12.1f} {o:>12.1f} {t:>12.1f}")

print("\nTop 5 compounds per cluster by relative odor index")
for cls in list(cfg["clusters"]):
    ranked = raw_matrix[cls].sort_values(ascending=False).head(5)
    print(f"\n{cls}:")
    for i, (cpd, val) in enumerate(ranked.items(), 1):
        cat  = all_odorants[cpd]['category']
        odor = all_odorants[cpd]['odor']
        print(f"   {i}. {clean_compound_name(cpd):<26} relative odor index={val:>10.1f}  [{cat}] {odor}")

print("\nOdor category contribution per cluster (sum of relative odor index)")
for cls in list(cfg["clusters"]):
    cat_sums = {}
    for compound in top_15_compounds:
        cat = all_odorants[compound]['category']
        cat_sums[cat] = cat_sums.get(cat, 0) + raw_matrix.loc[compound, cls]
    total = sum(cat_sums.values())
    ranked = sorted(cat_sums.items(), key=lambda x: -x[1])
    print(f"\n{cls}:")
    for cat, val in ranked[:5]:
        pct = val / total * 100 if total > 0 else 0
        print(f"   {cat:<14} relative odor index sum={val:>12.1f}  ({pct:>5.1f}%)")

n_vars  = len(top_15_compounds)
angles  = np.linspace(0, 2*np.pi, n_vars, endpoint=False).tolist()
angles += angles[:1]

fig, ax = plt.subplots(figsize=(13, 13), subplot_kw=dict(projection='polar'))

log_matrix = np.log10(raw_matrix + 1)

for cls in list(cfg["clusters"]):
    values  = [float(log_matrix.loc[c, cls]) for c in top_15_compounds]
    values += values[:1]
    ax.plot(angles, values, 'o-', linewidth=3, label=cls,
            color=colors_class[cls], markersize=9)
    ax.fill(angles, values, alpha=0.18, color=colors_class[cls])

ax.set_xticks(angles[:-1])
ax.set_xticklabels([clean_compound_name(c) for c in top_15_compounds],
                   fontsize=14, fontweight='bold')
ax.tick_params(axis='x', pad=35)
ax.set_ylabel('')
ax.set_title(OI('Chemotype Aroma Profiles — Absolute OAV\n(log₁₀ scale)'),
             fontsize=17, weight='bold', pad=45)
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.12),
          fontsize=13, title='Chemotype', title_fontsize=13)
ax.grid(True, alpha=0.4)
plt.tight_layout()
save_figure(out_path('plot_radar_odor_index.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("\n   Raw relative odor index spider saved: plot_radar_odor_index.png")

fig, ax = plt.subplots(figsize=(8, 11))
sns.heatmap(log_matrix.astype(float), annot=raw_matrix.astype(float),
            fmt='.0f', cmap='YlOrRd', ax=ax,
            cbar_kws={'label': OI('log₁₀(OAV + 1)'), 'shrink': 0.8},
            yticklabels=[clean_compound_name(c) for c in top_15_compounds],
            xticklabels=list(cfg["clusters"]),
            annot_kws={'size': 11, 'weight': 'bold'})
ax.set_title(OI('Absolute OAV per Chemotype'), fontsize=16, weight='bold', pad=15)
ax.set_xlabel('Chemotype', fontsize=14, weight='bold')
ax.set_ylabel('Compound', fontsize=14, weight='bold')
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=13)
plt.tight_layout()
save_figure(out_path('plot_heatmap_odor_index.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("   Raw relative odor index heatmap saved: plot_heatmap_odor_index.png")

write_table(raw_matrix, out_path('odor_index_per_cluster.xlsx'))
print("   Exported: odor_index_per_cluster.xlsx")


odor_index_df       = _odor_index_df
all_odorants = _all_odorants

_na_per_cpd = odor_index_df.isna().sum()
_complete_compounds = _na_per_cpd[_na_per_cpd == 0].index.tolist()
print(f"Complete-case compounds (0% missing): {len(_complete_compounds)}/{odor_index_df.shape[1]}")
print(f"  Excluded (had missing accessions): "
      f"{sorted(set(odor_index_df.columns) - set(_complete_compounds))}")

odor_index_complete = odor_index_df[_complete_compounds].dropna()

odor_index_log = np.log1p(odor_index_complete)
X = StandardScaler().fit_transform(odor_index_log)

print("INDEPENDENT AROMA-BASED CLUSTERING")
print(f"Input: {X.shape[0]} accessions x {X.shape[1]} odorants (log relative odor index, scaled, complete-case only)")

print("\nSilhouette scores:")
sil_scores = {}
for k in range(2, 8):
    km = KMeans(n_clusters=k, random_state=SEED, n_init=20)
    lbl = km.fit_predict(X)
    s = silhouette_score(X, lbl)
    sil_scores[k] = s
    print(f"   K={k}: {s:.3f}")

best_k = max(sil_scores, key=sil_scores.get)
print(f"\n   Optimal K = {best_k} (silhouette = {sil_scores[best_k]:.3f})")

km_final = KMeans(n_clusters=best_k, random_state=SEED, n_init=50)
aroma_labels = km_final.fit_predict(X)

aroma_clusters = pd.Series(aroma_labels, index=odor_index_complete.index, name='AromaCluster')
print(f"\nCluster sizes:")
for c in sorted(set(aroma_labels)):
    print(f"   Aroma-{c}: {(aroma_labels==c).sum()} accessions")

_lab2 = pd.read_csv(LABEL_FILE,
                    header=cfg["header_rows"], index_col=[0,1])
_classes_full2 = pd.Series(_lab2.index.get_level_values(1),
                           index=_lab2.index.get_level_values(0))

_classes_aligned = _classes_full2.reindex(odor_index_complete.index)
_valid = _classes_aligned.isin(list(cfg["clusters"]))

_classes_arr = _classes_aligned[_valid]
_aroma_aligned = pd.Series(aroma_labels, index=odor_index_complete.index)[_valid]

crosstab = pd.crosstab(pd.Series(_classes_arr.values, name='ML_Chemotype'),
                       pd.Series(_aroma_aligned.values, name='Aroma_Cluster'))
print("ML CLUSTER vs AROMA CLUSTER (cross-tabulation)")
print(crosstab.to_string())

ari = adjusted_rand_score(_classes_arr.values, _aroma_aligned.values)
print(f"\nAdjusted Rand Index: {ari:.3f}")

odor_index_with_ac = odor_index_complete.copy()
odor_index_with_ac['AromaCluster'] = aroma_labels
ac_means = odor_index_with_ac.groupby('AromaCluster').mean()

print("TOP 8 ODORANTS PER AROMA CLUSTER (raw relative odor index)")
for c in sorted(set(aroma_labels)):
    ranked = ac_means.loc[c].sort_values(ascending=False).head(8)
    n = (aroma_labels==c).sum()
    print(f"\nAroma-{c} (n={n}):")
    for i, (cpd, val) in enumerate(ranked.items(), 1):
        cat  = all_odorants[cpd]['category']
        odor = all_odorants[cpd]['odor']
        print(f"   {i}. {clean_compound_name(cpd):<26} relative odor index={val:>10.1f}  [{cat}] {odor}")

print("AROMA CATEGORY CONTRIBUTION PER CLUSTER (% of total relative odor index)")
for c in sorted(set(aroma_labels)):
    cat_sums = {}
    for cpd in ac_means.columns:
        cat = all_odorants[cpd]['category']
        cat_sums[cat] = cat_sums.get(cat, 0) + ac_means.loc[c, cpd]
    total = sum(cat_sums.values())
    ranked = sorted(cat_sums.items(), key=lambda x: -x[1])
    print(f"\nAroma-{c}:")
    for cat, val in ranked[:6]:
        pct = val/total*100 if total>0 else 0
        print(f"   {cat:<14} {pct:>5.1f}%   (relative odor index sum = {val:>10.1f})")

pca = PCA(n_components=2)
pcs = pca.fit_transform(X)

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

palette_ac = plt.cm.Set2(np.linspace(0, 1, best_k))
for c in sorted(set(aroma_labels)):
    m = aroma_labels == c
    axes[0].scatter(pcs[m,0], pcs[m,1], s=90, alpha=0.85,
                    color=palette_ac[c], edgecolors='black', linewidth=0.8,
                    label=f'Aroma-{c} (n={m.sum()})')
axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=13, weight='bold')
axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=13, weight='bold')
axes[0].set_title('Aroma-based Clusters (independent)', fontsize=15, weight='bold')
axes[0].legend(fontsize=11)
axes[0].grid(alpha=0.2)

for cls in list(cfg["clusters"]):
    m = np.array(_classes_arr) == cls
    axes[1].scatter(pcs[m,0], pcs[m,1], s=90, alpha=0.85,
                    color=CONFIG['colors_class'][cls], edgecolors='black', linewidth=0.8,
                    label=f'{cls} (n={m.sum()})')
axes[1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=13, weight='bold')
axes[1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=13, weight='bold')
axes[1].set_title('ML Chemotypes (for comparison)', fontsize=15, weight='bold')
axes[1].legend(fontsize=11)
axes[1].grid(alpha=0.2)

plt.tight_layout()
save_figure(out_path('aroma_clustering_comparison.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("\n   Saved: aroma_clustering_comparison.png")

out = pd.DataFrame({
    T('Genotype'): odor_index_df.index,
    'ML_Chemotype': list(_classes_arr),
    'Aroma_Cluster': aroma_labels
})
write_table(out, out_path('aroma_clusters.xlsx'), index=False)
write_table(ac_means, out_path('aroma_cluster_odor_index_means.xlsx'))
print("   Saved: aroma_clusters.xlsx, aroma_cluster_odor_index_means.xlsx")


odor_index_df       = _odor_index_df
all_odorants = _all_odorants
colors_class = CONFIG['colors_class']

_lab_x = pd.read_csv(LABEL_FILE, header=cfg["header_rows"], index_col=[0,1])
_classes_full_x = pd.Series(_lab_x.index.get_level_values(1),
                            index=_lab_x.index.get_level_values(0))
classes_s = _classes_full_x.loc[_classes_full_x.index.intersection(odor_index_df.index)]
classes_s = classes_s[classes_s.isin(list(cfg["clusters"]))]
classes_s.name = 'Class'

odor_index_log = np.log1p(odor_index_df)
z_samples = (odor_index_log - odor_index_log.mean()) / odor_index_log.std()
z_samples = z_samples.clip(-step["clip_limit"], step["clip_limit"])

print(f"Z-score per accession ({z_samples.shape[0]} accessions x {z_samples.shape[1]} compounds)")
_n_nan = int(z_samples.isna().sum().sum())
print(f"   Range: [{np.nanmin(z_samples.values):.2f}, {np.nanmax(z_samples.values):.2f}]  "
      f"({_n_nan} NaN cells: real non-detects, not imputed)")

cpd_order = sorted(odor_index_df.columns,
                   key=lambda c: all_odorants[c]['mean_oav'], reverse=True)
top_n     = 15
cpd_plot  = cpd_order[:top_n]

n_vars  = len(cpd_plot)
angles  = np.linspace(0, 2*np.pi, n_vars, endpoint=False).tolist()
angles += angles[:1]

fig, ax = plt.subplots(figsize=(13, 13), subplot_kw=dict(projection='polar'))

for cls in list(cfg["clusters"]):
    m   = classes_s == cls
    sub = z_samples.loc[m, cpd_plot]

    med = [float(sub[c].median()) for c in cpd_plot]
    q25 = [float(sub[c].quantile(0.25)) for c in cpd_plot]
    q75 = [float(sub[c].quantile(0.75)) for c in cpd_plot]
    med += med[:1]; q25 += q25[:1]; q75 += q75[:1]

    ax.fill_between(angles, q25, q75, color=colors_class[cls],
                    alpha=0.18, zorder=2, linewidth=0)
    ax.plot(angles, med, 'o-', linewidth=3.5, markersize=8,
            color=colors_class[cls], zorder=5,
            label=f"{cls} (n={m.sum()})")

ax.set_xticks(angles[:-1])
ax.set_xticklabels([clean_compound_name(c) for c in cpd_plot],
                   fontsize=14, fontweight='bold')
ax.tick_params(axis='x', pad=35)
ax.set_ylim(-3, 3)
ax.set_yticks([-3,-2,-1,0,1,2,3])
ax.set_yticklabels(['-3','-2','-1','0','1','2','3'], fontsize=11)
ax.set_title(OI(TT('Odour-Activity Profiles of Individual Accessions\n'
                   f'(Z-score across all {len(odor_index_df)} genotypes, clipped to ±3)')),
             fontsize=16, weight='bold', pad=45)
ax.legend(loc='upper right', bbox_to_anchor=(1.28, 1.12),
          fontsize=13, title='Chemotype', title_fontsize=13)
ax.grid(True, alpha=0.4)
plt.tight_layout()
save_figure(out_path('plot_spider_SAMPLE_zscore.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("\n   Saved: plot_spider_SAMPLE_zscore.png")

z_hm = z_samples[cpd_order]
row_order = classes_s.sort_values().index
z_hm = z_hm.loc[row_order]
row_colors = classes_s.loc[row_order].map(colors_class)

g = sns.clustermap(z_hm,
                   row_cluster=False, col_cluster=False,
                   row_colors=row_colors,
                   cmap='coolwarm', center=0, vmin=-3, vmax=3,
                   linewidths=0,
                   cbar_kws={'label': 'Z-score'},
                   figsize=(16, 20),
                   dendrogram_ratio=(0.10, 0.08),
                   colors_ratio=0.02,
                   yticklabels=False,
                   xticklabels=[clean_compound_name(c) for c in cpd_order])

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                             rotation=45, ha='right', fontsize=13, fontweight='bold')
g.ax_heatmap.set_xlabel(OI('Odour-active compound'), fontsize=16, weight='bold')
g.ax_heatmap.set_ylabel(f'Accession (n = {len(odor_index_df)})', fontsize=16, weight='bold')
g.ax_cbar.set_ylabel('Z-score', fontsize=14, weight='bold')
g.ax_cbar.tick_params(labelsize=12)

leg = [MPatch(facecolor=colors_class[c], edgecolor='black', label=f'{c} (n={(classes_s==c).sum()})')
       for c in list(cfg["clusters"])]
place_cluster_legend(g, leg)

save_figure(out_path('plot_heatmap_SAMPLE_zscore.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("   Saved: plot_heatmap_SAMPLE_zscore.png")

print("\nCluster median z-score")
print(f"{'Compound':<28} {'ZERO':>9} {'ONE':>9} {'TWO':>9}")
med = pd.DataFrame({cls: z_samples[classes_s==cls].median() for cls in list(cfg["clusters"])})
for c in cpd_order:
    print(f"{clean_compound_name(c)[:27]:<28} {med.loc[c,'ZERO']:>9.2f} "
          f"{med.loc[c,'ONE']:>9.2f} {med.loc[c,'TWO']:>9.2f}")

write_table(z_samples, out_path('odor_index_zscore_per_accession.xlsx'))
write_table(med, out_path('odor_index_zscore_cluster_medians.xlsx'))
print("\n   Saved: odor_index_zscore_per_accession.xlsx, odor_index_zscore_cluster_medians.xlsx")


odor_index_df       = _odor_index_df
all_odorants = _all_odorants
colors_class = CONFIG['colors_class']

_lab_x = pd.read_csv(LABEL_FILE, header=cfg["header_rows"], index_col=[0,1])
_classes_full_x = pd.Series(_lab_x.index.get_level_values(1),
                            index=_lab_x.index.get_level_values(0))
classes_s = _classes_full_x.loc[_classes_full_x.index.intersection(odor_index_df.index)]
classes_s = classes_s[classes_s.isin(list(cfg["clusters"]))]
classes_s.name = 'Class'


def pareto_scale(df):
    """(x - mean) / sqrt(SD): Pareto scaling"""
    mu    = df.mean(axis=0)
    sd    = df.std(axis=0)
    return (df - mu) / np.sqrt(sd.replace(0, np.nan))

odor_index_log   = np.log1p(odor_index_df)
par_full  = pareto_scale(odor_index_log)
par_full  = par_full.fillna(0)

print(f"Pareto scaling, (x - mean) / sqrt(SD), across {len(odor_index_df)} accessions")
print(f"   Range: [{par_full.values.min():.2f}, {par_full.values.max():.2f}]")

cpd_order = sorted(odor_index_df.columns,
                   key=lambda c: all_odorants[c]['mean_oav'], reverse=True)
cpd_plot  = cpd_order[:15]

med_par = pd.DataFrame({cls: par_full[classes_s == cls].median()
                        for cls in list(cfg["clusters"])})

print("CLASS MEDIAN: PARETO-SCALED (top 15 by relative odor index)")
print(f"{'Compound':<28} {'Category':<12} {'ZERO':>9} {'ONE':>9} {'TWO':>9}")
for c in cpd_plot:
    cat = all_odorants[c]['category']
    print(f"{clean_compound_name(c)[:27]:<28} {cat:<12} "
          f"{med_par.loc[c,'ZERO']:>9.2f} {med_par.loc[c,'ONE']:>9.2f} {med_par.loc[c,'TWO']:>9.2f}")

n_vars  = len(cpd_plot)
angles  = np.linspace(0, 2*np.pi, n_vars, endpoint=False).tolist()
angles += angles[:1]

fig, ax = plt.subplots(figsize=(13, 13), subplot_kw=dict(projection='polar'))

for cls in list(cfg["clusters"]):
    m   = classes_s == cls
    sub = par_full.loc[m, cpd_plot]
    med = [float(sub[c].median()) for c in cpd_plot]
    q25 = [float(sub[c].quantile(0.25)) for c in cpd_plot]
    q75 = [float(sub[c].quantile(0.75)) for c in cpd_plot]
    med += med[:1]; q25 += q25[:1]; q75 += q75[:1]
    ax.fill_between(angles, q25, q75, color=colors_class[cls],
                    alpha=0.16, zorder=2, linewidth=0)
    ax.plot(angles, med, 'o-', linewidth=3.5, markersize=8,
            color=colors_class[cls], zorder=5, label=f"{cls} (n={m.sum()})")

ax.set_xticks(angles[:-1])
ax.set_xticklabels([clean_compound_name(c) for c in cpd_plot],
                   fontsize=14, fontweight='bold')
ax.tick_params(axis='x', pad=35)
ax.set_title(OI('Odour-Activity Profiles by Chemotype\n'
                '(Pareto-scaled log OAV; median ± IQR)'),
             fontsize=16, weight='bold', pad=45)
ax.legend(loc='upper right', bbox_to_anchor=(1.28, 1.12),
          fontsize=13, title='Chemotype', title_fontsize=13)
ax.grid(True, alpha=0.4)
plt.tight_layout()
save_figure(out_path('plot_spider_PARETO.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("\n   Saved: plot_spider_PARETO.png")


par_hm     = par_full[cpd_order]
row_colors = classes_s.map(colors_class)
lim        = np.nanpercentile(np.abs(par_hm.values), 99)

row_link = scipy_linkage(pdist(par_hm.values, metric='euclidean'), method='ward')
col_link = scipy_linkage(pdist(par_hm.values.T, metric='euclidean'), method='ward')

g = sns.clustermap(par_hm,
                   row_cluster=True, col_cluster=True,
                   row_linkage=row_link, col_linkage=col_link,
                   row_colors=row_colors,
                   cmap='coolwarm', center=0, vmin=-lim, vmax=lim,
                   linewidths=0,
                   cbar_kws={'label': OI('Pareto-scaled log OAV')},
                   figsize=(16, 20),
                   dendrogram_ratio=(0.10, 0.08),
                   colors_ratio=0.02,
                   yticklabels=False,
                   xticklabels=[clean_compound_name(c) for c in cpd_order])

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                             rotation=45, ha='right', fontsize=13, fontweight='bold')
g.ax_heatmap.set_xlabel(OI('Odour-active compound'), fontsize=16, weight='bold')
g.ax_heatmap.set_ylabel(f'Accession (n = {len(odor_index_df)})', fontsize=16, weight='bold')
g.ax_cbar.set_ylabel('Pareto-scaled', fontsize=13, weight='bold')
g.ax_cbar.tick_params(labelsize=11)

leg = [MPatch(facecolor=colors_class[c], edgecolor='black',
              label=f'{c} (n={(classes_s==c).sum()})') for c in list(cfg["clusters"])]
place_cluster_legend(g, leg)

save_figure(out_path('plot_heatmap_PARETO.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("   Saved: plot_heatmap_PARETO.png")


z_full = pd.DataFrame(StandardScaler().fit_transform(odor_index_log),
                      index=odor_index_df.index, columns=odor_index_df.columns).clip(-step["clip_limit"], step["clip_limit"])

scalings = {
    OI('No scaling\n(log OAV)'): odor_index_log[cpd_plot],
    'Pareto\n(x-μ)/√SD':          par_full[cpd_plot],
    'Autoscale / Z-score\n(x-μ)/SD': z_full[cpd_plot],
}

fig, axes = plt.subplots(1, 3, figsize=(24, 9), subplot_kw=dict(projection='polar'))
for ax, (name, mat) in zip(axes, scalings.items()):
    for cls in list(cfg["clusters"]):
        m   = classes_s == cls
        med = [float(mat.loc[m, c].median()) for c in cpd_plot]
        med += med[:1]
        ax.plot(angles, med, 'o-', linewidth=2.5, markersize=6,
                color=colors_class[cls], label=cls)
        ax.fill(angles, med, alpha=0.12, color=colors_class[cls])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([clean_compound_name(c)[:14] for c in cpd_plot], fontsize=8)
    ax.tick_params(axis='x', pad=12)
    ax.set_title(name, fontsize=14, weight='bold', pad=28)
    ax.grid(True, alpha=0.35)
axes[-1].legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), fontsize=11)
plt.suptitle('Effect of scaling method on apparent chemotype separation',
             fontsize=17, weight='bold', y=1.02)
plt.tight_layout()
save_figure(out_path('plot_scaling_comparison.png'), panel=True, dpi=DPI, facecolor='white', bbox_inches='tight')
plt.show()
print("   Saved: plot_scaling_comparison.png")

with pd.ExcelWriter(out_path('odor_index_scalings.xlsx')) as w:
    write_table(odor_index_log, w, sheet_name='log_odor_index')
    write_table(par_full, w, sheet_name='Pareto')
    write_table(z_full, w, sheet_name='Zscore')
    write_table(med_par, w, sheet_name='Pareto_class_medians')
print("   Saved: odor_index_scalings.xlsx")


odor_index_df       = _odor_index_df
all_odorants = _all_odorants
K_ORDER      = list(cfg["clusters"])

_lab_x = pd.read_csv(LABEL_FILE, header=cfg["header_rows"], index_col=[0,1])
_classes_full_x = pd.Series(_lab_x.index.get_level_values(1),
                            index=_lab_x.index.get_level_values(0))
classes_s = _classes_full_x.loc[_classes_full_x.index.intersection(odor_index_df.index)]
classes_s = classes_s[classes_s.isin(K_ORDER)]
classes_s.name = 'Class'

X = np.log1p(odor_index_df)
class_means = odor_index_df.copy()
class_means['Class'] = classes_s.values
class_means = class_means.groupby('Class').mean()


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = (a[:, None] >  b[None, :]).sum()
    lt = (a[:, None] <  b[None, :]).sum()
    return (gt - lt) / (len(a) * len(b))


def bh_fdr(pvals):
    pv = np.asarray(pvals); n = len(pv)
    order = np.argsort(pv); q = np.empty(n); prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        val = min(prev, pv[idx] * n / i); q[idx] = val; prev = val
    return q

rows = []
for c in odor_index_df.columns:
    for k in K_ORDER:
        a = X.loc[classes_s == k, c].dropna().values
        b = X.loc[classes_s != k, c].dropna().values
        if len(a) < 3 or len(b) < 3:
            rows.append({
                'Compound': clean_compound_name(c), 'Class': k,
                'Category': all_odorants[c]['category'],
                'p': np.nan, 'delta': np.nan, 'n_a': len(a), 'n_b': len(b),
            })
            continue
        _, pv = stats.mannwhitneyu(a, b, alternative='two-sided')
        rows.append({
            'Compound': clean_compound_name(c),
            'Class': k,
            'Category': all_odorants[c]['category'],
            'Odor': all_odorants[c]['odor'],
            'Mean_OAV': float(class_means.loc[k, c]),
            'p': pv,
            'delta': cliffs_delta(a, b),
        })
res = pd.DataFrame(rows)
res['q'] = bh_fdr(res['p'].values)
res['Effect'] = pd.cut(res['delta'].abs(), [-.01, .147, .33, .474, 1.01],
                       labels=['negligible','small','medium','large'])
res['Direction'] = np.where(res['delta'] > 0, 'elevated', 'reduced')
res = res.sort_values(['Class','delta'], ascending=[True, False])

print("ONE-vs-REST  |  Mann-Whitney U on log relative odor index, BH-corrected  |  Cliff's delta")
print("thresholds: |d| 0.147 small, 0.33 medium, 0.474 large")
for k in K_ORDER:
    sub = res[res.Class == k]
    up = sub[(sub.q < 0.05) & (sub.delta >  0.33)]
    dn = sub[(sub.q < 0.05) & (sub.delta < -0.33)].sort_values('delta')
    print(f"\n### {k}  (n={int((classes_s==k).sum())}) : ELEVATED vs other two:")
    if len(up) == 0:
        print("    (none)")
    for _, r in up.iterrows():
        print(f"    {r.Compound:<20} d={r.delta:+.2f} ({r.Effect})  q={r.q:.1e}  "
              f"relative odor index={r.Mean_OAV:>10.0f}  [{r.Category}] {r.Odor}")
    print(f"### {k} : REDUCED vs other two:")
    if len(dn) == 0:
        print("    (none)")
    for _, r in dn.iterrows():
        print(f"    {r.Compound:<20} d={r.delta:+.2f} ({r.Effect})  q={r.q:.1e}  "
              f"relative odor index={r.Mean_OAV:>10.0f}  [{r.Category}] {r.Odor}")

print("ODOR-CATEGORY COMPOSITION PER CLUSTER (% of summed relative odor index)")
cat_of = {c: all_odorants[c]['category'] for c in odor_index_df.columns}
share = {}
for k in K_ORDER:
    s = class_means.loc[k]
    by_cat = s.groupby(pd.Series(cat_of)).sum()
    share[k] = by_cat / by_cat.sum() * 100
share = pd.DataFrame(share).sort_values('ZERO', ascending=False)
print(share.round(2).to_string())

print("MODULATION INDEX (non-sulfur fraction of summed relative odor index)")
for k in K_ORDER:
    tot  = class_means.loc[k].sum()
    sulf = class_means.loc[k][[c for c in odor_index_df.columns if cat_of[c]=='Sulfur']].sum()
    n    = int((classes_s == k).sum())
    print(f"  {k:<5} n={n:<4} sulfur {sulf/tot*100:5.1f}%  |  non-sulfur {100-sulf/tot*100:5.1f}%"
          f"  |  total relative odor index {tot:>10,.0f}")

with pd.ExcelWriter(out_path('odor_index_one_vs_rest.xlsx')) as w:
    write_table(res, w, sheet_name='one_vs_rest', index=False)
    write_table(res[(res.q < 0.05) & (res.delta.abs() > 0.33)], w, sheet_name='significant', index=False)
    write_table(share, w, sheet_name='category_share')
print("\n   Saved: odor_index_one_vs_rest.xlsx")


detect = _odor_index_df.copy()
detect['Class'] = classes_s.values
out = detect.groupby('Class')[_odor_index_df.columns].apply(lambda g: g.notna().mean()*100)
write_table(out.T.round(1), out_path('detection_per_cluster.xlsx'))
print(out.T.round(1).to_string())

finish_step(cfg, TREATMENT)

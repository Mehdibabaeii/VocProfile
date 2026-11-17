#!/usr/bin/env python3
#
# VocProfile - 17_plot_odor_radar.py
#
# Radar plots of the relative odor index for the 15 highest-ranked compounds under
# four scalings, drawn with cluster values only and with one thin line per accession
# (Fig. 4).
#
# in   REPLICATE_FILE, MAIN_DATASET, BLUP_CLUSTER_FILE (cluster of each accession)
#      ODOR_THRESHOLD_TABLE, CLUSTER_TABLE, SEED
# out  figures/17_plot_odor_radar/radar/radar_<scaling>_<mode>.png
#      results/17_plot_odor_radar/radar/radar_data_all_versions.xlsx
#      results/17_plot_odor_radar/17_qc.txt
#
# Usage:
#   python scripts/17_plot_odor_radar.py <config.txt>
#
# Notes: docs/NOTES.md, step 17
#
import os
import warnings
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, save_figure, write_table, check_same_accessions, use_font, read_own_table

warnings.filterwarnings('ignore')

cfg, step, args = load_config("17")
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
OUT_DIR = output_dir(cfg, "17_plot_odor_radar")


def out_path(name):
    return os.path.join(OUT_DIR, name)


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
    print("Calculating CV...")

    cv_results = {}

    for compound in odor_index_df.columns:
        class_values = [float(class_means.loc[cls, compound]) for cls in list(cfg["clusters"])]
        cv_between_class = np.std(class_values) / (np.mean(class_values) + 1e-10)

        all_sample_values = odor_index_df[compound].values
        cv_within_sample = np.std(all_sample_values) / (np.mean(all_sample_values) + 1e-10)

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
_stats_results  = results.get('stats_results', results.get('stats', {}))
_category_counts = Counter([info['category'] for info in results['all_odorants'].values()])
_class_means    = read_own_table(out_path('odor_index_cluster_means.xlsx'), index_col=0)


odor_index_df          = _odor_index_df
all_odorants    = _all_odorants
profiles        = _profiles
class_means_odor_index = _class_means
colors_class    = CONFIG['colors_class']
K_ORDER         = list(cfg["clusters"])

_lab_r = pd.read_csv(LABEL_FILE,
                     header=cfg["header_rows"], index_col=[0,1])
_classes_full_r = pd.Series(_lab_r.index.get_level_values(1),
                            index=_lab_r.index.get_level_values(0))
classes_s = _classes_full_r.reindex(odor_index_df.index)
classes_s = classes_s[classes_s.isin(K_ORDER)]
classes_s.name = 'Class'

CPD = list(all_odorants.keys())[:15]

OUTDIR = out_path(step["output_folder"])
os.makedirs(OUTDIR, exist_ok=True)

FIGSIZE   = (13, 13)
FS_LABEL  = 24
FS_TICK   = 18
PAD_LABEL = 43
LW_MAIN   = 5
MS_MAIN   = 11
ALPHA_FILL= 0.20
DPI       = step["dpi"]

label_names = [clean_compound_name(c) for c in CPD]
n_vars = len(CPD)
angles = np.linspace(0, 2*np.pi, n_vars, endpoint=False).tolist()
angles += angles[:1]

print(f"{len(odor_index_df)} accessions, {len(CPD)} compounds")
print(f"   Class sizes: " + ", ".join(f"{c}={int((classes_s==c).sum())}" for c in K_ORDER))


log_odor_index = np.log1p(odor_index_df[CPD])


def build_v1():
    dm = pd.DataFrame(index=CPD, columns=K_ORDER, dtype=float)
    for c in CPD:
        for k in K_ORDER:
            try:    dm.loc[c, k] = class_means_odor_index.loc[k, c]
            except: dm.loc[c, k] = np.nan
    z = pd.DataFrame(index=CPD, columns=K_ORDER, dtype=float)
    for c in CPD:
        row = dm.loc[c].values.astype(float)
        mu, sd = row.mean(), row.std()
        z.loc[c, :] = (row - mu) / sd if sd > 1e-10 else 0.0
    zmin, zmax = z.min().min(), z.max().max()
    cls_mat = ((z - zmin) / (zmax - zmin + 1e-10)) * 100
    zs = (log_odor_index - log_odor_index.mean()) / log_odor_index.std()
    smin, smax = np.nanmin(zs.values), np.nanmax(zs.values)
    smp_mat = ((zs - smin) / (smax - smin + 1e-10)) * 100
    return smp_mat, cls_mat.T[CPD], (0, 100), ['0','25','50','75','100'], [0,25,50,75,100]


def build_v2():
    smp = np.log10(odor_index_df[CPD] + 1)
    cls = pd.DataFrame({k: smp[classes_s == k].median() for k in K_ORDER}).T[CPD]
    hi  = float(np.ceil(np.nanmax(smp.values)))
    ticks = list(range(0, int(hi)+1))
    return smp, cls, (0, hi), [str(t) for t in ticks], ticks


def build_v3():
    smp = ((log_odor_index - log_odor_index.mean()) / log_odor_index.std()).clip(-step["clip_limit"], step["clip_limit"])
    cls = pd.DataFrame({k: smp[classes_s == k].median() for k in K_ORDER}).T[CPD]
    return smp, cls, (-3, 3), ['-3','-2','-1','0','1','2','3'], [-3,-2,-1,0,1,2,3]


def build_v4():
    mu, sd = log_odor_index.mean(axis=0), log_odor_index.std(axis=0)
    smp = ((log_odor_index - mu) / np.sqrt(sd.replace(0, np.nan))).fillna(0)
    cls = pd.DataFrame({k: smp[classes_s == k].median() for k in K_ORDER}).T[CPD]
    lim = float(np.ceil(np.nanpercentile(np.abs(smp.values), 99) * 2) / 2)
    ticks = [-lim, -lim/2, 0, lim/2, lim]
    return smp, cls, (-lim, lim), [f'{t:g}' for t in ticks], ticks

VERSIONS = {
    'v1_original_zscore3_0to100': ('Z-score across 3 class means, rescaled 0–100', build_v1),
    'v2_no_scaling_log10':        ('No scaling — log10(OAV+1)',                     build_v2),
    'v3_zscore_accessions_clip3': (f'Z-score across {len(log_odor_index)} accessions, clipped ±3',     build_v3),
    'v4_pareto_accessions':       (f'Pareto (x−μ)/√SD across {len(log_odor_index)} accessions',        build_v4),
}

BUILT = {}
for tag, (desc, fn) in VERSIONS.items():
    smp, cls, ylim, ytl, yt = fn()
    BUILT[tag] = dict(desc=desc, smp=smp, cls=cls, ylim=ylim, ytl=ytl, yt=yt)
    print(f"{tag:<30} range=[{np.nanmin(smp.values):.2f}, {np.nanmax(smp.values):.2f}]  "
          f"({int(np.isnan(smp.values).sum())} NaN cells: real non-detects, not imputed)")


def _style_axis(ax, ylim, yticks, yticklabels):
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(label_names, fontsize=FS_LABEL, fontweight='bold')
    ax.tick_params(axis='x', pad=PAD_LABEL)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=FS_TICK)
    ax.grid(True, alpha=0.4)


def radar_class(tag, show_legend=False):
    """Mode A: one line per cluster (class means/medians only)."""
    B = BUILT[tag]
    fig, ax = plt.subplots(figsize=FIGSIZE, subplot_kw=dict(projection='polar'))
    for k in K_ORDER:
        vals  = [float(B['cls'].loc[k, c]) for c in CPD]
        vals += vals[:1]
        lbl = f"{k}: {profiles[k]['label_short']}" if k in profiles else k
        ax.plot(angles, vals, 'o-', linewidth=LW_MAIN, markersize=MS_MAIN,
                color=colors_class[k], label=lbl, zorder=5)
        ax.fill(angles, vals, alpha=ALPHA_FILL, color=colors_class[k], zorder=3)
    _style_axis(ax, B['ylim'], B['yt'], B['ytl'])
    if show_legend:
        ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.15),
                  fontsize=13, title='Chemotype', title_fontsize=13)
    plt.tight_layout()
    p = f"{OUTDIR}/radar_{tag}_CLASS.png"
    save_figure(p, panel=True, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   {p}")


def radar_samples(tag, show_legend=False, line_alpha=0.3, line_w=1):
    """Mode B: every accession as a thin line coloured by cluster,"""
    B = BUILT[tag]
    fig, ax = plt.subplots(figsize=FIGSIZE, subplot_kw=dict(projection='polar'))
    for geno in B['smp'].index:
        k    = classes_s.loc[geno]
        vals = [float(B['smp'].loc[geno, c]) for c in CPD]
        vals += vals[:1]
        ax.plot(angles, vals, '-', linewidth=line_w, alpha=line_alpha,
                color=colors_class[k], zorder=2)
    for k in K_ORDER:
        m     = classes_s == k
        vals  = [float(B['smp'].loc[m, c].mean()) for c in CPD]
        vals += vals[:1]
        lbl = f"{k}: {profiles[k]['label_short']}" if k in profiles else k
        ax.plot(angles, vals, 'o-', linewidth=LW_MAIN + 1, markersize=MS_MAIN,
                color=colors_class[k], label=f"{lbl}  (n={int(m.sum())})", zorder=6,
                path_effects=None)
    _style_axis(ax, B['ylim'], B['yt'], B['ytl'])
    if show_legend:
        ax.legend(loc='upper right', bbox_to_anchor=(1.40, 1.15),
                  fontsize=13, title='Chemotype', title_fontsize=13)
    plt.tight_layout()
    p = f"{OUTDIR}/radar_{tag}_SAMPLES.png"
    save_figure(p, panel=True, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   {p}")


SHOW_LEGEND = step["show_legend"]

for tag in VERSIONS:
    print(f"\n{BUILT[tag]['desc']}")
    radar_class(tag,   show_legend=SHOW_LEGEND)
    radar_samples(tag, show_legend=SHOW_LEGEND)

print(f"\nAll radars written to {OUTDIR}")


xls = f"{OUTDIR}/radar_data_all_versions.xlsx"
with pd.ExcelWriter(xls, engine='openpyxl') as w:

    ref = pd.DataFrame({
        'Compound': [clean_compound_name(c) for c in CPD],
        'Category': [all_odorants[c]['category'] for c in CPD],
        'Odor':     [all_odorants[c]['odor'] for c in CPD],
        'ODT':      [all_odorants[c]['odt'] for c in CPD],
        'Mean_OAV': [all_odorants[c]['mean_oav'] for c in CPD],
        'ZERO_OAV': [class_means_odor_index.loc['ZERO', c] for c in CPD],
        'ONE_OAV':  [class_means_odor_index.loc['ONE',  c] for c in CPD],
        'TWO_OAV':  [class_means_odor_index.loc['TWO',  c] for c in CPD],
    })
    write_table(ref, w, sheet_name='00_compound_reference', index=False)

    write_table(pd.DataFrame({
        'Version':     list(VERSIONS.keys()),
        'Description': [BUILT[t]['desc'] for t in VERSIONS],
        'Axis_min':    [BUILT[t]['ylim'][0] for t in VERSIONS],
        'Axis_max':    [BUILT[t]['ylim'][1] for t in VERSIONS],
    }), w, sheet_name='00_versions', index=False)

    for tag in VERSIONS:
        B = BUILT[tag]
        short = tag.split('_')[0]

        cls_out = B['cls'].copy()
        cls_out.columns = [clean_compound_name(c) for c in cls_out.columns]
        write_table(cls_out, w, sheet_name=f'{short}_class')

        smp_out = B['smp'].copy()
        smp_out.columns = [clean_compound_name(c) for c in smp_out.columns]
        smp_out.insert(0, 'Class', classes_s.values)
        write_table(smp_out, w, sheet_name=f'{short}_samples')

        rows = []
        for k in K_ORDER:
            m = classes_s == k
            for c in CPD:
                v = B['smp'].loc[m, c]
                rows.append({'Class': k, 'Compound': clean_compound_name(c),
                             'n': int(m.sum()),
                             'Mean': v.mean(), 'Median': v.median(), 'SD': v.std(),
                             'Q25': v.quantile(0.25), 'Q75': v.quantile(0.75),
                             'Min': v.min(), 'Max': v.max()})
        write_table(pd.DataFrame(rows), w, sheet_name=f'{short}_stats', index=False)

print(f"Excel written: {xls}")

finish_step(cfg, TREATMENT)

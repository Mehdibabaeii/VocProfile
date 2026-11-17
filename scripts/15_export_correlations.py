#!/usr/bin/env python3
#
# VocProfile - 15_export_correlations.py
#
# Numbers behind the correlation matrix (Fig. 3b) and the accession heat map: Pearson
# correlations between compounds and clipped z-scores, for the accession means as they
# are, filled, and filled then transformed as ln(1 + x) (the version of the paper).
#
# in   MEANS_CLUSTER_FILE, MAIN_DATASET, CLIP_LIMIT, DISPLAY_NAMES_TABLE
# out  results/15_export_correlations/correlation_exports/<version>/CORRELATION_ANALYSIS.xlsx
#      results/15_export_correlations/correlation_exports/<version>/HEATMAP_ANALYSIS.xlsx
#      results/15_export_correlations/15_qc.txt
#
# Usage:
#   python scripts/15_export_correlations.py <config.txt>
#
# Notes: docs/NOTES.md, step 15
#
import gc
import os
import warnings

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list
from scipy.cluster.hierarchy import linkage as scipy_linkage
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

from lib_config import load_config, input_path, require_sheets, output_dir, display_names
from lib_run import start_step, finish_step, write_table

warnings.filterwarnings('ignore')

cfg, step, args = load_config("15")
start_step(cfg, 'means_clusters')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    '1_raw: values without detection are left out of every statistic; nothing is filled (pairwise in the correlations; set to zero before the z-scores)',
    '2_lod2 and 3_lod2_log: accession means without detection: set to NONDETECT_FRACTION of the smallest positive accession mean of the compound',
    '3_lod2_log: values transformed as ln(1 + x)',
)
FILE_PATH = input_path(cfg, step["input"])
SHEET_NAME = step["sheet"]
require_sheets(FILE_PATH, [cfg["sheet_of"][SHEET_NAME]])
# names shown in the exported tables; the data keep their own column names
SHOWN = display_names(cfg) if step["use_display_names"] else {}
OUTPUT_BASE = os.path.join(output_dir(cfg, "15_export_correlations"), step["output_folder"])


print(f"Loading {SHEET_NAME}...")
df_raw = pd.read_excel(FILE_PATH, sheet_name=cfg["sheet_of"][SHEET_NAME],
                        header=cfg["header_rows"], index_col=[0, 1], engine='openpyxl')
df_raw.index.names   = ['Acc', 'K']
df_raw.columns.names = ['Main_Class', 'Subclass', 'Compound']
df_raw = df_raw.dropna(axis=0, how='all')
print(f"  Shape: {df_raw.shape}")

k_series = df_raw.index.get_level_values('K').astype(str)
df_base  = df_raw.reset_index(level='K', drop=True)
df_base  = df_base.apply(pd.to_numeric, errors='coerce')


df_v1 = df_base.copy()

df_v2 = df_base.copy()
col_mins = df_v2.min(skipna=True)
fill_vals = (col_mins * cfg["settings"].value("NONDETECT_FRACTION")).where(col_mins > 0, 1e-9)
df_v2 = df_v2.fillna(fill_vals)

df_v3 = np.log1p(df_v2)

versions = [
    ("1_raw",        df_v1, "Raw (no imputation)"),
    ("2_lod2",       df_v2, "LOD/2 imputed"),
    ("3_lod2_log",   df_v3, "LOD/2 imputed + log1p"),
]


def with_k(df):
    out = df.copy()
    out.columns = out.columns.get_level_values('Compound')
    out['K'] = k_series.values
    return out


def export_correlation(df_plot, df_multiindex, out_folder, label):
    os.makedirs(out_folder, exist_ok=True)

    df_corr_calc = df_plot.drop(columns='K', errors='ignore')
    corr_matrix  = df_corr_calc.corr(method='pearson')

    cols_to_use = [col for col in df_multiindex.columns if col[2] in corr_matrix.columns]
    meta_df = pd.DataFrame(cols_to_use, columns=['Main_Class', 'Subclass', 'Compound'])

    def get_meta(compound, field):
        row = meta_df[meta_df['Compound'] == compound]
        return row[field].values[0] if len(row) > 0 else 'Unknown'

    pairs = []
    cols = corr_matrix.columns.tolist()
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            c1, c2 = cols[i], cols[j]
            r = corr_matrix.iloc[i, j]
            mc1 = get_meta(c1, 'Main_Class')
            mc2 = get_meta(c2, 'Main_Class')
            sg1 = get_meta(c1, 'Subclass')
            sg2 = get_meta(c2, 'Subclass')
            pairs.append({
                'Compound_1': c1, 'Compound_2': c2,
                'Pearson_r': r, 'Abs_r': abs(r), 'r_squared': r**2,
                'Main_Class_1': mc1, 'Main_Class_2': mc2,
                'Subclass_1': sg1,   'Subclass_2': sg2,
                'Same_Main_Class': 'Yes' if mc1 == mc2 else 'No',
                'Same_Subclass':   'Yes' if sg1 == sg2 else 'No',
            })

    df_pairs = pd.DataFrame(pairs).sort_values('Abs_r', ascending=False)
    within   = df_pairs[df_pairs['Same_Main_Class'] == 'Yes']
    cross    = df_pairs[df_pairs['Same_Main_Class'] == 'No']

    summary = pd.DataFrame({
        'Category':   ['All Pairs', 'Within Same Class', 'Cross Class'],
        'Count':      [len(df_pairs), len(within), len(cross)],
        'Mean_r':     [df_pairs['Pearson_r'].mean(), within['Pearson_r'].mean(), cross['Pearson_r'].mean()],
        'Mean_Abs_r': [df_pairs['Abs_r'].mean(),     within['Abs_r'].mean(),     cross['Abs_r'].mean()],
        'Median_r':   [df_pairs['Pearson_r'].median(), within['Pearson_r'].median(), cross['Pearson_r'].median()],
        'Max_r':      [df_pairs['Pearson_r'].max(),   within['Pearson_r'].max(),   cross['Pearson_r'].max()],
        'Min_r':      [df_pairs['Pearson_r'].min(),   within['Pearson_r'].min(),   cross['Pearson_r'].min()],
    })

    class_summary = []
    for mc in meta_df['Main_Class'].dropna().unique():
        cp = df_pairs[(df_pairs['Main_Class_1'] == mc) & (df_pairs['Main_Class_2'] == mc)]
        if len(cp) > 0:
            class_summary.append({
                'Main_Class': mc, 'N_Pairs': len(cp),
                'Mean_r': cp['Pearson_r'].mean(), 'Median_r': cp['Pearson_r'].median(),
                'Max_r':  cp['Pearson_r'].max(),  'Min_r':    cp['Pearson_r'].min(),
            })

    excel_path = os.path.join(out_folder, 'CORRELATION_ANALYSIS.xlsx')
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df_pairs = df_pairs.replace({'Compound_1': SHOWN, 'Compound_2': SHOWN})
        write_table(corr_matrix.rename(index=SHOWN, columns=SHOWN), writer, sheet_name='Correlation_Matrix')
        write_table(df_pairs, writer, sheet_name='All_Pairs', index=False)
        write_table(df_pairs.head(50), writer, sheet_name='Top_50_Correlations', index=False)
        write_table(summary, writer, sheet_name='Summary_Statistics', index=False)
        write_table(pd.DataFrame(class_summary), writer, sheet_name='By_Main_Class', index=False)

    print(f"  Correlation saved: {excel_path}")
    return excel_path


def export_heatmap(df_plot, out_folder, label):
    os.makedirs(out_folder, exist_ok=True)

    compound_traits = [c for c in df_plot.columns if c != 'K']
    temp_df = df_plot[compound_traits + ['K']].copy()

    X = StandardScaler().fit_transform(temp_df[compound_traits].fillna(0))
    scaled_df = pd.DataFrame(X, columns=compound_traits, index=temp_df.index)
    scaled_df = scaled_df.clip(-step["clip_limit"], step["clip_limit"])
    scaled_df['Cluster'] = temp_df['K'].values

    row_link = scipy_linkage(pdist(scaled_df[compound_traits].values, metric='euclidean'), method='ward')
    col_link = scipy_linkage(pdist(scaled_df[compound_traits].values.T, metric='euclidean'), method='ward')

    row_order = leaves_list(row_link)
    col_order = leaves_list(col_link)

    scaled_ordered = scaled_df[compound_traits].iloc[row_order, :].iloc[:, col_order]
    scaled_ordered['Cluster'] = scaled_df['Cluster'].iloc[row_order].values

    cluster_means = scaled_df.groupby('Cluster')[compound_traits].mean()
    cluster_std   = scaled_df.groupby('Cluster')[compound_traits].std()

    top_compounds = []
    for cluster in sorted(scaled_df['Cluster'].unique()):
        cd = scaled_df[scaled_df['Cluster'] == cluster]
        means = cd[compound_traits].mean().sort_values(ascending=False)
        for rank, (cpd, zscore) in enumerate(means.head(10).items(), 1):
            top_compounds.append({'Cluster': cluster, 'Rank': rank, 'Compound': cpd,
                                   'Mean_Z': zscore, 'Std_Z': cd[cpd].std()})

    cluster_summary = []
    for cluster in sorted(scaled_df['Cluster'].unique()):
        cd = scaled_df[scaled_df['Cluster'] == cluster]
        cluster_summary.append({
            'Cluster': cluster, 'N_Samples': len(cd),
            'Percent': f"{len(cd)/len(scaled_df)*100:.1f}%",
            'Mean_Z_All': cd[compound_traits].mean().mean(),
            'Std_Z_All':  cd[compound_traits].std().mean(),
            'Max_Z': cd[compound_traits].max().max(),
            'Min_Z': cd[compound_traits].min().min(),
        })

    excel_path = os.path.join(out_folder, 'HEATMAP_ANALYSIS.xlsx')
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        write_table(scaled_ordered.rename(columns=SHOWN), writer, sheet_name='Scaled_Data_Clustered')
        write_table(cluster_means.rename(columns=SHOWN), writer, sheet_name='Cluster_Means')
        write_table(cluster_std.rename(columns=SHOWN), writer, sheet_name='Cluster_StdDev')
        write_table(pd.DataFrame(top_compounds).replace({'Compound': SHOWN}), writer, sheet_name='Top_10_Per_Cluster', index=False)
        write_table(pd.DataFrame(cluster_summary), writer, sheet_name='Cluster_Summary', index=False)

    print(f"  Heatmap saved: {excel_path}")
    return excel_path


for tag, df_version, label in versions:
    print(f"\nVersion: {label}")
    out_folder = os.path.join(OUTPUT_BASE, tag)
    df_plot    = with_k(df_version)

    export_correlation(df_plot, df_raw, out_folder, label)
    export_heatmap(df_plot, out_folder, label)
    gc.collect()

print("\nAll exports complete.")

finish_step(cfg, TREATMENT)

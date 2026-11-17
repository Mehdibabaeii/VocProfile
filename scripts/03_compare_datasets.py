#!/usr/bin/env python3
#
# VocProfile - 03_compare_datasets.py
#
# Comparison of the quantification datasets from the replicate-level data: within-
# accession CV, between-accession CV, Spearman correlation of accession means between
# datasets, Shapiro-Wilk test, and the number of accession means without detection
# (Table S9).
#
# in   REPLICATE_FILE, DATASETS
# out  results/03_compare_datasets/dataset_comparison.xlsx
#      results/03_compare_datasets/03_qc.txt
#
# Usage:
#   python scripts/03_compare_datasets.py <config.txt>
#
# Notes: docs/NOTES.md, step 03
#
import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, shapiro

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, write_table

cfg, step, args = load_config("03")
start_step(cfg, 'replicates')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'values without detection are left out of every statistic; nothing is filled',
    'an accession mean is missing only when no replicate has a detection',
    'unequal numbers of replicates per accession (one to five) are used as they are; the within-accession CV needs at least two',
)
INPUT_FILE = input_path(cfg, step["input"])
SHEETS = step["sheets"]
SHORT = step["short_names"]
require_sheets(INPUT_FILE, [cfg["sheet_of"][name] for name in SHEETS])
OUTPUT_FILE = os.path.join(output_dir(cfg, "03_compare_datasets"), step["output_file"])


raw = {}
accession_means = {}

for sh in SHEETS:
    print(f"Loading: {sh}")
    df = pd.read_excel(INPUT_FILE, sheet_name=cfg["sheet_of"][sh], na_values=cfg["nondetect_code"], engine="openpyxl")
    compound_cols = [c for c in df.columns if c not in (cfg["accession_column"], cfg["replicate_column"])]
    df[compound_cols] = df[compound_cols].apply(pd.to_numeric, errors="coerce")
    raw[sh] = df

    means = df.groupby(cfg["accession_column"])[compound_cols].mean()
    accession_means[sh] = means
    print(f"  {df.shape[0]} replicate rows -> {means.shape[0]} accessions, "
          f"{len(compound_cols)} compounds")

compound_cols = [c for c in raw[SHEETS[0]].columns if c not in (cfg["accession_column"], cfg["replicate_column"])]
print(f"\nTotal compounds: {len(compound_cols)}")


missing_summary = {}
for sh in SHEETS:
    means = accession_means[sh]
    n_missing = int(means.isna().sum().sum())
    n_total = means.shape[0] * means.shape[1]
    missing_summary[sh] = n_missing
    print(f"{sh:<28} {n_missing:>5} / {n_total} missing "
          f"({n_missing/n_total*100:.1f}%)  [accession level]")


within_cv_summary = {}
per_compound_within = {}

for sh in SHEETS:
    df = raw[sh]
    cv_rows = []
    for cpd in compound_cols:
        g = df.groupby(cfg["accession_column"])[cpd].agg(["mean", "std", "count"])
        g = g[g["count"] >= 2]
        g["cv"] = g["std"] / g["mean"]
        cv_rows.append({
            "Compound": cpd,
            "Mean_within_CV": g["cv"].mean(skipna=True),
            "Median_within_CV": g["cv"].median(skipna=True),
            "N_accessions_used": len(g),
        })
    cv_df = pd.DataFrame(cv_rows)
    per_compound_within[sh] = cv_df
    within_cv_summary[sh] = cv_df["Median_within_CV"].median()

print("Within-accession CV (technical noise) computed.")


between_cv_summary = {}
per_compound_between = {}

for sh in SHEETS:
    means = accession_means[sh]
    rows = []
    for cpd in compound_cols:
        v = means[cpd].dropna()
        cv = v.std() / v.mean() if v.mean() != 0 else np.nan
        rows.append({"Compound": cpd, "Between_accession_CV": cv, "N": len(v)})
    bdf = pd.DataFrame(rows)
    per_compound_between[sh] = bdf
    between_cv_summary[sh] = bdf["Between_accession_CV"].median()

print("Between-accession CV (biological dispersion) computed.")


pairs = list(combinations(SHEETS, 2))
spearman_rows = []
for cpd in compound_cols:
    row = {"Compound": cpd}
    for a, b in pairs:
        va, vb = accession_means[a][cpd], accession_means[b][cpd]
        common = va.dropna().index.intersection(vb.dropna().index)
        if len(common) < 10:
            row[f"rho_{SHORT[a]}_vs_{SHORT[b]}"] = np.nan
            continue
        r, _ = spearmanr(va[common], vb[common])
        row[f"rho_{SHORT[a]}_vs_{SHORT[b]}"] = r
    spearman_rows.append(row)
spearman_df = pd.DataFrame(spearman_rows)

spearman_summary = []
for a, b in pairs:
    col = f"rho_{SHORT[a]}_vs_{SHORT[b]}"
    spearman_summary.append({
        "Pair": f"{SHORT[a]} vs {SHORT[b]}",
        "Mean_rho": spearman_df[col].mean(),
        "Median_rho": spearman_df[col].median(),
    })
spearman_summary_df = pd.DataFrame(spearman_summary)
print("Spearman correlation (accession means) computed.")


normality_rows = []
for sh in SHEETS:
    means = accession_means[sh]
    n_tested = n_nonnormal = 0
    for cpd in compound_cols:
        v = means[cpd].dropna()
        if len(v) < 8:
            continue
        n_tested += 1
        try:
            _, p = shapiro(v)
            if p < 0.05:
                n_nonnormal += 1
        except Exception:
            pass
    normality_rows.append({
        "Dataset": sh, "N_tested": n_tested, "N_nonnormal": n_nonnormal,
        "Pct_nonnormal": n_nonnormal / n_tested * 100 if n_tested else np.nan,
    })
normality_df = pd.DataFrame(normality_rows)
print("Shapiro-Wilk normality (accession means) computed.")


print("\nSummary")
header = f"{'Dataset':<28}{'within-CV':>14}{'between-CV':>16}{'missing (acc-level)':>22}"
print(header)
for sh in SHEETS:
    print(f"{sh:<28}{within_cv_summary[sh]:>14.3f}{between_cv_summary[sh]:>16.3f}"
          f"{missing_summary[sh]:>22}")

print("\nSpearman (accession means):")
print(spearman_summary_df.to_string(index=False))

print("\nShapiro-Wilk (accession means):")
print(normality_df.to_string(index=False))


with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    summary_df = pd.DataFrame({
        "Dataset": SHEETS,
        "Within_accession_CV_median": [within_cv_summary[s] for s in SHEETS],
        "Between_accession_CV_median": [between_cv_summary[s] for s in SHEETS],
        "Missing_accession_level": [missing_summary[s] for s in SHEETS],
    })
    write_table(summary_df, writer, sheet_name="00_summary", index=False)
    write_table(spearman_summary_df, writer, sheet_name="01_spearman_summary", index=False)
    write_table(normality_df, writer, sheet_name="02_shapiro_summary", index=False)
    write_table(spearman_df, writer, sheet_name="03_spearman_per_compound", index=False)
    for sh in SHEETS:
        write_table(per_compound_within[sh], writer, sheet_name=f"within_{SHORT[sh]}"[:31], index=False)
        write_table(per_compound_between[sh], writer, sheet_name=f"between_{SHORT[sh]}"[:31], index=False)

print(f"\nSaved: {OUTPUT_FILE}")

finish_step(cfg, TREATMENT)

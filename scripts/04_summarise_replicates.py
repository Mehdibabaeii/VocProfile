#!/usr/bin/env python3
#
# VocProfile - 04_summarise_replicates.py
#
# Descriptive statistics and coefficient of variation of each compound across all
# replicate runs, one sheet per dataset (Table S9).
#
# in   REPLICATE_FILE, DATASETS
# out  results/04_summarise_replicates/replicate_statistics.xlsx
#      results/04_summarise_replicates/04_qc.txt
#
# Usage:
#   python scripts/04_summarise_replicates.py <config.txt>
#
# Notes: docs/NOTES.md, step 04
#
import os

import numpy as np
import pandas as pd

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, write_table

cfg, step, args = load_config("04")
start_step(cfg, 'replicates')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'values without detection are left out of every statistic; nothing is filled',
    'unequal numbers of replicates per accession (one to five) are used as they are',
)
INPUT_FILE = input_path(cfg, step["input"])
sheets = step["sheets"]
require_sheets(INPUT_FILE, [cfg["sheet_of"][name] for name in sheets])
OUTPUT_FILE = os.path.join(output_dir(cfg, "04_summarise_replicates"), step["output_file"])

results = {}

for sh in sheets:
    print(f"Processing: {sh}")
    df = pd.read_excel(INPUT_FILE, sheet_name=cfg["sheet_of"][sh], na_values=cfg["nondetect_code"], engine="openpyxl")

    compound_cols = [c for c in df.columns if c not in (cfg["accession_column"], cfg["replicate_column"])]
    df[compound_cols] = df[compound_cols].apply(pd.to_numeric, errors="coerce")

    rows = []
    for cpd in compound_cols:
        data = df[cpd].dropna()
        n     = len(data)
        mean  = data.mean()
        std   = data.std()
        cv    = (std / mean) if mean != 0 else np.nan
        var   = data.var()
        med   = data.median()
        se    = data.sem()
        q1    = data.quantile(0.25)
        q3    = data.quantile(0.75)
        iqr   = q3 - q1
        mad   = (data - med).abs().median()
        rows.append({
            "Compound":    cpd,
            "Count":       n,
            "Missing":     len(df) - n,
            "Mean":        round(mean, 6),
            "StDev":       round(std,  6),
            "Mean ± SD":   f"{mean:.3f} ± {std:.3f}",
            "SE_Mean":     round(se,   6),
            "Median":      round(med,  6),
            "Variance":    round(var,  6),
            "CoefVar":     round(cv,   6),
            "Min":         round(data.min(), 6),
            "Q1":          round(q1,   6),
            "Q3":          round(q3,   6),
            "Max":         round(data.max(), 6),
            "Range":       round(data.max() - data.min(), 6),
            "IQR":         round(iqr,  6),
            "MAD":         round(mad,  6),
            "Skewness":    round(data.skew(), 6),
            "Kurtosis":    round(data.kurt(), 6),
        })

    result_df = pd.DataFrame(rows)
    clean_name = sh.replace("_(4)", "")
    results[clean_name] = result_df
    print(f"  Compounds: {len(rows)}, Total reps used: {len(df)}")

with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    for sheet_name, df_out in results.items():
        write_table(df_out, writer, sheet_name=sheet_name, index=False)

print(f"\nSaved: {OUTPUT_FILE}")

finish_step(cfg, TREATMENT)

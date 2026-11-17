#!/usr/bin/env python3
#
# VocProfile - 05_fill_nondetects.py
#
# Accession-mean workbook with the means without detection filled, one sheet per
# dataset. Kept for reference; no later step reads it.
#
# in   MEANS_FILE, DATASETS
# out  results/05_fill_nondetects/accession_means_filled.xlsx
#      results/05_fill_nondetects/05_qc.txt
#
# Usage:
#   python scripts/05_fill_nondetects.py <config.txt>
#
# Notes: docs/NOTES.md, step 05
#
import os

import pandas as pd

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, write_table, keep_data_names

cfg, step, args = load_config("05")
start_step(cfg, 'means')
keep_data_names()
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'accession means without detection: set to NONDETECT_FRACTION of the smallest positive accession mean of the compound',
)
file_path = input_path(cfg, step["input"])
sheets = step["sheets"]
require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets])
output_path = os.path.join(output_dir(cfg, "05_fill_nondetects"), step["output_file"])

filled_sheets = {}

for sheet in sheets:
    df = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet], header=cfg["header_rows"], index_col=[0,1], engine='openpyxl')
    df_num = df.apply(pd.to_numeric, errors='coerce')

    n_missing = int(df_num.isna().sum().sum())

    for col in df_num.columns:
        vals = df_num[col]
        if vals.isna().any():
            pos_min = vals[vals > 0].min()
            fill_val = pos_min * cfg["settings"].value("NONDETECT_FRACTION") if pd.notna(pos_min) else 1e-9
            df_num[col] = vals.fillna(fill_val)

    filled_sheets[sheet] = df_num
    print(f"{sheet}: {n_missing} missing values filled with LOD/2")

with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
    for sheet, df in filled_sheets.items():
        write_table(df, writer, sheet_name=sheet)

print(f"\nSaved: {output_path}")

finish_step(cfg, TREATMENT)

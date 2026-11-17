#!/usr/bin/env python3
#
# VocProfile - 08_build_cluster_inputs.py
#
# The two cluster-labelled files read by the later steps, built from the accession
# means and from the BLUP workbook of step 06: accession means with the cluster name
# in the second column, and BLUP estimates with cluster name and class header rows.
#
# in   MEANS_FILE, results of step 06 for MAIN_DATASET
#      CLUSTER_TABLE, CLUSTER_SIZES_EXPECTED
# out  results/08_build_cluster_inputs/accession_means_clusters.xlsx
#      results/08_build_cluster_inputs/blup_clusters.csv
#      results/08_build_cluster_inputs/08_qc.txt
#
# Usage:
#   python scripts/08_build_cluster_inputs.py <config.txt>
#
# Notes: docs/NOTES.md, step 08
#
import csv
import os
import sys

import pandas as pd

from lib_config import load_config, input_path, require_sheets, output_dir, die
from lib_run import start_step, finish_step, keep_data_names

cfg, step, args = load_config("08")
start_step(cfg, 'means')
keep_data_names()
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'clusters are named in the order of CLUSTER_TABLE, largest first',
    'rows of the BLUP file are in sorted accession order',
)
SHEET = cfg["main_sheet"]
HEADER = len(cfg["header_rows"])
MEANS_XLSX = input_path(cfg, step["input"])
require_sheets(MEANS_XLSX, [cfg["sheet_of"][SHEET]])
BLUP_XLSX = os.path.join(output_dir(cfg, "06_fit_blup_clusters"), step["blup_workbook"])
if not os.path.isfile(BLUP_XLSX):
    die(f"BLUP workbook not found: {BLUP_XLSX} (run step 06 first)")
BLUP_SHEET = step["blup_sheet"]
NAMES_BY_SIZE = cfg["clusters"]                   # largest cluster first
CLUSTER_NAMES = step["cluster_names"]             # null, or an explicit map, e.g. {0: ZERO, 2: ONE, 1: TWO}
EXPECTED_SIZES = step["expected_sizes"]           # null to skip the check
OUT = output_dir(cfg, "08_build_cluster_inputs")
OUT_XLSX = os.path.join(OUT, "accession_means_clusters.xlsx")
OUT_CSV = os.path.join(OUT, "blup_clusters.csv")


def cluster_id(v):
    s = str(v).replace("Cluster_", "").strip()
    return int(float(s))


def main():
    means = pd.read_excel(MEANS_XLSX, sheet_name=cfg["sheet_of"][SHEET], header=None)
    head, body = means.iloc[:HEADER].copy(), means.iloc[HEADER:].copy()
    acc = body.iloc[:, 0].astype(str).str.strip()

    blup = pd.read_excel(BLUP_XLSX, sheet_name=BLUP_SHEET, index_col=0)
    blup.index = blup.index.astype(str).str.strip()
    if "Cluster" not in blup.columns:
        sys.exit(f"No 'Cluster' column in sheet {BLUP_SHEET}")
    ids = blup.pop("Cluster").map(cluster_id)
    names = CLUSTER_NAMES or dict(zip(ids.value_counts().index, NAMES_BY_SIZE))
    labels = ids.map(names)
    if labels.isna().any():
        sys.exit("Cluster ids not covered by the cluster names")
    print("cluster id -> name:", names)

    sizes = labels.value_counts().to_dict()
    print("cluster sizes:", sizes)
    if EXPECTED_SIZES and sizes != EXPECTED_SIZES:
        sys.exit(f"Cluster sizes differ from expected {EXPECTED_SIZES}; check the cluster names")

    missing = sorted(set(acc) - set(labels.index))
    if missing:
        sys.exit(f"Accessions without a BLUP row: {missing}")

    body.iloc[:, 1] = acc.map(labels).values
    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as w:
        pd.concat([head, body]).to_excel(w, sheet_name=cfg["sheet_of"][SHEET], header=False, index=False)

    compounds = [str(c).strip() for c in head.iloc[2, 2:]]
    absent = [c for c in compounds if c not in blup.columns]
    if absent:
        sys.exit(f"Compounds missing from the BLUP sheet: {absent}")
    blup = blup[compounds].sort_index()
    with open(OUT_CSV, "w", newline="", encoding="utf-8-sig") as fh:
        wr = csv.writer(fh)
        wr.writerow(["", "Main Class"] + list(head.iloc[0, 2:]))
        wr.writerow(["", "Subclass"] + list(head.iloc[1, 2:]))
        wr.writerow(["Acc", "K"] + compounds)
        for a, row in blup.iterrows():
            wr.writerow([a, labels[a]] + [repr(float(x)) for x in row.values])
    print("written:", OUT_XLSX, OUT_CSV)


main()

finish_step(cfg, TREATMENT)

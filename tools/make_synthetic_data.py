#!/usr/bin/env python3
#
# VocProfile - make_synthetic_data.py
#
# Writes the small invented dataset of tests/data: a replicate workbook and a workbook of
# accession means with the layout docs/INPUTS.md describes. Accession names and all
# values are made up. Only the compound names and classes are real (config/compounds.tsv).
#
# in   config/compounds.tsv, config/compound_display_names.tsv
# out  <folder>/replicates.xlsx, <folder>/accession_means.xlsx
#
# Usage:
#   python tools/make_synthetic_data.py <output folder>
#
# Notes: docs/NOTES.md, test data
#
import os
import sys

import numpy as np
import pandas as pd

CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 7
N_ACCESSIONS = 48
GROUP_SIZES = (24, 15, 9)                    # three invented groups of accessions
SHEETS = ("Sheet1", "Sheet2", "Sheet3", "Sheet4")
REPLICATE_NOISE = 0.10                       # coefficient of variation among runs


def main(folder):
    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(SEED)
    compounds = pd.read_csv(os.path.join(CODE_DIR, "config", "compounds.tsv"), sep="\t")
    shown = pd.read_csv(os.path.join(CODE_DIR, "config", "compound_display_names.tsv"), sep="\t", comment="#")
    data_name = dict(zip(shown["display_name"], shown["data_name"]))
    names = [data_name.get(name, name) for name in compounds["Compound"]]
    n = len(names)

    accessions = [f"ACC-{i:03d}" for i in range(1, N_ACCESSIONS + 1)]
    group = np.repeat(np.arange(len(GROUP_SIZES)), GROUP_SIZES)
    rng.shuffle(group)
    k_group = rng.choice(["K1", "K2", "K3"], size=N_ACCESSIONS, p=[0.2, 0.3, 0.5])

    level = rng.normal(0.0, 1.2, n)                              # ln of the typical value of a compound
    shift = np.zeros((len(GROUP_SIZES), n))
    shift[1, rng.choice(n, 20, replace=False)] = rng.normal(1.2, 0.3, 20)
    shift[2, rng.choice(n, 20, replace=False)] = rng.normal(1.5, 0.3, 20)
    accession_effect = rng.normal(0.0, 0.45, (N_ACCESSIONS, n))
    truth = np.exp(level + shift[group] + accession_effect)

    # compounds found in part of the accessions only, more often in the second and third group
    partly = rng.choice(n, 14, replace=False)
    present = np.ones((N_ACCESSIONS, n), dtype=bool)
    for j in partly:
        chance = np.array([0.35, 0.75, 0.85])[group] * rng.uniform(0.6, 1.0)
        present[:, j] = rng.random(N_ACCESSIONS) < chance

    replicates = rng.choice([1, 2, 3, 4, 5], size=N_ACCESSIONS, p=[0.06, 0.2, 0.5, 0.16, 0.08])
    rows = [(accession, rep, i) for i, accession in enumerate(accessions) for rep in range(1, replicates[i] + 1)]
    index = np.array([i for _, _, i in rows])
    base = truth[index] * np.exp(rng.normal(0.0, REPLICATE_NOISE, (len(rows), n)))
    missed = np.zeros(n, dtype=bool)                             # compounds a run can miss
    missed[rng.choice(n, 18, replace=False)] = True
    detected = present[index] & ~(missed & (rng.random((len(rows), n)) < 0.08))

    scale = {                                                    # the four invented quantification datasets
        "Sheet1": np.exp(rng.normal(2.0, 0.5, n)),               # relative abundance
        "Sheet2": np.ones(n),                                    # external calibration
        "Sheet3": np.exp(rng.normal(0.0, 0.15, n)),              # internal calibration
        "Sheet4": np.full(n, 1.12),                              # spike corrected
    }
    main_class = compounds["Main_Class"].tolist()
    subclass = compounds["Subclass"].tolist()

    with pd.ExcelWriter(os.path.join(folder, "replicates.xlsx"), engine="openpyxl") as runs_book, \
            pd.ExcelWriter(os.path.join(folder, "accession_means.xlsx"), engine="openpyxl") as means_book:
        for sheet in SHEETS:
            values = base * scale[sheet] * np.exp(rng.normal(0.0, 0.02, base.shape))
            values = np.where(detected, np.round(values, 6), np.nan)
            runs = pd.DataFrame(values, columns=names)
            runs.insert(0, "Rep", [rep for _, rep, _ in rows])
            runs.insert(0, "Sample", [accession for accession, _, _ in rows])
            means = runs.drop(columns="Rep").groupby("Sample", sort=False).mean()
            runs.astype(object).where(runs.notna(), "ND").to_excel(runs_book, sheet_name=sheet, index=False)

            head = pd.DataFrame([[None, "Main Class"] + main_class, [None, "Subclass"] + subclass, [None, "K"] + names])
            body = pd.DataFrame(np.column_stack([np.array(accessions, dtype=object), k_group, means.loc[accessions].values.astype(object)]))
            body.iloc[:, 2:] = body.iloc[:, 2:].astype(float)
            pd.concat([head, body], ignore_index=True).to_excel(means_book, sheet_name=sheet, header=False, index=False)
    print(f"written: {folder}/replicates.xlsx ({len(rows)} runs, {N_ACCESSIONS} accessions, {n} compounds), "
          f"{folder}/accession_means.xlsx")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python tools/make_synthetic_data.py <output folder>")
    main(sys.argv[1])

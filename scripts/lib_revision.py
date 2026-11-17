#!/usr/bin/env python3
#
# VocProfile - lib_revision.py
#
# Shared by steps 18 to 21: settings, reading of the replicate workbook, filling of
# values without detection, BLUP estimates and the published cluster of each accession.
#
# in   REPLICATE_FILE, MEANS_CLUSTER_FILE, MARKER_TABLE, ODOR_THRESHOLD_TABLE
# out  results/<step>/ of the step that imports it
#
# Usage:
#   imported by steps 18 to 21
#
# Notes: docs/NOTES.md, steps 18 to 21
#
import os
import sys
import warnings

import numpy as np
import pandas as pd
from numpy.linalg import pinv
from sklearn.preprocessing import StandardScaler

from lib_config import load_config, input_path, output_dir
from lib_run import start_step, check_same_accessions, note, show_names_in_text_files
from lib_run import finish_step as write_qc

warnings.filterwarnings('ignore')

_script = os.path.splitext(os.path.basename(sys.argv[0]))[0]
cfg, step, args = load_config(_script[:2])
start_step(cfg, 'replicates', 'means_clusters')
check_same_accessions(cfg, 'replicates', 'means_clusters')
_settings = cfg["settings"]

RAW_FILE = input_path(cfg, 'replicates')
LABEL_FILE = input_path(cfg, 'means_clusters')
MAIN_SHEET = _settings.sheet(_settings.flat["MAIN_DATASET"])
DATASETS = {_settings.dataset_label(name): _settings.sheet(name) for name in _settings.datasets}
OUT_DIR = output_dir(cfg, _script)
SEED = cfg["seed"]
N_RESAMPLE = cfg["revision"]["N_RESAMPLE"]
N_SEEDS = cfg["revision"]["N_SEEDS"]
N_MONTE_CARLO = cfg["revision"]["N_MONTE_CARLO"]
CLUSTERS = cfg["clusters"]
PANEL = cfg["markers"]
THRESHOLDS = cfg["odor_thresholds"]
ACCESSION = cfg["accession_column"]
REPLICATE = cfg["replicate_column"]
HEADER_ROWS = len(cfg["header_rows"])
CLIP = _settings.value("CLIP_LIMIT")
FLOOR = _settings.value("ACCESSION_VARIANCE_FLOOR")
FRACTION = _settings.value("NONDETECT_FRACTION")
TREATMENT = (
    'values without detection in a run: ' + f'{FRACTION} x the smallest positive value of the compound over all runs, '
    'unless the step compares other rules (smallest value, smallest value / sqrt(2), zero)',
    'replicates per accession are unequal; the mean number enters the moment estimator of the accession variance',
    'a compound detected in part of the runs is used in full after the filling',
    'clusters compared with: second column of MEANS_CLUSTER_FILE',
)


def out(name):
    return os.path.join(OUT_DIR, name)


def from_step_18(name):
    """A file written by step 18, which must have run in the same WORK_DIR."""
    path = os.path.join(_settings.work_dir, "results", "18_test_cluster_stability", name)
    if not os.path.isfile(path):
        raise SystemExit(f"ERROR: {path} not found. Run step 18 first.")
    return path


def read_raw(sheet):
    df = pd.read_excel(RAW_FILE, sheet_name=sheet, na_values=cfg["nondetect_code"], engine='openpyxl')
    feat = [c for c in df.columns if c not in (ACCESSION, REPLICATE)]
    return df, feat


def fill_nd(df, feat, rule='half_min'):
    d = df.copy()
    for c in feat:
        pos = d[c][d[c] > 0]
        if not d[c].isna().any() or not len(pos):
            continue
        m = pos.min()
        d[c] = d[c].fillna({'half_min': m * FRACTION, 'min': m, 'min_sqrt2': m / np.sqrt(2), 'zero': 0.0}[rule])
    return d


def blup_matrix(df, feat):
    # accessions in sorted order, as in step 06
    geno = df[ACCESSION].values
    genos = np.array(sorted(pd.unique(geno)))
    n, k = len(geno), len(genos)
    idx = {g: i for i, g in enumerate(genos)}
    X = np.ones((n, 1))
    Z = np.zeros((n, k))
    for i, g in enumerate(geno):
        Z[i, idx[g]] = 1
    nrep = pd.Series(geno).value_counts().mean()
    est, vc = {}, {}
    for c in feat:
        y = df[c].values.astype(float)
        err = [np.var(y[geno == g], ddof=1) for g in genos if (geno == g).sum() > 1]
        se = np.mean(err) if err else np.var(y, ddof=1) * 0.5
        gm = np.array([y[geno == g].mean() for g in genos])
        sg = max(max(0, np.var(gm, ddof=1) - se / nrep), se * FLOOR)
        C = np.block([[X.T @ X, X.T @ Z], [Z.T @ X, Z.T @ Z + (se / sg) * np.eye(k)]])
        r = np.vstack([X.T @ y.reshape(-1, 1), Z.T @ y.reshape(-1, 1)])
        try:
            s = np.linalg.solve(C, r)
        except np.linalg.LinAlgError:
            s = pinv(C) @ r
        est[c] = (s[0] + s[1:k + 1].ravel()).astype(float)
        vc[c] = (sg, se)
    return pd.DataFrame(est, index=genos), vc


def zclip(M, clip=True):
    Z = pd.DataFrame(StandardScaler().fit_transform(M), index=M.index, columns=M.columns)
    return Z.clip(-CLIP, CLIP) if clip else Z


def published_labels():
    g = pd.read_excel(LABEL_FILE, sheet_name=MAIN_SHEET, header=None)
    return pd.Series(g.iloc[HEADER_ROWS:, 1].values, index=g.iloc[HEADER_ROWS:, 0].astype(str).str.strip().values)


def main_blup(rule='half_min'):
    df, feat = read_raw(MAIN_SHEET)
    return blup_matrix(fill_nd(df, feat, rule), feat)


def retention(ref, lab):
    """Largest number of accessions of each reference cluster that fall in one cluster of `lab`."""
    ct = pd.crosstab(ref, lab)
    return {c: int(ct.loc[c].max()) for c in CLUSTERS}


def finish_step(cfg, treatment=()):
    """Names shown to a reader in the tables of the step, then the QC file. blup.csv and
    ref_labels.csv are read by step 21 and keep the names of the data."""
    show_names_in_text_files(OUT_DIR, keep=("blup.csv", "ref_labels.csv"))
    write_qc(cfg, treatment)

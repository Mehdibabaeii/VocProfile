#!/usr/bin/env python3
#
# VocProfile - 21_run_nested_cv.py
#
# Nested cross-validation of the marker selection and of the classification of the
# clusters. Each outer fold is kept on disk; a second run continues from there (Table
# S17).
#
# in   REPLICATE_FILE, MEANS_CLUSTER_FILE (published cluster of each accession)
#      SEED, N_RESAMPLE, N_SEEDS, N_MONTE_CARLO, MARKER_TABLE, ODOR_THRESHOLD_TABLE
# out  results/21_run_nested_cv/S17_nested_cv_*.csv, S17_confusion_matrix.csv, S17_selection_frequency.csv
#      results/21_run_nested_cv/nested_cv_folds/ (one file per outer fold)
#      results/21_run_nested_cv/21_qc.txt
#
# Usage:
#   python scripts/21_run_nested_cv.py <config.txt>
#
# Notes: docs/NOTES.md, step 21
#
import pickle, sys
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif, RFE
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, recall_score
from lib_revision import *

K_GRID = list(range(2, 31))
TREES = [5, 10, 15, 20, 25, 30, 40, 50, 60, 100, 200]
FOLD_DIR = out('nested_cv_folds')
os.makedirs(FOLD_DIR, exist_ok=True)

B = pd.read_csv(from_step_18('blup.csv'), index_col=0)
y = pd.read_csv(from_step_18('ref_labels.csv'), index_col=0).iloc[:, 0].loc[B.index].values
X, feat = B.values, np.array(B.columns)
panel_idx = [list(feat).index(c) for c in PANEL]


def svc():
    return LinearSVC(random_state=SEED)


def mi_rank(A, t):
    return np.argsort(mutual_info_classif(A, t, discrete_features='auto', n_neighbors=3, random_state=SEED))[::-1]


def rfe_rank(A, t):
    return RFE(svc(), n_features_to_select=1, step=1).fit(A, t).ranking_


def select(Xtr, ytr):
    inner = StratifiedKFold(5, shuffle=True, random_state=SEED)
    a_mi, a_rfe, a_rf = np.zeros(len(K_GRID)), np.zeros(len(K_GRID)), np.zeros(len(TREES))
    for tr, va in inner.split(Xtr, ytr):
        sc = StandardScaler().fit(Xtr[tr]); A, V = sc.transform(Xtr[tr]), sc.transform(Xtr[va])
        r1, r2 = mi_rank(A, ytr[tr]), rfe_rank(A, ytr[tr])
        for j, k in enumerate(K_GRID):
            c = r1[:k]; a_mi[j] += svc().fit(A[:, c], ytr[tr]).score(V[:, c], ytr[va])
            c = np.where(r2 <= k)[0]; a_rfe[j] += svc().fit(A[:, c], ytr[tr]).score(V[:, c], ytr[va])
        for j, t in enumerate(TREES):
            a_rf[j] += RandomForestClassifier(n_estimators=t, random_state=SEED).fit(A, ytr[tr]).score(V, ytr[va])
    k_mi, k_rfe, n_t = K_GRID[int(np.argmax(a_mi))], K_GRID[int(np.argmax(a_rfe))], TREES[int(np.argmax(a_rf))]
    sc = StandardScaler().fit(Xtr); A = sc.transform(Xtr)
    s_mi = set(mi_rank(A, ytr)[:k_mi])
    s_rfe = set(np.where(rfe_rank(A, ytr) <= k_rfe)[0])
    imp = RandomForestClassifier(n_estimators=n_t, random_state=SEED).fit(A, ytr).feature_importances_
    s_rf = set(np.where(np.round(imp, 4) != 0)[0])
    cons = sorted(s_mi & s_rfe & s_rf) or sorted(s_mi)
    return sc, cons, dict(k_mi=k_mi, k_rfe=k_rfe, n_trees=n_t, n_consensus=len(cons))


splits = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED).split(X, y))
lo, hi = 0, len(splits)
stored = [f for f in range(lo, hi) if os.path.exists(os.path.join(FOLD_DIR, f'fold_{f:02d}.pkl'))]
if stored:
    print(f'{len(stored)} of {len(splits)} outer folds are read from {FOLD_DIR} and not computed again; '
          'delete that folder to compute all folds', flush=True)
    note(f'{len(stored)} of {len(splits)} outer folds reused from nested_cv_folds/')
for f in range(lo, hi):
    path = os.path.join(FOLD_DIR, f'fold_{f:02d}.pkl')
    if os.path.exists(path):
        continue
    tr, te = splits[f]
    sc, cons, info = select(X[tr], y[tr])
    A, T = sc.transform(X[tr]), sc.transform(X[te])
    pred = {'consensus': svc().fit(A[:, cons], y[tr]).predict(T[:, cons]),
            'all61': svc().fit(A, y[tr]).predict(T),
            'panel12': svc().fit(A[:, panel_idx], y[tr]).predict(T[:, panel_idx])}
    pickle.dump(dict(fold=f, test=te, pred=pred, cons=list(feat[cons]), **info), open(path, 'wb'))
    print('fold', f, 'consensus', len(cons), 'accuracy', round(accuracy_score(y[te], pred['consensus']), 3), flush=True)

done = sorted(os.listdir(FOLD_DIR))
if len(done) < len(splits):
    print(f'{len(done)}/{len(splits)} folds finished; run again to continue')
    sys.exit(0)

rows, cm, freq = [], {k: np.zeros((3, 3)) for k in ('consensus', 'all61', 'panel12')}, pd.Series(0, index=feat)
for name in done:
    d = pickle.load(open(os.path.join(FOLD_DIR, name), 'rb'))
    yt = y[d['test']]
    r = dict(fold=d['fold'], repeat=d['fold'] // 5, k_mi=d['k_mi'], k_rfe=d['k_rfe'], n_trees=d['n_trees'],
             n_consensus=d['n_consensus'], panel_overlap=len(set(d['cons']) & set(PANEL)))
    for k, p in d['pred'].items():
        r[f'accuracy_{k}'] = accuracy_score(yt, p)
        cm[k] += confusion_matrix(yt, p, labels=CLUSTERS)
    r['balanced_accuracy'] = balanced_accuracy_score(yt, d['pred']['consensus'])
    for c, v in zip(CLUSTERS, recall_score(yt, d['pred']['consensus'], labels=CLUSTERS, average=None)):
        r[f'recall_{c}'] = v
    rows.append(r)
    freq[d['cons']] += 1
R = pd.DataFrame(rows)
R.to_csv(out('S17_nested_cv_folds.csv'), index=False)


def summ(v):
    return dict(mean=v.mean(), SD=v.std(ddof=1), p2_5=np.percentile(v, 2.5), p97_5=np.percentile(v, 97.5))


summary = pd.DataFrame({k: summ(R[k]) for k in ['accuracy_consensus', 'balanced_accuracy', 'recall_ZERO', 'recall_ONE',
                                                 'recall_TWO', 'accuracy_all61', 'accuracy_panel12']}).T.round(3)
summary.to_csv(out('S17_nested_cv_summary.csv'))
c = cm['consensus']
cls = pd.DataFrame({'precision': np.diag(c) / c.sum(0), 'recall': np.diag(c) / c.sum(1),
                    'specificity': [(c.sum() - c[i].sum() - c[:, i].sum() + c[i, i]) / (c.sum() - c[i].sum()) for i in range(3)]},
                   index=CLUSTERS)
cls['F1'] = 2 * cls.precision * cls.recall / (cls.precision + cls.recall)
cls.round(3).to_csv(out('S17_nested_cv_by_cluster.csv'))
pd.DataFrame(c / 10, index=[f'true_{k}' for k in CLUSTERS], columns=[f'pred_{k}' for k in CLUSTERS]).round(1).to_csv(out('S17_confusion_matrix.csv'))
fq = (freq / len(R) * 100).sort_values(ascending=False).rename('selected_pct_of_folds')
fq.rename_axis('compound').to_csv(out('S17_selection_frequency.csv'))
print(summary.to_string()); print(cls.round(3).to_string()); print((c / 10).round(1))
print('consensus size: median %d, range %d-%d; repeat means %.3f-%.3f'
      % (R.n_consensus.median(), R.n_consensus.min(), R.n_consensus.max(),
         R.groupby('repeat').accuracy_consensus.mean().min(), R.groupby('repeat').accuracy_consensus.mean().max()))
print(fq.head(16).round(0).to_string())

finish_step(cfg, TREATMENT)

#!/usr/bin/env python3
#
# VocProfile - 18_test_cluster_stability.py
#
# Stability of the clusters: resampling of accessions and compounds, k-means
# initialisations, other algorithms, scalings, transformations, treatments of values
# without detection and compound subsets, agreement among the quantification datasets,
# silhouette per cluster (Table S8).
#
# in   REPLICATE_FILE, MEANS_CLUSTER_FILE (published cluster of each accession)
#      SEED, N_RESAMPLE, N_SEEDS, N_MONTE_CARLO, MARKER_TABLE, ODOR_THRESHOLD_TABLE
# out  results/18_test_cluster_stability/S8_*.csv, blup.csv, ref_labels.csv
#      results/18_test_cluster_stability/18_qc.txt
#
# Usage:
#   python scripts/18_test_cluster_stability.py <config.txt>
#
# Notes: docs/NOTES.md, step 18
#
from itertools import combinations
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score as ARI, silhouette_score, silhouette_samples
from lib_revision import *

rng = np.random.default_rng(SEED)
pub = published_labels()
B, _ = main_blup()
pub = pub.loc[B.index]
Z = zclip(B)


def km(X, k=3, seed=SEED):
    return KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(X)


ref = km(Z.values)
print('control: sizes', sorted(pd.Series(ref).value_counts(), reverse=True),
      'ARI vs published', round(ARI(pub, ref), 4), 'silhouette', round(silhouette_score(Z.values, ref), 4))
pub.rename('cluster').rename_axis('accession').to_csv(out('ref_labels.csv'))
B.to_csv(out('blup.csv'))

# resampling of accessions
rows, jac_rows = [], []
for k in (2, 3, 4):
    full = km(Z.values, k)
    aris, jac = [], {c: [] for c in range(k)}
    for b in range(N_RESAMPLE):
        idx = rng.choice(len(Z), int(0.8 * len(Z)), replace=False)
        Xs = np.clip(StandardScaler().fit_transform(B.values[idx]), -3, 3)
        lab = km(Xs, k, b)
        aris.append(ARI(full[idx], lab))
        for c in range(k):
            A = set(np.where(full[idx] == c)[0])
            jac[c].append(max(len(A & set(np.where(lab == d)[0])) / len(A | set(np.where(lab == d)[0])) for d in range(k)))
    a = np.array(aris)
    rows.append(dict(test='80% of accessions', k=k, mean_ARI=a.mean(), p2_5=np.percentile(a, 2.5),
                     p97_5=np.percentile(a, 97.5), share_ARI_above_0_8=(a > 0.8).mean()))
    if k == 3:
        for c in range(3):
            name = pub[full == c].mode()[0]
            jac_rows.append(dict(cluster=name, n=int((full == c).sum()), mean_jaccard=np.mean(jac[c]),
                                 p2_5=np.percentile(jac[c], 2.5)))

a = []
for b in range(N_RESAMPLE):
    cols = rng.choice(Z.shape[1], int(0.8 * Z.shape[1]), replace=False)
    a.append(ARI(ref, km(Z.values[:, cols], 3, b)))
a = np.array(a)
rows.append(dict(test='80% of compounds', k=3, mean_ARI=a.mean(), p2_5=np.percentile(a, 2.5),
                 p97_5=np.percentile(a, 97.5), share_ARI_above_0_8=(a > 0.8).mean()))
pd.DataFrame(rows).round(3).to_csv(out('S8_resampling.csv'), index=False)
pd.DataFrame(jac_rows).round(3).to_csv(out('S8_jaccard_by_cluster.csv'), index=False)
print(pd.DataFrame(rows).round(3).to_string(index=False))
print(pd.DataFrame(jac_rows).round(3).to_string(index=False))

# k-means initialization
res = []
for s in range(N_SEEDS):
    m = KMeans(n_clusters=3, random_state=s, n_init=10).fit(Z.values)
    res.append((ARI(pub, m.labels_), m.inertia_))
res = np.array(res)
seed_tab = pd.DataFrame([dict(n_seeds=N_SEEDS, share_ARI_1=(res[:, 0] > 0.999).mean(), share_ARI_ge_0_95=(res[:, 0] >= 0.95).mean(),
                              min_ARI=res[:, 0].min(), inertia_range_pct=(res[:, 1].max() - res[:, 1].min()) / res[:, 1].min() * 100)])
seed_tab.round(4).to_csv(out('S8_kmeans_seeds.csv'), index=False)
print(seed_tab.round(4).to_string(index=False))

# alternative choices
alt = []


def add(label, M=None, clip=True, labels=None, X=None):
    if labels is None:
        X = zclip(M, clip).values
        labels = km(X)
    r = retention(pub, pd.Series(labels, index=pub.index))
    alt.append(dict(analysis=label, ARI=ARI(pub, labels), silhouette=silhouette_score(X, labels),
                    sizes='/'.join(map(str, sorted(pd.Series(labels).value_counts(), reverse=True))),
                    kept_ZERO=r['ZERO'], kept_ONE=r['ONE'], kept_TWO=r['TWO']))


add('Published (BLUP, z-score, clip ±3, k-means)', B)
add('Ward linkage', labels=AgglomerativeClustering(n_clusters=3, linkage='ward').fit_predict(Z.values), X=Z.values)
add('No clipping', B, clip=False)
df_raw, feat = read_raw(MAIN_SHEET)
MEAN = fill_nd(df_raw, feat).groupby(ACCESSION)[feat].mean().loc[B.index]
add('Accession means instead of BLUP', MEAN)
add('log10 before scaling', np.log10(B.clip(lower=1e-4)))
add('Proportions of total', B.div(B.sum(axis=1), axis=0))
L = np.log(B.clip(lower=1e-4))
add('Centered log-ratio', L.sub(L.mean(axis=1), axis=0))
for rule, label in (('min_sqrt2', 'Non-detects = min/√2'), ('zero', 'Non-detects = 0'), ('min', 'Non-detects = min')):
    add(label, main_blup(rule)[0].loc[B.index])
det = df_raw[feat].notna().mean()
for th, label in ((0.8, 'Compounds detected in ≥80% of runs'), (1.0, 'Compounds detected in all runs')):
    add(f'{label} (n={(det >= th).sum()})', B[[c for c in feat if det[c] >= th]])
alt = pd.DataFrame(alt).round(3)
alt.to_csv(out('S8_sensitivity.csv'), index=False)
print(alt.to_string(index=False))

# quantification datasets
labs = {}
for name, sheet in DATASETS.items():
    d, f = read_raw(sheet)
    Bd, _ = blup_matrix(fill_nd(d, f), f)
    labs[name] = pd.Series(km(zclip(Bd).values), index=Bd.index).loc[B.index]
pair = pd.DataFrame([dict(dataset_1=a, dataset_2=b, ARI=round(ARI(labs[a], labs[b]), 3)) for a, b in combinations(labs, 2)])
pair.to_csv(out('S8_datasets_ARI.csv'), index=False)
keep = pd.DataFrame([dict(dataset=n, sizes='/'.join(map(str, sorted(labs[n].value_counts(), reverse=True))),
                          **{f'kept_{c}': v for c, v in retention(pub, labs[n]).items()}) for n in labs])
keep.to_csv(out('S8_datasets_retention.csv'), index=False)
print(pair.to_string(index=False)); print(keep.to_string(index=False))

ss = silhouette_samples(Z.values, pub.values)
sil = pd.DataFrame([dict(cluster=c, n=int((pub == c).sum()), mean_silhouette=ss[(pub == c).values].mean(),
                         share_negative=(ss[(pub == c).values] < 0).mean()) for c in CLUSTERS]).round(3)
sil.to_csv(out('S8_silhouette_by_cluster.csv'), index=False)
print(sil.to_string(index=False))

finish_step(cfg, TREATMENT)

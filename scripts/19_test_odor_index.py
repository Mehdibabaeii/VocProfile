#!/usr/bin/env python3
#
# VocProfile - 19_test_odor_index.py
#
# Relative odor index: Kruskal-Wallis across clusters with Benjamini-Hochberg and
# Bonferroni adjustment, category shares under three treatments of values without
# detection, and a Monte Carlo test of the ranking under compound-specific error
# (Table S18).
#
# in   REPLICATE_FILE, MEANS_CLUSTER_FILE (published cluster of each accession)
#      SEED, N_RESAMPLE, N_SEEDS, N_MONTE_CARLO, MARKER_TABLE, ODOR_THRESHOLD_TABLE
# out  results/19_test_odor_index/S18_*.csv
#      results/19_test_odor_index/19_qc.txt
#
# Usage:
#   python scripts/19_test_odor_index.py <config.txt>
#
# Notes: docs/NOTES.md, step 19
#
from scipy import stats
from lib_revision import *

ODT = {name: row["ODT"] for name, row in THRESHOLDS.items()}
CAT = {name: row["category"] for name, row in THRESHOLDS.items()}

rng = np.random.default_rng(SEED)
g = pd.read_excel(LABEL_FILE, sheet_name=MAIN_SHEET, header=None)
R = g.iloc[HEADER_ROWS:, 2:].apply(pd.to_numeric, errors='coerce')
R.columns = [str(c).strip() for c in g.iloc[2, 2:]]
R.index = g.iloc[HEADER_ROWS:, 0].astype(str).str.strip().values
lab = pd.Series(g.iloc[HEADER_ROWS:, 1].values, index=R.index)
cols = list(ODT)
odt = pd.Series(ODT)
I = R[cols].div(odt)
print('check: prenylthiol %.0f (284,159), α-pinene %.0f (58,143), camphene %.3f (0.078)'
      % (I['Prenylthiol'].mean(), I['α-Pinene'].mean(), I['Camphene'].mean()))


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    res = np.empty(n); res[o] = np.minimum(q, 1)
    return res


p = np.array([stats.kruskal(*[I.loc[lab == k, c].dropna() for k in CLUSTERS])[1] for c in cols])
kw = pd.DataFrame({'compound': cols, 'category': [CAT[c] for c in cols], 'n_accessions_detected': I.notna().sum().values,
                   'p_KW': p, 'q_BH': bh(p), 'p_Bonferroni': np.minimum(p * len(p), 1)})
kw.to_csv(out('S18_KW_adjusted.csv'), index=False)
for col in ('p_KW', 'q_BH', 'p_Bonferroni'):
    v = kw[col]
    print(col, 'significant:', int((v < 0.05).sum()), '| <0.001:', int((v < 0.001).sum()),
          '<0.01:', int(((v >= 0.001) & (v < 0.01)).sum()), '<0.05:', int(((v >= 0.01) & (v < 0.05)).sum()))

half_min = R[cols].apply(lambda s: s[s > 0].min() / 2)
shares, comp = [], []
for label, M in (('Detected values only', I), ('Non-detects = LOD/2', R[cols].fillna(half_min).div(odt)),
                 ('Non-detects = 0', R[cols].fillna(0).div(odt))):
    s = M.groupby(lab).mean().T[CLUSTERS].groupby(pd.Series(CAT)).sum()
    sh = s / s.sum() * 100
    shares.append(dict(treatment=label, **{f'sulfur_{c}': sh.loc['Sulfur', c] for c in CLUSTERS},
                       **{f'nonsulfur_{c}': 100 - sh.loc['Sulfur', c] for c in CLUSTERS}))
    ns = s.drop('Sulfur')
    comp.append((ns / ns.sum() * 100).assign(treatment=label))
pd.DataFrame(shares).round(1).to_csv(out('S18_sulfur_share.csv'), index=False)
pd.concat(comp).round(1).rename_axis('category').to_csv(out('S18_nonsulfur_composition.csv'))
print(pd.DataFrame(shares).round(1).to_string(index=False))
print(comp[0].round(1).to_string())

det = R['Prenylthiol']
print('prenylthiol detection (%):', (det.notna().groupby(lab).mean() * 100).round(1).to_dict(),
      '| median %.3f, range %.3f-%.3f' % (det.median(), det.min(), det.max()))

mean_all = I.mean()
cm = I.groupby(lab).mean().T[CLUSTERS]
mc = []
for fold in (3, 10):
    top = order = 0
    for _ in range(N_MONTE_CARLO):
        e = np.exp(rng.uniform(-np.log(fold), np.log(fold), len(cols)))
        top += cols[int(np.argmax(mean_all.values * e))] == 'Prenylthiol'
        s = cm.mul(e, axis=0).groupby(pd.Series(CAT)).sum()
        ns = 1 - s.loc['Sulfur'] / s.sum()
        order += ns['ZERO'] < ns['ONE'] < ns['TWO']
    mc.append(dict(max_fold_error=fold, prenylthiol_ranked_first_pct=top / N_MONTE_CARLO * 100,
                   nonsulfur_order_kept_pct=order / N_MONTE_CARLO * 100))
pd.DataFrame(mc).round(1).to_csv(out('S18_monte_carlo.csv'), index=False)
print(pd.DataFrame(mc).round(1).to_string(index=False))
print('prenylthiol / second-highest index: %.2f' % (mean_all['Prenylthiol'] / mean_all.drop('Prenylthiol').max()))

finish_step(cfg, TREATMENT)

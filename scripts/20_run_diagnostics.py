#!/usr/bin/env python3
#
# VocProfile - 20_run_diagnostics.py
#
# Runs and accessions in which each compound was detected, Kruskal-Wallis test of the
# accession effect, Shapiro-Wilk test of residuals, and correlation among the 12
# markers (Tables S9, S10, S17).
#
# in   REPLICATE_FILE, MEANS_CLUSTER_FILE (published cluster of each accession)
#      SEED, N_RESAMPLE, N_SEEDS, N_MONTE_CARLO, MARKER_TABLE, ODOR_THRESHOLD_TABLE
# out  results/20_run_diagnostics/S9_effective_n.csv, S10_diagnostics.csv, S17_marker_correlation.csv
#      results/20_run_diagnostics/20_qc.txt
#
# Usage:
#   python scripts/20_run_diagnostics.py <config.txt>
#
# Notes: docs/NOTES.md, step 20
#
from scipy import stats
from lib_revision import *

df, feat = read_raw(MAIN_SHEET)
eff = pd.DataFrame({'n_runs_detected': df[feat].notna().sum(),
                    'n_accessions_detected': df.groupby(ACCESSION)[feat].apply(lambda d: d.notna().any()).sum()})
eff['pct_runs'] = (eff.n_runs_detected / len(df) * 100).round(1)
eff.rename_axis('compound').to_csv(out('S9_effective_n.csv'))
print('runs per accession:', df.groupby(ACCESSION).size().value_counts().sort_index().to_dict())
print('compounds detected in all runs: %d; in all accessions: %d' % ((eff.n_runs_detected == len(df)).sum(),
                                                                     (eff.n_accessions_detected == df.Sample.nunique()).sum()))
print(eff.sort_values('n_runs_detected').head(6).to_string())

d = fill_nd(df, feat)
multi = d.groupby(ACCESSION)[feat[0]].transform('size') > 1
rows = []
for c in feat:
    p_kw = stats.kruskal(*[g[c].values for _, g in d.groupby(ACCESSION) if len(g) > 1])[1]
    y = np.log10(d[c])
    resid = (y - y.groupby(d[ACCESSION]).transform('mean'))[multi]
    rows.append(dict(compound=c, p_KW_accession=p_kw, p_Shapiro_residuals=stats.shapiro(resid)[1]))
diag = pd.DataFrame(rows)
diag.to_csv(out('S10_diagnostics.csv'), index=False)
print('KW accession effect at Bonferroni level: %d/%d' % ((diag.p_KW_accession < 0.05 / len(feat)).sum(), len(feat)))
print('Shapiro normality not rejected: %d/%d' % ((diag.p_Shapiro_residuals > 0.05).sum(), len(feat)))

B, vc = blup_matrix(d, feat)
C = np.log10(B[PANEL].clip(lower=1e-4)).corr()
C.round(3).to_csv(out('S17_marker_correlation.csv'))
iu = np.triu_indices(len(PANEL), 1)
a = C.abs().values[iu]
print('marker pairs |r|>0.9: %d/%d; median |r| %.2f' % ((a > 0.9).sum(), len(a), np.median(a)))
rep = pd.Series({c: vc[c][0] / (vc[c][0] + vc[c][1]) for c in feat})
print('repeatability mean %.3f, range %.3f-%.3f' % (rep.mean(), rep.min(), rep.max()))

finish_step(cfg, TREATMENT)

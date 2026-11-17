#!/usr/bin/env python3
#
# VocProfile - 07_evaluate_cluster_number.py
#
# Number of clusters: silhouette coefficient, Calinski-Harabasz index and within-
# cluster sum of squares (elbow criterion) for k = K_MIN to K_MAX on the unclipped
# z-scores of the BLUP matrix, with k-means and Ward compared at each k (Fig. S2,
# Table S8).
#
# in   REPLICATE_FILE, DATASETS
#      SEED, KMEANS_N_INIT, K_MIN, K_MAX, NONDETECT_FRACTION, ACCESSION_VARIANCE_FLOOR
# out  results/07_evaluate_cluster_number/cluster_number_<dataset>/ (criteria, BLUP matrix, ANOVA)
#      figures/07_evaluate_cluster_number/cluster_number_<dataset>/cluster_criteria_<dataset>.png
#      results/07_evaluate_cluster_number/07_qc.txt
#
# Usage:
#   python scripts/07_evaluate_cluster_number.py <config.txt>
#
# Notes: docs/NOTES.md, step 07
#
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import pinv
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.preprocessing import StandardScaler

from lib_config import load_config, input_path, require_sheets, output_dir
from lib_run import start_step, finish_step, save_figure, write_table, keep_data_names

cfg, step, args = load_config("07")
start_step(cfg, 'replicates')
keep_data_names()
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'runs without detection: set to NONDETECT_FRACTION of the smallest positive value of the compound, at replicate level',
    'z-score per compound without clipping',
    'the BLUP matrix is built with the accessions in sorted order',
)
SEED = cfg["seed"]
DPI = step["dpi"]

warnings.filterwarnings("ignore")
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
    "axes.edgecolor": "black",
    "axes.linewidth": 1.2,
})


class AccessionBLUPClusterAnalyzer:
    def __init__(self):
        self.raw_data = None
        self.accession_info = None
        self.feature_names = None
        self.blup_data = None
        self.normalized_data = None
        self.variance_components = None
        self.anova_results = None

    def load_data(self, df_raw, sample_col=cfg["accession_column"], replicate_col=cfg["replicate_column"]):
        self.feature_names = [c for c in df_raw.columns if c not in [sample_col, replicate_col]]
        self.accession_info = pd.DataFrame({
            "genotype": df_raw[sample_col].astype(str).values,
            "replicate": df_raw[replicate_col].values,
        })
        self.raw_data = df_raw[self.feature_names].reset_index(drop=True)

        n_geno = self.accession_info["genotype"].nunique()
        rep_counts = self.accession_info.groupby("genotype")["replicate"].count()
        print(f"Loaded: {self.raw_data.shape[0]} samples x {len(self.feature_names)} compounds")
        print(f"   Accessions: {n_geno} | reps/accession: min={rep_counts.min()}, "
              f"max={rep_counts.max()}, mean={rep_counts.mean():.1f}")

    def impute_lod2(self):
        """Set non-detects to half the smallest positive value of the compound, in place."""
        n_missing_total = int(self.raw_data.isna().sum().sum())
        for col in self.feature_names:
            vals = self.raw_data[col]
            if vals.isna().any():
                pos_min = vals[vals > 0].min()
                fill_val = pos_min * step["nondetect_fraction"] if pd.notna(pos_min) else 1e-9
                self.raw_data[col] = vals.fillna(fill_val)
        print(f"LOD/2 imputation done. Missing values filled: {n_missing_total}")

    def estimate_variance_components(self, trait):
        accessions = self.accession_info["genotype"].unique()
        accession_means, error_vars = [], []
        for g in accessions:
            mask = self.accession_info["genotype"] == g
            gd = trait[mask.values]
            accession_means.append(np.mean(gd))
            if len(gd) > 1:
                error_vars.append(np.var(gd, ddof=1))

        sigma_e2 = np.mean(error_vars) if error_vars else np.var(trait, ddof=1) * 0.5
        n_reps = self.accession_info.groupby("genotype").size().mean()
        sigma_g2 = max(0.0, np.var(accession_means, ddof=1) - sigma_e2 / n_reps)
        sigma_g2 = max(sigma_g2, sigma_e2 * step["variance_floor"])

        return {
            "genotype_variance": sigma_g2,
            "error_variance": sigma_e2,
            "heritability": sigma_g2 / (sigma_g2 + sigma_e2),
            "n_genotypes": len(accessions),
            "avg_reps_per_genotype": n_reps,
        }

    def compute_blup_estimates(self, trait, var_comp):
        accessions = self.accession_info["genotype"].unique()
        n_accessions = len(accessions)
        n_obs = len(trait)
        gidx = {g: i for i, g in enumerate(accessions)}

        X = np.ones((n_obs, 1))
        Z = np.zeros((n_obs, n_accessions))
        for i, g in enumerate(self.accession_info["genotype"]):
            Z[i, gidx[g]] = 1
        y = np.asarray(trait).reshape(-1, 1)

        lam = var_comp["error_variance"] / var_comp["genotype_variance"]
        C11 = X.T @ X
        C12 = X.T @ Z
        C21 = Z.T @ X
        C22 = Z.T @ Z + lam * np.eye(n_accessions)
        C = np.block([[C11, C12], [C21, C22]])
        rhs = np.vstack([X.T @ y, Z.T @ y])

        try:
            sol = np.linalg.solve(C, rhs)
        except np.linalg.LinAlgError:
            sol = pinv(C) @ rhs

        beta = sol[0]
        u = sol[1:n_accessions + 1]
        return {g: float(beta + u[i]) for i, g in enumerate(accessions)}

    def apply_blup_fitting(self):
        print(f"Fitting BLUP for {len(self.feature_names)} compounds...")
        blup_results, var_results = {}, {}
        for compound in self.feature_names:
            trait = self.raw_data[compound].values
            vc = self.estimate_variance_components(trait)
            var_results[compound] = vc
            blup_results[compound] = self.compute_blup_estimates(trait, vc)

        accessions = sorted(self.accession_info["genotype"].unique())
        matrix = [[blup_results[c][g] for c in self.feature_names] for g in accessions]
        self.blup_data = pd.DataFrame(matrix, index=accessions, columns=self.feature_names)
        self.variance_components = var_results

        h2_vals = [v["heritability"] for v in var_results.values()]
        print(f"BLUP fitting complete for {len(accessions)} accessions. "
              f"Mean repeatability={np.mean(h2_vals):.3f}, range=[{min(h2_vals):.3f}, {max(h2_vals):.3f}]")

    def normalize_blup(self):
        scaler = StandardScaler()
        arr = scaler.fit_transform(self.blup_data)
        self.normalized_data = pd.DataFrame(arr, index=self.blup_data.index, columns=self.blup_data.columns)
        print("BLUP matrix normalized (z-score) for clustering.")

    def perform_anova(self):
        scaler = StandardScaler()
        norm_raw = pd.DataFrame(
            scaler.fit_transform(self.raw_data[self.feature_names]),
            columns=self.feature_names,
        )
        norm_raw["genotype"] = self.accession_info["genotype"].values

        results = {}
        for compound in self.feature_names:
            groups = [g[compound].values for _, g in norm_raw.groupby("genotype")]
            try:
                f_stat, p_val = stats.f_oneway(*groups)
                grand_mean = norm_raw[compound].mean()
                ss_total = np.sum((norm_raw[compound] - grand_mean) ** 2)
                ss_between = sum(len(g) * (g[compound].mean() - grand_mean) ** 2
                                  for _, g in norm_raw.groupby("genotype"))
                eta_sq = ss_between / ss_total if ss_total > 0 else 0
                sig = "***" if p_val < LEVELS[0] else "**" if p_val < LEVELS[1] else "*" if p_val < LEVELS[2] else "ns"
                results[compound] = {"F_statistic": f_stat, "p_value": p_val,
                                      "eta_squared": eta_sq, "significance": sig}
            except Exception as e:
                results[compound] = {"F_statistic": np.nan, "p_value": np.nan,
                                      "eta_squared": np.nan, "significance": "error", "error": str(e)}

        self.anova_results = results
        n_sig = sum(1 for r in results.values() if r["significance"] not in ("ns", "error"))
        print(f"ANOVA complete: {n_sig}/{len(results)} compounds significant (p<0.05)")
        return results

    def evaluate_optimal_k(self, k_range=range(2, 8), methods=("kmeans", "hierarchical")):
        X = self.normalized_data.values
        all_results = {}

        for k in k_range:
            res = {}
            for method in methods:
                clusterer = (KMeans(n_clusters=k, random_state=SEED, n_init=step["n_init"]) if method == "kmeans"
                             else AgglomerativeClustering(n_clusters=k, linkage="ward"))
                labels = clusterer.fit_predict(X)
                if len(set(labels)) > 1:
                    res[method] = {
                        "silhouette": silhouette_score(X, labels),
                        "calinski_harabasz": calinski_harabasz_score(X, labels),
                        "labels": labels,
                    }
            if not res:
                continue
            # k-means or Ward, whichever has the higher silhouette at this k
            best_method = max(res.keys(), key=lambda m: res[m]["silhouette"])
            labels = res[best_method]["labels"]
            sil = res[best_method]["silhouette"]
            cal = res[best_method]["calinski_harabasz"]

            wcss = sum(np.sum((X[labels == i] - X[labels == i].mean(axis=0)) ** 2)
                       for i in range(k) if np.sum(labels == i) > 0)

            all_results[k] = {
                "best_method": best_method, "silhouette": sil, "calinski_harabasz": cal,
                "wcss": wcss, "labels": labels,
            }
            print(f"  k={k}: silhouette={sil:.3f}, calinski-harabasz={cal:.1f}, wcss={wcss:.1f}")

        # elbow: the k that follows the largest drop of the within-cluster sum of squares
        ks = sorted(all_results.keys())
        wcss_vals = [all_results[k]["wcss"] for k in ks]
        diffs = [wcss_vals[i] - wcss_vals[i + 1] for i in range(len(wcss_vals) - 1)]
        elbow_k = ks[diffs.index(max(diffs))] + 1 if diffs else ks[0]
        return all_results, elbow_k


METRIC_SPECS = [
    ("silhouette", "Silhouette Index", "#3B7DD8", max, "Score"),
    ("calinski_harabasz", "Calinski\u2013Harabasz Index", "#2E9E5B", max, "Index"),
    ("wcss", "Within-Cluster SS (Elbow)", "#7B4FA3", None, "WCSS"),
]


def plot_criteria(all_results, elbow_k, output_folder, sheet_name):
    k_vals = sorted(all_results.keys())
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.4))
    axes = axes.flatten()

    for idx, (metric, title, color, func, ylabel) in enumerate(METRIC_SPECS):
        ax = axes[idx]
        vals = [all_results[k][metric] for k in k_vals]

        if metric == "wcss":
            ax.plot(k_vals, vals, "o-", linewidth=2.2, markersize=9, color=color,
                     markeredgecolor="black", markeredgewidth=0.8)
            ax.axvline(x=elbow_k, color="black", linestyle="--", alpha=0.6, linewidth=1.5)
        else:
            ax.plot(k_vals, vals, "o-", linewidth=2.2, markersize=9, color=color,
                     markeredgecolor="black", markeredgewidth=0.8)
            if func is not None:
                best_k = k_vals[vals.index(func(vals))]
                ax.axvline(x=best_k, color="black", linestyle="--", alpha=0.5, linewidth=1.2)

        ax.set_xlabel("Number of clusters (k)", fontsize=11, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.set_title(f"{idx + 1}. {title}", fontsize=12, fontweight="bold")
        ax.set_xticks(k_vals)
        ax.grid(alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("black")
        ax.spines["bottom"].set_color("black")

    fig.suptitle(f"Cluster Validity Criteria Across k: {sheet_name}", fontsize=16, fontweight="bold", y=1.01)
    plt.tight_layout()
    save_path = os.path.join(output_folder, f"cluster_criteria_{sheet_name}.png")
    save_figure(save_path, dpi=DPI, bbox_inches="tight")
    plt.close("all")
    print(f"Saved: {save_path}")


def save_criteria_table(all_results, output_folder, sheet_name):
    df = pd.DataFrame([{
        "k": k,
        "Best_Method": all_results[k]["best_method"],
        "Silhouette": round(all_results[k]["silhouette"], 4),
        "Calinski_Harabasz": round(all_results[k]["calinski_harabasz"], 2),
        "WCSS": round(all_results[k]["wcss"], 2),
    } for k in sorted(all_results.keys())])

    save_path = os.path.join(output_folder, f"cluster_criteria_metrics_{sheet_name}.xlsx")
    write_table(df, save_path, index=False)
    print(f"Saved: {save_path}")
    return df


def main():
    file_path = input_path(cfg, step["input"])
    sheets = step["sheets"]
    require_sheets(file_path, [cfg["sheet_of"][name] for name in sheets])
    output_folder_base = os.path.join(output_dir(cfg, "07_evaluate_cluster_number"), step["output_folder"])
    k_range = range(step["k_min"], step["k_max"] + 1)

    summary_rows = []

    for sheet_name in sheets:
        print(f"\nDataset: {sheet_name}")
        output_folder = f"{output_folder_base}_{sheet_name}"
        os.makedirs(output_folder, exist_ok=True)

        df_raw = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], na_values=cfg["nondetect_code"], engine="openpyxl")

        analyzer = AccessionBLUPClusterAnalyzer()
        analyzer.load_data(df_raw, sample_col=cfg["accession_column"], replicate_col=cfg["replicate_column"])
        analyzer.impute_lod2()
        analyzer.apply_blup_fitting()
        analyzer.normalize_blup()
        analyzer.perform_anova()

        all_results, elbow_k = analyzer.evaluate_optimal_k(k_range=k_range)

        plot_criteria(all_results, elbow_k, output_folder, sheet_name)
        save_criteria_table(all_results, output_folder, sheet_name)

        write_table(analyzer.blup_data, os.path.join(output_folder, f"BLUP_matrix_{sheet_name}.xlsx"))
        anova_df = pd.DataFrame(analyzer.anova_results).T
        anova_df.index.name = "Compound"
        write_table(anova_df, os.path.join(output_folder, f"ANOVA_results_{sheet_name}.xlsx"))

        summary_rows.append({"Dataset": sheet_name, "Elbow_k": elbow_k})

    summary_df = pd.DataFrame(summary_rows)
    write_table(summary_df, f"{output_folder_base}_summary_across_datasets.xlsx", index=False)
    print("\nSummary across datasets")
    print(summary_df.to_string(index=False))
    print(f"\nSaved: {output_folder_base}_summary_across_datasets.xlsx")


main()

finish_step(cfg, TREATMENT)

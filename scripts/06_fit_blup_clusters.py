#!/usr/bin/env python3
#
# VocProfile - 06_fit_blup_clusters.py
#
# Variance components by the method of moments, BLUP per compound from the mixed-model
# equations, repeatability, ANOVA across accessions, and the clusters: z-score per
# compound, clipping, k-means and Ward clustering (Fig. 3a, Tables S10-S12).
#
# in   REPLICATE_FILE, DATASETS
#      SEED, N_CLUSTERS, KMEANS_N_INIT, CLIP_LIMIT, NONDETECT_FRACTION, ACCESSION_VARIANCE_FLOOR, CLUSTER_TABLE
# out  results/06_fit_blup_clusters/blup_<dataset>/blup_results_<dataset>.xlsx
#      figures/06_fit_blup_clusters/blup_<dataset>/cluster_figure_<dataset>.png
#      results/06_fit_blup_clusters/06_qc.txt
#
# Usage:
#   python scripts/06_fit_blup_clusters.py <config.txt>
#
# Notes: docs/NOTES.md, step 06
#
import gc
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.linalg import pinv
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler

from lib_config import load_config, input_path, require_sheets, output_dir, term_function
from lib_run import start_step, finish_step, save_figure, write_table, keep_data_names, save_legend, later, deferred

warnings.filterwarnings('ignore')

cfg, step, args = load_config("06")
start_step(cfg, 'replicates')
keep_data_names()
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'runs without detection: set to NONDETECT_FRACTION of the smallest positive value of the compound, at replicate level, before the BLUPs',
    'unequal numbers of replicates per accession (one to five) are used as they are; the mean number of replicates enters the variance among accessions',
    'the BLUP matrix is built with the accessions in sorted order',
    'z-score per compound, clipped at CLIP_LIMIT, then k-means and Ward; the one with the higher silhouette is kept',
    'clusters are numbered by size, largest first, unless CLUSTER_REMAP is set',
)
SEED = cfg["seed"]
DPI = step["dpi"]
OUT_ROOT = output_dir(cfg, "06_fit_blup_clusters")
T = term_function(cfg)                                  # axis labels and table headings
TT = term_function(cfg, "revised_terms_in_titles")      # titles and legend entries
CLUSTER_NAMES = cfg["labels"]["cluster_names"] or {}


def cluster_label(cluster, short=False):
    """Name of a cluster on the figure: its name in the label settings, else 'Cluster 0' / 'C0'."""
    cluster = int(str(cluster).replace('Cluster_', ''))
    submitted = f'C{cluster}' if short else f'Cluster {cluster}'
    return later(submitted, CLUSTER_NAMES.get(cluster, submitted))


class CannabisTerpeneBLUPAnalyzer:
    def __init__(self):
        self.raw_data = None
        self.blup_data = None
        self.normalized_data = None
        self.cluster_labels = None
        self.feature_names = None
        self.accession_info = None
        self.variance_components = None

    def load_data(self, data_path=None, data=None, sample_col=cfg["accession_column"], replicate_col=cfg["replicate_column"]):
        """Load GC-MS data with proper accession and replicate identification"""
        if data is not None:
            self.raw_data = data.copy()
        elif data_path:
            if data_path.endswith('.csv'):
                self.raw_data = pd.read_csv(data_path)
            elif data_path.endswith('.xlsx'):
                self.raw_data = pd.read_excel(data_path)

        if self.raw_data is not None:
            if sample_col not in self.raw_data.columns:
                print(f"Column '{sample_col}' not found. Available columns: {list(self.raw_data.columns)}")
                return

            if replicate_col not in self.raw_data.columns:
                print(f"Column '{replicate_col}' not found. Available columns: {list(self.raw_data.columns)}")
                return

            self.feature_names = [col for col in self.raw_data.columns
                                if col not in [sample_col, replicate_col]]

            self.accession_info = pd.DataFrame({
                'genotype': self.raw_data[sample_col],
                'replicate': self.raw_data[replicate_col],
                'sample_id': range(len(self.raw_data))
            })

            sample_names = [f"{sample}_{rep}" for sample, rep in
                          zip(self.raw_data[sample_col], self.raw_data[replicate_col])]

            compound_data = self.raw_data[self.feature_names].copy()
            compound_data.index = sample_names
            self.raw_data = compound_data

            print(f"Data loaded: {self.raw_data.shape[0]} samples, {len(self.feature_names)} compounds")
            print(f"   Accessions detected: {len(self.accession_info['genotype'].unique())}")

            rep_counts = self.accession_info.groupby('genotype')['replicate'].count()
            print(f"   Replicates per accession: min={rep_counts.min()}, max={rep_counts.max()}, mean={rep_counts.mean():.1f}")
            print(f"   Sample format example: {sample_names[0]}")
            print(f"   First 5 compounds: {self.feature_names[:5]}")

    def estimate_variance_components(self, trait):
        """Estimate variance components for BLUP using REML approach"""
        accessions = self.accession_info['genotype'].unique()
        n_accessions = len(accessions)
        n_obs = len(trait)

        X = np.ones((n_obs, 1))

        Z = np.zeros((n_obs, n_accessions))
        accession_to_idx = {g: i for i, g in enumerate(accessions)}

        for i, accession in enumerate(self.accession_info['genotype']):
            Z[i, accession_to_idx[accession]] = 1

        y = np.array(trait).reshape(-1, 1)

        accession_means = []
        for accession in accessions:
            mask = self.accession_info['genotype'] == accession
            accession_data = trait[mask]
            accession_means.append(np.mean(accession_data))

        error_vars = []
        for accession in accessions:
            mask = self.accession_info['genotype'] == accession
            accession_data = trait[mask]
            if len(accession_data) > 1:
                error_vars.append(np.var(accession_data, ddof=1))

        sigma_e2 = np.mean(error_vars) if error_vars else np.var(trait, ddof=1) * 0.5

        n_reps = self.accession_info.groupby('genotype').size().mean()
        sigma_g2 = max(0, (np.var(accession_means, ddof=1) - sigma_e2/n_reps))

        # accession variance not below 1 % of the error variance
        sigma_g2 = max(sigma_g2, sigma_e2 * step["variance_floor"])

        return {
            'genotype_variance': sigma_g2,
            'error_variance': sigma_e2,
            'heritability': sigma_g2 / (sigma_g2 + sigma_e2),
            'n_genotypes': n_accessions,
            'avg_reps_per_genotype': n_reps
        }

    def compute_blup_estimates(self, trait_data, variance_components):
        """Compute BLUP estimates for accession effects"""
        accessions = self.accession_info['genotype'].unique()
        n_accessions = len(accessions)
        n_obs = len(trait_data)

        X = np.ones((n_obs, 1))
        Z = np.zeros((n_obs, n_accessions))

        accession_to_idx = {g: i for i, g in enumerate(accessions)}
        for i, accession in enumerate(self.accession_info['genotype']):
            Z[i, accession_to_idx[accession]] = 1

        y = np.array(trait_data).reshape(-1, 1)

        sigma_g2 = variance_components['genotype_variance']
        sigma_e2 = variance_components['error_variance']
        lambda_val = sigma_e2 / sigma_g2

        XTX = X.T @ X
        XTZ = X.T @ Z
        ZTX = Z.T @ X
        ZTZ = Z.T @ Z
        I = np.eye(n_accessions)

        XTy = X.T @ y
        ZTy = Z.T @ y

        C11 = XTX
        C12 = XTZ
        C21 = ZTX
        C22 = ZTZ + lambda_val * I

        C = np.block([[C11, C12], [C21, C22]])
        rhs = np.vstack([XTy, ZTy])

        try:
            solution = np.linalg.solve(C, rhs)
        except np.linalg.LinAlgError:
            solution = pinv(C) @ rhs

        beta = solution[0]
        u = solution[1:n_accessions+1]

        blup_estimates = {}
        for i, accession in enumerate(accessions):
            blup_estimates[accession] = float(beta + u[i])

        return blup_estimates

    def apply_blup_fitting(self):
        """Apply BLUP fitting to all traits (compounds)"""
        if self.raw_data is None or self.accession_info is None:
            print("Please load data first")
            return

        print("Applying BLUP fitting with random accession and replicate effects...")

        compound_data = self.raw_data[self.feature_names].copy()

        blup_results = {}
        variance_results = {}

        print(f"Processing {len(self.feature_names)} compounds...")

        for compound in self.feature_names:
            trait_data = compound_data[compound].values

            var_components = self.estimate_variance_components(trait_data)
            variance_results[compound] = var_components

            blup_estimates = self.compute_blup_estimates(trait_data, var_components)
            blup_results[compound] = blup_estimates

            if var_components['heritability'] < 0.1:
                print(f"Warning: {compound}: Low repeatability ({var_components['heritability']:.3f})")

        accessions = sorted(self.accession_info['genotype'].unique())
        blup_matrix = []

        for accession in accessions:
            accession_row = []
            for compound in self.feature_names:
                accession_row.append(blup_results[compound][accession])
            blup_matrix.append(accession_row)

        self.blup_data = pd.DataFrame(blup_matrix,
                                     index=accessions,
                                     columns=self.feature_names)
        self.variance_components = variance_results

        print(f"BLUP fitting completed for {len(accessions)} accessions")

        avg_repeatability = np.mean([vc['heritability'] for vc in variance_results.values()])
        print(f"   Average repeatability: {avg_repeatability:.3f}")
        print(f"   Range of repeatabilities: {min(vc['heritability'] for vc in variance_results.values()):.3f} - {max(vc['heritability'] for vc in variance_results.values()):.3f}")

    def normalize_data(self, method='standard'):
        """Normalize the BLUP-fitted data"""
        if self.blup_data is None:
            print("Please apply BLUP fitting first")
            return

        print(f"Normalizing data using {method} scaling...")

        if method == 'standard':
            scaler = StandardScaler()
        elif method == 'minmax':
            scaler = MinMaxScaler()
        elif method == 'robust':
            scaler = RobustScaler()
        else:
            scaler = StandardScaler()

        normalized_array = scaler.fit_transform(self.blup_data)
        self.normalized_data = pd.DataFrame(normalized_array,
                                          index=self.blup_data.index,
                                          columns=self.blup_data.columns)
        # clip so that one compound does not dominate the clustering
        self.normalized_data = self.normalized_data.clip(lower=-step["clip_limit"], upper=step["clip_limit"])
        print("Data normalization completed (z-score, clipped to [-3, 3])")

    def perform_clustering(self, n_clusters=3, methods=['kmeans', 'hierarchical']):
        """Perform clustering using different methods"""
        if self.normalized_data is None:
            print("Please normalize data first")
            return

        print(f"Performing clustering with {len(methods)} methods...")

        results = {}
        X = self.normalized_data.values

        for method in methods:
            print(f"  - {method.title()}...")

            if method == 'kmeans':
                clusterer = KMeans(n_clusters=n_clusters, random_state=SEED, n_init=step["n_init"])
                labels = clusterer.fit_predict(X)

            elif method == 'hierarchical':
                clusterer = AgglomerativeClustering(n_clusters=n_clusters)
                labels = clusterer.fit_predict(X)

            elif method == 'dbscan':
                clusterer = DBSCAN(eps=0.5, min_samples=max(2, len(X)//10))
                labels = clusterer.fit_predict(X)
                n_clusters_found = len(set(labels)) - (1 if -1 in labels else 0)

            if len(set(labels)) > 1:
                silhouette = silhouette_score(X, labels)
                calinski = calinski_harabasz_score(X, labels)
                davies_bouldin = davies_bouldin_score(X, labels)
            else:
                silhouette = calinski = davies_bouldin = 0

            results[method] = {
                'labels': labels,
                'silhouette_score': silhouette,
                'calinski_harabasz_score': calinski,
                'davies_bouldin_score': davies_bouldin,
                'n_clusters': len(set(labels)) - (1 if -1 in labels else 0)
            }

        best_method = max(results.keys(), key=lambda k: results[k]['silhouette_score'])
        self.cluster_labels = results[best_method]['labels']

        print(f"Best clustering method: {best_method.title()}")
        print(f"   Number of clusters: {results[best_method]['n_clusters']}")
        print(f"   Silhouette Score: {results[best_method]['silhouette_score']:.3f}")
        print(f"   Calinski-Harabasz Score: {results[best_method]['calinski_harabasz_score']:.1f}")
        print(f"   Davies-Bouldin Score: {results[best_method]['davies_bouldin_score']:.3f}")

        return results

    def perform_anova_analysis(self):
        """Perform ANOVA for all compounds using raw data with replicates (normalized)"""
        if self.raw_data is None or self.accession_info is None:
            print("Please load data first")
            return

        print("Performing ANOVA analysis on all compounds...")

        scaler = StandardScaler()
        normalized_raw = pd.DataFrame(
            scaler.fit_transform(self.raw_data[self.feature_names]),
            index=self.raw_data.index,
            columns=self.feature_names
        )

        normalized_raw['genotype'] = self.accession_info['genotype'].values

        anova_results = {}

        for compound in self.feature_names:
            compound_data = normalized_raw[['genotype', compound]].copy()

            accession_groups = [group[compound].values for name, group in compound_data.groupby('genotype')]

            try:
                f_stat, p_value = stats.f_oneway(*accession_groups)

                if p_value < LEVELS[0]:
                    significance = '***'
                elif p_value < LEVELS[1]:
                    significance = '**'
                elif p_value < LEVELS[2]:
                    significance = '*'
                else:
                    significance = 'ns'

                grand_mean = compound_data[compound].mean()
                ss_total = np.sum((compound_data[compound] - grand_mean) ** 2)

                ss_between = 0
                for name, group in compound_data.groupby('genotype'):
                    group_mean = group[compound].mean()
                    group_size = len(group)
                    ss_between += group_size * (group_mean - grand_mean) ** 2

                eta_squared = ss_between / ss_total if ss_total > 0 else 0

                anova_results[compound] = {
                    'F_statistic': f_stat,
                    'p_value': p_value,
                    'significance': significance,
                    'eta_squared': eta_squared,
                    'n_genotypes': len(accession_groups),
                    'total_samples': len(compound_data)
                }

            except Exception as e:
                anova_results[compound] = {
                    'F_statistic': np.nan,
                    'p_value': np.nan,
                    'significance': 'error',
                    'eta_squared': np.nan,
                    'n_genotypes': len(accession_groups) if 'accession_groups' in locals() else 0,
                    'total_samples': len(compound_data) if 'compound_data' in locals() else 0,
                    'error': str(e)
                }

        self.anova_results = anova_results

        significant_compounds = sum(1 for r in anova_results.values()
                                  if r['significance'] not in ['ns', 'error'])
        total_compounds = len(anova_results)

        print(f"ANOVA completed for {total_compounds} compounds")
        print(f"   Significant compounds: {significant_compounds}/{total_compounds}")

        sig_counts = {'***': 0, '**': 0, '*': 0, 'ns': 0, 'error': 0}
        for result in anova_results.values():
            sig_counts[result['significance']] += 1

        print(f"   Significance levels:")
        print(f"     p < 0.001 (***): {sig_counts['***']}")
        print(f"     p < 0.01  (**):   {sig_counts['**']}")
        print(f"     p < 0.05  (*):     {sig_counts['*']}")
        print(f"     p >= 0.05 (ns):    {sig_counts['ns']}")
        if sig_counts['error'] > 0:
            print(f"     Errors:           {sig_counts['error']}")

        return anova_results

    def analyze_clusters(self):
        """Analyze cluster characteristics and create chemical profiles"""
        if self.cluster_labels is None:
            print("Please perform clustering first")
            return

        print("Analyzing cluster characteristics...")

        if not hasattr(self, 'anova_results'):
            self.perform_anova_analysis()

        analysis_data = self.normalized_data.copy()
        analysis_data['Cluster'] = self.cluster_labels

        blup_analysis = self.blup_data.copy()
        blup_analysis['Cluster'] = self.cluster_labels

        cluster_profiles = {}
        unique_clusters = sorted([c for c in set(self.cluster_labels) if c != -1])

        for cluster in unique_clusters:
            print(f"  Analyzing Cluster {cluster}...")

            cluster_data = analysis_data[analysis_data['Cluster'] == cluster]
            cluster_blup = blup_analysis[blup_analysis['Cluster'] == cluster]

            cluster_compounds = cluster_data.drop('Cluster', axis=1)
            cluster_blup_compounds = cluster_blup.drop('Cluster', axis=1)

            profile = {}
            profile['size'] = len(cluster_data)
            profile['percentage'] = (len(cluster_data) / len(analysis_data)) * 100
            profile['genotypes'] = cluster_data.index.tolist()

            profile['normalized_means'] = cluster_compounds.mean()
            profile['normalized_stds'] = cluster_compounds.std()

            profile['blup_means'] = cluster_blup_compounds.mean()
            profile['blup_stds'] = cluster_blup_compounds.std()

            profile['characteristic_compounds'] = profile['blup_means'].nlargest(10)

            discriminating_compounds = []

            for compound in self.feature_names:
                cluster_accessions = set(profile['genotypes'])

                raw_compound_data = self.raw_data[compound].copy()
                accession_labels = self.accession_info['genotype'].values

                cluster_values = []
                other_values = []

                for i, accession in enumerate(accession_labels):
                    if accession in cluster_accessions:
                        cluster_values.append(raw_compound_data.iloc[i])
                    else:
                        other_values.append(raw_compound_data.iloc[i])

                if len(cluster_values) > 1 and len(other_values) > 1:
                    try:
                        all_values = cluster_values + other_values
                        scaler = StandardScaler()
                        normalized_all = scaler.fit_transform(np.array(all_values).reshape(-1, 1)).flatten()

                        normalized_cluster = normalized_all[:len(cluster_values)]
                        normalized_other = normalized_all[len(cluster_values):]

                        f_stat, p_value = stats.f_oneway(normalized_cluster, normalized_other)
                        effect_size = (np.mean(normalized_cluster) - np.mean(normalized_other)) / np.sqrt(
                            ((np.var(normalized_cluster, ddof=1) + np.var(normalized_other, ddof=1)) / 2))

                        if p_value < LEVELS[0]:
                            significance = '***'
                        elif p_value < LEVELS[1]:
                            significance = '**'
                        elif p_value < LEVELS[2]:
                            significance = '*'
                        else:
                            significance = 'ns'

                        discriminating_compounds.append((compound, f_stat, p_value, effect_size, significance))

                    except Exception as e:
                        continue

            discriminating_compounds.sort(key=lambda x: x[1], reverse=True)
            profile['discriminating_compounds'] = discriminating_compounds

            cluster_heritabilities = {}
            for compound in self.feature_names:
                if compound in self.variance_components:
                    cluster_heritabilities[compound] = self.variance_components[compound]['heritability']

            profile['heritability_stats'] = {
                'mean_heritability': np.mean(list(cluster_heritabilities.values())),
                'compound_heritabilities': cluster_heritabilities
            }

            cluster_profiles[f'Cluster_{cluster}'] = profile

        self.cluster_profiles = cluster_profiles
        print(f"Cluster analysis completed for {len(cluster_profiles)} clusters")
        return cluster_profiles

    def create_visualizations(self):
        """Draw the nine-panel figure of the step."""
        if self.cluster_labels is None:
            print("Please perform clustering first")
            return


        fig = plt.figure(figsize=(20, 16))

        unique_clusters = np.unique(self.cluster_labels)
        n_clusters = len(unique_clusters)
        colors = plt.cm.Set1(np.linspace(0, 1, n_clusters))

        cluster_color_map = dict(cfg["cluster_colours_by_id"])

        ax1 = plt.subplot(3, 3, 1)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(self.normalized_data)

        for i, cluster in enumerate(unique_clusters):
            mask = self.cluster_labels == cluster
            cluster_color = cluster_color_map.get(cluster, colors[i])
            ax1.scatter(pca_result[mask, 0], pca_result[mask, 1],
                       c=[cluster_color], label=cluster_label(cluster),
                       alpha=0.7, s=60, edgecolors='black', linewidth=0.5)

        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=14, fontweight='bold')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=14, fontweight='bold')
        self.scatter_legends = [ax1.legend(loc='best', framealpha=0.9, fontsize=9)]
        ax1.grid(True, alpha=0.3)

        ax2 = plt.subplot(3, 3, 2)
        if len(self.normalized_data) > 5:
            perplexity = min(30, max(5, len(self.normalized_data)//4))
            tsne = TSNE(
                n_components=2,
                random_state=SEED,
                perplexity=perplexity,
                max_iter=1000,
                verbose=0,
                init='pca'
        )

            tsne_result = tsne.fit_transform(self.normalized_data)

            for i, cluster in enumerate(unique_clusters):
                mask = self.cluster_labels == cluster
                cluster_color = cluster_color_map.get(cluster, colors[i])
                ax2.scatter(tsne_result[mask, 0], tsne_result[mask, 1],
                           c=[cluster_color], label=cluster_label(cluster),
                           alpha=0.7, s=60, edgecolors='black', linewidth=0.5)

            self.scatter_legends.append(ax2.legend(loc='best', framealpha=0.9, fontsize=9))
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, 'Insufficient samples\nfor t-SNE',
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('t-SNE - Insufficient Data')

        ax3 = plt.subplot(3, 3, 3)
        unique, counts = np.unique(self.cluster_labels, return_counts=True)
        bar_colors = [cluster_color_map.get(c, colors[i]) for i, c in enumerate(unique)]
        bars = ax3.bar(range(len(unique)), counts, color=bar_colors, alpha=0.7, edgecolor='black', linewidth=1)
        ax3.set_xlabel('Cluster', fontsize=10)
        ax3.set_ylabel(deferred(T)('Number of Genotypes'), fontsize=16, fontweight='bold')
        ax3.set_title('Cluster Size Distribution', fontsize=12, fontweight='bold')
        ax3.set_xticks(range(len(unique)))
        ax3.set_xticklabels([cluster_label(i, short=True) for i in unique])
        ax3.grid(True, alpha=0.3, axis='y')

        for bar, count in zip(bars, counts):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    str(count), ha='center', va='bottom', fontsize=9, fontweight='bold')

        ax4 = plt.subplot(3, 3, 4)
        if self.variance_components:
            heritabilities = [vc['heritability'] for vc in self.variance_components.values()]
            ax4.hist(heritabilities, bins=15, alpha=0.7, color='skyblue', edgecolor='black')
            ax4.axvline(np.mean(heritabilities), color='red', linestyle='--', linewidth=2,
                       label=f'Mean: {np.mean(heritabilities):.3f}')
            ax4.set_xlabel(deferred(T)('Heritability'), fontsize=14, fontweight='bold')
            ax4.set_ylabel('Number of Compounds', fontsize=14, fontweight='bold')
            ax4.legend(fontsize=9)
            ax4.grid(True, alpha=0.3)

        if hasattr(self, 'cluster_profiles'):
            ax5 = plt.subplot(3, 3, 5)

            all_top_compounds = set()
            for profile in self.cluster_profiles.values():
                all_top_compounds.update(profile['characteristic_compounds'].head(8).index)

            heatmap_data = []
            cluster_names = []
            for cluster_name, profile in sorted(self.cluster_profiles.items()):
                cluster_means = []
                for compound in sorted(all_top_compounds):
                    cluster_means.append(profile['blup_means'][compound])
                heatmap_data.append(cluster_means)
                cluster_names.append(cluster_label(cluster_name, short=True))

            compound_labels = [c[:12] + '...' if len(c) > 12 else c for c in sorted(all_top_compounds)]

            heatmap_df = pd.DataFrame(heatmap_data,
                                    index=cluster_names,
                                    columns=compound_labels)

            sns.heatmap(heatmap_df, annot=True, cmap='RdYlBu_r',
                       center=heatmap_df.values.mean(), ax=ax5,
                       cbar_kws={'shrink': 0.8},
                       fmt='.2f',
                       annot_kws={'size': 7},
                       linewidths=0.5, linecolor='gray')

            ax5.set_title('Top Compounds by Cluster (BLUP)', fontsize=11, fontweight='bold')
            ax5.set_xlabel('Compounds', fontsize=9)
            ax5.set_ylabel('Clusters', fontsize=9)
            plt.setp(ax5.get_xticklabels(), rotation=45, ha='right', fontsize=7)
            plt.setp(ax5.get_yticklabels(), rotation=0, fontsize=9)

        ax6 = plt.subplot(3, 3, 6)
        if len(set(self.cluster_labels)) > 1:
            rf = RandomForestClassifier(n_estimators=100, random_state=SEED)
            rf.fit(self.normalized_data, self.cluster_labels)

            importances = rf.feature_importances_
            indices = np.argsort(importances)[::-1][:12]

            ax6.barh(range(len(indices)), importances[indices], alpha=0.7,
                    color='lightgreen', edgecolor='black', linewidth=0.8)
            ax6.set_yticks(range(len(indices)))
            compound_names = [self.feature_names[i][:20] for i in indices]
            ax6.set_yticklabels(compound_names, fontsize=8)
            ax6.set_xlabel('Feature Importance', fontsize=10)
            ax6.set_title('Most Discriminating Compounds', fontsize=12, fontweight='bold')
            ax6.grid(True, alpha=0.3, axis='x')
        ax7 = plt.subplot(3, 3, 7)
        try:

            linkage_matrix = linkage(self.normalized_data, method='ward')
            n_leaves = len(self.normalized_data)

            idx_to_cluster = {i: self.cluster_labels[i] for i in range(n_leaves)}

            node_majority = idx_to_cluster.copy()
            node_leaves = {i: [i] for i in range(n_leaves)}

            for i, row in enumerate(linkage_matrix):
                node_id = n_leaves + i
                left_child = int(row[0])
                right_child = int(row[1])

                combined_leaves = node_leaves[left_child] + node_leaves[right_child]
                node_leaves[node_id] = combined_leaves

                clusters = [idx_to_cluster[leaf] for leaf in combined_leaves]
                majority = max(set(clusters), key=clusters.count)
                node_majority[node_id] = majority

            def link_color_func_majority(link_id):
                majority = node_majority.get(link_id, -1)
                return cluster_color_map.get(majority, '#A8A9AE')

            dendro = dendrogram(
                linkage_matrix,
                ax=ax7,
                labels=self.normalized_data.index,
                leaf_font_size=3,
                color_threshold=0,
                above_threshold_color='#A8A9AE',
                link_color_func=link_color_func_majority,
                orientation='top',
                distance_sort='ascending'
            )

            accession_to_color = {g: cluster_color_map.get(c, 'grey')
                                 for g, c in zip(self.normalized_data.index, self.cluster_labels)}

            for lbl in ax7.get_xmajorticklabels():
                accession_name = lbl.get_text()
                if accession_name in accession_to_color:
                    lbl.set_color(accession_to_color[accession_name])
                    lbl.set_fontweight('bold')
                    lbl.set_fontsize(3)

            ax7.set_xlabel(deferred(T)('Genotypes'), fontsize=14, fontweight='bold')
            ax7.set_ylabel('Euclidean Distance', fontsize=14, fontweight='bold')
            ax7.tick_params(axis='x', which='major', pad=2)
            plt.setp(ax7.get_xticklabels(), rotation=90, ha='center')

            unique_clusters = np.unique(self.cluster_labels)
            legend_elements = [
                Patch(facecolor=cluster_color_map[c], edgecolor='black', linewidth=0.5, label=cluster_label(c))
                for c in sorted(unique_clusters) if c in cluster_color_map
            ]

            self.cluster_legend = legend_elements
            self.dendrogram_legend = ax7.legend(handles=legend_elements, loc='upper right', fontsize=7,
                                                framealpha=0.9, title='Clusters', title_fontsize=8)
            ax7.margins(x=0.01)

        except Exception as e:
            ax7.text(0.5, 0.5, f'Dendrogram Error:\n{str(e)[:50]}...',
                     ha='center', va='center', transform=ax7.transAxes, fontsize=9)
            ax7.set_title('Hierarchical Dendrogram - Error')

        ax8 = plt.subplot(3, 3, 8)
        if self.variance_components:
            compounds = list(self.variance_components.keys())[:15]
            gen_vars = [self.variance_components[c]['genotype_variance'] for c in compounds]
            err_vars = [self.variance_components[c]['error_variance'] for c in compounds]

            x = np.arange(len(compounds))
            width = 0.35

            ax8.bar(x - width/2, gen_vars, width, label=deferred(TT)('Genotype Variance'),
                   alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.8)
            ax8.bar(x + width/2, err_vars, width, label='Error Variance',
                   alpha=0.7, color='coral', edgecolor='black', linewidth=0.8)

            ax8.set_xlabel('Compounds', fontsize=9)
            ax8.set_ylabel('Variance', fontsize=10)
            ax8.set_title('Variance Components (First 15 Compounds)', fontsize=11, fontweight='bold')
            ax8.set_xticks(x)
            compound_labels = [c[:10] + '...' if len(c) > 10 else c for c in compounds]
            ax8.set_xticklabels(compound_labels, rotation=45, ha='right', fontsize=7)
            ax8.legend(fontsize=8)
            ax8.grid(True, alpha=0.3, axis='y')

        ax9 = plt.subplot(3, 3, 9)
        if hasattr(self, 'cluster_profiles'):
            cluster_sizes = [profile['size'] for profile in self.cluster_profiles.values()]
            cluster_names = [cluster_label(name, short=True) for name in self.cluster_profiles.keys()]

            wcss_values = []
            for i, (cluster_name, profile) in enumerate(self.cluster_profiles.items()):
                cluster_data = self.normalized_data[self.cluster_labels == i]
                if len(cluster_data) > 1:
                    centroid = cluster_data.mean()
                    wcss = np.sum((cluster_data - centroid) ** 2).sum()
                    wcss_values.append(wcss)
                else:
                    wcss_values.append(0)

            bars = ax9.bar(cluster_names, wcss_values, alpha=0.7, color='orange',
                          edgecolor='black', linewidth=1)
            ax9.set_xlabel('Clusters', fontsize=10)
            ax9.set_ylabel('Within-Cluster Sum of Squares', fontsize=9)
            ax9.set_title('Cluster Compactness', fontsize=12, fontweight='bold')
            ax9.tick_params(axis='x', rotation=0)
            ax9.grid(True, alpha=0.3, axis='y')

            for bar, value in zip(bars, wcss_values):
                ax9.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(wcss_values)*0.02,
                        f'{value:.1f}', ha='center', va='bottom', fontsize=8)

        plt.tight_layout()
        save_figure(os.path.join(OUT_ROOT, 'cluster_figure.png'), panel=True, dpi=DPI, bbox_inches='tight',
                    after_layout=self.revised_legend_titles)
        plt.show()

        print("Visualizations created and saved successfully")

    def revised_legend_titles(self, figure, box):
        """Legend titles of the label state, set once the layout is fixed."""
        self.dendrogram_legend.set_title(cfg["labels"]["dendrogram_legend_title"])
        for legend in self.scatter_legends:
            legend.set_title(cfg["labels"]["scatter_legend_title"])

    def generate_report(self):
        """Generate detailed analysis report"""
        if not hasattr(self, 'cluster_profiles'):
            print("Please analyze clusters first")
            return

        print(f"\nSummary:")
        print(f"   Total accessions analyzed: {len(self.blup_data)}")
        print(f"   Total compounds detected: {len(self.feature_names)}")
        print(f"   Original samples (with replicates): {len(self.raw_data)}")
        print(f"   Clusters identified: {len(self.cluster_profiles)}")

        if self.variance_components:
            avg_repeatability = np.mean([vc['heritability'] for vc in self.variance_components.values()])
            print(f"   Average repeatability: {avg_repeatability:.3f}")

            high_h2_compounds = [comp for comp, vc in self.variance_components.items()
                               if vc['heritability'] > 0.6]
            print(f"   High repeatability compounds (>0.6): {len(high_h2_compounds)}")

        print(f"\nCHEMOTYPE PROFILES:")

        for cluster_name, profile in self.cluster_profiles.items():
            print(f"\n{cluster_name.upper()}:")
            print(f"   Size: {profile['size']} accessions ({profile['percentage']:.1f}%)")
            print(f"   Accessions: {', '.join(profile['genotypes'][:5])}")
            if len(profile['genotypes']) > 5:
                print(f"        ... and {len(profile['genotypes'])-5} more")

            print(f"\n   TOP CHARACTERISTIC COMPOUNDS (BLUP estimates):")
            for i, (compound, concentration) in enumerate(profile['characteristic_compounds'].head(7).items(), 1):
                repeatability = ""
                if compound in self.variance_components:
                    h2 = self.variance_components[compound]['heritability']
                    repeatability = f" (h2={h2:.3f})"
                print(f"      {i}. {compound}: {concentration:.3f}{repeatability}")

            print(f"\n   MOST DISCRIMINATING COMPOUNDS (vs other clusters):")
            for i, (compound, f_stat, p_value, effect_size, significance) in enumerate(profile['discriminating_compounds'][:10], 1):
                print(f"      {i}. {compound}: F={f_stat:.2f}, p={p_value:.2e}, Effect={effect_size:.2f} {significance}")

            if 'heritability_stats' in profile:
                print(f"\n   REPEATABILITY STATISTICS:")
                print(f"      Average repeatability: {profile['heritability_stats']['mean_heritability']:.3f}")

                sorted_h2 = sorted(profile['heritability_stats']['compound_heritabilities'].items(),
                                 key=lambda x: x[1], reverse=True)
                print(f"      Highest repeatability compounds:")
                for compound, h2 in sorted_h2[:3]:
                    print(f"        - {compound}: {h2:.3f}")

        if hasattr(self, 'anova_results'):
            print(f"\nANOVA SUMMARY FOR ALL COMPOUNDS:")
            print(f"   (Based on normalized raw data with replicates)")

            sorted_anova = sorted(self.anova_results.items(),
                                key=lambda x: x[1]['F_statistic'] if not np.isnan(x[1]['F_statistic']) else 0,
                                reverse=True)

            print(f"\n   TOP 20 MOST SIGNIFICANT COMPOUNDS:")
            print(f"   {'Compound':<30} {'F-stat':<10} {'p-value':<12} {'Eta2':<8} {'Sig':<6}")
            print(f"   {'-'*66}")

            for i, (compound, results) in enumerate(sorted_anova[:20], 1):
                if not np.isnan(results['F_statistic']):
                    print(f"   {i:2d}. {compound:<26} {results['F_statistic']:8.2f} {results['p_value']:10.2e} {results['eta_squared']:6.3f} {results['significance']:>4s}")

        if hasattr(self, 'anova_results'):
            sig_counts = {'***': 0, '**': 0, '*': 0, 'ns': 0}
            for result in self.anova_results.values():
                if result['significance'] in sig_counts:
                    sig_counts[result['significance']] += 1

            print(f"\n   ANOVA significance distribution:")
            print(f"     - Highly significant (p<0.001): {sig_counts['***']} compounds")
            print(f"     - Very significant (p<0.01): {sig_counts['**']} compounds")
            print(f"     - Significant (p<0.05): {sig_counts['*']} compounds")
            print(f"     - Non-significant (p>=0.05): {sig_counts['ns']} compounds")

        return self.cluster_profiles

    def export_results(self, filename='cannabis_terpene_analysis_results'):
        """Export all results to Excel file with multiple sheets"""
        if not hasattr(self, 'cluster_profiles'):
            print("Please complete analysis first")
            return

        print(f"Exporting results to {filename}.xlsx...")

        with pd.ExcelWriter(f'{filename}.xlsx', engine='openpyxl') as writer:
            if self.raw_data is not None:
                export_raw = self.raw_data.copy()
                export_raw = export_raw.join(self.accession_info.set_index('sample_id')[['genotype', 'replicate']])
                write_table(export_raw.rename(columns=T), writer, sheet_name='1_Raw_Data')

            if self.blup_data is not None:
                export_blup = self.blup_data.copy()
                if hasattr(self, 'cluster_labels'):
                    cluster_mapping = dict(zip(self.blup_data.index, self.cluster_labels))
                    export_blup['Cluster'] = [cluster_mapping.get(idx, -1) for idx in export_blup.index]
                write_table(export_blup, writer, sheet_name='2_BLUP_Estimates')

            if self.normalized_data is not None:
                export_normalized = self.normalized_data.copy()
                export_normalized['Cluster'] = self.cluster_labels
                write_table(export_normalized, writer, sheet_name='3_Normalized_Data')

            if self.variance_components:
                variance_data = []
                for compound, components in self.variance_components.items():
                    variance_data.append({
                        'Compound': compound,
                        'Genotype_Variance': components['genotype_variance'],
                        'Error_Variance': components['error_variance'],
                        'Heritability': components['heritability'],
                        'N_Genotypes': components['n_genotypes'],
                        'Avg_Reps_Per_Genotype': components['avg_reps_per_genotype']
                    })

                variance_df = pd.DataFrame(variance_data)
                write_table(variance_df.rename(columns=T), writer, sheet_name='4_Variance_Components', index=False)

            profile_data = []
            for cluster_name, profile in self.cluster_profiles.items():
                for compound in self.feature_names:
                    profile_data.append({
                        'Cluster': cluster_name,
                        'Compound': compound,
                        'BLUP_Mean': profile['blup_means'].get(compound, 0),
                        'BLUP_Std': profile['blup_stds'].get(compound, 0),
                        'Normalized_Mean': profile['normalized_means'].get(compound, 0),
                        'Normalized_Std': profile['normalized_stds'].get(compound, 0),
                        'Heritability': self.variance_components.get(compound, {}).get('heritability', 0)
                    })

            profile_df = pd.DataFrame(profile_data)
            write_table(profile_df.rename(columns=T), writer, sheet_name='5_Cluster_Profiles', index=False)

            discriminating_data = []
            for cluster_name, profile in self.cluster_profiles.items():
                for item in profile.get('discriminating_compounds', []):
                    data_dict = {'Cluster': cluster_name}

                    col_names = ['Compound', 'F_Statistic', 'P_Value', 'Effect_Size', 'Significance']
                    for i in range(min(len(item), len(col_names))):
                        data_dict[col_names[i]] = item[i]

                    discriminating_data.append(data_dict)

            discriminating_df = pd.DataFrame(discriminating_data)
            write_table(discriminating_df, writer, sheet_name='6_Discriminating_Compounds', index=False)

            accession_assignments = []
            for i, accession in enumerate(self.blup_data.index):
                cluster_label = self.cluster_labels[i]
                accession_assignments.append({
                    'Genotype': accession,
                    'Cluster': cluster_label,
                    'Cluster_Name': f'Cluster_{cluster_label}' if cluster_label != -1 else 'Outlier'
                })

            assignment_df = pd.DataFrame(accession_assignments)
            write_table(assignment_df.rename(columns=T), writer, sheet_name='7_Accession_Assignments', index=False)

            if not hasattr(self, 'anova_results'):
                self.perform_anova_analysis()

            anova_df = pd.DataFrame.from_dict(self.anova_results, orient='index')
            write_table(anova_df.rename(columns=T), writer, sheet_name='8_ANOVA_Accessions', index=True)

        print("Results exported.")


def load_and_impute_sheet(file_path, sheet_name, na_values=cfg["nondetect_code"], engine="openpyxl"):
    """Load one sheet and set non-detects to half the smallest positive value of the compound."""
    df = pd.read_excel(file_path, sheet_name=cfg["sheet_of"][sheet_name], na_values=na_values, engine=engine)
    feature_cols = [c for c in df.columns if c not in [cfg["accession_column"], cfg["replicate_column"]]]
    n_missing = int(df[feature_cols].isna().sum().sum())
    for col in feature_cols:
        vals = df[col]
        if vals.isna().any():
            pos_min = vals[vals > 0].min()
            fill_val = pos_min * step["nondetect_fraction"] if pd.notna(pos_min) else 1e-9
            df[col] = vals.fillna(fill_val)
    print(f"  LOD/2 imputation: {n_missing} missing values filled.")
    gc.collect()
    return df


FILE_PATH = input_path(cfg, step["input"])
SHEETS = step["sheets"]
require_sheets(FILE_PATH, [cfg["sheet_of"][name] for name in SHEETS])
N_CLUSTERS = cfg["n_clusters"]
OUTPUT_BASE = os.path.join(OUT_ROOT, step["output_folder"])
# None: clusters are renumbered by size, largest = 0. Otherwise a map {old label: new label}.
MANUAL_REMAP = step["cluster_remap"]

for sheet_name in SHEETS:
    print(f"\nDataset: {sheet_name}")
    output_folder = f"{OUTPUT_BASE}_{sheet_name}"
    os.makedirs(output_folder, exist_ok=True)

    df = load_and_impute_sheet(FILE_PATH, sheet_name)
    print(f"  Loaded: {df.shape[0]} rows x {df.shape[1]-2} compounds, "
          f"{df[cfg["accession_column"]].nunique()} genotypes")

    analyzer = CannabisTerpeneBLUPAnalyzer()
    analyzer.load_data(data=df, sample_col=cfg["accession_column"], replicate_col=cfg["replicate_column"])
    analyzer.apply_blup_fitting()
    analyzer.normalize_data(method="standard")
    analyzer.perform_clustering(n_clusters=N_CLUSTERS, methods=["kmeans", "hierarchical"])

    # renumber the clusters by size, largest = 0
    _labels = analyzer.cluster_labels
    _sizes = {c: np.sum(_labels == c) for c in np.unique(_labels)}
    _sorted = sorted(_sizes.keys(), key=lambda c: -_sizes[c])
    _remap_auto = {old_c: new_c for new_c, old_c in enumerate(_sorted)}

    _remap = MANUAL_REMAP if MANUAL_REMAP is not None else _remap_auto
    analyzer.cluster_labels = np.array([_remap[l] for l in _labels])
    print(f"  Cluster sizes: { {c: int(np.sum(_labels==c)) for c in np.unique(_labels)} }")
    print(f"  Remapping applied: {_remap}")
    analyzer.perform_anova_analysis()
    analyzer.analyze_clusters()

    analyzer.create_visualizations()
    save_figure(f"{output_folder}/cluster_figure_{sheet_name}.png", panel=True, dpi=DPI, bbox_inches="tight")
    save_legend(f"{output_folder}/legend_clusters_{sheet_name}.png", analyzer.cluster_legend,
                [patch.get_label() for patch in analyzer.cluster_legend],
                cfg["labels"]["dendrogram_legend_title"], DPI, fontsize=10, title_fontsize=11)
    plt.close("all")

    analyzer.export_results(
        f"{output_folder}/blup_results_{sheet_name}"
    )

    analyzer.generate_report()

    del analyzer, df
    gc.collect()
    print(f"\n  Done: {sheet_name}")

print("\nAll datasets processed.")

finish_step(cfg, TREATMENT)

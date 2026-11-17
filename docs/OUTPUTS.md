# Outputs

All paths are relative to `WORK_DIR`. `<dataset>` is the label of a quantification dataset
(`LABEL_<DATASET>`), `<set>` an input set of `config/input_sets.tsv`. Every step also writes
`results/<step>/<id>_qc.txt`.

## 01_report_compounds

| File | Content |
|---|---|
| `figures/.../compound_report/all_compounds_structures_A4_high_quality.png` | structures of all compounds with a SMILES as a four-column A4 grid (Fig. 1) |
| `figures/.../compound_report/all_compounds_structures.png` | the same structures grouped by class |
| `figures/.../compound_report/stacked_bar_classes_with_counts.png` | number of compounds per class and subclass |
| `figures/.../compound_report/compound_structures/*.png` | one drawing per compound |
| `figures/.../compound_report/complete_compound_table.html` | compound table with structures |
| `results/.../compound_report/complete_compound_table_with_images.xlsx` | the same table as a workbook (Table S6) |

## 02_describe_accession_means

| File | Content |
|---|---|
| `results/.../statistics_<dataset>/unified_statistical_report_original_<dataset>.xlsx` | statistics per compound, subclass and main class on the accession means as they are |
| `results/.../statistics_<dataset>/unified_statistical_report_processed_<dataset>.xlsx` | the same after filling |
| `results/.../statistics_<dataset>/processed_and_filled_data_<dataset>.xlsx` | accession means after filling |
| `results/.../statistics_<dataset>/log_transformed_data_<dataset>.xlsx` | natural logarithm of the filled means |
| `results/.../statistics_by_K_<dataset>/` | the report per K group and a consolidated workbook |
| `figures/.../statistics_<dataset>/`, `statistics_by_K_<dataset>/` | share of accessions without detection per compound |

## 03_compare_datasets, 04_summarise_replicates, 05_fill_nondetects

| File | Content |
|---|---|
| `results/03_compare_datasets/dataset_comparison.xlsx` | coefficient of variation within and among accessions per dataset, rank correlation between datasets, Shapiro-Wilk tests |
| `results/04_summarise_replicates/replicate_statistics.xlsx` | statistics of each compound over runs, one sheet per dataset |
| `results/05_fill_nondetects/accession_means_filled.xlsx` | accession means of every dataset after filling |

## 06_fit_blup_clusters

| File | Content |
|---|---|
| `results/.../blup_<dataset>/blup_results_<dataset>.xlsx` | sheets: runs after filling, BLUP estimates with the cluster number, scaled matrix, variance components and repeatability, cluster profiles, compounds that distinguish each cluster, accession assignments, ANOVA across accessions |
| `figures/.../blup_<dataset>/cluster_figure_<dataset>.png` | nine panels: dendrogram (Fig. 3a), PCA, t-SNE, cluster sizes, repeatability, cluster heat map and others |
| `figures/.../blup_<dataset>/legend_clusters_<dataset>.png` | the cluster legend alone, for the assembled Fig. 3 |
| `figures/.../cluster_figure.png` | the figure of the last dataset processed |

## 07_evaluate_cluster_number

| File | Content |
|---|---|
| `results/.../cluster_number_<dataset>/cluster_criteria_metrics_<dataset>.xlsx` | per k: method with the higher silhouette, silhouette, Calinski-Harabasz index, within-cluster sum of squares |
| `results/.../cluster_number_<dataset>/BLUP_matrix_<dataset>.xlsx`, `ANOVA_results_<dataset>.xlsx` | the matrix evaluated and the ANOVA across accessions |
| `results/.../cluster_number_summary_across_datasets.xlsx` | k at the elbow per dataset |
| `figures/.../cluster_number_<dataset>/cluster_criteria_<dataset>.png` | the three criteria against k |

## 08_build_cluster_inputs

| File | Content |
|---|---|
| `results/.../accession_means_clusters.xlsx` | accession means of the main dataset with the cluster name in the second column |
| `results/.../blup_clusters.csv` | BLUP estimates in the same layout |

## 09_summarise_classes

`figures/.../<set>/class_summaries_<label>/`, where `<label>` is the dataset, or
`blup_clusters` for the two cluster sets:

| File | Content |
|---|---|
| `1_<label>_main_class_distribution.png` | share of each main class |
| `2_<label>_subclass_analysis.png`, `donut_<label>.png` | shares of the subclasses |
| `combined_subclass_mainclass_<label>.png` | number of compounds per main class and per subclass (Fig. 2a, 2b) |
| `3_<label>_detailed_compound_view.png` | compounds by mean value |
| `4_` to `9_<label>_accession_*.png` | profiles of the accessions, stacked by compound, main class and all compounds, as values and as percentages (plot 8 of `clusters_all`: Fig. 2d; plot 9 of `clusters_markers`: Fig. 3c) |
| `legend_compounds_<label>.png` | the compound legend of plot 9 alone |

## 10_plot_detection, 11_plot_distributions, 12_plot_filled_distributions

| File | Content |
|---|---|
| `figures/10_plot_detection/detection/3_<dataset>_detailed_compound_view.png` | number of accessions with a detection per compound (Fig. 2c) |
| `figures/11_plot_distributions/distributions_<dataset>/violin_box_<dataset>.png` | distributions of the compounds, values as they are |
| `figures/12_plot_filled_distributions/distributions_filled_<dataset>/violin_box_filled_<dataset>.png` | the same after filling |

## 13_plot_group_profiles

`figures/.../<set>/group_profiles_<label>/` (`group_profiles_markers_<label>/` for the marker
set): plots 4 to 9, group means stacked by main class, subclass and compound, as values and
as percentages. `K` in a file name stands for the group of the input set.

## 14_plot_group_distributions

`figures/.../<set>/plots_<dataset>/`:

| File | Content |
|---|---|
| `<class>_distributions.png` | violin and box plots of the compounds of a class by group, with Tukey letters |
| `compound_heatmap.png` | clustered heat map, accessions by compounds, clipped z-scores |
| `compound_pca.png` | PCA of the accessions |
| `compound_correlation_heatmap.png` | clustered Pearson correlation matrix between compounds (set `clusters_log1p_fig3b`: Fig. 3b) |
| `legend_main_class.png`, `legend_subclass.png` | the two legends of Fig. 3b alone (set `clusters_log1p_fig3b`) |

## 15_export_correlations

`results/.../correlation_exports/<version>/`, with the versions `1_raw` (values as they are),
`2_lod2` (filled) and `3_lod2_log` (filled, ln(1 + x)):

| File | Content |
|---|---|
| `CORRELATION_ANALYSIS.xlsx` | correlation matrix, all pairs, the 50 strongest pairs, summary, correlations by main class |
| `HEATMAP_ANALYSIS.xlsx` | clipped z-scores in clustered order, means and standard deviations per cluster, summary |

## 16_compute_odor_index

| File | Content |
|---|---|
| `results/.../odor_index_accessions.xlsx` | relative odor index per accession and compound |
| `results/.../odor_index_cluster_means.xlsx` | mean per cluster and compound |
| `results/.../odor_index_compounds.xlsx` | per compound: threshold, value over all accessions, cluster means, test across clusters |
| `results/.../odor_index_per_cluster.xlsx`, `detection_per_cluster.xlsx` | cluster values and detection counts |
| `results/.../odor_index_zscore_per_accession.xlsx`, `odor_index_zscore_cluster_medians.xlsx`, `odor_index_scalings.xlsx` | scaled versions |
| `results/.../odor_index_one_vs_rest.xlsx` | each cluster against the other two, per compound, with adjusted p and Cliff's delta |
| `results/.../aroma_clusters.xlsx`, `aroma_cluster_odor_index_means.xlsx` | clusters of accessions formed on the index, for comparison |
| `figures/.../plot1_concentration_vs_odor_index.png` | concentration against index (Fig. 4a) |
| `figures/.../plot2_odor_index_ranking.png` | ranking of the compounds (Fig. S5) |
| `figures/.../plot3_category_distribution.png` | number of compounds per odor category (Fig. 4b) |
| `figures/.../plot_radar_odor_index.png`, `plot4_spider.png`, `plot_spider_*.png` | radar figures per cluster |
| `figures/.../plot_heatmap_odor_index.png`, `plot5_heatmap.png`, `plot_heatmap_*.png` | heat maps |
| `figures/.../plot_scaling_comparison.png`, `aroma_clustering_comparison.png` | comparisons of scalings and of clusterings |

## 17_plot_odor_radar

| File | Content |
|---|---|
| `figures/.../radar/radar_<version>_CLASS.png` | one line per cluster |
| `figures/.../radar/radar_<version>_SAMPLES.png` | one line per accession (`v4_pareto_accessions`: Fig. 4c) |
| `results/.../radar/radar_data_all_versions.xlsx` | values behind the figures |
| `results/.../odor_index_*.xlsx` | the index tables as computed by this step |

Versions: `v1_original_zscore3_0to100`, `v2_no_scaling_log10`, `v3_zscore_accessions_clip3`,
`v4_pareto_accessions`.

## 18 to 21

| File | Content |
|---|---|
| `results/18_test_cluster_stability/S8_resampling.csv`, `S8_jaccard_by_cluster.csv` | agreement with the clusters after resampling accessions and compounds |
| `.../S8_kmeans_seeds.csv` | agreement over k-means seeds |
| `.../S8_sensitivity.csv` | other algorithms, scalings, transformations, fill rules and compound subsets |
| `.../S8_datasets_ARI.csv`, `S8_datasets_retention.csv` | agreement among the quantification datasets |
| `.../S8_silhouette_by_cluster.csv` | silhouette per cluster |
| `.../blup.csv`, `ref_labels.csv` | BLUP matrix and cluster of each accession, read by step 21 |
| `results/19_test_odor_index/S18_KW_adjusted.csv` | Kruskal-Wallis across clusters with adjusted p |
| `.../S18_sulfur_share.csv`, `S18_nonsulfur_composition.csv` | category shares under three treatments of values without detection |
| `.../S18_monte_carlo.csv` | ranking under compound-specific error |
| `results/20_run_diagnostics/S9_effective_n.csv` | runs and accessions with a detection per compound |
| `.../S10_diagnostics.csv` | test of the accession effect and of residuals per compound |
| `.../S17_marker_correlation.csv` | correlation among the markers |
| `results/21_run_nested_cv/S17_nested_cv_*.csv`, `S17_confusion_matrix.csv`, `S17_selection_frequency.csv` | nested cross-validation |
| `.../nested_cv_folds/` | one file per outer fold, used to continue an interrupted run |

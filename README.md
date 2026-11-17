# VocProfile

volatile profiling and odor screening of cannabis accessions

## Citation

Babaei M, Goulet C, Torkamaneh D (2026). Volatile profiling and estimated odor-activity
analysis of commercial drug-type cannabis accessions. *bioRxiv*.
[https://doi.org/10.64898/2026.08.27.747640](https://doi.org/10.64898/2026.08.27.747640)

- Data: Figshare, [https://doi.org/10.6084/m9.figshare.33359469](https://doi.org/10.6084/m9.figshare.33359469)
- Feature selection: CannaFeatML, release v1.0.0, [https://github.com/Mehdibabaeii/CannaFeatML](https://github.com/Mehdibabaeii/CannaFeatML)

Cite the version you ran. It is the content of the file `VERSION` and the first line of
every QC file a run writes.

VocProfile takes GC measurements of volatile compounds in a panel of accessions and
produces the tables and figures of the paper: descriptive statistics,
variance components with BLUPs and repeatability, clusters of accessions with the criteria
used to choose their number, profiles by biosynthetic class, correlations among compounds,
and a relative odor index that weights concentrations by published odor thresholds. Every
choice of the analysis is a setting in one configuration file, and every step leaves a
record of what it read and how it treated values without detection.

## Overview

| Step | Script | What it computes | In the paper |
|---|---|---|---|
| 01 | `01_report_compounds.py` | compound table with structures, structure grid, class counts | Fig. 1, Table S6 |
| 02 | `02_describe_accession_means.py` | descriptive statistics per compound, class and K group | supplementary statistics |
| 03 | `03_compare_datasets.py` | agreement among the four quantification datasets | |
| 04 | `04_summarise_replicates.py` | statistics over runs per compound | |
| 05 | `05_fill_nondetects.py` | accession means after filling | |
| 06 | `06_fit_blup_clusters.py` | variance components, BLUPs, repeatability, k-means clusters, nine-panel figure | Fig. 3a, Fig. S2 |
| 07 | `07_evaluate_cluster_number.py` | Calinski-Harabasz index, elbow criterion and silhouette for k = 2 to 7 | choice of k |
| 08 | `08_build_cluster_inputs.py` | accession means and BLUPs with the cluster of each accession | |
| 09 | `09_summarise_classes.py` | profiles by main class, subclass and compound, per input set | Fig. 2a, 2b, 2d, Fig. 3c |
| 10 | `10_plot_detection.py` | number of accessions in which each compound was detected | Fig. 2c |
| 11 | `11_plot_distributions.py` | distributions of the compounds by class | |
| 12 | `12_plot_filled_distributions.py` | the same after filling | |
| 13 | `13_plot_group_profiles.py` | stacked profiles of K groups and clusters | |
| 14 | `14_plot_group_distributions.py` | box plots by group, heat map, PCA, correlation matrix, per input set | Fig. 3b |
| 15 | `15_export_correlations.py` | tables behind the heat map and the correlation matrix | |
| 16 | `16_compute_odor_index.py` | relative odor index per accession, cluster and compound, with figures | Fig. 4a, 4b, Figs. S5, S6 |
| 17 | `17_plot_odor_radar.py` | radar figures of the relative odor index in four scalings | Fig. 4c |
| 18 | `18_test_cluster_stability.py` | stability of the clusters under resampling and other analysis choices | Table S8 |
| 19 | `19_test_odor_index.py` | adjusted tests and sensitivity of the relative odor index | Table S18 |
| 20 | `20_run_diagnostics.py` | detection counts, tests of residuals, correlation among markers | Tables S9, S10, S17 |
| 21 | `21_run_nested_cv.py` | nested cross-validation of the marker panel | Table S17 |

Step 07 reports the three criteria named in the Methods of the paper (Calinski-Harabasz,
elbow, silhouette) for k = 2 to 7.

### Feature selection with CannaFeatML

The 12 marker compounds of `config/markers.tsv` were selected with
[CannaFeatML](https://github.com/Mehdibabaeii/CannaFeatML) (release v1.0.0), a separate pipeline
that is not part of this repository. The target of the selection is the cluster of each accession
written by step 08 (`BLUP_CLUSTER_FILE`: ZERO, ONE, TWO), and the features are the 61 BLUP
estimates of the same file. CannaFeatML ran three selectors:

| Selector | Model | Setting tuned |
|---|---|---|
| MI-SVM | mutual information ranking, linear support vector machine | number of features, 2 to 30 |
| RFE-SVM | recursive feature elimination with a linear support vector machine | number of features, 2 to 30 |
| Random forest | importance ranking | number of trees, 5 to 200 |

The accessions were split into a training part (80%) and a test part (20%), stratified by
cluster. Scaling and the three selectors were fitted on the training part, and the settings
were chosen by the accuracy on the test part. A compound is a marker when all three
selectors keep it; the 12 compounds kept by all three form the panel. Every model used
`random_state = 42`.

Because the settings were chosen on the test part, the panel is evaluated again here
without that dependence. Step 21 repeats the whole selection inside a nested
cross-validation (outer loop: stratified five-fold, repeated ten times; inner loop:
stratified five-fold for the settings), fits a linear support vector classifier on the
consensus of each outer training part and tests it on the outer test part. It reports how
often each compound is selected and compares the accuracy with the full 61-compound
profile. Step 20 gives the correlations among the 12 markers.

To use another panel, run CannaFeatML on your data with the clusters of step 08 as
classes, write the compounds it keeps in a copy of `markers.tsv` (one name per line under
the header `compound`, spelled as in the data), and point `MARKER_TABLE` at the copy.

## Dependencies

Python 3.12 and the packages of `requirements.txt`, at the versions the pipeline was run
with:

| Package | Version |
|---|---|
| numpy | 2.2.6 |
| pandas | 2.3.3 |
| scipy | 1.18.1 |
| scikit-learn | 1.9.1 |
| statsmodels | 0.15.0 |
| matplotlib | 3.11.2 |
| seaborn | 0.13.2 |
| openpyxl | 3.1.5 |
| rdkit | 2026.09.2 |
| Pillow | 13.0.0 |
| adjustText | 1.4.0 |

```bash
python3 -m venv venv
```

```bash
venv/bin/pip install -r requirements.txt
```

`run.sh` uses the `python3` on the path, or the interpreter named in the variable `PYTHON`.
No step needs network access.

## Input files

`docs/INPUTS.md` describes each file with an invented example. The data of the paper are
in the [Figshare record](https://doi.org/10.6084/m9.figshare.33359469) and the supplementary
tables of the article; they are not in
this repository.

### Compound list

`COMPOUND_FILE`. One row per compound: name, SMILES, main class, subclass. The list of the
paper is shipped as `config/compounds.tsv` and is the default.

### Replicate workbook

`REPLICATE_FILE`. One sheet per quantification dataset, one row per run, the accession name
in `ACCESSION_COLUMN`, the run number in `REPLICATE_COLUMN`, one column per compound, `ND`
or an empty cell for a run without detection. Take the replicate-level tables of the
published record, one sheet per dataset, and give the sheet names in `SHEET_<DATASET>`.

### Accession means

`MEANS_FILE`. One sheet per dataset with the same sheet names. Three header rows (main
class, subclass, compound), the accession name in the first column, the K group of the
accession in the second, an empty cell where the compound was not detected in any run of
the accession.

### Cluster files

`MEANS_CLUSTER_FILE` and `BLUP_CLUSTER_FILE` are written by step 08 and read from
`WORK_DIR` by the later steps. Nothing has to be prepared.

## Running it

```bash
cp config.txt my_config.txt
```

Set `WORK_DIR`, `INPUT_DIR`, the file names and the sheet names in the copy, then:

```bash
bash run.sh my_config.txt
```

| Command | Effect |
|---|---|
| `bash run.sh my_config.txt` | all steps |
| `bash run.sh my_config.txt 06` | one step |
| `bash run.sh my_config.txt from 09` | from a step to the end |
| `bash run.sh my_config.txt --dry-run` | print the plan, run nothing |
| `bash run.sh my_config.txt --set LABELS=submitted` | change a setting for this run |
| `bash run.sh --list` | list the steps |

A run prints its title, the config file, the folders it reads and writes, the settings
changed on the command line and the value of `LABELS`. Each step is announced with its
script and its log. The run closes with one row per step (status, seconds, start time); the
same rows are in `logs/run.tsv`. A step that fails stops the run with a non-zero exit, and
the steps after it are reported as not run.

A step can be run alone, with the same arguments:

```bash
python scripts/14_plot_group_distributions.py my_config.txt --input-set clusters_log1p_fig3b
```

`--input-set` limits steps 09, 13 and 14 to one row of `config/input_sets.tsv`.

Steps 18 to 21 run k-means many times and are much faster with one thread:

```bash
OMP_NUM_THREADS=1 bash run.sh my_config.txt from 18
```

Step 21 stores each outer fold of the cross-validation in
`results/21_run_nested_cv/nested_cv_folds/`, so that an interrupted run can go on. A later
call reads the folds found there instead of computing them, finishes in a few seconds, and
prints how many folds were reused. To compute everything again, for example after the data
or a setting changed, delete that folder first:

```bash
rm -r /path/to/work/results/21_run_nested_cv/nested_cv_folds
```

On a cluster with a scheduler, submit the same command. A template for SLURM, with
placeholders:

```bash
#!/bin/bash
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=03:00:00
export OMP_NUM_THREADS=1
bash run.sh /path/to/my_config.txt
```

## Configuration

`config.txt` is the only file to edit. One setting per line, `KEY = VALUE`; `#` starts a
comment; lists are separated by commas. Each setting has a comment that says what it does
and what the paper used. A key the pipeline does not know, a key that is missing or a file
that does not exist stops the run at the start and names the key.

| Group | Settings |
|---|---|
| Run | `WORK_DIR`, `INPUT_DIR`, `BANNER`, `LABELS` |
| Input files and their layout | `COMPOUND_FILE`, `REPLICATE_FILE`, `MEANS_FILE`, `MEANS_CLUSTER_FILE`, `BLUP_CLUSTER_FILE`, `ACCESSION_COLUMN`, `REPLICATE_COLUMN`, `NONDETECT_CODE`, `MEANS_HEADER_ROWS`, `DATASETS`, `SHEET_*`, `LABEL_*`, `SHORT_*`, `MAIN_DATASET`, `ABUNDANCE_DATASET` |
| Tables | `*_TABLE`: the files of `config/`. Point a setting at a copy to change a table |
| Method | `SEED`, `JITTER_SEED`, `NONDETECT_FRACTION`, `ACCESSION_VARIANCE_FLOOR`, `CLIP_LIMIT`, `N_CLUSTERS`, `KMEANS_N_INIT`, `K_MIN`, `K_MAX`, `CLUSTER_REMAP`, `CLUSTER_SIZES_EXPECTED`, `SIGNIFICANCE_LEVELS`, `TUKEY_ALPHA`, `N_RESAMPLE`, `N_SEEDS`, `N_MONTE_CARLO` |
| Figures | `FONT_FAMILY`, `STRUCTURE_LAYOUT_FILE`, `STRUCTURE_*`, `TABLE_STRUCTURE_PX`, `DPI_*`, `RADAR_LEGEND` |

The shipped values of the method settings are those of the paper. `CLUSTER_SIZES_EXPECTED`
holds the cluster sizes of the paper (90, 53, 22) and stops step 08 when a run gives other
sizes; leave it empty for other data.

`LABELS = revised` (default) writes the figures and table headings as in the revised paper.
`LABELS = submitted` writes those of the first submission. The numbers are the same in
both; `docs/NOTES.md` lists what differs. A figure has the same canvas and the same position
of its axes in both states: the layout is computed with the submitted texts, and the revised
texts are put in when the figure is saved. With `LABELS = revised` the panels of the paper
are saved at 600 dpi or more (`MINIMUM_PANEL_DPI` in `config/labels.tsv`), and compound names
are shown as listed in `config/compound_display_names.tsv` on every figure and in every
table except those a later step reads back.

In file settings, `CODE_DIR` at the start of a value stands for the folder of the pipeline
and `WORK_DIR` for the folder of the run. Other file names are taken from `INPUT_DIR`.

## Outputs

Everything is written under `WORK_DIR`; the folder of the code is never written to.

```
WORK_DIR/results/<step>/    tables and <id>_qc.txt
WORK_DIR/figures/<step>/    figures
WORK_DIR/logs/              step_<id>.log, run.tsv
```

`docs/OUTPUTS.md` lists the files of each step. The ones behind the figures and tables of
the paper:

| Item | File |
|---|---|
| Fig. 1 | `figures/01_report_compounds/compound_report/all_compounds_structures_A4_high_quality.png` |
| Table S6 | `results/01_report_compounds/compound_report/complete_compound_table_with_images.xlsx` |
| Fig. 2a, 2b | `figures/09_summarise_classes/<set>/class_summaries_<label>/combined_subclass_mainclass_<label>.png` |
| Fig. 2c | `figures/10_plot_detection/detection/3_<dataset>_detailed_compound_view.png` |
| Fig. 2d | `figures/09_summarise_classes/clusters_all/class_summaries_blup_clusters/8_blup_clusters_accession_all_compounds_profiles.png` |
| Fig. 3a, Fig. S2 | `figures/06_fit_blup_clusters/blup_<dataset>/cluster_figure_<dataset>.png` (dendrogram; PCA and t-SNE) |
| Fig. 3b | `figures/14_plot_group_distributions/clusters_log1p_fig3b/plots_<dataset>/compound_correlation_heatmap.png` |
| Fig. 3c | `figures/09_summarise_classes/clusters_markers/class_summaries_blup_clusters/9_blup_clusters_accession_all_compounds_percentages.png` |
| Legends of Fig. 3 | `legend_clusters_<dataset>.png` (step 06), `legend_main_class.png`, `legend_subclass.png` (step 14, Fig. 3b set), `legend_compounds_blup_clusters.png` (step 09, marker set) |
| Fig. 4a | `figures/16_compute_odor_index/plot1_concentration_vs_odor_index.png` |
| Fig. 4b | `figures/16_compute_odor_index/plot3_category_distribution.png` |
| Fig. 4c | `figures/17_plot_odor_radar/radar/radar_v4_pareto_accessions_SAMPLES.png` |
| Fig. S5 | `figures/16_compute_odor_index/plot2_odor_index_ranking.png` |
| BLUPs, variance components, repeatability, clusters | `results/06_fit_blup_clusters/blup_<dataset>/blup_results_<dataset>.xlsx` |
| Criteria for the number of clusters | `results/07_evaluate_cluster_number/cluster_number_<dataset>/cluster_criteria_metrics_<dataset>.xlsx` |
| Relative odor index per compound | `results/16_compute_odor_index/odor_index_compounds.xlsx` |
| Tables S8, S18, S9, S10, S17 | `results/18_...`, `results/19_...`, `results/20_...`, `results/21_...` |

The figures of the paper are assembled from these panels. The legends of Fig. 3 stand
between the panels there, so each is also written as an image of its own.

Each QC file starts with the pipeline version and states the inputs with their sizes, the
treatment of values without detection, the seeds, all settings and the package versions.

## Tested example

`tests/data/` holds an invented dataset of the same shape as the real one (48 accessions,
61 compounds, 135 runs, one to five runs per accession, values without detection). One
command runs all steps on it and compares the tables with `tests/expected/`:

```bash
bash tests/run_test.sh /path/to/test_work
```

The run takes about 15 minutes. Its last lines:

```
  step   status     seconds  started
  01     ok              44  2026-11-17 01:00:57
  ...
  21     ok             213  2026-11-17 01:15:18
  table      /path/to/test_work/logs/run.tsv
tables compared: 92, equal to 6 significant digits: 92, different: 0, tables present on one side only: 0
machine-dependent tables (21_run_nested_cv), not part of the result: 5 equal, 0 different
RESULT: all tables are equal
```

What is compared, and how:

| Tables | Comparison |
|---|---|
| steps 01 to 20 (92 tables) | every cell: numbers to six significant digits, text exactly. A difference fails the test |
| step 21 (5 tables: `S17_confusion_matrix.csv`, `S17_nested_cv_by_cluster.csv`, `S17_nested_cv_folds.csv`, `S17_nested_cv_summary.csv`, `S17_selection_frequency.csv`) | compared the same way, reported in a block of their own, and not part of the result |
| figures, QC files | not compared |

Step 21 gives the same tables in repeated runs on one machine, but not on every machine. In
each fold it selects compounds from rankings (mutual information, random forest importance)
and tunes a classifier; where two candidates are nearly tied, differences in the last digits
between processors and numerical libraries decide, and a different selection changes the
accuracies that follow. On the invented data, with few accessions per fold, such ties are
common. A line "differs on this machine" for a table of step 21 is therefore expected on
another machine and does not mean an error.

The run of the paper: 165 accessions, 61 compounds, 490 runs, four quantification datasets,
one CPU.

| Step | Seconds | Peak memory (GB) |
|---|---|---|
| 01 | 40 | 1.9 |
| 02 | 18 | 0.2 |
| 03 | 4 | 0.2 |
| 04 | 3 | 0.2 |
| 05 | 3 | 0.2 |
| 06 | 65 | 1.2 |
| 07 | 14 | 0.3 |
| 08 | 2 | 0.2 |
| 09 | 510 | 1.9 |
| 10 | 20 | 0.9 |
| 11 | 44 | 0.9 |
| 12 | 44 | 0.9 |
| 13 | 24 | 0.3 |
| 14 | 1483 | 8.8 |
| 15 | 6 | 0.2 |
| 16 | 33 | 1.3 |
| 17 | 31 | 0.7 |
| 18 | 40 | 0.2 |
| 19 | 29 | 0.2 |
| 20 | 5 | 0.2 |
| 21 | 269 | 0.2 |

All steps: about 44 minutes. Steps 18 to 21 were run with `OMP_NUM_THREADS=1`.

## Running it on other data

The tables in `config/` hold the compounds, classes, colours, thresholds and markers of the
paper, so that the shipped configuration reproduces it. They are inputs read through
`config.txt`, not values written in the scripts. For another compound set, copy the tables
that change, edit the copies, and point the settings at them; the scripts stay as they are.

| What changes | Table | Setting |
|---|---|---|
| compounds, their main class, subclass and structure | `compounds.tsv` (columns read: `Compound`, `Name`, `Main_Class`, `Subclass`, `SMILES`; the PubChem columns are descriptive) | `COMPOUND_FILE` |
| colour of each compound | `compound_colours.tsv` | `COMPOUND_COLOUR_TABLE` |
| colours of main classes and K groups | `colours.tsv` | `COLOUR_TABLE` |
| order of the subclasses | `subclass_order.tsv` | `SUBCLASS_ORDER_TABLE` |
| names, colours and descriptors of the clusters | `clusters.tsv` | `CLUSTER_TABLE` |
| odor thresholds, descriptors and their sources | `odor_thresholds.tsv` | `ODOR_THRESHOLD_TABLE` |
| marker panel evaluated in steps 20 and 21 | `markers.tsv` | `MARKER_TABLE` |
| names shown on figures and in tables where they differ from the data | `compound_display_names.tsv` | `DISPLAY_NAMES_TABLE` |
| 2D layout of the structure drawings | `structure_layouts.sdf` | `STRUCTURE_LAYOUT_FILE` |

`structure_layouts.sdf` holds one record per compound: the molecule with 2D coordinates,
explicit hydrogens on stereocentres and wedge bonds, and the compound name as record title.
Step 01 draws a compound from its record, so the drawings are the same with every RDKit
version. A compound without a record is laid out from its SMILES by the installed RDKit, so
for another compound set the file can be left as it is, or the setting can be left empty. To
fix the layouts of other compounds, write them with RDKit (`Chem.SDWriter`, after
`rdDepictor.Compute2DCoords` or your own coordinates) under the names of the compound list.

Compound names must be written the same way in the measurement files and in these tables.
With other data, also set `CLUSTER_SIZES_EXPECTED` empty, set `DATASETS` and the sheet names to
the datasets you have, and choose `N_CLUSTERS` from step 07. The limits that remain are listed
under "Scope and limitations".

## Scope and limitations

- The pipeline reproduces one analysis. It expects the layout of `docs/INPUTS.md`, four
  quantification datasets of the same runs, and compound names equal to those of the tables
  in `config/`.
- Colours of main classes, clusters and K groups, label sizes of the heat maps, resolutions
  and the wording of the two label states are settings. Other sizes, line widths and
  positions of the figures are written in the plotting code.
- Three clusters are named in `config/clusters.tsv`. A run with another `N_CLUSTERS` needs
  further rows there, and the figures of steps 13, 16 and 17 are laid out for three.
- Concentrations are nonyl acetate equivalents (one response factor for all compounds). The
  relative odor index is a screening
  quantity: thresholds come from the literature and were measured in other matrices.
- The relative odor index leaves out runs and accessions without detection
  (`docs/NOTES.md`, step 16).
- Step 14 needs about 9 GB of memory at 300 dpi.
- The feature selection itself (CannaFeatML) is outside this repository; steps 20 and 21
  only evaluate the panel.

## Common issues

| Message or symptom | Cause and remedy |
|---|---|
| `unknown key(s) in ...` | a key of your config is not in `config.txt` of this version; compare the two files |
| `key(s) missing from ...` | your config was copied from another version; add the keys named |
| `REPLICATE_FILE: file not found` | `INPUT_DIR` or the file name is wrong; file names are relative to `INPUT_DIR` |
| `sheet '...' not found` | set `SHEET_<DATASET>` to the sheet names of your workbooks |
| `accession names differ between ...` | the two files do not hold the same accessions; the message lists the names |
| `Cluster sizes differ from expected` | other data or other package versions than the paper; empty `CLUSTER_SIZES_EXPECTED` to go on |
| `BLUP workbook not found ... (run step 06 first)` | step 08 was started alone in a new `WORK_DIR` |
| step 14 is killed without a Python error | not enough memory; give the job 12 GB or lower `DPI_GROUP_DISTRIBUTIONS` |
| step 21 finishes in a few seconds and prints that folds were read from `nested_cv_folds/` | folds of an earlier run were reused; delete `results/21_run_nested_cv/nested_cv_folds/` to compute them again |
| the test lists a table of step 21 under "differs on this machine" | expected on a machine other than the one that wrote `tests/expected`; see Tested example |
| steps 18 to 21 take hours | several threads per process; set `OMP_NUM_THREADS=1` |
| names in the A4 structure grid look different | `STRUCTURE_GRID_FONT_FILE` was not found and the default font was used; the QC file of the step names the font |
| the points of the box plots sit elsewhere than in another run | `JITTER_SEED` is empty or differs |

## Licence

MIT. See `LICENSE`.

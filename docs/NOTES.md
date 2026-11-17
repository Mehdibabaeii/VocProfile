# Notes on the implementation

What each step does where the choice is not obvious from the output, and what to keep in
mind when the tables and figures are read or the code is reused. Numbers quoted here were
obtained on the dataset of the paper (165 accessions, 61 compounds, 490 runs).

## What every step writes

Each step writes `results/<step>/<id>_qc.txt`. Its first line is the pipeline version. It
lists the input files with their numbers of runs, accessions and compounds, how values
without detection were treated and at which level, the notes of the step, the seeds, every
setting of the run, and the versions of Python and of the packages.

Before it computes, a step checks its inputs: the sheets and columns named in the config are
present, measured columns are numeric, no accession name occurs twice in a file of accession
means, and the accession names are equal across the sheets of a workbook and across the
files a step combines. A mismatch stops the step and prints the names that differ.

Figures (`.png`, `.html`) go to `figures/<step>/`, tables to `results/<step>/`.

## Label states

`LABELS = revised` gives the wording, sizes and resolutions of the revised figures and
tables. `LABELS = submitted` gives those of the first submission. No number differs between
the two. The settings that differ are the rows of `config/labels.tsv`, the columns
`*_submitted` and `*_revised` of `config/input_sets.tsv` and `config/figure_styles.tsv`, and
the word pairs of `config/label_wording.tsv`, which are applied only when `LABELS = revised`.

A figure has the same canvas and the same position of its axes in both states. The scripts
draw every figure with the submitted texts; when the figure is saved with `LABELS = revised`,
its canvas is fixed first and the texts are replaced afterwards (terms, cluster names, legend
titles, compound names). Font sizes are the same in both states. Two things change the
layout on purpose, in the revised state only: the cluster legend of the accession profiles
of step 09 is left out because it covers bars, and the cluster legend of the clustered heat
maps of step 16 is moved (see step 16).

With `LABELS = revised` the panels of the paper drawn by steps 06, 16 and 17 and the
correlation matrix of the Fig. 3b input set are saved at `MINIMUM_PANEL_DPI` (600) or more;
they are rendered at that resolution, not enlarged. Steps 09 and 10 are at 600 dpi in both
states, and the structure grid of step 01 has four times the pixels of the submitted one.

Three things are the same in both states:

- the jitter of the points in the box plots of step 14 is seeded with `JITTER_SEED`
  (see Seeds);
- the colour of cluster ONE is the one of `config/clusters.tsv` in every step;
- the descriptors of the clusters shown in legends are those of `config/clusters.tsv`.

"Germacrene β" is the name of the compound in every input file and table of `config/`. It is
never renamed in the data. With `LABELS = revised` the name of `config/compound_display_names.tsv`
("Germacrene B") is shown on every figure, in the tables of steps 02, 03, 04 and 15 to 21 and
in the QC files. Tables that a later step reads back keep the name of the data, so that
names match between files: the filled and log-transformed means of step 02, the outputs of
steps 05 to 08, and `blup.csv` and `ref_labels.csv` of step 18.

## Seeds

`SEED` (42) is used by every k-means, t-SNE and random forest call, by the resampling of
steps 18 to 21 and by the numpy generator of steps 16 and 17. The same seed gives the same
tables on the same package versions.

`JITTER_SEED` seeds the sideways jitter of the points in the box plots of step 14 before
each figure is drawn. The jitter moves points horizontally and nothing else: values, boxes,
violins and letters do not depend on it. With an empty `JITTER_SEED` the jitter is unseeded
and two runs give figures whose points sit at different horizontal positions.
`LABELS = submitted` does not change the seed.

## Values without detection

A run without detection is `ND` or an empty cell in the replicate workbook; an accession
without any detection of a compound has an empty accession mean. Two rules are used, at two
levels, and the QC file of each step says which one applied.

- Replicate level (steps 06, 07, 18 to 21): a run without detection is set to
  `NONDETECT_FRACTION` (0.5) times the smallest positive value of the compound over all
  runs. This is the matrix from which variance components, BLUPs and clusters are computed.
- Accession-mean level (steps 02, 05, 12, 14, 15): an empty accession mean is set to 0.5
  times the smallest accession mean of the compound. An accession with a detection in part
  of its runs keeps the mean of the runs with a detection.

For 30 of the 42 compounds with values without detection the two fill values differ (ratio
up to 1.30, median 1.00). The correlation matrix and the heat map of step 14 are therefore
not computed from the matrix that was clustered in step 06.

Descriptive statistics on the data as they are (step 02), detection counts (step 10), the
distributions of step 11 and the relative odor index (steps 16, 17) use the values present
and leave the others out.

## Steps

### Step 01

- Compounds without SMILES are listed in the table and are not drawn. In the compound list
  of the paper these are two of 61, so the structure figure shows 59 structures.
- Every structure is drawn from its record in `STRUCTURE_LAYOUT_FILE` when the compound has one:
  stored 2D coordinates, explicit hydrogens on stereocentres and wedges. The drawings then do
  not depend on the RDKit version. Only the side on which the second line of a double bond is
  drawn can differ between versions. A compound without a record is laid out from its SMILES.
  The QC file says how many structures came from the file. Molecular weight, logP and the
  number of rotatable bonds are computed from the SMILES.
- The structures of the table (`compound_structures/`, workbook, HTML) are drawn at
  `TABLE_STRUCTURE_PX` (150) in both label states.
- Fig. 1 is the A4 grid, `all_compounds_structures_A4_high_quality.png`: four columns on
  2480 x 3508 px. With `LABELS = revised` the grid has four times as many pixels in each
  direction (`STRUCTURE_GRID_SCALE`), every structure is drawn at that size with line width and
  atom labels scaled, and the names are 15 % larger (`STRUCTURE_GRID_NAME_FACTOR`). At 20 % the
  two longest neighbouring names come within 5 px of each other at the submitted scale.
- The grid asks for the font file `STRUCTURE_GRID_FONT_FILE` and uses the default PIL font when
  the file is not found; the QC file says so.
- `all_compounds_structures.png` shows the same structures grouped by class. Its names have the
  same size in both states; with `LABELS = revised` it is saved at 600 dpi with structures drawn
  at 600 px. A name longer than 20 characters is split at hyphens and the hyphen at the break
  is not drawn. The stacked bar of the class counts is saved at the default resolution of
  matplotlib.

### Step 02

- Two reports per dataset: on the accession means as they are, and after filling at
  accession-mean level followed by the natural logarithm for the exported transformed table.
- In the report on the data as they are, an empty accession mean enters the geometric mean
  as 1e-9. For a compound with empty means the column `GeometricMean` is far below every
  measured value and should not be used; the report after filling is not affected.
- Statistics of a main class or subclass are computed over all accession-by-compound cells
  of the class pooled together (`Count` is the number of cells), not over sums per accession.
- Standard deviations and coefficients of variation use ddof = 1.

### Step 03

Computes accession means from the replicate workbook for each quantification dataset (runs
with a detection), the coefficient of variation within and among accessions per compound, Shapiro-Wilk tests, and
the Spearman correlation between each pair of datasets over the accessions with a value in both.
A pair with fewer than 10 such accessions gives no value.

### Step 04

Statistics of each compound over all runs with a detection, per dataset: mean, standard
deviation (ddof = 1), coefficient of variation, median, quartiles.

### Step 05

Writes the accession means of every dataset after filling at accession-mean level.

### Step 06

- The matrix of runs is filled at replicate level, then each compound is analysed alone.
- The error variance is the mean of the within-accession variances (ddof = 1) of the
  accessions with more than one run. The accession variance is the variance of the accession
  means (ddof = 1) minus the error variance divided by the mean number of runs per accession
  (2.97 in the paper), not below `ACCESSION_VARIANCE_FLOOR` times the error variance. This is
  the method-of-moments estimator. The floor is never reached in the main dataset.
- The BLUPs are the solution of the mixed-model equations with these two variances. The
  repeatability is that of a single run, not of an accession mean.
- Accessions enter the BLUP matrix in sorted order of their names. The QC file says so.
- Clustering: z-score per compound (ddof = 0), clipped at plus and minus `CLIP_LIMIT`,
  k-means with `N_CLUSTERS`, `SEED` and `KMEANS_N_INIT`. Ward linkage on the same matrix is
  computed as well, and the step keeps the method with the higher silhouette. k-means has
  the higher silhouette in all four datasets, so the clusters are k-means clusters. The
  dendrogram is drawn with Ward linkage.
- Clusters are numbered by size, largest first, unless `CLUSTER_REMAP` says otherwise.
- The ANOVA across accessions and the one-against-rest tests of each cluster use the runs as
  observations and are not adjusted for the number of compounds.

### Step 07

Calinski-Harabasz index, within-cluster sum of squares (elbow) and silhouette for
k = `K_MIN` to `K_MAX`, per dataset, on the z-scored BLUP matrix without clipping. At each k
both k-means and Ward are run and the one with the higher silhouette is reported
(`Best_Method`). The elbow is the k that follows the largest drop of the
within-cluster sum of squares.
The figure marks the elbow with a dashed line; its value is in the summary table and in the log.

### Step 08

Names the clusters of step 06 in the order of `config/clusters.tsv`, largest first, and
writes the two cluster files the later steps read: the accession means of the main dataset
and the BLUP estimates, each with the cluster name in the second column. The step stops if
the cluster sizes differ from `CLUSTER_SIZES_EXPECTED`. The rows of the BLUP file are in
sorted accession order.

### Step 09

Group profiles by main class, subclass and compound, once per input set of
`config/input_sets.tsv`. Group means use the accessions with a value. In the stacked bars a
compound without any value in a group counts as zero.

In the cluster input sets the accession names are coloured by cluster. With
`LABELS = submitted` a cluster legend is drawn inside the accession profiles, over the bars;
with `LABELS = revised` it is left out, and the legend written by step 06
(`legend_clusters_<dataset>.png`) is used in its place. The compound legend at the right of
plot 9 is also written as `legend_compounds_<label>.png`.

### Step 10

Number of accessions in which each compound was detected, from the accession means as they
are.

### Steps 11 and 12

Distributions of the compounds by class: step 11 on the accession means as they are, step 12
after filling at accession-mean level.

### Step 13

Three layouts in one step, chosen by the column `layout` of the input set: K groups on the
accession means of every dataset, clusters on the BLUP estimates for all compounds, clusters
for the markers of `config/markers.tsv`. The layouts differ in group order, palette and
figure sizes. Empty group means count as zero in the stacks.

### Step 14

- One run per input set. The column `transform` decides the treatment: `log1p` (filled at
  accession-mean level, then ln(1 + x)), `filled`, or `none`.
- Box plots: Tukey test at `TUKEY_ALPHA` on the groups; the letters are assigned by a rule
  that is correct for three groups and can miss shared letters with more.
- Heat map: z-score per compound (ddof = 0), clipped at plus and minus `CLIP_LIMIT`, Ward
  linkage on rows and columns. With `transform = none` empty means are set to zero before the
  z-scores, so that heat map is a zero-filled one. The PCA uses the accessions without empty
  means and is not written when none is left.
- Correlation matrix: Pearson correlation between compounds over the accession means after
  the treatment of the input set. The matrix of `clusters_log1p_fig3b` is Fig. 3b. With
  `LABELS = revised` it is saved at 600 dpi on the canvas of the submitted figure, its legend
  has the titles "Main class (b)" and "Subclass (b)" without prefixes in the entries, and the
  two legends are also written as `legend_main_class.png` and `legend_subclass.png`.
- Compound names of Fig. 3b: the setting `compound_label_fit_max` of `config/figure_styles.tsv`
  asks for the largest size up to 40 pt at which every name stays on the submitted canvas,
  clear of its neighbours, of the axis title and of the colour bar. On the data of the paper
  no larger size passes: at the submitted 30 pt the longest row name already reaches the
  colour bar. The names therefore stay at 30 pt, and the QC file says so. Larger names need the
  colour bar further right (`colourbar_left`), which widens the canvas.
- In the heat map of `clusters_log1p_fig3b` the accession names are drawn larger than the
  rows are high and overlap.
- Value axis of the box plots: `LABELS = submitted` writes "log(...)" for all three
  treatments, as in the first submission, although only `log1p` transforms the values.
  `LABELS = revised` writes "ln(1 + ...)" for `log1p` and no logarithm for the other two.
- The step needs about 9 GB of memory at 300 dpi because of the size of the heat maps.

### Step 15

Tables behind the heat map and the correlation matrix of step 14 for the three treatments
(as they are, filled, filled and ln(1 + x)). In the version without filling, the sheet of
scaled values is computed after empty means are set to zero, as in the heat map.

### Step 16

- Clustered heat maps (`plot_heatmap_PARETO.png`, `plot_heatmap_SAMPLE_zscore.png`): with
  `LABELS = submitted` the cluster legend sits in the upper left corner over the colour bar and
  the strip is named "Class" below its lower end. With `LABELS = revised` the legend is one row
  above the column dendrogram and the strip is named "Cluster" above its upper end. The heat map
  keeps its size; the canvas is 0.8 inch higher.
- Accession value of a compound: mean of the runs with a detection, divided by the odor
  threshold of `config/odor_thresholds.tsv`. Cluster mean: mean over the accessions with a
  value. Value over all accessions: the cluster means weighted by the full cluster sizes.
- Runs and accessions without detection are left out at both levels, so the index describes
  the accessions in which the compound was found. For a compound found in part of the
  accessions the value depends on this rule. Prenylthiol, detected in 106 of 165 accessions:
  284,159 as computed; 290,593 as the simple mean of the accessions with a value; 186,684
  with accessions without detection set to zero; 151,992 with runs without detection set to
  zero; 209,135 with the replicate-level rule of step 06. It ranks first under every rule.
  For compounds detected everywhere the rules give the same value. Step 19 reports category
  shares under three treatments.
- In the figure of concentration against index, the concentration is the unweighted mean of
  the three cluster means and the index is weighted by cluster size.
- Pareto-scaled values set empty accession values to zero, that is to the mean of the
  accessions with a value, before medians and the clustered heat map are computed.
- z-scores over the three cluster means use ddof = 0; z-scores over accessions use ddof = 1.
- Cluster names are attached to two tables by position. This is correct while the cluster
  file and the index table are both in sorted accession order, as written by step 08.
- The correction of the one-against-rest tests counts tests without a p-value in the number
  of tests; all tests have one in the paper.

### Step 17

Radar figures in several scalings from the same index values as step 16. The statistics
sheets report n as the cluster size, not the number of accessions with a value. In step 17
the column `CV_within_sample` is empty for a compound that lacks a value in some accession;
step 16 computes it over the accessions with a value.

### Steps 18 to 21

- All four start from the replicate workbook, fill at replicate level and compute the BLUP
  matrix as step 06 does, and compare with the cluster of each accession in
  `MEANS_CLUSTER_FILE`.
- Step 18 writes `blup.csv` and `ref_labels.csv`; step 21 reads them, so step 18 runs first.
- Step 21 keeps each outer fold in `nested_cv_folds/` and continues from the folds present. It
  prints how many folds were reused and writes the number to its QC file. Delete that folder
  to compute all folds again.
- The tables of step 21 are the same in repeated runs on one machine. Between machines they
  can differ: the selection of compounds inside each fold rests on rankings (mutual
  information, random forest importance) in which near ties are broken by differences in the
  last digits, which depend on the processor and the numerical libraries. The test therefore
  reports these tables apart (see Comparing two runs).
- These steps run k-means many times. With several threads per process they can be much
  slower; set `OMP_NUM_THREADS=1`.

## Not written by any step

The total volatile output per cluster is the median, per cluster, of the row sums of the
BLUP estimates over all compounds in `blup_clusters.csv` (step 08). No step writes it. On the
dataset of the paper the medians are 462, 858 and 1,292 µg/g; the same medians from the sums
of the accession means are 438, 855 and 1,302.

## Fonts

Steps 01, 16 and 17 set `FONT_FAMILY`. When the family is not installed the step uses
DejaVu Sans and writes both names to its QC file. The other steps use the default font of
matplotlib. Text extents depend on the font, so figures made with different fonts differ in
layout.

## Comparing two runs

`tools/compare_outputs.py` compares the tables of two run folders cell by cell: numbers to
six significant digits, text exactly. Figures are not compared. Tables under a folder given
with `--machine-dependent` are compared and listed in a block of their own, and a difference
there does not change the result. `tests/run_test.sh` uses this for step 21.
`tests/output_name_map.tsv` lists, for outputs that carried other names in earlier versions
of the scripts, the earlier and the present relative path; give it with `--map` when the
first folder uses the earlier names.

## Test data

`tests/data/` holds an invented dataset written by `tools/make_synthetic_data.py`: 48
accessions with invented names, one to five runs each, values without detection, compounds
found in part of the accessions, and the four sheets of the real layout. Only the compound
names and classes are real. `tests/expected/` holds the tables of a run on it.

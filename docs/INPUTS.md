# Input files

The pipeline reads three files of measurements and writes two more itself. Their names,
sheets and columns are settings of `config.txt`; nothing about them is fixed in the code.
The examples below are invented.

## Compound list (`COMPOUND_FILE`)

A tab-separated table (or a workbook with the same columns), one row per compound. The list
of the paper is shipped as `config/compounds.tsv`.

| Column | Content |
|---|---|
| `Compound` | name as shown on figures. It is the column heading of the measurement files, except for the names listed in `config/compound_display_names.tsv` |
| `SMILES` | structure; empty when no structure is assigned (the compound is listed, not drawn) |
| `Main_Class` | main biosynthetic class |
| `Subclass` | subclass |

## Replicate workbook (`REPLICATE_FILE`)

One sheet per quantification dataset (`SHEET_<DATASET>`). One row per run, one header row.

| Sample | Rep | Methyl butyrate | Hexanal | ... |
|---|---|---|---|---|
| ACC-001 | 1 | 0.512 | ND | ... |
| ACC-001 | 2 | 0.498 | 0.031 | ... |
| ACC-002 | 1 | 1.204 | 0.117 | ... |

- `ACCESSION_COLUMN` (here `Sample`) holds the accession name, `REPLICATE_COLUMN` (here
  `Rep`) the number of the run. All other columns are compounds.
- A run without detection is `NONDETECT_CODE` (here `ND`) or an empty cell.
- The number of runs per accession may differ (one to five in the paper).
- The accession names must be the same in every sheet.

## Accession means (`MEANS_FILE`)

One sheet per dataset, same sheet names as the replicate workbook. `MEANS_HEADER_ROWS` (3)
header rows: main class, subclass, compound. The first column holds the accession name, the
second the group of the accession (in the paper the genetic group K1, K2 or K3). An empty
cell is an accession without any detection of the compound.

|  | Main Class | Fatty Acid-Derived | Fatty Acid-Derived | ... |
|---|---|---|---|---|
|  | Subclass | Acyl-CoA ester | Lipoxygenase-derived | ... |
|  | K | Methyl butyrate | Hexanal | ... |
| ACC-001 | K1 | 0.505 | 0.031 | ... |
| ACC-002 | K3 | 1.204 | 0.117 | ... |

Each accession occurs once. The accession names must be those of the replicate workbook.

## Cluster files (`MEANS_CLUSTER_FILE`, `BLUP_CLUSTER_FILE`)

Written by step 08 under `WORK_DIR/results/08_build_cluster_inputs/` and read from there by
the later steps with the shipped settings.

- `accession_means_clusters.xlsx`: the sheet of `MAIN_DATASET` of the accession means, with
  the cluster name (ZERO, ONE, TWO) in the second column.
- `blup_clusters.csv`: the BLUP estimates of step 06 in the same layout (three header rows,
  accession, cluster, one column per compound), rows in sorted accession order.

To use cluster files made elsewhere, point the two settings at them. They must have this
layout and the accession names of the other files.

## Tables of `config/`

| File | Setting | Content |
|---|---|---|
| `compounds.tsv` | `COMPOUND_FILE` | the 61 compounds with class, subclass and SMILES |
| `structure_layouts.sdf` | `STRUCTURE_LAYOUT_FILE` | 2D layout of each structure drawing (record title = compound name) |
| `compound_colours.tsv` | `COMPOUND_COLOUR_TABLE` | colour of each compound in the detection figure |
| `compound_display_names.tsv` | `DISPLAY_NAMES_TABLE` | name shown where it differs from the name in the data |
| `colours.tsv` | `COLOUR_TABLE` | colours of the main classes and of the K groups |
| `clusters.tsv` | `CLUSTER_TABLE` | number, name, colour and descriptor of each cluster |
| `subclass_order.tsv` | `SUBCLASS_ORDER_TABLE` | order of the subclasses in the figures |
| `markers.tsv` | `MARKER_TABLE` | the 12 marker compounds |
| `odor_thresholds.tsv` | `ODOR_THRESHOLD_TABLE` | odor threshold, odor note and category of 33 compounds |
| `input_sets.tsv` | `INPUT_SETS_TABLE` | the runs of steps 09, 13 and 14 |
| `figure_styles.tsv` | `FIGURE_STYLES_TABLE` | label sizes and strip labels of the heat maps and correlation matrices |
| `labels.tsv` | `LABELS_TABLE` | settings that differ between the two label states |
| `label_wording.tsv` | `LABEL_WORDING_TABLE` | word pairs applied when `LABELS = revised` |

The 12 markers were selected with
[CannaFeatML](https://github.com/Mehdibabaeii/CannaFeatML), a separate pipeline that is not
part of this repository. The README section "Feature selection with CannaFeatML" describes
how.

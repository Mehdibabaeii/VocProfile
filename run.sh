#!/usr/bin/env bash
#
# VocProfile - run.sh
#
# Runs the steps in order, each with its own log, and closes with a table of what ran.
#
# Usage:
#   bash run.sh <config.txt>              all steps
#   bash run.sh <config.txt> 03           one step
#   bash run.sh <config.txt> from 05      from a step to the end
#   bash run.sh <config.txt> --dry-run    print the plan, run nothing
#   bash run.sh <config.txt> --set KEY=VALUE [--set ...]   change a setting for this run
#   bash run.sh --list                    list the steps
#
set -u
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="$(tr -d '[:space:]' < "$CODE_DIR/VERSION")"
NAME="VocProfile"
TAGLINE="volatile profiling and odor screening of cannabis accessions"
PYTHON="${PYTHON:-python3}"
export PYTHONDONTWRITEBYTECODE=1

STEPS=(
"01|01_report_compounds|compound table and structure grid (Fig. 1, Table S6)"
"02|02_describe_accession_means|descriptive statistics of the accession means"
"03|03_compare_datasets|agreement among the quantification datasets"
"04|04_summarise_replicates|variation among replicate runs"
"05|05_fill_nondetects|accession means with values without detection filled"
"06|06_fit_blup_clusters|variance components, BLUP, repeatability, clustering (Fig. 3a)"
"07|07_evaluate_cluster_number|Calinski-Harabasz, elbow and silhouette for k = 2 to 7"
"08|08_build_cluster_inputs|cluster files for the later steps"
"09|09_summarise_classes|profiles by class, subclass and compound (Fig. 2, Fig. 3c)"
"10|10_plot_detection|number of accessions in which each compound was detected (Fig. 2c)"
"11|11_plot_distributions|distributions of the compounds"
"12|12_plot_filled_distributions|distributions after filling"
"13|13_plot_group_profiles|stacked profiles of K groups and clusters"
"14|14_plot_group_distributions|box plots, heat map, PCA, correlation matrix (Fig. 3b)"
"15|15_export_correlations|tables behind the correlation matrix"
"16|16_compute_odor_index|relative odor index (Fig. 4, Figs. S5 and S6)"
"17|17_plot_odor_radar|radar figures of the relative odor index"
"18|18_test_cluster_stability|stability of the clusters under resampling (Table S8)"
"19|19_test_odor_index|tests of the relative odor index (Table S18)"
"20|20_run_diagnostics|detection counts, residual tests, marker correlation"
"21|21_run_nested_cv|nested cross-validation of the markers (Table S17)"
)
TOTAL=${#STEPS[@]}

list_steps() {
    local entry
    for entry in "${STEPS[@]}"; do
        IFS='|' read -r id script what <<< "$entry"
        printf '  %s  %-32s %s\n' "$id" "$script.py" "$what"
    done
}

usage() { sed -n '7,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

if [ "${1:-}" = "--list" ]; then list_steps; exit 0; fi
if [ $# -lt 1 ] || [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then usage; exit 1; fi

CONFIG="$1"; shift
if [ ! -f "$CONFIG" ]; then echo "ERROR: config file not found: $CONFIG" >&2; exit 1; fi
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

FIRST="01"; ONLY=""; DRY=0; SETS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --list) list_steps; exit 0 ;;
        --dry-run) DRY=1 ;;
        --set) shift; [ $# -gt 0 ] || { echo "ERROR: --set needs KEY=VALUE" >&2; exit 1; }; SETS+=(--set "$1") ;;
        --set=*) SETS+=(--set "${1#--set=}") ;;
        from) shift; FIRST="${1:-}" ;;
        [0-9][0-9]) ONLY="$1" ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
    shift
done

known() { local entry; for entry in "${STEPS[@]}"; do [ "${entry%%|*}" = "$1" ] && return 0; done; return 1; }
for id in "$FIRST" ${ONLY:+"$ONLY"}; do
    known "$id" || { echo "ERROR: no step '$id'. bash run.sh --list shows the steps." >&2; exit 1; }
done

setting() {
    "$PYTHON" "$CODE_DIR/scripts/lib_config.py" "$CONFIG" ${SETS[@]+"${SETS[@]}"} --print "$1"
}

BANNER="$(setting BANNER)" || exit 1
WORK_DIR="$(setting WORK_DIR)" || exit 1

if [ "$BANNER" = "0" ]; then
    echo "$NAME $VERSION"
    echo "$TAGLINE"
else
    line() { printf '*  %-72s  *\n' "$1"; }
    printf '%.0s*' $(seq 78); echo
    line ""
    while IFS= read -r art; do line "  $art"; done <<'ART'
__     __         ____             __ _ _
\ \   / /__   ___|  _ \ _ __ ___  / _(_) | ___
 \ \ / / _ \ / __| |_) | '__/ _ \| |_| | |/ _ \
  \ V / (_) | (__|  __/| | | (_) |  _| | |  __/
   \_/ \___/ \___|_|   |_|  \___/|_| |_|_|\___|
ART
    line ""
    line "$TAGLINE"
    line "version $VERSION"
    line ""
    printf '%.0s*' $(seq 78); echo
fi
"$PYTHON" "$CODE_DIR/scripts/lib_config.py" "$CONFIG" ${SETS[@]+"${SETS[@]}"} || exit 1

SELECTED=()
for entry in "${STEPS[@]}"; do
    id="${entry%%|*}"
    if [ -n "$ONLY" ]; then [ "$id" = "$ONLY" ] && SELECTED+=("$entry")
    elif [ "$id" \> "$FIRST" ] || [ "$id" = "$FIRST" ]; then SELECTED+=("$entry"); fi
done

if [ "$DRY" = "1" ]; then
    echo "  steps of this run (nothing is run with --dry-run):"
    for entry in "${SELECTED[@]}"; do
        IFS='|' read -r id script what <<< "$entry"
        printf '  %s  %-32s %s\n' "$id" "$script.py" "$what"
    done
    exit 0
fi

LOG_DIR="$WORK_DIR/logs"
mkdir -p "$LOG_DIR" || exit 1
RUN_TABLE="$LOG_DIR/run.tsv"
ROWS=()
FAILED=0
RULE="$(printf '%.0s=' $(seq 78))"

for entry in "${SELECTED[@]}"; do
    IFS='|' read -r id script what <<< "$entry"
    if [ "$FAILED" = "1" ]; then ROWS+=("$id	not run	0	-"); continue; fi
    log="$LOG_DIR/step_$id.log"
    started="$(date '+%Y-%m-%d %H:%M:%S')"
    echo "$RULE"
    printf '  [%d/%d]  step %s    %s\n' "$((10#$id))" "$TOTAL" "$id" "$what"
    echo "  started    $started"
    echo "  script     scripts/$script.py"
    echo "  log        $log"
    echo "$RULE"
    t0=$(date +%s)
    "$PYTHON" "$CODE_DIR/scripts/$script.py" "$CONFIG" ${SETS[@]+"${SETS[@]}"} 2>&1 | tee "$log"
    status=${PIPESTATUS[0]}
    seconds=$(( $(date +%s) - t0 ))
    if [ "$status" = "0" ]; then ROWS+=("$id	ok	$seconds	$started")
    else ROWS+=("$id	failed	$seconds	$started"); FAILED=1; fi
done

{
    printf 'step\tstatus\tseconds\tstarted\n'
    printf '%s\n' "${ROWS[@]}"
} > "$RUN_TABLE"
echo "$RULE"
awk -F'\t' '{ printf "  %-6s %-9s %8s  %s\n", $1, $2, $3, $4 }' "$RUN_TABLE"
echo "  table      $RUN_TABLE"
if [ "$FAILED" = "1" ]; then
    echo "  A step failed. Its log is in $LOG_DIR. The steps after it were not run."
    exit 1
fi

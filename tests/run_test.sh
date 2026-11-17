#!/usr/bin/env bash
#
# VocProfile - run_test.sh
#
# Runs all steps on the invented dataset of tests/data and compares the tables with
# tests/expected. The tables of step 21 depend on the machine and are reported apart.
#
# Usage:
#   bash tests/run_test.sh <work folder>
#
set -u
if [ $# -ne 1 ]; then echo "usage: bash tests/run_test.sh <work folder>" >&2; exit 1; fi
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python3}"
mkdir -p "$1" || exit 1
WORK="$(cd "$1" && pwd)"
bash "$CODE_DIR/run.sh" "$CODE_DIR/tests/config_test.txt" --set WORK_DIR="$WORK" || exit 1
"$PYTHON" "$CODE_DIR/tools/compare_outputs.py" "$CODE_DIR/tests/expected" "$WORK/results" \
    --machine-dependent 21_run_nested_cv

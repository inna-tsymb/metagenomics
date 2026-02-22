#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT_DIR/.venv/bin/python}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "Python interpreter not found or not executable: $PYTHON_BIN"
  echo "Set PYTHON_BIN explicitly, e.g.:"
  echo "  PYTHON_BIN=/path/to/python ./run_all_dual_mode.sh"
  exit 1
fi

run_step() {
  local cwd="$1"
  local script="$2"
  shift 2
  echo
  echo "=== $(basename "$cwd") :: $script ==="
  (
    cd "$cwd"
    env "$@" "$PYTHON_BIN" "$script"
  )
}

run_lecture03_mode() {
  local mode="$1"
  local normalize="$2"
  local l3="$ROOT_DIR/lecture_03"

  echo
  echo "########################################################################"
  echo "Lecture 03 mode: $mode (APPLY_SAMPLE_SIZE_NORMALIZATION=$normalize)"
  echo "########################################################################"

  run_step "$l3" scripts/filtering_clr_analysis.py \
    APPLY_SAMPLE_SIZE_NORMALIZATION="$normalize" \
    OUTPUT_DIR="output/filtering_clr_analysis/$mode"

  run_step "$l3" scripts/association_analysis.py \
    APPLY_SAMPLE_SIZE_NORMALIZATION="$normalize" \
    ABUNDANCE_FILE="output/filtering_clr_analysis/$mode/abundance_clr_aligned.csv" \
    METADATA_FILE="output/filtering_clr_analysis/$mode/metadata_aligned.csv" \
    OUTPUT_DIR="output/association_analysis/$mode"

  run_step "$l3" scripts/binary_logistic_analysis.py \
    APPLY_SAMPLE_SIZE_NORMALIZATION="$normalize" \
    ABUNDANCE_FILE="output/filtering_clr_analysis/$mode/abundance_clr_aligned.csv" \
    METADATA_FILE="output/filtering_clr_analysis/$mode/metadata_aligned.csv" \
    LINEAR_RESULTS_FILE="output/association_analysis/$mode/association_results_all.csv" \
    OUTPUT_DIR="output/binary_logistic_analysis/$mode"

  run_step "$l3" scripts/pairwise_comparisons.py \
    APPLY_SAMPLE_SIZE_NORMALIZATION="$normalize" \
    INPUT_ABUNDANCE="output/filtering_clr_analysis/$mode/abundance_clr_aligned.csv" \
    INPUT_METADATA="output/filtering_clr_analysis/$mode/metadata_aligned.csv" \
    OUTPUT_DIR="output/pairwise_comparisons/$mode"

  run_step "$l3" scripts/association_statements.py \
    ABUNDANCE_FILE="output/filtering_clr_analysis/$mode/abundance_clr_aligned.csv" \
    METADATA_FILE="output/filtering_clr_analysis/$mode/metadata_aligned.csv" \
    RESULTS_FILE="output/association_analysis/$mode/association_results_all.csv" \
    OUTPUT_DIR="output/association_analysis/$mode"

  run_step "$l3" scripts/comparison_statements.py \
    OUTPUT_DIR="output/binary_logistic_analysis/$mode" \
    LINEAR_RESULTS_FILE="output/association_analysis/$mode/association_results_all.csv" \
    LOGISTIC_RESULTS_FILE="output/binary_logistic_analysis/$mode/logistic_regression_results_all.csv" \
    COMPARISON_FILE="output/binary_logistic_analysis/$mode/comparison_linear_vs_logistic.csv"

  run_step "$l3" scripts/sensitivity_balanced.py \
    OUTPUT_DIR="output/sensitivity_balanced/$mode" \
    CLR_FILE="output/filtering_clr_analysis/$mode/abundance_clr_aligned.csv" \
    METADATA_FILE="output/filtering_clr_analysis/$mode/metadata_aligned.csv" \
    ORIGINAL_RESULTS="output/association_analysis/$mode/association_results_all.csv"
}

run_lecture04_mode() {
  local mode="$1"
  local normalize="$2"
  local l4="$ROOT_DIR/lecture_04"

  echo
  echo "########################################################################"
  echo "Lecture 04 mode: $mode (APPLY_SAMPLE_SIZE_NORMALIZATION=$normalize)"
  echo "########################################################################"

  run_step "$l4" scripts/meta_analysis.py \
    APPLY_SAMPLE_SIZE_NORMALIZATION="$normalize" \
    OUTPUT_DIR="output/meta_analysis/$mode"

  run_step "$l4" scripts/sensitivity_meta_balanced.py \
    OUTPUT_DIR="output/sensitivity_meta_balanced/$mode" \
    CLR_FILE="../lecture_03/output/sensitivity_balanced/$mode/abundance_clr_balanced.csv" \
    METADATA_FILE="../lecture_03/output/sensitivity_balanced/$mode/metadata_balanced.csv" \
    ORIGINAL_META_RESULTS="output/meta_analysis/$mode/meta_analysis_results_all.csv"
}

echo "Using Python: $PYTHON_BIN"
echo "Workspace: $ROOT_DIR"

run_lecture03_mode normalized true
run_lecture03_mode non_normalized false

run_lecture04_mode normalized true
run_lecture04_mode non_normalized false

echo
echo "✅ Dual-mode pipeline completed for lecture_03 and lecture_04"
echo "Outputs are available in */output/*/normalized and */output/*/non_normalized"

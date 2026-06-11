#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}"
IGORCALL="$IGORBIN -w $OUTDIR -j 1"

$IGORCALL init
$IGORCALL config set model.source custom
$IGORCALL config set model.parms "$TESTINPUT/TRB_model_parms.txt"
$IGORCALL config set model.marginals "$TESTINPUT/TRB_uniform_model_marginals.txt"
$IGORCALL config set infer.iterations 4

# Copy reference alignment files
mkdir -p "$OUTDIR/aligns"
cp "$TESTREF"/aligns/* "$OUTDIR/aligns/"

# Run the inference with the demo parameters
$IGORCALL config set infer.likelihood_threshold 1e-35
$IGORCALL config set infer.probability_ratio_threshold 0.0001
$IGORCALL config set output.scenarios 10
$IGORCALL config set output.Pgen true
$IGORCALL config set output.coverage VJ_gene
$IGORCALL -b demo infer

# Run the inference with the default parameters
$IGORCALL config set infer.likelihood_threshold 1e-60
$IGORCALL config set infer.probability_ratio_threshold 1e-5
$IGORCALL config set output.scenarios 0
$IGORCALL config set output.Pgen false
$IGORCALL config set output.coverage ""
$IGORCALL -b default infer
# Evaluate sequences to generate outputs
$IGORCALL config set model.source last_inferred
$IGORCALL config set evaluate.likelihood_threshold 1e-35
$IGORCALL config set evaluate.probability_ratio_threshold 0.0001
$IGORCALL config set output.scenarios 10
$IGORCALL config set output.Pgen true
$IGORCALL config set output.coverage VJ_gene
$IGORCALL -b default evaluate

# ------------------------------------------------------------------
# 2️⃣ Test output file regression
# ------------------------------------------------------------------
LOGFILE="$OUTDIR/infer_regression.log"

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/assert_regression.sh"

resolve_sort_columns() {
    local filename="$1"
    case "$filename" in
        # Don't sort model files
        *_parms.txt)           echo "None" ;;
        *_marginals.txt)       echo "None" ;;
        initial_model.txt)     echo "None" ;;
        iteration_*.txt)       echo "None" ;;
        likelihoods.out)       echo "None" ;;
        # Sort most counters by seq_id
        best_scenarios_counts.csv) echo "all" ;;
        *counts.csv)           echo "col1" ;;
        sequence_mutation_frequency.csv) echo "col1" ;;
        scenarios_background_and_errors.csv) echo "col1" ;;
        # Sort err and coverage counter by iteration and gene
        *_genes_cov_and_err.csv) echo "col1,col2" ;;
        # Sort inference logs by seq_id
        inference_logs.txt)    echo "col1" ;;
        # Catch undefined patterns
        *)                     echo "Undefined" ;;   # fallback
    esac
}

for batch in "demo" "default"
do

    # Drop non reproducible seq processing order and time elapsed per sequence
    tmp="$(mktemp)"                                   # create a safe temp name
    cut -d';' -f 1,3-10  "$OUTDIR/${batch}_inference/inference_logs.txt" >"$tmp"
    mv "$tmp" "$OUTDIR/${batch}_inference/inference_logs.txt"

    assert_regression "$TESTREF/${batch}_inference" "$OUTDIR/${batch}_inference" "$LOGFILE"

    # Drop not always reproducible scenario ranks (possible likelihood ties)
    tmp="$(mktemp)"                                   # create a safe temp name
    cut -d';' -f 1,3-  "$OUTDIR/${batch}_output/best_scenarios_counts.csv" >"$tmp"
    mv "$tmp" "$OUTDIR/${batch}_output/best_scenarios_counts.csv"

    assert_regression "$TESTREF/${batch}_output" "$OUTDIR/${batch}_output" "$LOGFILE"

done
# The script exits with the same status that run_regression returned

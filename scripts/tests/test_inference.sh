#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}"
IGORCALL="$IGORBIN -set_wd $OUTDIR"

# Copy reference alignment files
if [[ ! -d "$OUTDIR/aligns" ]]; then
    cp -r "$TESTREF/aligns" "$OUTDIR"
fi

# Run the inference with the demo parameters
$IGORCALL -batch demo -set_custom_model "$TESTINPUT/TRB_model_parms.txt" "$TESTINPUT/TRB_uniform_model_marginals.txt" -infer --N_iter 4  --L_thresh 1e-35 --P_ratio_thresh 0.0001 -output --scenarios 10 --Pgen --coverage VJ_gene

# Run the inference with the default parameters
$IGORCALL -batch default -set_custom_model "$TESTINPUT/TRB_model_parms.txt" "$TESTINPUT/TRB_uniform_model_marginals.txt" -infer --N_iter 4 
# Evaluate sequences to generate outputs
$IGORCALL -batch default -load_last_inferred -evaluate --L_thresh 1e-35 --P_ratio_thresh 0.0001 -output --scenarios 10 --Pgen --coverage VJ_gene

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
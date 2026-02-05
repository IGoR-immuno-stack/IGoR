#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}" # Create a temp dir if none passed
IGORCALL="$IGORBIN -set_wd $OUTDIR"

###################################################
# Generate sequences with fixed seeds
###################################################

# Generate random sequences with a single thread for reproducibility
$IGORCALL -batch seed42 -threads 1 -set_custom_model "$TESTREF/demo_inference/final_parms.txt" "$TESTREF/demo_inference/final_marginals.txt" -generate 100 --seed 42 
$IGORCALL -batch seedRd -threads 1 -set_custom_model "$TESTREF/demo_inference/final_parms.txt" "$TESTREF/demo_inference/final_marginals.txt" -generate 100 --seed 8557 
# Check generation without error
$IGORCALL -batch noerr -threads 1 -set_custom_model "$TESTREF/demo_inference/final_parms.txt" "$TESTREF/demo_inference/final_marginals.txt" -generate 100 --seed 35863 --noerr

# ------------------------------------------------------------------
# 2️⃣ Test output file regression
# ------------------------------------------------------------------
LOGFILE="$OUTDIR/generate_regression.log"

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/assert_regression.sh"

resolve_sort_columns() {
    local filename="$1"
    case "$filename" in
        # Sorting not necessary with fixed seed and single thread
        generated_realizations_*err.csv) echo "None" ;;
        generated_seqs_*err.csv) echo "None" ;;

        # Catch undefined patterns
        *)                echo "Undefined" ;;   # fallback
    esac
}

for batch in "seed42" "seedRd" "noerr"
do
    assert_regression "$TESTREF/${batch}_generated" "$OUTDIR/${batch}_generated" "$LOGFILE"
done
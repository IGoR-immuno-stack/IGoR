#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}" # Create a temp dir if none passed
IGORCALL="$IGORBIN -w $OUTDIR"

$IGORCALL init
$IGORCALL config set model.source custom
$IGORCALL config set model.parms "$TESTREF/demo_inference/final_parms.txt"
$IGORCALL config set model.marginals "$TESTREF/demo_inference/final_marginals.txt"
$IGORCALL config set generate.error true
$IGORCALL config set generate.fast false

###################################################
# Generate sequences with fixed seeds
###################################################

# Generate random sequences with a single thread for reproducibility
$IGORCALL config set generate.seed 42
$IGORCALL -b seed42 -j 1 generate 100
$IGORCALL config set generate.seed 8557
$IGORCALL -b seedRd -j 1 generate 100
# Check generation without error
$IGORCALL config set generate.seed 35863
$IGORCALL config set generate.error false
$IGORCALL -b noerr -j 1 generate 100

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

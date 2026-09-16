#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------
# Regression cover for the no_d_align exhaustive position scan.
#
# Gene_choice falls back to scanning every position a D could occupy when the
# aligner hands it no D alignment at all. The standard inference corpus reaches
# that branch, but only for a handful of scenarios per iteration, so a change
# that broke it for the other 99.8 % would pass unnoticed. Emptying the D
# alignment file forces the branch for roughly half of all D choices instead,
# which is the manipulation docs/ITERATE_GENERIC_REWRITE_PLAN.md section 6
# prescribes: it needs no aligner change, only an alignment set with a header
# and no rows.
#
# One evaluate pass over the reference TRB model and the usual 300 demo
# sequences, so what is pinned is the counters that pass produces: the ten best
# scenarios per sequence, Pgen, and the V and J coverage/error counters.
#
# Provenance of the golden data, which is not the obvious choice. It was NOT
# taken from before the refactoring: this fixture is what established that the
# pre-B8 behaviour of this path was an *uninitialised read*. `a58808b` (the
# overlap safety map's port to LayeredArray) turned that read into an abort, and
# `09d1b4c` replaced it with a written, conservative "not established safe"
# verdict. So the pre-refactor output differs -- 346 of 300 Pgen rows, among
# others -- because it was undefined, not because anything regressed. The golden
# data is therefore this branch's, from the first commit at which the path has
# defined behaviour.
#
# best_scenarios_counts.csv will need regenerating when fix/scenario_tie lands:
# that branch makes degenerate-likelihood ties explicit, which changes which
# members of a tied group are reported. Nothing else here moves with it -- the
# same comparison with the D alignments left intact is bitwise across the two
# branches for Pgen and both coverage counters.
# ------------------------------------------------------------------

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}"
IGORCALL="$IGORBIN -w $OUTDIR"

$IGORCALL init
$IGORCALL config set model.source custom
$IGORCALL config set model.parms "$TESTREF/demo_inference/final_parms.txt"
$IGORCALL config set model.marginals "$TESTREF/demo_inference/final_marginals.txt"

# Reference alignments, with the D set emptied down to its header row.
mkdir -p "$OUTDIR/aligns"
cp "$TESTREF/aligns/demo_indexed_sequences.csv" "$OUTDIR/aligns/"
cp "$TESTREF/aligns/demo_V_alignments.csv" "$OUTDIR/aligns/"
cp "$TESTREF/aligns/demo_J_alignments.csv" "$OUTDIR/aligns/"
head -n 1 "$TESTREF/aligns/demo_D_alignments.csv" > "$OUTDIR/aligns/demo_D_alignments.csv"

$IGORCALL config set evaluate.likelihood_threshold 1e-35
$IGORCALL config set evaluate.probability_ratio_threshold 0.0001
$IGORCALL config set output.scenarios 10
$IGORCALL config set output.Pgen true
$IGORCALL config set output.coverage VJ_gene
$IGORCALL -b demo evaluate

# ------------------------------------------------------------------
# Output regression
# ------------------------------------------------------------------
LOGFILE="$OUTDIR/no_d_align_regression.log"

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/assert_regression.sh"

resolve_sort_columns() {
    local filename="$1"
    case "$filename" in
        best_scenarios_counts.csv)       echo "all" ;;
        *counts.csv)                     echo "col1" ;;
        *_genes_cov_and_err.csv)         echo "col1,col2" ;;
        *)                               echo "Undefined" ;;
    esac
}

# Drop the scenario rank, as test_inference.sh does: equal likelihoods can be
# ranked either way round and the rank is not what this test is about.
tmp="$(mktemp)"
cut -d';' -f 1,3- "$OUTDIR/demo_output/best_scenarios_counts.csv" >"$tmp"
mv "$tmp" "$OUTDIR/demo_output/best_scenarios_counts.csv"

# assert_regression walks the reference directory, so a missing or empty one
# passes vacuously. Say so instead.
if [[ ! -d "$TESTREF/no_d_align_output" ]] || [[ -z "$(ls -A "$TESTREF/no_d_align_output")" ]]; then
    echo "\xe2\x9d\x8c No golden data in $TESTREF/no_d_align_output - nothing was compared" >&2
    exit 1
fi

assert_regression "$TESTREF/no_d_align_output" "$OUTDIR/demo_output" "$LOGFILE"

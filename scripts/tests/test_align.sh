#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}" # Create a temp dir if none passed
IGORCALL="$IGORBIN -w $OUTDIR"

$IGORCALL init
$IGORCALL config set genomic.V "$TESTINPUT/genomicVs_with_primers.fasta"
$IGORCALL config set genomic.D "$TESTINPUT/genomicDs.fasta"
$IGORCALL config set genomic.J "$TESTINPUT/genomicJs_all_curated.fasta"

###################################################
# Align sequences with default parameters
###################################################

# Read files
$IGORCALL -b default import-seqs "$TESTINPUT/murugan_naive1_noncoding_demo_seqs.txt"
# Run alignments with the demo parameters
$IGORCALL -b default align --gene V
$IGORCALL -b default align --gene D
$IGORCALL -b default align --gene J

###################################################
# Align sequences with demo parameters
###################################################

# Read files
$IGORCALL -b demo import-seqs "$TESTINPUT/murugan_naive1_noncoding_demo_seqs.txt"
# Run alignments with the demo parameters
$IGORCALL config set alignment.V.threshold 50
$IGORCALL config set alignment.V.left_offset -999
$IGORCALL config set alignment.V.right_offset -155
$IGORCALL config set alignment.D.threshold 0
$IGORCALL config set alignment.J.threshold 10
$IGORCALL config set alignment.J.left_offset 42
$IGORCALL config set alignment.J.right_offset 48
$IGORCALL -b demo align --gene V
$IGORCALL -b demo align --gene D
$IGORCALL -b demo align --gene J
# ------------------------------------------------------------------
# 2️⃣ Test output file regression
# ------------------------------------------------------------------
LOGFILE="$OUTDIR/align_regression.log"

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/assert_regression.sh"

resolve_sort_columns() {
    local filename="$1"
    case "$filename" in
        # sort indexed seqs by first column only
        *indexed_sequences.csv) echo "col1" ;;

        # Sort alignments by seq_index and gene name
        *_alignments.csv) echo "col1,col2" ;;

        # Catch undefined patterns
        *)                echo "Undefined" ;;   # fallback
    esac
}

assert_regression "$TESTREF/aligns" "$OUTDIR/aligns" "$LOGFILE"
# The script exits with the same status that run_regression returned

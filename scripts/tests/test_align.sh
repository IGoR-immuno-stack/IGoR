#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh
OUTDIR="${1:-$(mktemp -d)}" # Create a temp dir if none passed
IGORCALL="$IGORBIN -set_wd $OUTDIR"

###################################################
# Align sequences with default parameters
###################################################

# Read files
$IGORCALL -batch default -read_seqs "$TESTINPUT/murugan_naive1_noncoding_demo_seqs.txt"
# Run alignments with the demo parameters
$IGORCALL -batch default -align --V -set_genomic --V "$TESTINPUT/genomicVs_with_primers.fasta"
$IGORCALL -batch default -align --D -set_genomic --D "$TESTINPUT/genomicDs.fasta"
$IGORCALL -batch default -align --J -set_genomic --J "$TESTINPUT/genomicJs_all_curated.fasta" 

###################################################
# Align sequences with demo parameters
###################################################

# Read files
$IGORCALL -batch demo -read_seqs "$TESTINPUT/murugan_naive1_noncoding_demo_seqs.txt"
# Run alignments with the demo parameters
$IGORCALL -batch demo -align --V ---thresh 50 ---offset_bounds -999 -155 ---best_align_only true ---gap_penalty 50 -set_genomic --V "$TESTINPUT/genomicVs_with_primers.fasta"
$IGORCALL -batch demo -align --D ---thresh 0 ---gap_penalty 50 -set_genomic --D "$TESTINPUT/genomicDs.fasta"
$IGORCALL -batch demo -align --J ---thresh 10 ---offset_bounds 42 48 ---best_align_only true ---gap_penalty 50 -set_genomic --J "$TESTINPUT/genomicJs_all_curated.fasta" 
# ------------------------------------------------------------------
# 2️⃣ Test output file regression
# ------------------------------------------------------------------
LOGFILE="align_regression.log"

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/assert_regression.sh"

declare -A SORT_PATTERNS=(
    # sort indexed seqs by first column only
    ["*indexed_sequences.csv"]="col1"

    # Sort alignments by seq_index and gene name
    ["*_alignments.csv"]="col1,col2"
)

assert_regression "$TESTREF/aligns" "$OUTDIR/aligns" "$LOGFILE"
# The script exits with the same status that run_regression returned
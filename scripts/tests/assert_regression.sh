#!/usr/bin/env bash
# --------------------------------------------------------------
# assert_regression.sh – reusable diff/check routine
# --------------------------------------------------------------
# Usage (from another script):
#   source "$(dirname "${BASH_SOURCE[0]}")/assert_regression.sh"
#   assert_regression "<ref_dir>" "<new_dir>" "<logfile>"
# --------------------------------------------------------------

set -euo pipefail

# ------------------------------------------------------------------
# Helper: compare two files and report a mismatch
# ------------------------------------------------------------------
compare_file() {
    # Arguments
    #   $1 – reference CSV file (golden)
    #   $2 – current CSV file to test
    #   $3 – optional sort mode:
    #        "col1"      – sort by the first column only (default)
    #        "col1,col2" – sort by the first column, then the second

    local ref="$1"
    local cur="$2"
    local mode="${3:-None}"   # default = sort by first column only

    # Create temporary sorted versions
    local ref_sorted=$(mktemp)
    local cur_sorted=$(mktemp)

    # Choose the appropriate sort key(s)
    case "$mode" in
        col1)
            sort -t, -k1,1 "$ref" >"$ref_sorted"
            sort -t, -k1,1 "$cur" >"$cur_sorted"
            ;;
        col1,col2)
            sort -t, -k1,1 -k2,2 "$ref" >"$ref_sorted"
            sort -t, -k1,1 -k2,2 "$cur" >"$cur_sorted"
            ;;
        None)
            cp "$ref" "$ref_sorted"
            cp "$cur" "$cur_sorted"
            ;;
        *)
            echo "⚠️ Unknown sort mode '$mode'. Falling back to unsorted comparison."
            cp "$ref" "$ref_sorted"
            cp "$cur" "$cur_sorted"
            ;;
    esac

    # Run the diff – diff returns 0 when files are identical,
    # 1 when they differ, and >1 on error.
    if diff -q "$ref_sorted" "$cur_sorted" >/dev/null 2>&1; then
        # No differences
        echo "✅ OK: $(basename "$ref")"
        rm -f "$ref_sorted" "$cur_sorted"
        return 0
    else
        # Differences found
        echo "❌ MISMATCH: $(basename "$ref")" | tee -a "$LOGFILE"
        diff -u "$ref_sorted" "$cur_sorted" | sed 's/^/    /' >>"$LOGFILE"
        rm -f "$ref_sorted" "$cur_sorted"
        return 1
    fi
}

# ------------------------------------------------------------
# Helper: resolve a filename to its column list
# ------------------------------------------------------------
# Expects a pre declared SORT_PATTERNS with following format:
# declare -A SORT_PATTERNS=(
#     # sort indexed seqs by first column only
#     ["indexed_sequences.csv"]="col1"

#     # Sort alignments by seq_index and gene name
#     ["*_alignments.csv"]="col1,col2"
# )
resolve_sort_columns() {
    local filename="$1"               # just the basename, e.g. "foo_v2.csv"
    local cols="Undefined"                    # default if nothing matches

    # Iterate over the associative array in the order it was declared.
    # The first pattern that matches wins – this mimics typical wildcard
    # precedence (more specific patterns should be placed earlier).
    for pat in "${!SORT_PATTERNS[@]}"; do
        if [[ "$filename" == $pat ]]; then
            cols="${SORT_PATTERNS[$pat]}"
            break
        fi
    done

    printf '%s' "$cols"
}

# ------------------------------------------------------------------
# Main routine – called by each test script
# ------------------------------------------------------------------
assert_regression() {
    local REF_DIR="$1"
    local NEW_DIR="$2"
    LOGFILE="${3:-regression.log}"   # default logfile if none supplied

    >"$LOGFILE"                       # start fresh for this test suite
    local suite_status=0

    shopt -s nullglob

    # --------------------------------------------------------------
    # 1️⃣  Walk through every reference file – ensure a counterpart exists
    # --------------------------------------------------------------
    for ref_path in "$REF_DIR"/*; do
        local base=$(basename "$ref_path")
        local new_path="$NEW_DIR/$base"

        if [[ ! -e "$new_path" ]]; then
            echo "❌  MISSING: $base (no corresponding output)" | tee -a "$LOGFILE"
            suite_status=1
            continue
        fi

        # Resolve the column spec for this particular file
        local col_spec="$(resolve_sort_columns "$base")"   # default to column 1

        if ! compare_file "$ref_path" "$new_path" "$col_spec"; then
            suite_status=1
        fi
    done

    # --------------------------------------------------------------
    # Final reporting
    # --------------------------------------------------------------
    if (( suite_status == 0 )); then
        echo "🎉 All checks passed for $(basename "$REF_DIR")"
        rm -f "$LOGFILE"
    else
        echo "🚨 Failures in $(basename "$REF_DIR") – see $LOGFILE"
    fi

    return $suite_status
}
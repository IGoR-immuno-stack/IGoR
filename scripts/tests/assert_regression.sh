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

    # Sorted versions are kept under $SORTED_DIR (set by assert_regression)
    # instead of throwaway tmp files, so they remain available for
    # diagnostics after a failing run.
    local base="$(basename "$ref")"
    local ref_sorted="$SORTED_DIR/${base}.ref.sorted"
    local cur_sorted="$SORTED_DIR/${base}.cur.sorted"

    # Choose the appropriate sort key(s). The header line (line 1) is kept
    # in place and never fed to sort.
    case "$mode" in
        col1)
            { head -n1 "$ref"; tail -n+2 "$ref" | sort --field-separator=";" -k1,1; } >"$ref_sorted"
            { head -n1 "$cur"; tail -n+2 "$cur" | sort --field-separator=";" -k1,1; } >"$cur_sorted"
            ;;
        col1,col2)
            { head -n1 "$ref"; tail -n+2 "$ref" | sort --field-separator=";" -k1,1 -k2,2; } >"$ref_sorted"
            { head -n1 "$cur"; tail -n+2 "$cur" | sort --field-separator=";" -k1,1 -k2,2; } >"$cur_sorted"
            ;;
        all)
            { head -n1 "$ref"; tail -n+2 "$ref" | sort --field-separator=";"; } >"$ref_sorted"
            { head -n1 "$cur"; tail -n+2 "$cur" | sort --field-separator=";"; } >"$cur_sorted"
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
        return 0
    else
        # Differences found
        echo "❌ MISMATCH: $(basename "$ref")" | tee -a "$LOGFILE"
        echo "" >>"$LOGFILE"
        echo "    ----------------------------------------" >>"$LOGFILE"
        echo "    Reference: $ref" >>"$LOGFILE"
        echo "    Generated: $cur" >>"$LOGFILE"
        echo "    ----------------------------------------" >>"$LOGFILE"
        diff -u "$ref_sorted" "$cur_sorted" | head -100 >>"$LOGFILE"
        if [[ $(diff -u "$ref_sorted" "$cur_sorted" | wc -l) -gt 100 ]]; then
            echo "    ... (diff truncated, showing first 100 lines)" >>"$LOGFILE"
        fi
        echo "    ----------------------------------------" >>"$LOGFILE"
        echo "" >>"$LOGFILE"
        return 1
    fi
}

# ------------------------------------------------------------------
# Main routine – called by each test script
# ------------------------------------------------------------------
assert_regression() {
    local REF_DIR="$1"
    local NEW_DIR="$2"
    LOGFILE="${3:-regression.log}"   # default logfile if none supplied

    # Ensure LOGFILE directory exists
    local LOGDIR=$(dirname "$LOGFILE")
    if [[ ! -d "$LOGDIR" ]]; then
        mkdir -p "$LOGDIR"
    fi

    # Where compare_file stores sorted copies of ref/cur files. Scoped under
    # NEW_DIR (rather than LOGDIR) because assert_regression may be called
    # several times against the same LOGFILE with different NEW_DIR/REF_DIR
    # pairs (e.g. one per batch/config) — files sharing a basename across
    # those calls would otherwise collide. Reset for every run, then kept
    # around (not deleted) whenever the suite fails, so the sorted files
    # remain available for diagnostics.
    SORTED_DIR="$NEW_DIR/sorted"
    rm -rf "$SORTED_DIR"
    mkdir -p "$SORTED_DIR"

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
        rm -rf "$SORTED_DIR"
    else
        echo "🚨 Failures in $(basename "$REF_DIR") – see $LOGFILE"
        echo "    Sorted ref/cur files for inspection: $SORTED_DIR"
    fi

    return $suite_status
}

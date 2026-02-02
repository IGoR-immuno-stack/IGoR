#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $SCRIPT_DIR/config.sh

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
REPORT_DIR="regression_reports"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
REPORT_BASE_DIR="$REPORT_DIR/regression_$TIMESTAMP"
mkdir -p "$REPORT_BASE_DIR"

# Test definitions with name, script, and description
declare -A TEST_NAMES=(
    [1]="align"
    [2]="inference"
    [3]="generate"
)

declare -A TEST_SCRIPTS=(
    [1]="$SCRIPT_DIR/test_align.sh"
    [2]="$SCRIPT_DIR/test_inference.sh"
    [3]="$SCRIPT_DIR/test_generate.sh"
)

declare -A TEST_DESCRIPTIONS=(
    [1]="Alignment tests"
    [2]="Inference tests"
    [3]="Generation tests"
)

# ------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------
usage() {
    cat << EOF
Usage: $0 [OPTIONS] [TEST_SPEC]

Run regression tests and generate detailed reports.

Arguments:
  TEST_SPEC           Test to run (index, name, or 'all')
                      Index: 1, 2, 3
                      Name: align, inference, generate
                      Default: all

Options:
  -h, --help          Show this help message
  -l, --list          List available tests
  -c, --clean         Clean old reports (before running)
  -k, --keep          Keep all output files for debugging

Examples (pixi):
  pixi run test_regression               Run all tests
  pixi run test_regression 1             Run test 1 (align)
  pixi run test_regression align         Run align test
  pixi run test_regression 1,2           Run tests 1 and 2
  pixi run test_regression align,inference Run align and inference tests
  pixi run test_regression --list        List available tests
  pixi run test_regression -k 1          Run test 1 and keep output

Direct execution:
  $0                  Run all tests
  $0 1                Run test 1 (align)
  $0 -l               List available tests

EOF
}

# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------
KEEP_OUTPUT=false
CLEAN_REPORTS=false
LIST_TESTS=false
TEST_SPEC="all"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -l|--list)
            LIST_TESTS=true
            shift
            ;;
        -c|--clean)
            CLEAN_REPORTS=true
            shift
            ;;
        -k|--keep)
            KEEP_OUTPUT=true
            shift
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            TEST_SPEC="$1"
            shift
            ;;
    esac
done

# ------------------------------------------------------------------
# List tests if requested
# ------------------------------------------------------------------
if $LIST_TESTS; then
    echo "Available regression tests:"
    echo "---------------------------"
    printf "%-3s | %-15s | %s\n" "Idx" "Name" "Description"
    echo "-----|-----------------|-------------------------------------"
    for i in $(seq 1 ${#TEST_NAMES[@]}); do
        printf "%-3s | %-15s | %s\n" "$i" "${TEST_NAMES[$i]}" "${TEST_DESCRIPTIONS[$i]}"
    done
    echo "-----|-----------------|-------------------------------------"
    exit 0
fi

# ------------------------------------------------------------------
# Clean old reports if requested
# ------------------------------------------------------------------
if $CLEAN_REPORTS && [[ -d "$REPORT_DIR" ]]; then
    echo "Cleaning old reports in $REPORT_DIR..."
    rm -rf "$REPORT_DIR"
fi

# ------------------------------------------------------------------
# Parse test specification
# ------------------------------------------------------------------
parse_test_spec() {
    local spec="$1"
    local tests=()

    if [[ "$spec" == "all" ]]; then
        for i in $(seq 1 ${#TEST_NAMES[@]}); do
            tests+=("$i")
        done
        echo "${tests[@]}"
        return
    fi

    # Handle comma-separated values
    IFS=',' read -ra ITEMS <<< "$spec"
    for item in "${ITEMS[@]}"; do
        item=$(echo "$item" | xargs)  # trim whitespace

        # Handle range (e.g., 1-3)
        if [[ "$item" =~ ^[0-9]+-[0-9]+$ ]]; then
            local start=$(echo "$item" | cut -d'-' -f1)
            local end=$(echo "$item" | cut -d'-' -f2)
            for ((i=start; i<=end; i++)); do
                if [[ -n "${TEST_NAMES[$i]:-}" ]]; then
                    tests+=("$i")
                else
                    echo "Error: Invalid test index '$i'" >&2
                    echo "Run '$0 --list' to see available tests." >&2
                    exit 1
                fi
            done
        # Check if it's an index
        elif [[ "$item" =~ ^[0-9]+$ ]]; then
            if [[ -n "${TEST_NAMES[$item]:-}" ]]; then
                tests+=("$item")
            else
                echo "Error: Invalid test index '$item'" >&2
                echo "Run '$0 --list' to see available tests." >&2
                exit 1
            fi
        # Check if it's a name
        else
            local found=false
            for i in $(seq 1 ${#TEST_NAMES[@]}); do
                if [[ "${TEST_NAMES[$i]}" == "$item" ]]; then
                    tests+=("$i")
                    found=true
                    break
                fi
            done
            if ! $found; then
                echo "Error: Invalid test name '$item'" >&2
                echo "Run '$0 --list' to see available tests." >&2
                exit 1
            fi
        fi
    done

    # Remove duplicates and sort
    printf '%s\n' "${tests[@]}" | sort -nu | tr '\n' ' '
}

TESTS_TO_RUN=($(parse_test_spec "$TEST_SPEC"))

# ------------------------------------------------------------------
# Run tests
# ------------------------------------------------------------------
overall_status=0
declare -a TEST_RESULTS

echo ""
echo "=============================================================================="
echo "REGRESSION TEST SUITE"
echo "=============================================================================="
echo "Report directory: $REPORT_BASE_DIR"
echo "Tests to run: ${TESTS_TO_RUN[*]}"
echo "=============================================================================="
echo ""

for test_idx in "${TESTS_TO_RUN[@]}"; do
    test_name="${TEST_NAMES[$test_idx]}"
    test_script="${TEST_SCRIPTS[$test_idx]}"
    test_desc="${TEST_DESCRIPTIONS[$test_idx]}"
    test_outdir="$REPORT_BASE_DIR/${test_name}_output"
    test_report="$REPORT_BASE_DIR/${test_name}_report.txt"

    printf "Running test %d: %s (%s)\n" "$test_idx" "$test_name" "$test_desc"

    # Create test-specific output directory
    mkdir -p "$test_outdir"

    # Initialize report file
    cat > "$test_report" << EOF
==============================================================================
REGRESSION TEST REPORT: $test_name
==============================================================================
Test:       $test_name (index $test_idx)
Description: $test_desc
Script:     $test_script
Started:    $(date)
Output dir: $test_outdir

------------------------------------------------------------------------------
EXECUTION LOG
------------------------------------------------------------------------------

EOF

    echo "" >> "$test_report"
    echo "------------------------------------------------------------------------------" >> "$test_report"
    echo "DIFF OUTPUT" >> "$test_report"
    echo "------------------------------------------------------------------------------" >> "$test_report"
    echo "" >> "$test_report"

    # Run the test and capture output
    test_status=0

    if "$test_script" "$test_outdir" >> "$test_report" 2>&1; then
        test_status=0
        echo "    ✅ PASSED"
        echo "" >> "$test_report"
        echo "STATUS: PASSED" >> "$test_report"
    else
        test_status=1
        overall_status=1
        echo "    ❌ FAILED"
        echo "" >> "$test_report"
        echo "STATUS: FAILED" >> "$test_report"

        # Extract and display diff summary
        echo "    Differences found. See $test_report for details."

        # Count mismatches
        mismatches=$(grep -c "^❌ MISMATCH:" "$test_report" || true)
        if [[ $mismatches -gt 0 ]]; then
            echo "    Files with mismatches: $mismatches"
        fi
    fi

    echo "    Finished: $(date)"
    echo "    Finalized: $(date)" >> "$test_report"
    echo "" >> "$test_report"

    TEST_RESULTS+=("$test_idx|$test_name|$test_status|$test_report")
done

# Clean up output directories if not keeping them
if ! $KEEP_OUTPUT; then
    for test_idx in "${TESTS_TO_RUN[@]}"; do
        test_name="${TEST_NAMES[$test_idx]}"
        test_outdir="$REPORT_BASE_DIR/${test_name}_output"

        # Preserve any .log files before cleanup
        if [[ -d "$test_outdir" ]]; then
            while IFS= read -r -d '' log_file; do
                # Extract relative path from test_outdir for better organization
                rel_path="${log_file#$test_outdir/}"
                # Replace '/' with '_' for filesystem safety
                safe_name=$(echo "$rel_path" | tr '/' '_')
                cp "$log_file" "$REPORT_BASE_DIR/${safe_name}"
            done < <(find "$test_outdir" -type f -name "*.log" -print0 2>/dev/null)

            rm -rf "$test_outdir"
            echo "Cleaned: $test_outdir"
        fi
    done
fi

# ------------------------------------------------------------------
# Generate summary
# ------------------------------------------------------------------
SUMMARY_FILE="$REPORT_BASE_DIR/SUMMARY.txt"

cat > "$SUMMARY_FILE" << EOF
==============================================================================
REGRESSION TEST SUMMARY
==============================================================================
Timestamp:  $(date)
Report dir: $REPORT_BASE_DIR
Tests run:   ${TESTS_TO_RUN[*]}

------------------------------------------------------------------------------
RESULTS
------------------------------------------------------------------------------
EOF

echo ""
echo "=============================================================================="
echo "REGRESSION TEST SUMMARY"
echo "=============================================================================="
printf "%-3s | %-15s | %-11s | %s\n" "Idx" "Test" "Status" "Report"
echo "-----|-----------------|-------------|------------------------------------"

for result in "${TEST_RESULTS[@]}"; do
    IFS='|' read -r idx name status report <<< "$result"
    status_str="✅ PASSED"
    report_rel="${report##$(pwd)/}"
    if [[ $status -eq 1 ]]; then
        status_str="❌ FAILED"
    fi
    printf "%-3s | %-15s | %-11s | %s\n" "$idx" "$name" "$status_str" "$report_rel"

    # Also write to summary file
    printf "%-3s | %-15s | %-11s | %s\n" "$idx" "$name" "$status_str" "$report_rel" >> "$SUMMARY_FILE"
done

echo "-----|-----------------|-------------|------------------------------------"

# Add detailed diff summary
printf "" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "------------------------------------------------------------------------------" >> "$SUMMARY_FILE"
echo "DIFF DETAILS" >> "$SUMMARY_FILE"
echo "------------------------------------------------------------------------------" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

had_diffs=false
for result in "${TEST_RESULTS[@]}"; do
    IFS='|' read -r idx name status report <<< "$result"
    if [[ $status -eq 1 ]]; then
        had_diffs=true
        echo "Test $idx ($name):" >> "$SUMMARY_FILE"
        # Check for mismatch patterns in the report file
        grep "^❌ MISMATCH:" "$report" 2>/dev/null | sed 's/❌ MISMATCH: /  - /' >> "$SUMMARY_FILE" || true
        echo "" >> "$SUMMARY_FILE"
    fi
done

if ! $had_diffs; then
    echo "No differences found (all tests passed)." >> "$SUMMARY_FILE"
fi

# ------------------------------------------------------------------
# Final report
# ------------------------------------------------------------------
echo ""
if [[ $overall_status -eq 0 ]]; then
    echo "🎉 All regression tests PASSED"
else
    echo "🚨 One or more tests FAILED"
    echo ""
    echo "To view differences:"
    for result in "${TEST_RESULTS[@]}"; do
        IFS='|' read -r idx name status report <<< "$result"
        if [[ $status -eq 1 ]]; then
            echo "  Test $idx ($name): cat $report"
        fi
    done

    # Show preserved log files for failed tests only
    echo ""
    log_files=()
    for result in "${TEST_RESULTS[@]}"; do
        IFS='|' read -r idx name status report <<< "$result"
        if [[ $status -eq 1 ]]; then
            for log in "$REPORT_BASE_DIR/${name}_"*.log; do
                if [[ -f "$log" ]]; then
                    log_files+=("$log")
                fi
            done
        fi
    done

    if [[ ${#log_files[@]} -gt 0 ]]; then
        echo "Detailed log files:"
        for log in "${log_files[@]}"; do
            echo "  - $log"
        done
    fi
fi

echo ""
echo "Summary report: $SUMMARY_FILE"
echo "=============================================================================="

exit $overall_status

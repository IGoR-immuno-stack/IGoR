#!/usr/bin/env bash
set -euo pipefail

source config.sh
OUTDIR=$(mktemp -d)

# ------------------------------------------------------------------
# Run each unit‑test script in sequence.
# If any script exits non‑zero, remember that fact but keep going
# so we get a full report from every suite.
# ------------------------------------------------------------------
overall_status=0

run_one() {
    local script="$1"
    echo -e "\n=== Running ${script} ==="
    if ./"$script" "$OUTDIR"; then
        echo "✅ ${script} succeeded"
    else
        echo "❌ ${script} FAILED"
        overall_status=1
    fi
}

run_one test_align.sh
run_one test_inference.sh
run_one test_evaluate.sh
run_one test_seqgen.sh

echo -e "\n===== REGRESSION SUMMARY ====="
if (( overall_status == 0 )); then
    echo "🎉 All test suites passed."
else
    echo "🚨 One or more suites failed – inspect the *.log files."
fi

exit $overall_status
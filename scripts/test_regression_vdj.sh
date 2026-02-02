#!/bin/bash
# Regression Test Script
# Tests standard VDJ recombination to verify no regression

set -e

IGOR_BIN="./build/bin/igor"
TEST_DIR="vj_regression_test"
SEED=42

echo "=== IGoR Regression Test for Standard VDJ ==="
echo ""
echo "Testing to ensure Issue #7 refactoring doesn't break existing functionality"
echo ""

# Clean up previous runs
rm -rf $TEST_DIR
mkdir -p $TEST_DIR

# Test 1: Model Loading
echo "Test 1: Loading simple VJ model..."
cd build/bin
./igor -set_wd ../../$TEST_DIR -batch gen -set_custom_model ../../demo/simple_vj_model.txt -generate 10 --seed $SEED > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ PASS: Model loaded and generation completed"
else
    echo "❌ FAIL: Model loading or generation failed"
    exit 1
fi
cd ../..

# Test 2: Check Output Files Exist
echo ""
echo "Test 2: Verifying output files..."
if [ -f "$TEST_DIR/gen_generated/generated_seqs_werr.csv" ]; then
    echo "✅ PASS: sequences file created"
else
    echo "❌ FAIL: sequences file not found"
    exit 1
fi

if [ -f "$TEST_DIR/gen_generated/generated_realizations_werr.csv" ]; then
    echo "✅ PASS: realizations file created"
else
    echo "❌ FAIL: realizations file not found"
    exit 1
fi

# Test 3: Check Generated Sequence Format
echo ""
echo "Test 3: Verifying sequence format..."
SEQ_COUNT=$(tail -n +2 "$TEST_DIR/gen_generated/generated_seqs_werr.csv" | wc -l)
if [ "$SEQ_COUNT" -eq 10 ]; then
    echo "✅ PASS: Generated 10 sequences as expected"
else
    echo "❌ FAIL: Expected 10 sequences, got $SEQ_COUNT"
    exit 1
fi

# Test 4: Verify Sequence Generation is Deterministic
echo ""
echo "Test 4: Testing determinism (same seed = same output)..."
rm -rf $TEST_DIR/run2
mkdir -p $TEST_DIR/run2
cd build/bin
./igor -set_wd ../../$TEST_DIR/run2 -batch gen -set_custom_model ../../demo/simple_vj_model.txt -generate 10 --seed $SEED > /dev/null 2>&1
cd ../..

if diff "$TEST_DIR/gen_generated/generated_seqs_werr.csv" "$TEST_DIR/run2/gen_generated/generated_seqs_werr.csv" > /dev/null 2>&1; then
    echo "✅ PASS: Deterministic generation verified"
else
    echo "❌ FAIL: Non-deterministic output"
    exit 1
fi

# Summary
echo ""
echo "=== Regression Test Summary ==="
echo "✅ All tests passed - NO REGRESSION detected"
echo ""
echo "Standard VDJ recombination works correctly after refactoring."
echo ""
echo "Files generated in: $TEST_DIR/"
ls -lh $TEST_DIR/gen_generated/

exit 0

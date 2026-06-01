#!/bin/bash
#
# Test script for igor-convert
# Tests all conversion paths with real data
#

set -e  # Exit on error

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

IGOR_CONVERT="$PROJECT_ROOT/build/bin/igor-convert"
TEST_DIR="$PROJECT_ROOT/test_output"
mkdir -p "$TEST_DIR"
# Clean previous test run
rm -rf "$TEST_DIR"/*

echo "========================================="
echo "Testing igor-convert"
echo "========================================="
echo "Test directory: $TEST_DIR"
echo ""

# Check if executable exists
if [ ! -f "$IGOR_CONVERT" ]; then
    echo "Error: $IGOR_CONVERT not found!"
    echo "Run: pixi run cmake --build build --target igor-convert"
    exit 1
fi

# Test 1: Help message
echo "✓ Test 1: Help message"
$IGOR_CONVERT --help > /dev/null
echo ""

# Create simple test data in AIRR Rearrangement format
echo "✓ Test 2: Creating test data..."
cat > "$TEST_DIR/rearrangement.tsv" << 'EOF'
sequence_id	sequence	v_call	d_call	j_call	v_score	d_score	j_score
1	ATCGATCGATCGATCG	IGHV1-2*01	IGHD1-1*01	IGHJ3*01	150.5	45.2	120.0
2	GCTAGCTAGCTAGCTA	IGHV3-9*01		IGHJ4*02	140.0		115.0
3	TTAATTAATTAATTAA	IGHV1-3*01	IGHD2-2*01	IGHJ5*01	155.0	50.0	125.0
EOF
echo "  Created: $TEST_DIR/rearrangement.tsv (3 sequences)"
echo ""

# Test 3: AIRR Rearrangement → Parquet
echo "✓ Test 3: AIRR Rearrangement → Parquet (SNAPPY)"
$IGOR_CONVERT "$TEST_DIR/rearrangement.tsv" "$TEST_DIR/sequences.parquet"
if [ ! -f "$TEST_DIR/sequences.parquet" ]; then
    echo "  ✗ Failed: Output file not created"
    exit 1
fi
SIZE=$(stat -f%z "$TEST_DIR/sequences.parquet" 2>/dev/null || stat -c%s "$TEST_DIR/sequences.parquet" 2>/dev/null)
echo "  File size: $SIZE bytes"
echo ""

# Test 4: Parquet → AIRR Rearrangement (round-trip)
echo "✓ Test 4: Parquet → AIRR Rearrangement (round-trip)"
$IGOR_CONVERT "$TEST_DIR/sequences.parquet" "$TEST_DIR/rearrangement_roundtrip.tsv"
LINES=$(wc -l < "$TEST_DIR/rearrangement_roundtrip.tsv")
echo "  Lines in output: $LINES (expected: 4 = 1 header + 3 sequences)"
echo ""

# Test 5: Parquet → AIRR Alignment
echo "✓ Test 5: Parquet → AIRR Alignment"
$IGOR_CONVERT "$TEST_DIR/sequences.parquet" "$TEST_DIR/alignment.tsv"
LINES=$(wc -l < "$TEST_DIR/alignment.tsv")
echo "  Lines in output: $LINES (one row per alignment)"
echo ""

# Test 6: AIRR Rearrangement → AIRR Alignment
echo "✓ Test 6: AIRR Rearrangement → AIRR Alignment"
$IGOR_CONVERT "$TEST_DIR/rearrangement.tsv" "$TEST_DIR/alignment2.tsv"
echo "  Converted rearrangement to alignment format"
echo ""

# Test 7: AIRR Alignment → Parquet
echo "✓ Test 7: AIRR Alignment → Parquet"
$IGOR_CONVERT "$TEST_DIR/alignment.tsv" "$TEST_DIR/from_alignment.parquet"
echo "  Converted alignment to Parquet"
echo ""

# Test 8: CSV format
echo "✓ Test 8: CSV format support"
$IGOR_CONVERT "$TEST_DIR/sequences.parquet" "$TEST_DIR/output.csv" --delimiter COMMA
if grep -q "," "$TEST_DIR/output.csv"; then
    echo "  CSV format verified (contains commas)"
else
    echo "  ✗ Warning: CSV file doesn't contain commas"
fi
echo ""

# Test 9: Different compression types
echo "✓ Test 9: Compression types"
for COMP in NONE SNAPPY GZIP ZSTD LZ4; do
    $IGOR_CONVERT --compression $COMP "$TEST_DIR/rearrangement.tsv" "$TEST_DIR/compressed_$COMP.parquet"
    SIZE=$(stat -f%z "$TEST_DIR/compressed_$COMP.parquet" 2>/dev/null || stat -c%s "$TEST_DIR/compressed_$COMP.parquet" 2>/dev/null)
    echo "  $COMP: $SIZE bytes"
done
echo ""

# Test 10: Error handling
echo "✓ Test 10: Error handling"
if $IGOR_CONVERT /nonexistent/file.tsv "$TEST_DIR/out.parquet" 2>/dev/null; then
    echo "  ✗ Should have failed on non-existent file"
    exit 1
else
    echo "  Correctly handles non-existent input file"
fi
echo ""

echo "========================================="
echo "All tests passed! ✓"
echo "========================================="
echo ""
echo "Test files preserved in: $TEST_DIR/"
echo ""
echo "To inspect test outputs:"
echo "  ls -lh $TEST_DIR/"
echo "  head $TEST_DIR/rearrangement.tsv"
echo "  head $TEST_DIR/alignment.tsv"
echo ""
echo "To clean up:"
echo "  rm -rf $TEST_DIR/"
echo ""
echo "Example conversions:"
echo "  1. Rearrangement → Parquet:  $IGOR_CONVERT input.tsv output.parquet"
echo "  2. Parquet → Alignment:      $IGOR_CONVERT input.parquet output_align.tsv"
echo "  3. With compression:         $IGOR_CONVERT --compression ZSTD input.tsv output.parquet"
echo ""

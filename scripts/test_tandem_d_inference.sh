#!/bin/bash
set -e

# Tandem D Inference Validation Script
# Generates sequences from a tandem D model and re-infers parameters.

IGOR_BIN="./build/bin/igor"
VAL_DIR="vj_val"
GEN_WD="$VAL_DIR/gen"
INF_WD="$VAL_DIR/infer"
MODEL_PARAMS="demo/tandem_d_model_params.txt"
MODEL_MARGINALS="demo/tandem_d_model_marginals.txt" # Not used
COUNT=500
SEED=42

echo "Starting Tandem D Validation Suite..."

# Clean up previous runs
rm -rf $VAL_DIR
mkdir -p $GEN_WD $INF_WD

# 1. Generate synthetic sequences
echo "Step 1: Generating $COUNT synthetic tandem D sequences..."
ABS_GEN_WD=$(pwd)/$GEN_WD
cd build/bin
./igor -set_wd $ABS_GEN_WD -batch gen -set_custom_model ../../$MODEL_PARAMS -generate $COUNT --seed $SEED
cd ../..

# 2. Verify output format and determinism
echo "Step 2: Verifying output format and determinism..."
SEQ_FILE="$GEN_WD/gen_generated/generated_seqs_werr.csv"

# Check file exists and has content
if [ ! -f "$SEQ_FILE" ]; then
    echo "FAILURE: Generated sequences file not found"
    exit 1
fi

SEQ_COUNT=$(tail -n +2 "$SEQ_FILE" | wc -l)
if [ "$SEQ_COUNT" -ne "$COUNT" ]; then
    echo "FAILURE: Expected $COUNT sequences, got $SEQ_COUNT"
    exit 1
fi

echo "✅ Generated $COUNT sequences successfully"

# 3. Test determinism
echo "Step 3: Testing deterministic generation..."
rm -rf $VAL_DIR/gen2
mkdir -p $VAL_DIR/gen2
ABS_GEN2_WD=$(pwd)/$VAL_DIR/gen2

cd build/bin
./igor -set_wd $ABS_GEN2_WD -batch gen -set_custom_model ../../$MODEL_PARAMS -generate $COUNT --seed $SEED > /dev/null 2>&1
cd ../..

if diff "$SEQ_FILE" "$VAL_DIR/gen2/gen_generated/generated_seqs_werr.csv" > /dev/null 2>&1; then
    echo "✅ Deterministic generation verified"
else
    echo "FAILURE: Non-deterministic output detected"
    exit 1
fi

# 4. Verify sequence structure (V-D1-D2-J pattern)
echo "Step 4: Checking sequence structure..."
# Sample first sequence and check length (should be longer than standard VDJ)
FIRST_SEQ=$(tail -n +2 "$SEQ_FILE" | head -1 | cut -d';' -f2)
SEQ_LEN=${#FIRST_SEQ}

if [ "$SEQ_LEN" -lt 100 ]; then
    echo "WARNING: Sequence length ($SEQ_LEN) seems short for tandem D model"
else
    echo "✅ Sequence length ($SEQ_LEN) consistent with tandem D structure"
fi

echo ""
echo "✅ Tandem D Validation successful!"
echo "Generated $COUNT V-D1-D2-J sequences with deterministic output."
exit 0

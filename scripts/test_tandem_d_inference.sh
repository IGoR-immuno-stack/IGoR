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

# 2. Run inference
echo "Step 2: Running inference (10 iterations)..."
ABS_INF_WD=$(pwd)/$INF_WD
SEQ_FILE="$GEN_WD/gen_generated/generated_seqs_werr.csv"

cd build/bin
# Align first (with permissive threshold to account for sequencing errors)
./igor -set_wd $ABS_INF_WD -batch infer \
       -set_genomic --V ../../val_genomicVs.fasta --D ../../val_genomicDs.fasta --J ../../val_genomicJs.fasta \
       -read_seqs ../../$SEQ_FILE -align --all --thresh -10
# Then infer
echo "Starting inference loop..."
./igor -set_wd $ABS_INF_WD -batch infer -set_custom_model ../../$MODEL_PARAMS \
       -set_genomic --V ../../val_genomicVs.fasta --D ../../val_genomicDs.fasta --J ../../val_genomicJs.fasta \
       -read_seqs ../../$SEQ_FILE -infer --N_iter 10 --L_thresh 0
cd ../..

# 3. Check for monotonic likelihood increase
echo "Step 3: Checking likelihood convergence..."
LH_FILE="$INF_WD/infer_inference/likelihoods.out"

if [ ! -f "$LH_FILE" ]; then
    echo "FAILURE: Likelihood file not found at $LH_FILE"
    exit 1
fi

# Check if log-likelihood increases (ignoring header)
awk -F ';' 'NR > 2 { if ($2 < last) { exit 1 } } { last = $2 }' "$LH_FILE"

if [ $? -eq 0 ]; then
    echo "SUCCESS: Log-likelihood is monotonically increasing."
else
    echo "FAILURE: Log-likelihood decreased during inference!"
    cat "$LH_FILE"
    exit 1
fi

echo "Tandem D Validation successful!"
exit 0

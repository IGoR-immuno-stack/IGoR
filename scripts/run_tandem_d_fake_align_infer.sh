#!/usr/bin/env bash
set -euo pipefail

# Run tandem D generation and fake alignment workflow
# This script generates sequences from the tandem D model and creates fake alignments.
# 
# NOTE: Inference step is currently disabled because tandem D inference requires
# D gene alignments which are not yet implemented in the fake alignment generator.
# To run inference, you would need to:
# 1. Generate real D alignments using the aligner
# 2. Or extend make_fake_alignments.py to create D alignments
#
# Usage:
#   ./scripts/run_tandem_d_fake_align_infer.sh
#
# Environment variables:
#   IGOR_BIN - Path to igor binary (default: build/bin/igor)
#   WORKDIR - Working directory (default: /tmp/tandem_d_fake)
#   BATCH_NAME - Batch name (default: tandem_d)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IGOR_BIN="${IGOR_BIN:-$ROOT_DIR/build/bin/igor}"
MODEL_PARAMS="$ROOT_DIR/demo/tandem_d_model_params.txt"
WORKDIR="${WORKDIR:-/tmp/tandem_d_fake}"
BATCH_NAME="${BATCH_NAME:-tandem_d}"

# Note: Marginals file is optional. If not provided, uniform marginals will be used.
# MODEL_MARGINALS="$ROOT_DIR/demo/tandem_d_model_marginals.txt"

mkdir -p "$WORKDIR"

echo "=== Step 1: Generating sequences from tandem D model ==="
"$IGOR_BIN" \
  -set_wd "$WORKDIR" \
  -batch "$BATCH_NAME" \
  -set_custom_model "$MODEL_PARAMS" \
  -generate 1000 --noerr

GENERATED_SEQS="$WORKDIR/${BATCH_NAME}_generated/generated_seqs_noerr.csv"

echo ""
echo "=== Step 2: Indexing generated sequences ==="
"$IGOR_BIN" \
  -set_wd "$WORKDIR" \
  -batch "$BATCH_NAME" \
  -read_seqs "$GENERATED_SEQS"

ALIGN_DIR="$WORKDIR/aligns"
INDEXED_SEQS="$ALIGN_DIR/${BATCH_NAME}_indexed_sequences.csv"

echo ""
echo "=== Step 3: Creating fake alignments ==="
python3 "$ROOT_DIR/scripts/make_fake_alignments.py" \
  --model-params "$MODEL_PARAMS" \
  --indexed-seqs "$INDEXED_SEQS" \
  --out-v "$ALIGN_DIR/${BATCH_NAME}_V_alignments.csv" \
  --out-j "$ALIGN_DIR/${BATCH_NAME}_J_alignments.csv"

echo ""
echo "=== Workflow complete ==="
echo "Generated files:"
echo "  Sequences: $GENERATED_SEQS"
echo "  V alignments: $ALIGN_DIR/${BATCH_NAME}_V_alignments.csv"
echo "  J alignments: $ALIGN_DIR/${BATCH_NAME}_J_alignments.csv"
echo ""
echo "NOTE: Inference step is skipped. To run inference, you need D gene alignments."
echo "You can either:"
echo "  1. Use real alignments: $IGOR_BIN -align --all ..."
echo "  2. Extend make_fake_alignments.py to create D alignments"
echo "  3. Run inference without D alignments (VJ only):"
echo "     $IGOR_BIN -set_wd $WORKDIR -batch $BATCH_NAME -set_custom_model $MODEL_PARAMS -infer --N_iter 5"

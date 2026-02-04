#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IGOR_BIN="${IGOR_BIN:-$ROOT_DIR/build/bin/igor}"
MODEL_PARAMS="$ROOT_DIR/demo/tandem_d_model_params.txt"
MODEL_MARGINALS="$ROOT_DIR/demo/tandem_d_model_marginals.txt"
WORKDIR="${WORKDIR:-/tmp/tandem_d_fake}"
BATCH_NAME="${BATCH_NAME:-tandem_d}"

mkdir -p "$WORKDIR"

"$IGOR_BIN" \
  -set_wd "$WORKDIR" \
  -batch "$BATCH_NAME" \
  -set_custom_model "$MODEL_PARAMS" "$MODEL_MARGINALS" \
  -generate 1000 --noerr

GENERATED_SEQS="$WORKDIR/${BATCH_NAME}_generated/generated_seqs_noerr.csv"

"$IGOR_BIN" \
  -set_wd "$WORKDIR" \
  -batch "$BATCH_NAME" \
  -read_seqs "$GENERATED_SEQS"

ALIGN_DIR="$WORKDIR/aligns"
INDEXED_SEQS="$ALIGN_DIR/${BATCH_NAME}_indexed_sequences.csv"

python3 "$ROOT_DIR/scripts/make_fake_alignments.py" \
  --model-params "$MODEL_PARAMS" \
  --indexed-seqs "$INDEXED_SEQS" \
  --out-v "$ALIGN_DIR/${BATCH_NAME}_V_alignments.csv" \
  --out-j "$ALIGN_DIR/${BATCH_NAME}_J_alignments.csv"

"$IGOR_BIN" \
  -set_wd "$WORKDIR" \
  -batch "$BATCH_NAME" \
  -set_custom_model "$MODEL_PARAMS" "$MODEL_MARGINALS" \
  -infer --N_iter 5

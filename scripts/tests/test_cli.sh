#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/config.sh"

OUTDIR="${1:-$(mktemp -d)}"
IGORCALL="$IGORBIN -w $OUTDIR"

"$IGORBIN" --help >/dev/null
"$IGORBIN" --version >/dev/null
"$IGORBIN" datadir >/dev/null
"$IGORBIN" datadir --help >/dev/null
"$IGORBIN" init --help >/dev/null
"$IGORBIN" config --help >/dev/null
"$IGORBIN" import-seqs --help >/dev/null
"$IGORBIN" align --help >/dev/null
"$IGORBIN" infer --help >/dev/null
"$IGORBIN" evaluate --help >/dev/null
"$IGORBIN" generate --help >/dev/null
"$IGORBIN" run --help >/dev/null
"$IGORBIN" config list | grep -q '^model.source'

$IGORCALL init
test -f "$OUTDIR/.igor/config.toml"
test -d "$OUTDIR/.igor/runs"
test -d "$OUTDIR/aligns"
test -d "$OUTDIR/output"

$IGORCALL config get model.source >/dev/null
$IGORCALL config set model.source custom
$IGORCALL config set model.parms "$TESTREF/demo_inference/final_parms.txt"
$IGORCALL config set model.marginals "$TESTREF/demo_inference/final_marginals.txt"
$IGORCALL config set generate.seed 42
$IGORCALL config set generate.error false

if $IGORCALL config get does.not.exist >/dev/null 2>&1; then
    echo "expected invalid config key to fail" >&2
    exit 1
fi

if $IGORCALL config set generate.threads -1 >/dev/null 2>&1; then
    echo "expected invalid config value to fail" >&2
    exit 1
fi

if "$IGORBIN" -set_wd "$OUTDIR" datadir >/dev/null 2>&1; then
    echo "expected legacy -set_wd flag to fail" >&2
    exit 1
fi

if "$IGORBIN" -w "$OUTDIR/missing" infer >/dev/null 2>&1; then
    echo "expected missing workdir config to fail" >&2
    exit 1
fi

$IGORCALL -b smoke -j 1 generate 5 >/dev/null
MANIFEST=$(ls -t "$OUTDIR/.igor/runs/"*generate.toml | head -1)
test -f "$MANIFEST"

$IGORCALL run --replay "$MANIFEST" >/dev/null

MODEL_COPY="$OUTDIR/model_parms_mutated.txt"
cp "$TESTREF/demo_inference/final_parms.txt" "$MODEL_COPY"
$IGORCALL config set model.parms "$MODEL_COPY"
$IGORCALL -b mutated -j 1 generate 1 >/dev/null
BAD_MANIFEST=$(ls -t "$OUTDIR/.igor/runs/"*generate.toml | head -1)
printf '\n# mutate for replay hash mismatch\n' >> "$MODEL_COPY"
if $IGORCALL run --replay "$BAD_MANIFEST" >/dev/null 2>&1; then
    echo "expected replay hash mismatch to fail" >&2
    exit 1
fi

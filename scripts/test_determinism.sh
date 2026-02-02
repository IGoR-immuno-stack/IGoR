#!/bin/bash
set -e

# Determinism test for IGoR
# Generates sequences twice with the same seed and checks if they are identical.

IGOR_BIN="./build/bin/igor"
WD1="val_1"
WD2="val_2"
SEED=12345
COUNT=1000

echo "Running determinism test with seed $SEED..."

# Clean up
rm -rf $WD1 $WD2

# Get absolute path for working directories
mkdir -p $WD1 $WD2
ABS_WD1=$(pwd)/$WD1
ABS_WD2=$(pwd)/$WD2

# Run from build/bin to satisfy hardcoded relative paths for models
cd build/bin

# First run
./igor -set_wd $ABS_WD1 -batch test -species human -chain beta -generate $COUNT --seed $SEED > /dev/null 2>&1

# Second run
./igor -set_wd $ABS_WD2 -batch test -species human -chain beta -generate $COUNT --seed $SEED > /dev/null 2>&1

cd ../..

# Compare outputs
FILE1="$WD1/test_generated/generated_seqs_werr.csv"
FILE2="$WD2/test_generated/generated_seqs_werr.csv"

if diff "$FILE1" "$FILE2" > /dev/null; then
    echo "SUCCESS: Outputs are identical."
    rm -rf $WD1 $WD2
    exit 0
else
    echo "FAILURE: Outputs differ!"
    diff "$FILE1" "$FILE2" | head -n 20
    exit 1
fi

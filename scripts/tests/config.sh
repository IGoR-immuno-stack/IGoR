# The binary the tracks run. Set IGORBIN to run them with another one, as the A/B oracle does
# (ab_oracle.sh runs every track once with the base binary and once with the candidate).
IGORBIN="${IGORBIN:-$SCRIPT_DIR/../../build/bin/igor}"
TESTINPUT="$SCRIPT_DIR/data/input"
TESTREF="$SCRIPT_DIR/data/reference"

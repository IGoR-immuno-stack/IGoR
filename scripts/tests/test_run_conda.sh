#!/usr/bin/env bash
# Runs the regression tracks (test_run.sh) with the igor of a built conda package instead of
# build/bin/igor: the binary that ships, built with the package's compiler and flags, and run
# from its install prefix, where it finds its data directory and its libraries.
#
#     pixi build --output-dir=bundle
#     pixi run test_regression_conda [test_run.sh options] [TEST_SPEC]
#
# BUNDLE names the directory that holds the igor_vdj .conda file (default: bundle/). The package
# is installed with pixi exec, from a local channel plus the workspace's channels, in an
# environment of its own. The tracks run in the calling environment: only the binary changes.
# Linux and macOS.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUNDLE="${BUNDLE:-$ROOT/bundle}"

shopt -s nullglob
packages=("$BUNDLE"/igor_vdj-*.conda)
shopt -u nullglob
if [ ${#packages[@]} -ne 1 ]; then
    echo "Expected one igor_vdj .conda file in $BUNDLE, found ${#packages[@]}." >&2
    echo "Run: pixi build --output-dir=bundle" >&2
    exit 1
fi

# The channel subdirectory of this platform, which the package was built for.
case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) subdir=linux-64 ;;
    Darwin-arm64) subdir=osx-arm64 ;;
    Darwin-x86_64) subdir=osx-64 ;;
    *) echo "No conda platform known for $(uname -s)-$(uname -m)" >&2; exit 1 ;;
esac

# A fixed path: pixi exec keys its cached environments on the channel URL, so a fresh temporary
# directory would leave one environment behind per run.
channel="${TMPDIR:-/tmp}/igor_test_run_conda_channel"
rm -rf "$channel"
trap 'rm -rf "$channel"' EXIT
mkdir -p "$channel/$subdir" "$channel/noarch"
cp "${packages[0]}" "$channel/$subdir/"
pixi exec rattler-index fs "$channel" >/dev/null

# The channels after the local one are the workspace's (pixi.toml), for the run dependencies.
# --force-reinstall: a rebuilt package keeps its file name, so the cached environment would
# still hold the package of an earlier build.
IGORBIN="$(pixi exec --force-reinstall -c "file://$channel" -c igor-immuno -c conda-forge \
    -s igor_vdj -- sh -c 'command -v igor')"
export IGORBIN
echo "Regression tracks with $IGORBIN ($("$IGORBIN" --version))"

"$ROOT/scripts/tests/test_run.sh" "$@"

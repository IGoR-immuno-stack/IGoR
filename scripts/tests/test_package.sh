#!/usr/bin/env bash
# Checks the installed CMake package (cmake/igorConfig.cmake.in) from a consumer project.
#
# Installs the existing build into a scratch prefix, then configures, builds and runs
# tst/package against it once per scenario. Run inside the pixi environment:
#
#     pixi run test_package [build-dir]
#
# The build directory must already be built; this script does not build IGoR.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${1:-$ROOT/build}"
WORK="$BUILD/package_test"
PREFIX="$WORK/prefix"

if [ ! -f "$BUILD/CMakeCache.txt" ]; then
    echo "No build found in $BUILD. Run: pixi run build" >&2
    exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK"
cmake --install "$BUILD" --prefix "$PREFIX" > "$WORK/install.log" 2>&1 \
    || { cat "$WORK/install.log"; echo "FAILED: install" >&2; exit 1; }

# The consumer is built like IGoR was, so that the runtime libraries match on Windows.
BUILD_TYPE="$(grep '^CMAKE_BUILD_TYPE:' "$BUILD/CMakeCache.txt" | cut -d= -f2)"
BUILD_TYPE_ARG=()
if [ -n "$BUILD_TYPE" ]; then
    BUILD_TYPE_ARG=("-DCMAKE_BUILD_TYPE=$BUILD_TYPE")
fi

# Windows finds the IGoR DLLs through PATH; elsewhere the build rpath is enough.
export PATH="$PREFIX/bin:$PREFIX/lib:$PATH"

failed=0
for scenario in all model core streaming umbrella unknown; do
    dir="$WORK/$scenario"
    log="$WORK/$scenario.log"
    if cmake -S "$ROOT/tst/package" -B "$dir" -G Ninja "${BUILD_TYPE_ARG[@]}" \
             -DCMAKE_PREFIX_PATH="$PREFIX" -DSCENARIO="$scenario" > "$log" 2>&1 \
       && cmake --build "$dir" >> "$log" 2>&1 \
       && { [ "$scenario" = unknown ] || "$dir/consumer" >> "$log" 2>&1; }; then
        echo "ok      $scenario"
    else
        echo "FAILED  $scenario (log: $log)"
        tail -n 25 "$log"
        failed=1
    fi
done

exit $failed

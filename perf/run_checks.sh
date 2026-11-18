#!/bin/sh
# Builds and runs, in one go:
#   1. the repository's existing unit tests (unittest/*.cpp, unmodified)
#   2. the ParallelCDT benchmark / consistency checker
# No build-system files are touched; everything is compiled directly.
set -e

REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${TMPDIR:-/tmp}/p2t_checks
mkdir -p "$OUT"

LIB_SRC=$(ls "$REPO"/poly2tri/common/*.cc "$REPO"/poly2tri/sweep/*.cc)

echo "=== [1/2] existing repository unit tests ==="
# Boost is not installed and may not be added as a dependency, so the
# repo's own test sources are compiled against the local compatibility
# shim in perf/shim (test sources themselves are unmodified).
g++ -std=c++14 -O2 -Wall \
    -I"$REPO" -I"$REPO/perf/shim" \
    -DP2T_BASE_DIR="\"$REPO\"" \
    "$REPO"/unittest/main.cpp "$REPO"/unittest/TriangleTest.cpp \
    $LIB_SRC \
    -o "$OUT/test_poly2tri"
"$OUT/test_poly2tri"

echo
echo "=== [2/2] ParallelCDT benchmark / consistency checks ==="
g++ -std=c++14 -O2 -Wall -pthread \
    -I"$REPO" \
    "$REPO"/perf/parallel_cdt_bench.cpp \
    $LIB_SRC \
    -o "$OUT/parallel_cdt_bench"
"$OUT/parallel_cdt_bench"

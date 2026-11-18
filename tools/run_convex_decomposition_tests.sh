#!/usr/bin/env bash
# Builds and runs everything needed to validate the convex decomposition:
#   1. the existing Boost.Test suite plus the new decomposition tests,
#   2. the four-scene demo (pieces, vertices, area error, timing),
#   3. a fingerprint check proving the CDT triangulation output is unchanged
#      relative to the pristine HEAD sources,
#   4. the repository's own CMake/ctest suite.
# Build scripts are not modified; sources are compiled directly.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CXX="${CXX:-g++}"
CMAKE="${CMAKE:-cmake}"
BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT

BOOST_PREFIX="${BOOST_PREFIX:-/home/ys/boost-prefix/usr}"
BOOST_INC="$BOOST_PREFIX/include"
BOOST_LIB="$BOOST_PREFIX/lib/x86_64-linux-gnu"
COMMON_FLAGS=(-std=c++14 -O2 -Wall -Wextra -I"$ROOT")
if [ -d "$BOOST_INC" ]; then
  COMMON_FLAGS+=(-I"$BOOST_INC")
fi

echo "== 1. building unit tests (existing + convex decomposition) =="
"$CXX" "${COMMON_FLAGS[@]}" \
  "$ROOT/unittest/main.cpp" \
  "$ROOT/unittest/TriangleTest.cpp" \
  "$ROOT/unittest/ConvexDecompositionTest.cpp" \
  "$ROOT/poly2tri/common/shapes.cc" \
  "$ROOT"/poly2tri/sweep/*.cc \
  -isystem "$BOOST_INC" -L"$BOOST_LIB" -lboost_unit_test_framework -lboost_filesystem \
  -o "$BUILD/test_poly2tri"

echo "== running unit tests =="
LD_LIBRARY_PATH="$BOOST_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  "$BUILD/test_poly2tri" --log_level=test_suite 2>&1 | grep -E 'Entering test case|No errors|failure|error'

echo
echo "== 2. building and running the demo =="
"$CXX" "${COMMON_FLAGS[@]}" \
  "$ROOT/tools/convex_decomposition_demo.cpp" \
  "$ROOT/poly2tri/common/shapes.cc" \
  "$ROOT"/poly2tri/sweep/*.cc \
  -o "$BUILD/convex_decomposition_demo"
"$BUILD/convex_decomposition_demo" --with-triangles

echo
echo "== 3. triangulation fingerprint: pristine HEAD vs current =="
PRISTINE="$BUILD/pristine"
mkdir -p "$PRISTINE"
git -C "$ROOT" archive HEAD | tar -x -C "$PRISTINE"
"$CXX" -std=c++14 -O2 -DP2T_PRISTINE_FINGERPRINT -I"$PRISTINE" \
  "$ROOT/tools/convex_decomposition_demo.cpp" \
  "$PRISTINE/poly2tri/common/shapes.cc" \
  "$PRISTINE"/poly2tri/sweep/*.cc \
  -o "$BUILD/demo_pristine"
"$CXX" "${COMMON_FLAGS[@]}" \
  "$ROOT/tools/convex_decomposition_demo.cpp" \
  "$ROOT/poly2tri/common/shapes.cc" \
  "$ROOT"/poly2tri/sweep/*.cc \
  -o "$BUILD/demo_current"
"$BUILD/demo_pristine" --fingerprint > "$BUILD/fp_pristine.txt"
"$BUILD/demo_current" --fingerprint > "$BUILD/fp_current.txt"
cat "$BUILD/fp_current.txt"
if ! diff -u "$BUILD/fp_pristine.txt" "$BUILD/fp_current.txt"; then
  echo "TRIANGULATION OUTPUT CHANGED" >&2
  exit 1
fi
echo "TRIANGULATION OUTPUT UNCHANGED"

echo
echo "== 4. repository CMake/ctest suite =="
"$CMAKE" -S "$ROOT" -B "$BUILD/cmake" -GNinja -DP2T_BUILD_TESTS=ON \
  -DBOOST_ROOT="$BOOST_PREFIX" -DBoost_NO_SYSTEM_PATHS=ON -Wno-dev > /dev/null
"$CMAKE" --build "$BUILD/cmake" > /dev/null
(cd "$BUILD/cmake" && LD_LIBRARY_PATH="$BOOST_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  ctest --output-on-failure)

echo
echo "ALL CHECKS PASSED"

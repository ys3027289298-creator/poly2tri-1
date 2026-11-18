#!/usr/bin/env bash
# Builds and runs the size-field refinement demo/test.
#
# The demo lives outside the library build on purpose: it uses only the
# headers and sources already in the repository (no third party libraries,
# no network) and does not touch CMakeLists.txt or meson.build.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
build="${TMPDIR:-/tmp}/poly2tri-size-field-demo"

mkdir -p "$build"
"${CXX:-c++}" -std=c++14 -O2 -Wall -Wextra \
  -I "$root" \
  -DP2T_BASE_DIR="\"$root\"" \
  "$here/main.cpp" \
  "$root"/poly2tri/common/*.cc \
  "$root"/poly2tri/sweep/*.cc \
  -o "$build/size_field_demo"

exec "$build/size_field_demo"

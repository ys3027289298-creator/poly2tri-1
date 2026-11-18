#!/usr/bin/env bash
# Build and run the quality-driven remeshing demo in one shot.
# Uses only the repo sources and the STL; no third-party dependency is added.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/remesh_demo/build"
mkdir -p "${BUILD_DIR}"

LIB_SOURCES=()
for f in "${ROOT}"/poly2tri/common/*.cc "${ROOT}"/poly2tri/sweep/*.cc; do
  LIB_SOURCES+=("${f}")
done

g++ -std=c++14 -O2 -Wall -Wextra -DP2T_STATIC_EXPORTS \
  -I"${ROOT}" \
  "${ROOT}/remesh_demo/remesh_demo.cpp" \
  "${LIB_SOURCES[@]}" \
  -o "${BUILD_DIR}/remesh_demo"

"${BUILD_DIR}/remesh_demo"

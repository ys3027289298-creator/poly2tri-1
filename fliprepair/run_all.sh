#!/usr/bin/env bash
# One-shot runner: builds and runs (1) the repo's existing unit tests and
# (2) the flip-repair scenarios. No network, no third-party downloads.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$ROOT/fliprepair/build"
mkdir -p "$BUILD"

LIB_SOURCES=("$ROOT"/poly2tri/common/*.cc "$ROOT"/poly2tri/sweep/*.cc)

echo "==================================================================="
echo "1/2 仓库现有单元测试 (unittest/, 经自带 shim 编译, 测试代码未改动)"
echo "==================================================================="
g++ -std=c++14 -Wall -Wextra \
    -I"$ROOT" -I"$ROOT/fliprepair/boost_shim" \
    -DP2T_BASE_DIR="\"$ROOT\"" \
    "${LIB_SOURCES[@]}" \
    "$ROOT/unittest/main.cpp" "$ROOT/unittest/TriangleTest.cpp" \
    -o "$BUILD/test_poly2tri"
"$BUILD/test_poly2tri"

echo
echo "==================================================================="
echo "2/2 翻边修复四组用例 (fliprepair/flip_repair_scenarios.cpp)"
echo "==================================================================="
g++ -std=c++14 -Wall -Wextra \
    -I"$ROOT" -I"$ROOT/fliprepair" \
    "${LIB_SOURCES[@]}" \
    "$ROOT/fliprepair/flip_repair_scenarios.cpp" \
    -o "$BUILD/flip_repair_scenarios"
"$BUILD/flip_repair_scenarios"

echo
echo "全部完成: 现有测试与翻边修复用例均通过。"

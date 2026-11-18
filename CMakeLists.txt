cmake_minimum_required(VERSION 3.12)

project(poly2tri LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)

option(P2T_BUILD_TESTS "Build tests" OFF)
option(P2T_BUILD_TESTBED "Build the testbed application" OFF)
option(P2T_BUILD_BENCH "Build the winding normalization benchmark" OFF)

file(GLOB SOURCES poly2tri/common/*.cc poly2tri/sweep/*.cc)
file(GLOB HEADERS poly2tri/*.h poly2tri/common/*.h poly2tri/sweep/*.h)
add_library(poly2tri ${SOURCES} ${HEADERS})
target_include_directories(poly2tri INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

get_target_property(poly2tri_target_type poly2tri TYPE)
if(poly2tri_target_type STREQUAL SHARED_LIBRARY)
  target_compile_definitions(poly2tri PRIVATE P2T_SHARED_EXPORTS)
else()
  target_compile_definitions(poly2tri PUBLIC P2T_STATIC_EXPORTS)
endif()

if(P2T_BUILD_TESTS)
    enable_testing()
    add_subdirectory(unittest)
endif()

if(P2T_BUILD_BENCH)
    add_executable(winding_bench bench/winding_bench.cc)
    target_link_libraries(winding_bench PRIVATE poly2tri)
    # The benchmark is meaningless unoptimized; force -O2 in non-debug builds
    # without changing the build type of the library itself.
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(winding_bench PRIVATE -O2)
    endif()
endif()

if(P2T_BUILD_TESTBED)
    add_subdirectory(testbed)
endif()

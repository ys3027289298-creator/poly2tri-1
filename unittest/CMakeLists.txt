# Dependencies
if (WIN32)
    set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.69 REQUIRED COMPONENTS
    filesystem
    unit_test_framework
)

# Build Unit Tests
add_executable(test_poly2tri
    main.cpp
    TriangleTest.cpp
)

target_include_directories(test_poly2tri
    PRIVATE
    ${Boost_INCLUDE_DIRS}
)

target_compile_definitions(test_poly2tri
    PRIVATE
    P2T_BASE_DIR="${PROJECT_SOURCE_DIR}"
)

target_link_libraries(test_poly2tri
    PRIVATE
    poly2tri
    ${Boost_LIBRARIES}
)

add_test(NAME poly2tri COMMAND test_poly2tri)

# Zero-dependency winding order tests: no Boost required, so they can run in
# minimal build environments.
add_executable(test_winding winding_test.cpp)
target_link_libraries(test_winding PRIVATE poly2tri)
add_test(NAME winding COMMAND test_winding)

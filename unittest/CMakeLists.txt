# Build Unit Tests
find_package(Boost 1.69 COMPONENTS
    filesystem
    unit_test_framework
)

if(Boost_FOUND)
    add_executable(test_poly2tri
        main.cpp
        TriangleTest.cpp
        TaubinSmoothingTest.cpp
    )

    target_include_directories(test_poly2tri
        PRIVATE
        ${Boost_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}
    )

    target_compile_definitions(test_poly2tri
        PRIVATE
        P2T_BASE_DIR="${PROJECT_SOURCE_DIR}"
        MESH_HAVE_BOOST_TEST
    )

    target_link_libraries(test_poly2tri
        PRIVATE
        poly2tri
        ${Boost_LIBRARIES}
    )

    add_test(NAME poly2tri COMMAND test_poly2tri)
else()
    # No Boost available: build the self-contained Taubin test with its
    # minimal Boost.Test-compatible shim instead of skipping tests entirely.
    add_executable(test_taubin
        TaubinSmoothingTest.cpp
    )

    target_include_directories(test_taubin
        PRIVATE
        ${PROJECT_SOURCE_DIR}
    )

    add_test(NAME taubin COMMAND test_taubin)
endif()

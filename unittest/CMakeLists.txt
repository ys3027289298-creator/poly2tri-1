# Dependency-free exercise of self-intersection detection and resolution.
add_executable(test_self_intersection self_intersection_demo.cpp)
target_link_libraries(test_self_intersection PRIVATE poly2tri)
add_test(NAME self_intersection COMMAND test_self_intersection)

# Boost based regression tests are only enabled when Boost is available,
# keeping the default build runnable on minimal toolchains.
find_package(Boost 1.69 COMPONENTS
    filesystem
    unit_test_framework
)

if (Boost_FOUND)
    if (WIN32)
        set(Boost_USE_STATIC_LIBS ON)
    endif()

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
endif()

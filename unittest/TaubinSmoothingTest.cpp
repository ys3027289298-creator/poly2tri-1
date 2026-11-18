#ifdef MESH_HAVE_BOOST_TEST
#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE TaubinSmoothingTest
#include <boost/test/unit_test.hpp>
#else
#include "boost_test_shim.h"
#endif

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "mesh/generators.h"
#include "mesh/taubin.h"
#include "mesh/trimesh.h"
#include "mesh/vec3.h"

using namespace meshproc;

namespace {

const double PI = 3.14159265358979323846;
const double VOLUME_TOLERANCE = 0.01;   // required: volume change below 1%
const double FEATURE_THRESHOLD = 0.5;   // radians; folds over ~30 deg are feature edges
const int SMOOTH_ITERATIONS = 10;

void report(const std::string& name, double rmseNoisy, double rmseSmoothed,
            double volChange, double featureAngleChange,
            double minFeatureAngle, int featureEdges) {
    std::cout << "\n=== " << name << " ===\n" << std::setprecision(4)
              << "  RMSE noisy / smoothed : " << rmseNoisy << " / " << rmseSmoothed << "\n"
              << "  volume change         : " << (volChange * 100.0) << " %\n";
    if (featureEdges > 0) {
        std::cout << "  feature edges         : " << featureEdges << "\n"
                  << "  mean feature angle d  : " << (featureAngleChange * 180.0 / PI)
                  << " deg\n"
                  << "  min feature angle aft : " << (minFeatureAngle * 180.0 / PI)
                  << " deg (reference 90 deg)\n";
    }
}

double boundaryDisplacement(const TriMesh& mesh, const std::vector<Vec3>& reference) {
    double maxDisp = 0.0;
    for (int i = 0; i < mesh.vertexCount(); ++i) {
        if (mesh.isBoundary[i]) {
            maxDisp = std::max(maxDisp, (mesh.positions[i] - reference[i]).length());
        }
    }
    return maxDisp;
}

} // namespace

BOOST_AUTO_TEST_CASE(SphereVolumePreservationTest)
{
    // Closed mesh: sphere with high-frequency surface noise.
    TriMesh mesh = makeSphere(1.0, 40, 80);
    std::vector<Vec3> clean = mesh.positions;
    double volumeBefore = mesh.volume();
    addNoise(mesh, 0.02, 42);

    double rmseNoisy = rmse(mesh, clean);
    double volumeNoisy = mesh.volume();

    taubinSmooth(mesh, {SMOOTH_ITERATIONS, 0.5, -0.53});

    double rmseSmoothed = rmse(mesh, clean);
    double volChange = relativeVolumeChange(volumeBefore, mesh.volume());

    // Laplacian-only counterpart on the same noisy input.
    TriMesh shrunk = makeSphere(1.0, 40, 80);
    addNoise(shrunk, 0.02, 42);
    laplacianSmooth(shrunk, SMOOTH_ITERATIONS, 0.5);
    double laplacianVolChange = relativeVolumeChange(volumeBefore, shrunk.volume());

    report("Sphere (closed)", rmseNoisy, rmseSmoothed, volChange, 0.0, 0.0, 0);
    std::cout << "  Laplacian vol change  : " << (laplacianVolChange * 100.0) << " %\n";

    BOOST_CHECK(volChange < VOLUME_TOLERANCE);
    BOOST_CHECK(rmseSmoothed < rmseNoisy);
    // Taubin must decisively outperform plain Laplacian on volume retention.
    BOOST_CHECK(volChange < laplacianVolChange);
    (void)volumeNoisy;
}

BOOST_AUTO_TEST_CASE(PlaneBoundaryPinnedTest)
{
    // Open surface: a planar grid whose boundary must remain exactly fixed.
    TriMesh mesh = makePlane(1.0, 40, 40);
    std::vector<Vec3> clean = mesh.positions;
    BOOST_REQUIRE(!mesh.boundaryEdges.empty());
    addNoise(mesh, 0.05, 7);

    double rmseNoisy = rmse(mesh, clean);
    taubinSmooth(mesh, {SMOOTH_ITERATIONS, 0.5, -0.53});
    double rmseSmoothed = rmse(mesh, clean);

    report("Plane (open boundary)", rmseNoisy, rmseSmoothed, 0.0, 0.0, 0.0, 0);
    std::cout << "  max boundary move     : "
              << boundaryDisplacement(mesh, clean) << "\n";

    BOOST_CHECK_SMALL(boundaryDisplacement(mesh, clean), 1e-12);
    BOOST_CHECK(rmseSmoothed < rmseNoisy);
}

BOOST_AUTO_TEST_CASE(HemisphereBoundaryTest)
{
    // Open curved surface: hemisphere with an open circular boundary.
    TriMesh mesh = makeHemisphere(1.0, 30, 60);
    std::vector<Vec3> clean = mesh.positions;
    BOOST_REQUIRE(!mesh.boundaryEdges.empty());
    double volumeBefore = mesh.volume();
    addNoise(mesh, 0.02, 99);

    double rmseNoisy = rmse(mesh, clean);
    taubinSmooth(mesh, {SMOOTH_ITERATIONS, 0.5, -0.53});
    double rmseSmoothed = rmse(mesh, clean);
    double volChange = relativeVolumeChange(volumeBefore, mesh.volume());

    report("Hemisphere (open curved)", rmseNoisy, rmseSmoothed, volChange, 0.0, 0.0, 0);
    std::cout << "  max boundary move     : "
              << boundaryDisplacement(mesh, clean) << "\n";

    BOOST_CHECK_SMALL(boundaryDisplacement(mesh, clean), 1e-12);
    BOOST_CHECK(rmseSmoothed < rmseNoisy);
    BOOST_CHECK(volChange < VOLUME_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(SharpFeaturePreservationTest)
{
    // Subdivided cube: true 90 degree dihedral edges must survive smoothing
    // while the surface noise on the flat faces is removed.
    TriMesh mesh = makeSubdividedCube(1.0, 12);
    std::vector<Vec3> clean = mesh.positions;
    double volumeBefore = mesh.volume();
    auto anglesClean = dihedralAngles(mesh);

    addNoise(mesh, 0.01, 123);
    auto anglesNoisy = dihedralAngles(mesh);
    double rmseNoisy = rmse(mesh, clean);

    SmoothParameters params;
    params.iterations = SMOOTH_ITERATIONS;
    params.preserveFeatures = true;
    params.featureAngleThreshold = FEATURE_THRESHOLD;
    // Detect features from the clean geometry so noise cannot hide edges,
    // and restore the exact clean feature-line positions.
    params.referenceAngles = &anglesClean;
    params.referencePositions = &clean;
    taubinSmooth(mesh, params);
    auto anglesSmoothed = dihedralAngles(mesh);
    double rmseSmoothed = rmse(mesh, clean);
    double volChange = relativeVolumeChange(volumeBefore, mesh.volume());

    auto feature = featureEdgePreservation(anglesClean, anglesSmoothed, FEATURE_THRESHOLD);
    report("Cube (sharp features)", rmseNoisy, rmseSmoothed, volChange,
           feature.meanAngleChangeRad, feature.minFeatureAngleRad,
           feature.featureEdgeCount);

    // The detector must find exactly the 12 cube edges, each subdivided
    // into n=12 triangle edges: 144 unique feature edges.
    BOOST_REQUIRE_EQUAL(static_cast<int>(mesh.featureEdges.size()), 12 * 12);
    BOOST_REQUIRE_EQUAL(feature.featureEdgeCount, 12 * 12);

    // Laplacian-only comparison on the same input.
    TriMesh rounded = makeSubdividedCube(1.0, 12);
    addNoise(rounded, 0.01, 123);
    laplacianSmooth(rounded, SMOOTH_ITERATIONS, 0.5);
    auto anglesLaplacian = dihedralAngles(rounded);
    auto lapFeature = featureEdgePreservation(anglesClean, anglesLaplacian, FEATURE_THRESHOLD);
    std::cout << "  Laplacian min feat ang: "
              << (lapFeature.minFeatureAngleRad * 180.0 / PI) << " deg\n";

    BOOST_REQUIRE(feature.featureEdgeCount > 0);
    // Pinned to clean reference positions, the folds keep their 90 degree
    // angle: mean drift below 2 degrees and no feature edge below 85 degrees.
    BOOST_CHECK(feature.meanAngleChangeRad < 2.0 * PI / 180.0);
    BOOST_CHECK(feature.minFeatureAngleRad > 85.0 * PI / 180.0);
    BOOST_CHECK(feature.meanAngleChangeRad < lapFeature.meanAngleChangeRad);
    BOOST_CHECK(volChange < VOLUME_TOLERANCE);
    BOOST_CHECK(rmseSmoothed < rmseNoisy);
    (void)anglesNoisy;
}

BOOST_AUTO_TEST_CASE(PerformanceLargeMeshTest)
{
    // ~100k vertices, 10 iterations must complete inside 5 seconds.
    TriMesh mesh = makeSphere(1.0, 224, 448);
    BOOST_REQUIRE_GE(mesh.vertexCount(), 100000);
    addNoise(mesh, 0.02, 2024);

    auto start = std::chrono::steady_clock::now();
    taubinSmooth(mesh, {SMOOTH_ITERATIONS, 0.5, -0.53});
    auto end = std::chrono::steady_clock::now();

    double seconds = std::chrono::duration<double>(end - start).count();
    std::cout << "\n=== Performance ===\n"
              << "  vertices / faces      : " << mesh.vertexCount()
              << " / " << mesh.faceCount() << "\n"
              << "  10 Taubin iterations  : " << seconds << " s\n";
    BOOST_CHECK(seconds < 5.0);
}

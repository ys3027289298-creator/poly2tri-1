#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif

#include <poly2tri/poly2tri.h>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace {

// Owns points and triangles so the raw pointers handed to
// p2t::ExtractBoundaryLoops stay valid for the duration of a test.
struct MeshBuilder {
  std::deque<p2t::Point> points;
  std::vector<std::unique_ptr<p2t::Triangle>> triangles;

  std::size_t AddPoint(double x, double y)
  {
    points.emplace_back(x, y);
    return points.size() - 1;
  }

  void AddTriangle(std::size_t a, std::size_t b, std::size_t c)
  {
    triangles.emplace_back(new p2t::Triangle(points[a], points[b], points[c]));
  }

  // Two CCW triangles filling the quad a-b-c-d.
  void AddQuad(std::size_t a, std::size_t b, std::size_t c, std::size_t d)
  {
    AddTriangle(a, b, c);
    AddTriangle(a, c, d);
  }

  std::vector<p2t::Triangle*> Triangles()
  {
    std::vector<p2t::Triangle*> result;
    result.reserve(triangles.size());
    for (const auto& t : triangles) {
      result.push_back(t.get());
    }
    return result;
  }
};

const p2t::BoundaryLoop& OnlyLoopAtDepth(const p2t::BoundaryLoops& loops,
                                         int depth)
{
  const p2t::BoundaryLoop* found = nullptr;
  for (const auto& loop : loops.loops) {
    if (loop.depth == depth) {
      BOOST_REQUIRE_MESSAGE(found == nullptr,
                            "more than one loop at depth " << depth);
      found = &loop;
    }
  }
  BOOST_REQUIRE_MESSAGE(found != nullptr, "no loop at depth " << depth);
  return *found;
}

} // namespace

// Mesh 1: a single quad. Hand-annotated: 1 loop, 4 points, depth 0.
BOOST_AUTO_TEST_CASE(BoundaryLoopsQuadTest)
{
  MeshBuilder mesh;
  const auto a = mesh.AddPoint(0, 0);
  const auto b = mesh.AddPoint(4, 0);
  const auto c = mesh.AddPoint(4, 3);
  const auto d = mesh.AddPoint(0, 3);
  mesh.AddQuad(a, b, c, d);

  const auto loops = p2t::ExtractBoundaryLoops(mesh.Triangles());
  BOOST_REQUIRE_EQUAL(loops.loops.size(), 1);
  BOOST_REQUIRE_EQUAL(loops.OuterLoops().size(), 1);
  const auto& loop = loops.loops[0];
  BOOST_CHECK_EQUAL(loop.points.size(), 4);
  BOOST_CHECK_EQUAL(loop.depth, 0);
  BOOST_CHECK_EQUAL(loop.parent, -1);
  BOOST_CHECK(!loop.is_hole);
  BOOST_CHECK_CLOSE(loop.signed_area, 12.0, 1e-9);
}

// Mesh 2: CDT of a square with a square hole.
// Hand-annotated: 2 loops of 4 points each, depths 0 and 1.
BOOST_AUTO_TEST_CASE(BoundaryLoopsCdtHoleTest)
{
  std::vector<p2t::Point*> polyline{
    new p2t::Point(0, 0),
    new p2t::Point(10, 0),
    new p2t::Point(10, 10),
    new p2t::Point(0, 10),
  };
  std::vector<p2t::Point*> hole{
    new p2t::Point(3, 3),
    new p2t::Point(7, 3),
    new p2t::Point(7, 7),
    new p2t::Point(3, 7),
  };
  p2t::CDT cdt{ polyline };
  cdt.AddHole(hole);
  cdt.Triangulate();

  const auto loops = p2t::ExtractBoundaryLoops(cdt.GetTriangles());
  BOOST_REQUIRE_EQUAL(loops.loops.size(), 2);

  const auto& outer = OnlyLoopAtDepth(loops, 0);
  BOOST_CHECK_EQUAL(outer.points.size(), 4);
  BOOST_CHECK(!outer.is_hole);
  BOOST_CHECK_EQUAL(outer.parent, -1);
  BOOST_CHECK_CLOSE(outer.signed_area, 100.0, 1e-9);

  const auto& inner = OnlyLoopAtDepth(loops, 1);
  BOOST_CHECK_EQUAL(inner.points.size(), 4);
  BOOST_CHECK(inner.is_hole);
  BOOST_CHECK_CLOSE(inner.signed_area, -16.0, 1e-9);
  BOOST_CHECK_EQUAL(loops.loops[inner.parent].depth, 0);

  for (const auto p : polyline) {
    delete p;
  }
  for (const auto p : hole) {
    delete p;
  }
}

// Mesh 3: square ring (outer + hole) plus an island inside the hole.
// Hand-annotated: 3 loops of 4 points each, depths 0, 1 and 2.
BOOST_AUTO_TEST_CASE(BoundaryLoopsNestedIslandTest)
{
  MeshBuilder mesh;
  const auto o0 = mesh.AddPoint(0, 0);
  const auto o1 = mesh.AddPoint(10, 0);
  const auto o2 = mesh.AddPoint(10, 10);
  const auto o3 = mesh.AddPoint(0, 10);
  const auto h0 = mesh.AddPoint(3, 2);
  const auto h1 = mesh.AddPoint(8, 3);
  const auto h2 = mesh.AddPoint(7, 8);
  const auto h3 = mesh.AddPoint(2, 7);
  // Ring between the outer square and the hole (all triangles CCW).
  mesh.AddTriangle(o0, o1, h1);
  mesh.AddTriangle(o0, h1, h0);
  mesh.AddTriangle(o1, o2, h2);
  mesh.AddTriangle(o1, h2, h1);
  mesh.AddTriangle(o2, o3, h3);
  mesh.AddTriangle(o2, h3, h2);
  mesh.AddTriangle(o3, o0, h0);
  mesh.AddTriangle(o3, h0, h3);
  // Island inside the hole.
  const auto i0 = mesh.AddPoint(4, 4);
  const auto i1 = mesh.AddPoint(6, 4);
  const auto i2 = mesh.AddPoint(6, 6);
  const auto i3 = mesh.AddPoint(4, 6);
  mesh.AddQuad(i0, i1, i2, i3);

  const auto loops = p2t::ExtractBoundaryLoops(mesh.Triangles());
  BOOST_REQUIRE_EQUAL(loops.loops.size(), 3);
  BOOST_REQUIRE_EQUAL(loops.OuterLoops().size(), 1);

  const auto& outer = OnlyLoopAtDepth(loops, 0);
  BOOST_CHECK_EQUAL(outer.points.size(), 4);
  BOOST_CHECK(!outer.is_hole);
  BOOST_CHECK_CLOSE(outer.signed_area, 100.0, 1e-9);

  const auto& hole = OnlyLoopAtDepth(loops, 1);
  BOOST_CHECK_EQUAL(hole.points.size(), 4);
  BOOST_CHECK(hole.is_hole);
  BOOST_CHECK_CLOSE(hole.signed_area, -26.0, 1e-9);
  BOOST_CHECK_EQUAL(loops.loops[hole.parent].depth, 0);

  const auto& island = OnlyLoopAtDepth(loops, 2);
  BOOST_CHECK_EQUAL(island.points.size(), 4);
  BOOST_CHECK(!island.is_hole);
  BOOST_CHECK_CLOSE(island.signed_area, 4.0, 1e-9);
  BOOST_CHECK_EQUAL(loops.loops[island.parent].depth, 1);
}

// Mesh 4: two disconnected components, a quad and a square ring.
// Hand-annotated: 3 loops of 4 points each, depths 0, 0 and 1.
BOOST_AUTO_TEST_CASE(BoundaryLoopsTwoComponentsTest)
{
  MeshBuilder mesh;
  // Component 1: plain quad.
  const auto a0 = mesh.AddPoint(0, 0);
  const auto a1 = mesh.AddPoint(4, 0);
  const auto a2 = mesh.AddPoint(4, 4);
  const auto a3 = mesh.AddPoint(0, 4);
  mesh.AddQuad(a0, a1, a2, a3);
  // Component 2: square ring (outer square + hole).
  const auto o0 = mesh.AddPoint(10, 0);
  const auto o1 = mesh.AddPoint(20, 0);
  const auto o2 = mesh.AddPoint(20, 10);
  const auto o3 = mesh.AddPoint(10, 10);
  const auto h0 = mesh.AddPoint(13, 2);
  const auto h1 = mesh.AddPoint(18, 3);
  const auto h2 = mesh.AddPoint(17, 8);
  const auto h3 = mesh.AddPoint(12, 7);
  mesh.AddTriangle(o0, o1, h1);
  mesh.AddTriangle(o0, h1, h0);
  mesh.AddTriangle(o1, o2, h2);
  mesh.AddTriangle(o1, h2, h1);
  mesh.AddTriangle(o2, o3, h3);
  mesh.AddTriangle(o2, h3, h2);
  mesh.AddTriangle(o3, o0, h0);
  mesh.AddTriangle(o3, h0, h3);

  const auto loops = p2t::ExtractBoundaryLoops(mesh.Triangles());
  BOOST_REQUIRE_EQUAL(loops.loops.size(), 3);
  BOOST_CHECK_EQUAL(loops.OuterLoops().size(), 2);

  int num_depth0 = 0, num_depth1 = 0;
  for (const auto& loop : loops.loops) {
    BOOST_CHECK_EQUAL(loop.points.size(), 4);
    if (loop.depth == 0) {
      ++num_depth0;
      BOOST_CHECK(!loop.is_hole);
      BOOST_CHECK_EQUAL(loop.parent, -1);
    } else if (loop.depth == 1) {
      ++num_depth1;
      BOOST_CHECK(loop.is_hole);
      BOOST_CHECK_CLOSE(loop.signed_area, -26.0, 1e-9);
      BOOST_CHECK_EQUAL(loops.loops[loop.parent].depth, 0);
    } else {
      BOOST_FAIL("unexpected loop depth " << loop.depth);
    }
  }
  BOOST_CHECK_EQUAL(num_depth0, 2);
  BOOST_CHECK_EQUAL(num_depth1, 1);
}

// Three triangles sharing one edge must be rejected, and the error must
// identify the offending half-edge (triangle and corner).
BOOST_AUTO_TEST_CASE(BoundaryLoopsNonManifoldEdgeTest)
{
  MeshBuilder mesh;
  mesh.AddPoint(0, 0);  // 0
  mesh.AddPoint(1, 0);  // 1
  mesh.AddPoint(0, 1);  // 2
  mesh.AddPoint(0, -1); // 3
  mesh.AddPoint(1, 1);  // 4
  mesh.AddTriangle(0, 1, 2); // half-edge 0 is 0->1
  mesh.AddTriangle(1, 0, 3); // half-edge 3 is 1->0
  mesh.AddTriangle(0, 1, 4); // half-edge 6 is 0->1: third user of edge (0,1)

  try {
    p2t::ExtractBoundaryLoops(mesh.Triangles());
    BOOST_FAIL("expected NonManifoldEdgeError");
  } catch (const p2t::NonManifoldEdgeError& e) {
    BOOST_CHECK_EQUAL(e.half_edge, 0);
    BOOST_CHECK_EQUAL(e.triangle, 0);
    BOOST_CHECK_EQUAL(e.corner, 0);
    const std::string what = e.what();
    BOOST_CHECK(what.find("half-edge 0") != std::string::npos);
    BOOST_CHECK(what.find("triangle 0") != std::string::npos);
  }
}

// Two triangles winding the same edge in the same direction are also
// non-manifold.
BOOST_AUTO_TEST_CASE(BoundaryLoopsSameWindingEdgeTest)
{
  MeshBuilder mesh;
  mesh.AddPoint(0, 0); // 0
  mesh.AddPoint(1, 0); // 1
  mesh.AddPoint(0, 1); // 2
  mesh.AddPoint(2, 1); // 3
  mesh.AddTriangle(0, 1, 2); // half-edge 0 is 0->1
  mesh.AddTriangle(0, 1, 3); // half-edge 3 is 0->1 again

  try {
    p2t::ExtractBoundaryLoops(mesh.Triangles());
    BOOST_FAIL("expected NonManifoldEdgeError");
  } catch (const p2t::NonManifoldEdgeError& e) {
    BOOST_CHECK_EQUAL(e.half_edge, 0);
    BOOST_CHECK_EQUAL(e.triangle, 0);
    BOOST_CHECK_EQUAL(e.corner, 0);
  }
}

// 224x224 cells = 100,352 faces: extraction must take less than 150 ms.
BOOST_AUTO_TEST_CASE(BoundaryLoopsPerformanceTest)
{
  constexpr int n = 224;
  MeshBuilder mesh;
  for (int j = 0; j <= n; ++j) {
    for (int i = 0; i <= n; ++i) {
      mesh.AddPoint(i, j);
    }
  }
  const auto vertex = [n](int i, int j) {
    return static_cast<std::size_t>(j) * (n + 1) + i;
  };
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      mesh.AddQuad(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1),
                   vertex(i, j + 1));
    }
  }
  const auto triangles = mesh.Triangles();
  BOOST_REQUIRE_EQUAL(triangles.size(), 100352);

  const auto start = std::chrono::steady_clock::now();
  const auto loops = p2t::ExtractBoundaryLoops(triangles);
  const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                           std::chrono::steady_clock::now() - start)
                           .count();
  BOOST_TEST_MESSAGE("ExtractBoundaryLoops on " << triangles.size()
                                                << " faces took " << elapsed
                                                << " ms");
  BOOST_CHECK_LT(elapsed, 150);

  BOOST_REQUIRE_EQUAL(loops.loops.size(), 1);
  const auto& loop = loops.loops[0];
  BOOST_CHECK_EQUAL(loop.points.size(), 4 * n);
  BOOST_CHECK_EQUAL(loop.depth, 0);
  BOOST_CHECK(!loop.is_hole);
  BOOST_CHECK_CLOSE(loop.signed_area,
                    static_cast<double>(n) * static_cast<double>(n), 1e-9);
}

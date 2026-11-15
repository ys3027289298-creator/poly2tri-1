#include <poly2tri/poly2tri.h>
#include <poly2tri/sweep/refinement.h>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

using Vertex = std::pair<double, double>;
using TriVerts = std::array<Vertex, 3>;

TriVerts SortedVerts(p2t::Triangle* t)
{
  TriVerts verts{{
      {t->GetPoint(0)->x, t->GetPoint(0)->y},
      {t->GetPoint(1)->x, t->GetPoint(1)->y},
      {t->GetPoint(2)->x, t->GetPoint(2)->y},
  }};
  std::sort(verts.begin(), verts.end());
  return verts;
}

std::vector<TriVerts> SortedMesh(const std::vector<p2t::Triangle*>& triangles)
{
  std::vector<TriVerts> mesh;
  mesh.reserve(triangles.size());
  for (p2t::Triangle* t : triangles) {
    mesh.push_back(SortedVerts(t));
  }
  std::sort(mesh.begin(), mesh.end());
  return mesh;
}

double MaxEdge(const std::vector<p2t::Triangle*>& triangles)
{
  double max_length = 0.0;
  for (p2t::Triangle* t : triangles) {
    for (int i = 0; i < 3; i++) {
      const p2t::Point* a = t->GetPoint(i);
      const p2t::Point* b = t->GetPoint((i + 1) % 3);
      max_length = std::max(max_length, std::hypot(a->x - b->x, a->y - b->y));
    }
  }
  return max_length;
}

void DeletePoints(std::vector<p2t::Point*>& points)
{
  for (p2t::Point* p : points) {
    delete p;
  }
  points.clear();
}

} // namespace

// With refinement disabled the output must be vertex-by-vertex identical
// to a plain CDT triangulation.
BOOST_AUTO_TEST_CASE(RefinementDisabledMatchesCdt)
{
  std::vector<p2t::Point*> polyline{
    new p2t::Point(0, 0),
    new p2t::Point(10, 0),
    new p2t::Point(10, 10),
    new p2t::Point(5, 4),
    new p2t::Point(0, 10),
  };
  p2t::Point* steiner = new p2t::Point(5, 8);

  p2t::CDT cdt{ polyline };
  cdt.AddPoint(steiner);
  cdt.Triangulate();
  const std::vector<TriVerts> expected = SortedMesh(cdt.GetTriangles());

  p2t::Refiner refiner{ polyline };
  refiner.AddPoint(steiner);
  p2t::RefinementCriteria criteria;
  criteria.enabled = false;
  criteria.max_edge_length = 0.5; // ignored while disabled
  criteria.min_area = 1000.0;     // ignored while disabled
  const p2t::RefinementReport& report = refiner.Triangulate(criteria);
  const std::vector<TriVerts> actual = SortedMesh(refiner.GetTriangles());

  BOOST_CHECK(report.converged);
  BOOST_CHECK_EQUAL(report.iterations, 0);
  BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
  BOOST_CHECK_EQUAL(report.triangle_count, expected.size());
  for (std::size_t i = 0; i < expected.size(); i++) {
    BOOST_CHECK(actual[i] == expected[i]);
  }
  BOOST_CHECK_EQUAL(report.max_edge_length, MaxEdge(refiner.GetTriangles()));

  DeletePoints(polyline);
  delete steiner;
}

// Max-edge-length constraint: a coarse 10x10 square is refined until every
// triangle edge is at most 1.0.
BOOST_AUTO_TEST_CASE(RefinementMaxEdgeLengthConverges)
{
  std::vector<p2t::Point*> polyline{
    new p2t::Point(0, 0),
    new p2t::Point(10, 0),
    new p2t::Point(10, 10),
    new p2t::Point(0, 10),
  };

  p2t::CDT cdt{ polyline };
  cdt.Triangulate();
  const std::size_t coarse_count = cdt.GetTriangles().size();
  const double coarse_max_edge = MaxEdge(cdt.GetTriangles());

  p2t::Refiner refiner{ polyline };
  p2t::RefinementCriteria criteria;
  criteria.max_edge_length = 1.0;
  criteria.max_iterations = 64;
  const p2t::RefinementReport& report = refiner.Triangulate(criteria);
  const std::vector<p2t::Triangle*> triangles = refiner.GetTriangles();

  BOOST_CHECK(report.converged);
  BOOST_CHECK_GT(report.iterations, 0);
  BOOST_CHECK_EQUAL(report.triangle_count, triangles.size());
  BOOST_CHECK_GT(triangles.size(), coarse_count);
  BOOST_CHECK_LT(report.max_edge_length, coarse_max_edge);
  BOOST_CHECK_LE(report.max_edge_length, 1.0 + 1e-12);
  BOOST_CHECK_LE(MaxEdge(triangles), 1.0 + 1e-12);
  BOOST_CHECK_EQUAL(report.max_edge_length, MaxEdge(triangles));

  DeletePoints(polyline);
}

// Min-area constraint: an already satisfied constraint needs no iteration,
// an impossible one stops at the iteration cap and reports non-convergence.
BOOST_AUTO_TEST_CASE(RefinementMinAreaConstraint)
{
  std::vector<p2t::Point*> polyline{
    new p2t::Point(0, 0),
    new p2t::Point(10, 0),
    new p2t::Point(10, 10),
    new p2t::Point(0, 10),
  };

  {
    p2t::Refiner refiner{ polyline };
    p2t::RefinementCriteria criteria;
    criteria.min_area = 1e-9; // already satisfied
    const p2t::RefinementReport& report = refiner.Triangulate(criteria);
    BOOST_CHECK(report.converged);
    BOOST_CHECK_EQUAL(report.iterations, 0);
    BOOST_CHECK_EQUAL(report.triangle_count, 2);
  }
  {
    p2t::Refiner refiner{ polyline };
    p2t::RefinementCriteria criteria;
    criteria.min_area = 200.0; // larger than the whole 10x10 domain
    criteria.max_iterations = 3;
    const p2t::RefinementReport& report = refiner.Triangulate(criteria);
    BOOST_CHECK(!report.converged);
    BOOST_CHECK_EQUAL(report.iterations, 3);
    BOOST_CHECK_EQUAL(report.triangle_count, refiner.GetTriangles().size());
  }

  DeletePoints(polyline);
}

// A constrained edge crossed by a triangle edge must be detected and
// reported, never silently accepted.
BOOST_AUTO_TEST_CASE(ConstraintCrossingIsDetected)
{
  std::vector<p2t::Point*> polyline{
    new p2t::Point(0, 0),
    new p2t::Point(1, 0),
    new p2t::Point(1, 1),
    new p2t::Point(0, 1),
  };
  p2t::CDT cdt{ polyline };
  cdt.Triangulate();
  const std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();
  BOOST_REQUIRE_EQUAL(triangles.size(), 2);

  // The genuine boundary is respected.
  std::vector<std::vector<p2t::Point*>> genuine{ polyline };
  BOOST_CHECK_NO_THROW(p2t::CheckConstraintsIntact(triangles, genuine));

  // A bogus "constrained" edge cutting through the interior diagonal of
  // the mesh must be reported.
  p2t::Point a(0.5, -1.0), b(0.5, 2.0), c(-1.0, 0.5);
  std::vector<std::vector<p2t::Point*>> bogus{ { &a, &b, &c } };
  BOOST_CHECK_THROW(p2t::CheckConstraintsIntact(triangles, bogus),
                    std::runtime_error);

  DeletePoints(polyline);
}

// Refinement with a hole: constraints stay intact and the hole boundary is
// refined together with the interior.
BOOST_AUTO_TEST_CASE(RefinementWithHoleKeepsConstraints)
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
  const std::size_t coarse_count = cdt.GetTriangles().size();

  p2t::Refiner refiner{ polyline };
  refiner.AddHole(hole);
  p2t::RefinementCriteria criteria;
  criteria.max_edge_length = 1.0;
  criteria.max_iterations = 64;
  const p2t::RefinementReport& report = refiner.Triangulate(criteria);
  const std::vector<p2t::Triangle*> triangles = refiner.GetTriangles();

  BOOST_CHECK(report.converged);
  BOOST_CHECK_GT(triangles.size(), coarse_count);
  BOOST_CHECK_LE(report.max_edge_length, 1.0 + 1e-12);
  // Boundary and hole edges remain uncrossed.
  std::vector<std::vector<p2t::Point*>> polylines{ polyline, hole };
  BOOST_CHECK_NO_THROW(p2t::CheckConstraintsIntact(triangles, polylines));
  // No triangle vertex lands inside the hole.
  for (p2t::Triangle* t : triangles) {
    for (int i = 0; i < 3; i++) {
      const p2t::Point* p = t->GetPoint(i);
      const bool inside_hole = p->x > 3.0 + 1e-12 && p->x < 7.0 - 1e-12 &&
                               p->y > 3.0 + 1e-12 && p->y < 7.0 - 1e-12;
      BOOST_CHECK(!inside_hole);
    }
  }

  DeletePoints(polyline);
  DeletePoints(hole);
}

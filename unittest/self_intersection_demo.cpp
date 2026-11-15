// Self-contained, dependency-free exercise for the self-intersection
// detection and even-odd resolution policy. It is intentionally small and
// deterministic so it can run together with the regular unit tests without
// Boost or any additional library.

#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

double TriangleArea(p2t::Triangle& triangle)
{
  const p2t::Point* a = triangle.GetPoint(0);
  const p2t::Point* b = triangle.GetPoint(1);
  const p2t::Point* c = triangle.GetPoint(2);
  return 0.5 * std::fabs((b->x - a->x) * (c->y - a->y)
                         - (b->y - a->y) * (c->x - a->x));
}

double TotalArea(const std::vector<p2t::Triangle*>& triangles)
{
  double total = 0.0;
  for (p2t::Triangle* triangle : triangles) {
    total += TriangleArea(*triangle);
  }
  return total;
}

double Shoelace(const std::vector<p2t::Point*>& points)
{
  double twice_area = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    const p2t::Point& a = *points[i];
    const p2t::Point& b = *points[(i + 1) % points.size()];
    twice_area += a.x * b.y - b.x * a.y;
  }
  return 0.5 * twice_area;
}

int g_failures = 0;

void Check(bool condition, const char* message)
{
  std::printf("  [%s] %s\n", condition ? "PASS" : "FAIL", message);
  if (!condition) {
    ++g_failures;
  }
}

std::vector<p2t::Point*> Bowtie()
{
  return {
    new p2t::Point(0.0, 0.0), new p2t::Point(1.0, 1.0),
    new p2t::Point(1.0, 0.0), new p2t::Point(0.0, 1.0)
  };
}

std::vector<p2t::Point*> PentagonStar()
{
  std::vector<p2t::Point*> points;
  for (int i = 0; i < 5; ++i) {
    const double angle = -M_PI / 2.0 + i * 4.0 * M_PI / 5.0;
    points.push_back(new p2t::Point(std::cos(angle), std::sin(angle)));
  }
  return points;
}

std::vector<p2t::Point*> Square()
{
  return {
    new p2t::Point(0.0, 0.0), new p2t::Point(1.0, 0.0),
    new p2t::Point(1.0, 1.0), new p2t::Point(0.0, 1.0)
  };
}

void FreePoints(std::vector<p2t::Point*>& points)
{
  for (p2t::Point* point : points) {
    delete point;
  }
}

void RejectPolicyTest()
{
  std::printf("Reject policy on self-intersecting inputs\n");
  auto bowtie = Bowtie();
  try {
    p2t::CDT cdt(bowtie);
    cdt.Triangulate();
    Check(false, "bowtie must be rejected");
  } catch (const p2t::SelfIntersectionError& error) {
    Check(error.code() == p2t::SI_PROPER_CROSSING,
          "bowtie reports SI_PROPER_CROSSING");
  }
  FreePoints(bowtie);

  auto star = PentagonStar();
  try {
    p2t::CDT cdt(star);
    cdt.Triangulate();
    Check(false, "pentagram must be rejected");
  } catch (const p2t::SelfIntersectionError& error) {
    Check(error.code() == p2t::SI_PROPER_CROSSING,
          "pentagram reports SI_PROPER_CROSSING");
  }
  FreePoints(star);

  std::vector<p2t::Point*> duplicate{
    new p2t::Point(0.0, 0.0), new p2t::Point(1.0, 0.0),
    new p2t::Point(1.0, 1.0), new p2t::Point(1.0, 1.0),
    new p2t::Point(0.0, 1.0)
  };
  try {
    p2t::CDT cdt(duplicate);
    cdt.Triangulate();
    Check(false, "duplicate point must be rejected");
  } catch (const p2t::SelfIntersectionError& error) {
    Check(error.code() == p2t::SI_DUPLICATE_POINT,
          "duplicate point reports SI_DUPLICATE_POINT");
  }
  FreePoints(duplicate);
}

void EvenOddPolicyTest()
{
  std::printf("Even-odd policy resolves self-intersecting inputs\n");
  auto bowtie = Bowtie();
  p2t::EvenOddResolver reference(bowtie);
  reference.Resolve();
  {
    p2t::CDT cdt(bowtie);
    cdt.SetSelfIntersectionPolicy(p2t::RESOLVE_EVEN_ODD);
    cdt.Triangulate();
    const auto triangles = cdt.GetTriangles();
    const double area = TotalArea(triangles);
    const double expected = reference.FilledArea();
    std::printf("  bowtie: triangles=%zu area=%.6f expected=%.6f\n",
                triangles.size(), area, expected);
    Check(triangles.size() == 2, "bowtie yields two triangles");
    Check(std::fabs(area - expected) <= 1e-9 * std::max(1.0, expected),
          "bowtie area conserved within tolerance");
  }
  FreePoints(bowtie);

  auto star = PentagonStar();
  p2t::EvenOddResolver star_reference(star);
  star_reference.Resolve();
  {
    p2t::CDT cdt(star);
    cdt.SetSelfIntersectionPolicy(p2t::RESOLVE_EVEN_ODD);
    cdt.Triangulate();
    const auto triangles = cdt.GetTriangles();
    const double area = TotalArea(triangles);
    const double expected = star_reference.FilledArea();
    std::printf("  pentagram: triangles=%zu area=%.6f expected=%.6f "
                "shoelace=%.6f\n",
                triangles.size(), area, expected,
                std::fabs(Shoelace(star)));
    Check(triangles.size() == 5, "pentagram yields five arm triangles");
    Check(std::fabs(area - expected) <= 1e-9 * std::max(1.0, expected),
          "pentagram area conserved within tolerance");
  }
  FreePoints(star);
}

void SimplePolygonUnchangedTest()
{
  std::printf("Simple polygons keep identical triangulation\n");
  auto square = Square();
  p2t::CDT cdt(square);
  cdt.Triangulate();
  const auto triangles = cdt.GetTriangles();
  Check(triangles.size() == 2, "square yields two triangles");
  Check(std::fabs(TotalArea(triangles) - 1.0) < 1e-12,
        "square area equals the input area");
  bool uses_original_points = true;
  for (p2t::Triangle* triangle : triangles) {
    for (int vertex = 0; vertex < 3; ++vertex) {
      const p2t::Point* point = triangle->GetPoint(vertex);
      const bool known = std::find(square.begin(), square.end(), point)
        != square.end();
      uses_original_points = uses_original_points && known;
    }
  }
  Check(uses_original_points,
        "triangle vertices are the original Point objects");
  FreePoints(square);
}

} // namespace

int main()
{
  RejectPolicyTest();
  EvenOddPolicyTest();
  SimplePolygonUnchangedTest();
  std::printf(g_failures == 0 ? "ALL CHECKS PASSED\n" : "CHECKS FAILED: %d\n",
              g_failures);
  return g_failures == 0 ? 0 : 1;
}

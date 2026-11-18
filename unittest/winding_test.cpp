// Zero-dependency tests for the winding order helpers in
// poly2tri/common/winding.h. A minimal assertion framework stands in for the
// Boost-based suite, which is not available in every build environment.

#include <poly2tri/common/winding.h>
#include <poly2tri/poly2tri.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

int g_failures = 0;
int g_checks = 0;
std::string g_current_case;

void Check(bool condition, const char* expression, const char* file, int line)
{
  ++g_checks;
  if (!condition) {
    ++g_failures;
    std::printf("  FAIL %s:%d: %s\n", file, line, expression);
  }
}

void CheckClose(double actual, double expected, double eps, const char* file, int line)
{
  ++g_checks;
  if (std::fabs(actual - expected) > eps) {
    ++g_failures;
    std::printf("  FAIL %s:%d: |%g - %g| > %g\n", file, line, actual, expected, eps);
  }
}

#define CHECK(cond) Check((cond), #cond, __FILE__, __LINE__)
#define CHECK_CLOSE(actual, expected, eps) CheckClose((actual), (expected), (eps), __FILE__, __LINE__)

// Square split into two triangles.
p2t::Point MakePoint(double x, double y)
{
  return p2t::Point(x, y);
}

void MixedWindingIsNormalizedToCCW()
{
  // Vertices of a unit square; tri 0 is CCW, tri 1 deliberately wound CW.
  std::vector<p2t::Point> points{
    MakePoint(0, 0), MakePoint(1, 0), MakePoint(1, 1), MakePoint(0, 1),
  };
  std::vector<int> triangles{
    0, 1, 2, // CCW: area +0.5
    0, 3, 2, // CW:  area -0.5
  };

  p2t::WindingStats before = p2t::AnalyzeWinding(triangles, points);
  CHECK(before.face_count == 2);
  CHECK(before.ccw_faces == 1);
  CHECK(before.cw_faces == 1);
  CHECK(before.OverallOrientation() == p2t::Winding::Degenerate);

  const std::size_t flipped = p2t::NormalizeWinding(triangles, points);
  CHECK(flipped == 1);

  p2t::WindingStats after = p2t::VerifyWinding(triangles, points);
  CHECK(after.ccw_faces == 2);
  CHECK(after.cw_faces == 0);
  CHECK(after.IsManifoldConsistent());
  CHECK(after.inconsistent_edges == 0);
  // Four boundary edges plus the shared diagonal.
  CHECK(after.boundary_edges == 4);
}

void NormalizeToClockwise()
{
  std::vector<p2t::Point> points{
    MakePoint(0, 0), MakePoint(1, 0), MakePoint(0, 1),
  };
  std::vector<int> triangles{ 0, 1, 2 }; // CCW

  const std::size_t flipped =
      p2t::NormalizeWinding(triangles, points, p2t::Winding::CW);
  CHECK(flipped == 1);

  p2t::WindingStats stats = p2t::VerifyWinding(triangles, points);
  CHECK(stats.cw_faces == 1);
  CHECK(stats.ccw_faces == 0);
  CHECK(stats.OverallOrientation() == p2t::Winding::CW);
  CHECK(stats.IsManifoldConsistent());
}

void DegenerateTrianglesAreLeftUntouched()
{
  std::vector<p2t::Point> points{
    MakePoint(0, 0), MakePoint(1, 0), MakePoint(2, 0), MakePoint(0, 1),
  };
  // First triangle is collinear (zero area), the second a normal CCW one.
  // The degenerate face traverses shared edge (0,1) backwards, matching the
  // opposite-direction rule of a consistent manifold.
  std::vector<int> triangles{
    1, 0, 2,
    0, 1, 3,
  };

  p2t::WindingStats before = p2t::AnalyzeWinding(triangles, points);
  CHECK(before.degenerate_faces == 1);
  CHECK(before.ccw_faces == 1);

  const std::size_t flipped = p2t::NormalizeWinding(triangles, points);
  CHECK(flipped == 0);

  p2t::WindingStats after = p2t::VerifyWinding(triangles, points);
  CHECK(after.degenerate_faces == 1);
  // The collinear edge (0,1) is shared, but its two directed occurrences
  // oppose each other, so it is not flagged as inconsistent.
  CHECK(after.inconsistent_edges == 0);
}

void InconsistentSharedEdgeIsDetected()
{
  // A valid manifold pair of CCW triangles sharing diagonal (0,2). Flipping
  // one face reverses its traversal of the diagonal and must be reported.
  std::vector<p2t::Point> points{
    MakePoint(0, 0), MakePoint(1, 0), MakePoint(1, 1), MakePoint(0, 1),
  };
  std::vector<int> triangles{
    0, 1, 2, // CCW, traverses diagonal 2 -> 0
    0, 2, 3, // CCW, traverses diagonal 0 -> 2
  };

  p2t::WindingStats healthy = p2t::VerifyWinding(triangles, points);
  CHECK(healthy.IsManifoldConsistent());

  // Break the winding by reversing the second face.
  p2t::FlipTriangle(&triangles[3]);
  p2t::WindingStats broken = p2t::VerifyWinding(triangles, points);
  CHECK(broken.inconsistent_edges == 1);
  CHECK(!broken.IsManifoldConsistent());

  // Per-face CCW normalization restores the shared-edge consistency.
  p2t::NormalizeWinding(triangles, points);
  p2t::WindingStats fixed = p2t::VerifyWinding(triangles, points);
  CHECK(fixed.IsManifoldConsistent());
  CHECK(fixed.ccw_faces == 2);
}

void NonManifoldEdgeIsReported()
{
  // Three triangles all sharing edge (0,1).
  std::vector<p2t::Point> points{
    MakePoint(0, 0), MakePoint(1, 0), MakePoint(1, 1),
    MakePoint(0, 1), MakePoint(-1, 1),
  };
  std::vector<int> triangles{
    0, 1, 2,
    0, 1, 3,
    4, 0, 1, // contains (0,1) as well
  };

  p2t::WindingStats stats = p2t::VerifyWinding(triangles, points);
  CHECK(stats.non_manifold_edges == 1);
  CHECK(!stats.IsManifoldConsistent());
}

void AreaIsPreservedWithinEpsilon()
{
  // Random-ish coordinates far from unit scale to exercise double precision.
  std::vector<p2t::Point> points{
    MakePoint(1234.5678, -8765.4321),
    MakePoint(2345.6789, -7654.3210),
    MakePoint(3456.7891, -6543.2109),
    MakePoint(4567.8912, -5432.1098),
  };
  std::vector<int> original{
    0, 1, 2,
    0, 3, 1, // deliberately CW
  };
  std::vector<int> triangles = original;

  const p2t::WindingStats before = p2t::AnalyzeWinding(triangles, points);
  p2t::NormalizeWinding(triangles, points);
  const p2t::WindingStats after = p2t::VerifyWinding(triangles, points);

  // Index swaps preserve each double area exactly; only its sign changes.
  CHECK_CLOSE(after.unsigned_double_area, before.unsigned_double_area, 0.0);
  CHECK_CLOSE(std::fabs(after.signed_double_area), before.unsigned_double_area, 1e-9);
  CHECK_CLOSE(after.signed_double_area - after.unsigned_double_area, 0.0, 1e-9);

  // Hand-computed unsigned double area of the two triangles.
  const double expected =
      std::fabs(p2t::SignedDoubleArea(points[0], points[1], points[2])) +
      std::fabs(p2t::SignedDoubleArea(points[0], points[3], points[1]));
  CHECK_CLOSE(after.unsigned_double_area, expected, 1e-9);
}

void GridMeshEndToEndThroughCDT()
{
  // Triangulate a polygon with CDT, deliberately scramble face winding, then
  // run the full pipeline and confirm every shared edge is consistent.
  const int n = 8;
  std::vector<p2t::Point*> polyline;
  std::vector<p2t::Point> storage;
  storage.reserve(4 * n);
  for (int i = 0; i <= n; ++i) {
    storage.emplace_back(static_cast<double>(i), 0.0);
  }
  for (int i = 1; i <= n; ++i) {
    storage.emplace_back(static_cast<double>(n), static_cast<double>(i));
  }
  for (int i = n - 1; i >= 0; --i) {
    storage.emplace_back(static_cast<double>(i), static_cast<double>(n));
  }
  for (int i = n - 1; i >= 1; --i) {
    storage.emplace_back(0.0, static_cast<double>(i));
  }
  for (auto& p : storage) {
    polyline.push_back(&p);
  }

  p2t::CDT cdt{ polyline };
  cdt.Triangulate();
  const std::vector<p2t::Triangle*> tris = cdt.GetTriangles();
  // Boundary-only polygon with V vertices: T = 2V - B - 2, and B == V here,
  // so T = V - 2 (32 boundary vertices -> 30 triangles).
  CHECK(tris.size() == polyline.size() - 2);

  std::vector<p2t::Point> points;
  std::vector<int> triangles;
  p2t::ExtractTriangles(tris, points, triangles);
  CHECK(points.size() == polyline.size());
  CHECK(triangles.size() == tris.size() * 3);

  // Scramble: reverse every second face.
  for (std::size_t f = 0; f < tris.size(); f += 2) {
    p2t::FlipTriangle(&triangles[3 * f]);
  }
  p2t::WindingStats scrambled = p2t::VerifyWinding(triangles, points);
  CHECK(scrambled.inconsistent_edges > 0);

  p2t::WindingStats healed =
      p2t::NormalizeAndVerifyWinding(triangles, points);
  CHECK(healed.IsManifoldConsistent());
  CHECK(healed.ccw_faces == tris.size());
  CHECK(healed.cw_faces == 0);
  // Square of side n has double area 2 * n^2.
  CHECK_CLOSE(healed.unsigned_double_area, 2.0 * n * n, 1e-9);
}

void RunCase(const char* name, void (*fn)())
{
  g_current_case = name;
  std::printf("[ RUN      ] %s\n", name);
  fn();
  std::printf("[       OK ] %s\n", name);
}

} // namespace

int main()
{
  RunCase("MixedWindingIsNormalizedToCCW", MixedWindingIsNormalizedToCCW);
  RunCase("NormalizeToClockwise", NormalizeToClockwise);
  RunCase("DegenerateTrianglesAreLeftUntouched", DegenerateTrianglesAreLeftUntouched);
  RunCase("InconsistentSharedEdgeIsDetected", InconsistentSharedEdgeIsDetected);
  RunCase("NonManifoldEdgeIsReported", NonManifoldEdgeIsReported);
  RunCase("AreaIsPreservedWithinEpsilon", AreaIsPreservedWithinEpsilon);
  RunCase("GridMeshEndToEndThroughCDT", GridMeshEndToEndThroughCDT);

  std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
  return g_failures == 0 ? 0 : 1;
}

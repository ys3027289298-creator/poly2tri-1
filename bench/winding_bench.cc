// Benchmark for the winding normalization pipeline:
//
//   AnalyzeWinding -> NormalizeWinding -> VerifyWinding
//
// on a structured triangle grid of roughly 100k faces with deliberately
// scrambled winding. Reports per-stage wall-clock time and checks the area
// error against the exact grid area. Build standalone with:
//
//   g++ -std=c++14 -O2 -I.. bench/winding_bench.cc -o winding_bench

#include <poly2tri/common/winding.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

double NowSeconds()
{
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

// Builds a gridx by gridy quad grid, two triangles per cell: 2*gridx*gridy
// faces and (gridx+1)*(gridy+1) vertices. Every third face is reversed to
// simulate an export with mixed winding.
void BuildScrambledGrid(int gridx,
                        int gridy,
                        std::vector<p2t::Point>& points,
                        std::vector<int>& triangles)
{
  points.clear();
  triangles.clear();
  points.reserve(static_cast<std::size_t>((gridx + 1) * (gridy + 1)));
  triangles.reserve(static_cast<std::size_t>(2 * gridx * gridy) * 3);

  for (int j = 0; j <= gridy; ++j) {
    for (int i = 0; i <= gridx; ++i) {
      points.emplace_back(static_cast<double>(i), static_cast<double>(j));
    }
  }
  auto vid = [gridx](int i, int j) { return j * (gridx + 1) + i; };

  for (int j = 0; j < gridy; ++j) {
    for (int i = 0; i < gridx; ++i) {
      const int a = vid(i, j);
      const int b = vid(i + 1, j);
      const int c = vid(i + 1, j + 1);
      const int d = vid(i, j + 1);
      triangles.insert(triangles.end(), { a, b, c });
      triangles.insert(triangles.end(), { a, c, d });
    }
  }

  // Deterministic scrambling: reverse every third face so a third of the
  // shared edges become inconsistent.
  for (std::size_t f = 0; f < triangles.size() / 3; f += 3) {
    p2t::FlipTriangle(&triangles[3 * f]);
  }
}

} // namespace

int main(int argc, char** argv)
{
  const int gridx = argc > 1 ? std::atoi(argv[1]) : 224;
  const int gridy = argc > 2 ? std::atoi(argv[2]) : 224;
  // Several iterations per run: a single scheduler hiccup must not decide
  // the timing result, so we report the median and the best observed time.
  const int iterations = argc > 3 ? std::atoi(argv[3]) : 10;
  const std::size_t face_count =
      static_cast<std::size_t>(2) * gridx * gridy;

  std::vector<p2t::Point> points;
  std::vector<int> reference;
  BuildScrambledGrid(gridx, gridy, points, reference);

  std::printf("grid %dx%d: %zu vertices, %zu faces, %.0f scrambled, %d iterations\n",
              gridx, gridy, points.size(), face_count, face_count / 3.0, iterations);

  std::vector<double> totals;
  totals.reserve(iterations);
  p2t::WindingStats verified;
  for (int it = 0; it < iterations; ++it) {
    // Restore the scrambled winding each round so every iteration does the
    // same amount of work.
    std::vector<int> triangles = reference;

    const double t0 = NowSeconds();
    p2t::AnalyzeWinding(triangles, points);
    p2t::NormalizeWinding(triangles, points);
    verified = p2t::VerifyWinding(triangles, points);
    totals.push_back((NowSeconds() - t0) * 1e3);
  }

  std::sort(totals.begin(), totals.end());
  const double best = totals.front();
  const double median = totals[totals.size() / 2];
  const double worst = totals.back();
  std::printf("  pipeline time: best=%6.2f ms  median=%6.2f ms  worst=%6.2f ms"
              "  (budget 100 ms)\n", best, median, worst);

  // One detailed run for the per-stage breakdown.
  std::vector<int> triangles = reference;
  double t0 = NowSeconds();
  p2t::WindingStats analysis = p2t::AnalyzeWinding(triangles, points);
  double t1 = NowSeconds();
  const std::size_t flipped = p2t::NormalizeWinding(triangles, points);
  double t2 = NowSeconds();
  verified = p2t::VerifyWinding(triangles, points);
  double t3 = NowSeconds();

  std::printf("  per-stage (last run): analyze=%.2f normalize=%.2f verify=%.2f ms\n",
              (t1 - t0) * 1e3, (t2 - t1) * 1e3, (t3 - t2) * 1e3);
  std::printf("pre-normalization: ccw=%zu cw=%zu deg=%zu overall=%s\n",
              analysis.ccw_faces, analysis.cw_faces, analysis.degenerate_faces,
              analysis.OverallOrientation() == p2t::Winding::CCW ? "CCW" :
              analysis.OverallOrientation() == p2t::Winding::CW ? "CW" : "degenerate");
  std::printf("  %zu faces flipped\n", flipped);

  std::printf("\nresults:\n");
  std::printf("  ccw=%zu cw=%zu deg=%zu\n", verified.ccw_faces, verified.cw_faces,
              verified.degenerate_faces);
  std::printf("  inconsistent=%zu non_manifold=%zu boundary=%zu\n",
              verified.inconsistent_edges, verified.non_manifold_edges,
              verified.boundary_edges);

  // Exact double area of the grid rectangle: 2 * gridx * gridy.
  const double expected_double_area = 2.0 * gridx * gridy;
  const double area_error =
      std::fabs(verified.unsigned_double_area - expected_double_area);
  std::printf("  unsigned double area=%.6f expected=%.6f error=%.3e (eps 1e-9)\n",
              verified.unsigned_double_area, expected_double_area, area_error);

  // The median is the honest figure on a loaded machine; best and worst are
  // reported alongside it so a single preemption stays visible.
  const bool timing_ok = median < 100.0;
  const bool ok = verified.IsManifoldConsistent() && verified.cw_faces == 0 &&
                  area_error <= 1e-9 && timing_ok;
  std::printf("\n%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}

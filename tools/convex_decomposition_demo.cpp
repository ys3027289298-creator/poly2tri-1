#include <poly2tri/poly2tri.h>
#include <chrono>
#include <cstring>
#include <exception>
#include <memory>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {
const double kPi = 3.14159265358979323846;

struct Scene {
  std::string name;
  std::vector<p2t::Point*> outer;
  std::vector<std::vector<p2t::Point*>> holes;
  std::vector<std::unique_ptr<p2t::Point>> keep;
  void add(std::vector<p2t::Point*>& dst, double x, double y)
  {
    keep.emplace_back(new p2t::Point(x, y));
    dst.push_back(keep.back().get());
  }
};

// Per-vertex fingerprint of a triangulation.  tools/run_convex_decomposition_tests.sh
// prints these for the pre-change library and for the current one and requires
// them to be identical, so any change to the triangulation output fails.
std::uint64_t mix(std::uint64_t h, double d)
{
  std::uint64_t bits;
  std::memcpy(&bits, &d, sizeof bits);
  return (h ^ (bits + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2)));
}

std::uint64_t fingerprint(const std::vector<p2t::Triangle*>& triangles)
{
  std::uint64_t h = 1469598103934665603ULL;
  h = mix(h, static_cast<double>(triangles.size()));
  for (auto* t : triangles)
    for (int i = 0; i < 3; ++i) {
      h = mix(h, t->GetPoint(i)->x);
      h = mix(h, t->GetPoint(i)->y);
    }
  return h;
}

double ring_area(const std::vector<p2t::Point*>& r)
{
  double a = 0;
  for (size_t i = 0; i < r.size(); ++i)
    a += r[i]->x * r[(i + 1) % r.size()]->y - r[(i + 1) % r.size()]->x * r[i]->y;
  return std::fabs(a * 0.5);
}

// The four groups required by the task:
//   convex, single reflex, multiple reflexes, one hole.
std::vector<Scene> make_scenes()
{
  std::vector<Scene> scenes;

  {
    Scene s;
    s.name = "1. convex (hexagon)";
    for (int i = 0; i < 6; ++i)
      s.add(s.outer, 2 * std::cos(i * kPi / 3), 2 * std::sin(i * kPi / 3));
    scenes.push_back(std::move(s));
  }
  {
    Scene s;
    s.name = "2. single reflex (L)";
    const double pts[][2] = { { 0, 0 }, { 4, 0 }, { 4, 1 }, { 1, 1 }, { 1, 4 }, { 0, 4 } };
    for (auto& p : pts)
      s.add(s.outer, p[0], p[1]);
    scenes.push_back(std::move(s));
  }
  {
    Scene s;
    s.name = "3. multiple reflexes (comb)";
    const double pts[][2] = { { 0, 0 },    { 4, 0 },    { 4, 3 },    { 3.2, 3 }, { 3.2, 1.4 },
                              { 2.4, 1.4 }, { 2.4, 3 }, { 1.6, 3 }, { 1.6, 1.4 }, { 0.8, 1.4 },
                              { 0.8, 3 },  { 0, 3 } };
    for (auto& p : pts)
      s.add(s.outer, p[0], p[1]);
    scenes.push_back(std::move(s));
  }
  {
    Scene s;
    s.name = "4. square with hole";
    const double outer[][2] = { { 0, 0 }, { 4, 0 }, { 4, 4 }, { 0, 4 } };
    for (auto& p : outer)
      s.add(s.outer, p[0], p[1]);
    s.holes.emplace_back();
    const double hole[][2] = { { 1.2, 1.2 }, { 2.8, 1.2 }, { 2.8, 2.8 }, { 1.2, 2.8 } };
    for (auto& p : hole)
      s.add(s.holes.back(), p[0], p[1]);
    scenes.push_back(std::move(s));
  }
  return scenes;
}
} // namespace

int main(int argc, char** argv)
{
  bool with_triangles = false;
  bool fingerprint_only = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--with-triangles")
      with_triangles = true;
    else if (arg == "--fingerprint")
      fingerprint_only = true;
  }
  const auto scenes = make_scenes();

  int failures = 0;
  if (!fingerprint_only)
    std::printf("%-27s %6s %6s %8s %7s %7s %8s %9s %10s\n", "polygon", "in_v", "reflex",
                "pieces", "out_v", "ears", "bridges", "area_err", "time_us");
  for (const auto& scene : scenes) {
#ifndef P2T_PRISTINE_FINGERPRINT
    if (fingerprint_only) {
      try {
        p2t::DecomposePolygon(scene.outer, scene.holes);
      } catch (const std::exception& e) {
        std::printf("%-27s ERROR: %s\n", scene.name.c_str(), e.what());
        ++failures;
      }
      continue;
    }
#endif
#ifdef P2T_PRISTINE_FINGERPRINT
    continue;
#else
    try {
      const auto start = std::chrono::steady_clock::now();
      const auto result = p2t::DecomposePolygon(scene.outer, scene.holes);
      const auto stop = std::chrono::steady_clock::now();
      const double us = std::chrono::duration<double, std::micro>(stop - start).count();

      double piece_sum = 0;
      for (const auto& piece : result.pieces)
        piece_sum += piece.Area();
      const double expected = ring_area(scene.outer) -
                              (scene.holes.empty() ? 0.0 : ring_area(scene.holes.front()));

      std::printf("%-27s %6zu %6zu %8zu %7zu %7zu %8zu %9.2e %10.2f\n", scene.name.c_str(),
                  result.input_vertices, result.reflex_vertices, result.pieces.size(),
                  result.output_vertices, result.ears_clipped, result.bridges, result.area_error,
                  us);
      if (std::fabs(piece_sum - expected) > 1e-9 * std::max(1.0, expected)) {
        std::printf("  ERROR: piece area %.12f != polygon area %.12f\n", piece_sum, expected);
        ++failures;
      }
    } catch (const std::exception& e) {
      std::printf("%-27s ERROR: %s\n", scene.name.c_str(), e.what());
      ++failures;
    }
#endif
  }

  if (with_triangles || fingerprint_only) {
    if (!fingerprint_only)
      std::printf("\nTriangulation output (must not change):\n");
    for (const auto& scene : scenes) {
      p2t::CDT cdt(scene.outer);
      for (const auto& hole : scene.holes)
        cdt.AddHole(hole);
      cdt.Triangulate();
      const auto triangles = cdt.GetTriangles();
      const std::uint64_t h = fingerprint(triangles);
      if (fingerprint_only) {
        std::printf("FP|%s|%zu|%016llx\n", scene.name.c_str(), triangles.size(),
                    static_cast<unsigned long long>(h));
      } else {
        std::printf("  %-27s triangles=%2zu fingerprint=%016llx\n", scene.name.c_str(),
                    triangles.size(), static_cast<unsigned long long>(h));
      }
    }
  }

  if (!fingerprint_only)
    std::printf("\n%s\n", failures == 0 ? "ALL FOUR GROUPS OK" : "FAILURES");
  return failures == 0 ? 0 : 1;
}

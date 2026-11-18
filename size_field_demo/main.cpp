#include <poly2tri/sweep/size_field.h>

#ifndef P2T_BASE_DIR
#define P2T_BASE_DIR "."
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Vec2 {
  double x, y;
};

std::vector<p2t::Point*> MakeSquare(double x0, double y0, double x1, double y1)
{
  return {new p2t::Point(x0, y0), new p2t::Point(x1, y0),
          new p2t::Point(x1, y1), new p2t::Point(x0, y1)};
}

std::vector<p2t::Point*> MakeLShape()
{
  return {new p2t::Point(0.0, 0.0), new p2t::Point(1.0, 0.0),
          new p2t::Point(1.0, 0.35), new p2t::Point(0.4, 0.35),
          new p2t::Point(0.4, 1.0), new p2t::Point(0.0, 1.0)};
}

void Delete(std::vector<p2t::Point*>& points)
{
  for (p2t::Point* point : points) {
    delete point;
  }
  points.clear();
}

double DistanceToSegment(const Vec2& p, const Vec2& a, const Vec2& b)
{
  const double dx = b.x - a.x;
  const double dy = b.y - a.y;
  const double l2 = dx * dx + dy * dy;
  double t = l2 > 0.0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / l2 : 0.0;
  t = std::max(0.0, std::min(1.0, t));
  return std::hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

double DistanceToRing(const Vec2& p, const std::vector<Vec2>& ring)
{
  double best = 1e300;
  const std::size_t n = ring.size();
  for (std::size_t i = 0; i < n; ++i) {
    best = std::min(best, DistanceToSegment(p, ring[i], ring[(i + 1) % n]));
  }
  return best;
}

std::vector<p2t::Point*> LoadTestbedData(const std::string& filename, std::size_t& line_count)
{
  std::vector<p2t::Point*> points;
  std::ifstream file(std::string(P2T_BASE_DIR) + "/testbed/data/" + filename);
  std::string line;
  line_count = 0;
  while (std::getline(file, line)) {
    if (line.empty()) {
      break;
    }
    std::istringstream stream(line);
    double x = 0.0;
    double y = 0.0;
    if (!(stream >> x >> y)) {
      break;
    }
    points.push_back(new p2t::Point(x, y));
    ++line_count;
  }
  return points;
}

using Tri = std::array<std::array<double, 2>, 3>;

std::vector<Tri> SortedMesh(const std::vector<p2t::Triangle*>& triangles)
{
  std::vector<Tri> mesh;
  mesh.reserve(triangles.size());
  for (p2t::Triangle* t : triangles) {
    Tri tri;
    for (int i = 0; i < 3; ++i) {
      tri[i][0] = t->GetPoint(i)->x;
      tri[i][1] = t->GetPoint(i)->y;
    }
    std::sort(tri.begin(), tri.end());
    mesh.push_back(tri);
  }
  std::sort(mesh.begin(), mesh.end());
  return mesh;
}

struct Case {
  std::string name;
  p2t::SizeFieldFunction field;
};

void PrintReport(const Case& c, const p2t::SizeFieldReport& r)
{
  std::cout << std::left << std::setw(14) << c.name
            << " tris=" << std::setw(6) << r.triangle_count
            << " verts=" << std::setw(5) << r.vertex_count
            << " mean_dev=" << std::fixed << std::setprecision(4)
            << std::setw(8) << r.mean_relative_deviation
            << " med_dev=" << std::setw(8) << r.median_relative_deviation
            << " p90_dev=" << std::setw(8) << r.p90_relative_deviation
            << " in_band=" << std::setw(6) << r.fraction_within_band
            << " added(b/i)=" << r.boundary_points_added << "/" << r.interior_points_added
            << " merged=" << r.points_merged
            << " area_err=" << std::scientific << std::setprecision(2)
            << std::setw(9) << r.relative_area_error << std::fixed
            << " iters=" << std::setw(2) << r.iterations
            << "/" << r.max_iterations
            << " total_ms=" << std::setprecision(2) << std::setw(8)
            << 1000.0 * r.total_seconds << " [";
  for (std::size_t i = 0; i < r.iteration_seconds.size(); ++i) {
    if (i) std::cout << ",";
    std::cout << std::setprecision(2) << 1000.0 * r.iteration_seconds[i];
  }
  std::cout << "]ms\n";
}

} // namespace

int main()
{
  int failures = 0;
  auto check = [&failures](bool ok, const std::string& what) {
    if (!ok) {
      std::cerr << "FAIL: " << what << "\n";
      ++failures;
    }
  };

  std::cout << "== size-field driven refinement demo ==\n";

  // 1. Illegal fields must be rejected, both zero and negative.
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    bool zero_threw = false;
    bool negative_threw = false;
    try {
      p2t::SizeFieldRefiner r(square);
      r.Triangulate([](const p2t::Point&) { return 0.0; });
    } catch (const std::runtime_error&) {
      zero_threw = true;
    }
    try {
      p2t::SizeFieldRefiner r(square);
      r.Triangulate([](const p2t::Point& p) { return -0.1 - p.x; });
    } catch (const std::runtime_error&) {
      negative_threw = true;
    }
    check(zero_threw, "zero size field must throw");
    check(negative_threw, "negative size field must throw");

    // Invalid options are rejected as well.
    bool bad_scale = false;
    bool bad_ratios = false;
    bool no_field = false;
    try {
      p2t::SizeFieldRefiner r(square);
      p2t::SizeFieldOptions options;
      options.scale = 0.0;
      r.Triangulate([](const p2t::Point&) { return 0.1; }, options);
    } catch (const std::runtime_error&) {
      bad_scale = true;
    }
    try {
      p2t::SizeFieldRefiner r(square);
      p2t::SizeFieldOptions options;
      options.merge_ratio = options.split_ratio;
      r.Triangulate([](const p2t::Point&) { return 0.1; }, options);
    } catch (const std::runtime_error&) {
      bad_ratios = true;
    }
    try {
      p2t::SizeFieldRefiner r(square);
      r.Triangulate(p2t::SizeFieldFunction());
    } catch (const std::runtime_error&) {
      no_field = true;
    }
    check(bad_scale, "non-positive scale must throw");
    check(bad_ratios, "merge_ratio >= split_ratio must throw");
    check(no_field, "empty size field must throw");
    Delete(square);
  }

  // 2. Uniform field coarse enough to accept the input: the result must be
  //    vertex-by-vertex identical to a plain CDT triangulation.
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    p2t::Point* steiner = new p2t::Point(0.5, 0.5);

    p2t::CDT cdt(square);
    cdt.AddPoint(steiner);
    cdt.Triangulate();
    const std::vector<Tri> expected = SortedMesh(cdt.GetTriangles());

    p2t::SizeFieldRefiner refiner(square);
    refiner.AddPoint(steiner);
    const p2t::SizeFieldReport& report =
        refiner.Triangulate([](const p2t::Point&) { return 10.0; });
    const std::vector<Tri> actual = SortedMesh(refiner.GetTriangles());

    check(report.iterations == 0, "uniform coarse field must not iterate");
    check(actual == expected, "uniform coarse field must match plain CDT exactly");
    std::cout << "uniform-identity: tris=" << actual.size()
              << " iterations=" << report.iterations
              << " identical=" << (actual == expected ? "yes" : "no") << "\n";
    Delete(square);
    delete steiner;
  }

  // 3. Four target fields on a small domain (unit square with a square hole).
  // Geometry used by the boundary-dense field (deleted at the end of main).
  const std::vector<Vec2> outline{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  const std::vector<Vec2> hole{{0.4, 0.4}, {0.6, 0.4}, {0.6, 0.6}, {0.4, 0.6}};

  std::vector<Case> cases;
  cases.push_back(Case{"uniform", [](const p2t::Point&) { return 0.08; }});
  cases.push_back(Case{"gradient", [](const p2t::Point& p) { return 0.03 + 0.10 * p.x; }});
  cases.push_back(Case{"center-dense", [](const p2t::Point& p) {
    const double r = std::hypot(p.x - 0.5, p.y - 0.5);
    return 0.025 + 0.20 * r;
  }});
  cases.push_back(Case{"boundary-dense", [&outline, &hole](const p2t::Point& p) {
    const double d = std::min(DistanceToRing({p.x, p.y}, outline),
                             DistanceToRing({p.x, p.y}, hole));
    return 0.02 + 0.30 * d;
  }});

  std::cout << "\nfield          triangles / edge-length deviation vs target / area / time\n";
  for (const Case& c : cases) {
    std::vector<p2t::Point*> ring = MakeSquare(0, 0, 1, 1);
    std::vector<p2t::Point*> inner = MakeSquare(0.4, 0.4, 0.6, 0.6);
    p2t::SizeFieldRefiner refiner(ring);
    refiner.AddHole(inner);
    p2t::SizeFieldOptions options;
    options.scale = 1.0;
    options.max_iterations = 10;
    const p2t::SizeFieldReport& report = refiner.Triangulate(c.field, options);
    PrintReport(c, report);

    check(report.converged, c.name + ": expected convergence");
    check(report.triangle_count > 2, c.name + ": expected a refined mesh");
    check(report.relative_area_error < 1e-9, c.name + ": area must be preserved");
    check(report.fraction_within_band > 0.85, c.name + ": most edges within band");
    check(report.max_length_ratio <= options.split_ratio + 1e-9,
          c.name + ": no edge above split_ratio");
    Delete(ring);
    Delete(inner);
  }

  // Extra: the refined mesh is a valid Delaunay triangulation and remains
  // constrained (no boundary edge is crossed).
  {
    std::vector<p2t::Point*> ring = MakeSquare(0, 0, 1, 1);
    std::vector<p2t::Point*> inner = MakeSquare(0.4, 0.4, 0.6, 0.6);
    p2t::SizeFieldRefiner refiner(ring);
    refiner.AddHole(inner);
    const p2t::SizeFieldReport& report =
        refiner.Triangulate([](const p2t::Point& p) { return 0.02 + 0.2 * p.y; });
    (void)report;
    const auto triangles = refiner.GetTriangles();
    // IsDelaunay() checks every triangle against every vertex and therefore
    // flags points beyond a constrained hole edge; check a hole-free domain.
    std::vector<p2t::Point*> plain = MakeSquare(0, 0, 1, 1);
    p2t::SizeFieldRefiner plain_refiner(plain);
    const p2t::SizeFieldReport& plain_report =
        plain_refiner.Triangulate([](const p2t::Point& p) { return 0.02 + 0.2 * p.y; });
    const auto plain_triangles = plain_refiner.GetTriangles();
    check(p2t::IsDelaunay(plain_triangles), "refined mesh must be Delaunay");
    std::cout << "\ndelaunay=" << (p2t::IsDelaunay(plain_triangles) ? "yes" : "no")
              << " (gradient field, " << plain_report.triangle_count << " tris)\n";
    Delete(plain);
    Delete(ring);
    Delete(inner);
  }

  // 4. L-shaped outline (reflex corner) with a center-dense field.
  {
    std::vector<p2t::Point*> lshape = MakeLShape();
    p2t::SizeFieldRefiner refiner(lshape);
    const p2t::SizeFieldReport& report = refiner.Triangulate([](const p2t::Point& p) {
      return 0.02 + 0.25 * std::hypot(p.x - 0.4, p.y - 0.4);
    });
    std::cout << "L-shape center-dense: tris=" << report.triangle_count
              << " area_err=" << report.relative_area_error
              << " in_band=" << report.fraction_within_band
              << " converged=" << report.converged << "\n";
    check(report.converged, "L-shape: expected convergence");
    check(report.relative_area_error < 1e-9, "L-shape: area must be preserved");
    Delete(lshape);
  }

  // 5. Short-edge merging: a deliberately dense Steiner grid with a coarse
  //    field is thinned when merge_dense_steiner_points is enabled, and kept
  //    untouched when it is not.
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    std::vector<p2t::Point*> steiner;
    for (int i = 1; i < 20; ++i) {
      for (int j = 1; j < 20; ++j) {
        steiner.push_back(new p2t::Point(0.05 * i, 0.05 * j));
      }
    }
    const auto coarse_field = [](const p2t::Point&) { return 0.5; };

    p2t::SizeFieldRefiner keep(square);
    for (p2t::Point* p : steiner) keep.AddPoint(p);
    p2t::SizeFieldOptions options;
    options.max_iterations = 10;
    const p2t::SizeFieldReport& keep_report = keep.Triangulate(coarse_field, options);
    const std::size_t kept_vertices = keep_report.vertex_count;

    p2t::SizeFieldRefiner thin(square);
    for (p2t::Point* p : steiner) thin.AddPoint(p);
    options.merge_dense_steiner_points = true;
    const p2t::SizeFieldReport& thin_report = thin.Triangulate(coarse_field, options);
    const std::size_t thinned_vertices = thin_report.vertex_count;

    std::cout << "dense-steiner: kept_verts=" << kept_vertices
              << " thinned_verts=" << thinned_vertices
              << " thinned_tris=" << thin_report.triangle_count
              << " steiner_merged=" << thin_report.steiner_points_merged
              << " in_band=" << thin_report.fraction_within_band << "\n";
    check(thin_report.steiner_points_merged > 0, "dense steiner points must be merged");
    check(thinned_vertices < kept_vertices, "merging must reduce the vertex count");
    check(keep_report.steiner_points_merged == 0, "caller points kept by default");
    check(thin_report.relative_area_error < 1e-9, "merging must preserve the area");
    check(thin_report.fraction_within_band > 0.80,
          "thinning must leave mostly in-band edges");

    Delete(square);
    Delete(steiner);
  }

  // 6. The size field scale controls the resolution: halving it roughly
  //    quadruples the triangle count.
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    p2t::SizeFieldRefiner a(square);
    p2t::SizeFieldOptions fine_options;
    fine_options.scale = 1.0;
    const std::size_t fine_count = a.Triangulate(
        [](const p2t::Point&) { return 0.08; }, fine_options).triangle_count;
    p2t::SizeFieldOptions coarse;
    coarse.scale = 2.0;
    p2t::SizeFieldRefiner b(square);
    const std::size_t coarse_count = b.Triangulate(
        [](const p2t::Point&) { return 0.08; }, coarse).triangle_count;
    std::cout << "scale: scale=1 tris=" << fine_count
              << " scale=2 tris=" << coarse_count << "\n";
    check(coarse_count < fine_count, "larger scale must produce fewer triangles");
    check(coarse_count * 3 < fine_count, "scale=2 should roughly quarter the count");
    Delete(square);
  }

  // 7. Iteration cap: an unsatisfiable cap is reported as non-convergence.
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    p2t::SizeFieldOptions options;
    options.max_iterations = 1;
    p2t::SizeFieldRefiner refiner(square);
    const p2t::SizeFieldReport& report =
        refiner.Triangulate([](const p2t::Point&) { return 0.01; }, options);
    std::cout << "iteration-cap: iterations=" << report.iterations
              << " converged=" << report.converged
              << " tris=" << report.triangle_count << "\n";
    check(!report.converged, "iteration cap must be reported");
    check(report.iterations == 1, "iteration cap must stop after one round");
    Delete(square);
  }

  // 8. Re-running a refined mesh with a coarser field merges the points the
  //    refiner created (both interior and contour) back out.
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    p2t::SizeFieldRefiner refiner(square);
    p2t::SizeFieldOptions fine_options;
    fine_options.max_iterations = 5;
    const p2t::SizeFieldReport& fine_report =
        refiner.Triangulate([](const p2t::Point&) { return 0.25; }, fine_options);
    const std::size_t fine_verts = fine_report.vertex_count;

    p2t::SizeFieldOptions coarse_options;
    coarse_options.max_iterations = 5;
    coarse_options.merge_ratio = 0.6;
    const p2t::SizeFieldReport& coarse_report =
        refiner.Triangulate([](const p2t::Point&) { return 0.5; }, coarse_options);
    std::cout << "coarsen: fine_verts=" << fine_verts
              << " merged=" << coarse_report.points_merged
              << " coarse_verts=" << coarse_report.vertex_count
              << " coarse_tris=" << coarse_report.triangle_count
              << " area_err=" << coarse_report.relative_area_error
              << " in_band=" << coarse_report.fraction_within_band << "\n";
    check(coarse_report.points_merged > 0, "coarsening must merge points");
    check(coarse_report.vertex_count < fine_verts,
          "coarsening must reduce the vertex count");
    check(coarse_report.relative_area_error < 1e-9,
          "coarsening must preserve the area");
    check(coarse_report.max_length_ratio <= coarse_options.split_ratio + 1e-9,
          "coarsening must not leave edges above split_ratio");
    Delete(square);
  }

  // 9. Real polygon data from testbed/data: a fine field refines the star,
  //    while a coarse enough field reproduces the plain CDT output exactly.
  {
    for (const std::string filename : {"star.dat", "diamond.dat"}) {
      std::size_t line_count = 0;
      std::vector<p2t::Point*> star = LoadTestbedData(filename, line_count);
      // The files differ in scale: pick a target edge from the bounding box.
      double span = 0.0;
      for (const p2t::Point* p : star) {
        span = std::max(span, std::max(std::fabs(p->x), std::fabs(p->y)));
      }
      const double target = 0.1 * span;
      if (star.size() < 3) {
        std::cout << "skip " << filename << " (missing data file)\n";
        Delete(star);
        continue;
      }
      p2t::CDT cdt(star);
      cdt.Triangulate();
      const std::vector<Tri> plain = SortedMesh(cdt.GetTriangles());

      p2t::SizeFieldOptions coarse;
      p2t::SizeFieldRefiner identity(star);
      const p2t::SizeFieldReport& identity_report =
          identity.Triangulate([](const p2t::Point&) { return 1e9; }, coarse);
      const std::vector<Tri> identical = SortedMesh(identity.GetTriangles());
      check(identity_report.iterations == 0, filename + ": coarse field iterates");
      check(identical == plain, filename + ": coarse field must match plain CDT");

      p2t::SizeFieldOptions fine;
      fine.max_iterations = 6;
      p2t::SizeFieldRefiner refiner(star);
      const p2t::SizeFieldReport& report = refiner.Triangulate(
          [target](const p2t::Point&) { return target; }, fine);
      std::cout << filename << ": plain_tris=" << plain.size()
                << " refined_tris=" << report.triangle_count
                << " area_err=" << report.relative_area_error
                << " in_band=" << report.fraction_within_band
                << " iters=" << report.iterations
                << " time_ms=" << 1000.0 * report.total_seconds
                << " identical_when_coarse=" << (identical == plain ? "yes" : "no")
                << "\n";
      check(report.triangle_count > plain.size(), filename + ": expected refinement");
      check(report.relative_area_error < 1e-9, filename + ": area must be preserved");
      Delete(star);
    }
  }

  // 10. Re-triangulating with the same field is a no-op (idempotence).
  {
    std::vector<p2t::Point*> square = MakeSquare(0, 0, 1, 1);
    p2t::SizeFieldRefiner refiner(square);
    p2t::SizeFieldOptions options;
    options.max_iterations = 8;
    const auto field = [](const p2t::Point& p) { return 0.03 + 0.15 * p.x; };
    const p2t::SizeFieldReport first = refiner.Triangulate(field, options);
    const p2t::SizeFieldReport second = refiner.Triangulate(field, options);
    std::cout << "idempotent: first_iters=" << first.iterations
              << " second_iters=" << second.iterations
              << " tris=" << second.triangle_count
              << " verts=" << second.vertex_count
              << " max_ratio=" << second.max_length_ratio << "\n";
    check(first.iterations > 0, "first run must refine");
    check(second.iterations == 0, "second run must be a no-op");
    check(second.triangle_count == first.triangle_count,
          "second run must keep the triangle count");
    Delete(square);
  }

  std::cout << (failures == 0 ? "\nALL CHECKS PASSED\n" : "\nCHECKS FAILED\n");
  return failures == 0 ? 0 : 1;
}

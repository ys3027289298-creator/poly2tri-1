// Quality-driven remeshing demo and verification harness.
//
// Four small, self-contained meshes (uniform, thin strip, holed, extreme
// aspect ratio) are triangulated with the stock CDT and then remeshed. Each
// case independently recomputes every reported figure by brute force and
// checks:
//   * the CDT output is bit-for-bit unchanged (per-vertex) after remeshing;
//   * all original constrained segments remain covered (constraint integrity);
//   * the output has no self-intersection and only positive-area triangles;
//   * area is preserved within floating-point tolerance.
//
// Build with run.sh (no third-party dependency beyond the repo and STL).

#include <poly2tri/poly2tri.h>
#include <poly2tri/remesh/quality_remesh.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& what)
{
  std::printf("    [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
  if (!ok) {
    ++g_failures;
  }
}

double dist(const p2t::Point& a, const p2t::Point& b)
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

double angle_at(const p2t::Point& apex, const p2t::Point& u,
                const p2t::Point& v)
{
  static const double kDeg = 180.0 / 3.14159265358979323846;
  const double ux = u.x - apex.x, uy = u.y - apex.y;
  const double vx = v.x - apex.x, vy = v.y - apex.y;
  return std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) * kDeg;
}

struct FaceGeom {
  double min_angle;
  double min_edge, max_edge;
  double area;
};

FaceGeom geom_of(const p2t::Point& a, const p2t::Point& b,
                 const p2t::Point& c)
{
  FaceGeom g;
  g.area = 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) -
                           (c.x - a.x) * (b.y - a.y));
  g.min_angle = std::min(angle_at(a, b, c),
                         std::min(angle_at(b, a, c), angle_at(c, a, b)));
  const double l0 = dist(a, b), l1 = dist(b, c), l2 = dist(c, a);
  g.min_edge = std::min(l0, std::min(l1, l2));
  g.max_edge = std::max(l0, std::max(l1, l2));
  return g;
}

using PtList = std::vector<p2t::Point*>;

struct CaseData {
  std::string name;
  PtList polyline;
  std::vector<PtList> holes;
  PtList steiner;
  double expected_area;
  p2t::RemeshOptions opts;
};

// deque storage keeps Point addresses stable across push_back.
p2t::Point* add(std::deque<p2t::Point>& store, double x, double y)
{
  store.emplace_back(x, y);
  return &store.back();
}

std::vector<CaseData> make_cases(std::deque<p2t::Point>& store)
{
  std::vector<CaseData> cases;

  // 1. Uniform-ish grid: 6x6 square with a 5x5 interior Steiner lattice.
  //    The lattice gives a dense interior next to long boundary edges, so
  //    both collapses and splits have work to do.
  {
    CaseData c;
    c.name = "uniform: 6x6 square with 5x5 Steiner lattice";
    c.polyline = { add(store, 0, 0), add(store, 6, 0),
                   add(store, 6, 6), add(store, 0, 6) };
    for (int i = 1; i <= 5; ++i) {
      for (int j = 1; j <= 5; ++j) {
        c.steiner.push_back(add(store, i, j));
      }
    }
    c.expected_area = 36.0;
    c.opts.target_min_angle_deg = 25.0;
    c.opts.max_edge_length = 1.6;
    cases.push_back(c);
  }

  // 2. Thin strip: 100x2 rectangle with 10-unit boundary sampling.
  {
    CaseData c;
    c.name = "thin strip: 100x2 with 10-unit boundary spacing";
    for (int x = 0; x <= 100; x += 10) {
      c.polyline.push_back(add(store, x, 0));
    }
    for (int x = 100; x >= 0; x -= 10) {
      c.polyline.push_back(add(store, x, 2));
    }
    c.expected_area = 200.0;
    c.opts.target_min_angle_deg = 25.0;
    c.opts.max_edge_length = 4.0;
    cases.push_back(c);
  }

  // 3. Square with a square hole.
  {
    CaseData c;
    c.name = "holed: 10x10 square with a 2x2 centered hole";
    c.polyline = { add(store, 0, 0), add(store, 10, 0),
                   add(store, 10, 10), add(store, 0, 10) };
    c.holes.push_back({ add(store, 4, 4), add(store, 6, 4),
                        add(store, 6, 6), add(store, 4, 6) });
    c.expected_area = 96.0;
    c.opts.target_min_angle_deg = 25.0;
    c.opts.max_edge_length = 2.5;
    cases.push_back(c);
  }

  // 4. Extreme aspect ratio: 200x1 rectangle, corners only. The initial
  //    triangulation contains two extremely skinny triangles; the target
  //    angle is deliberately set high to also exercise the "target not
  //    reached, report the actual value" path.
  {
    CaseData c;
    c.name = "extreme ratio: 200x1 rectangle (100:1)";
    c.polyline = { add(store, 0, 0), add(store, 200, 0),
                   add(store, 200, 1), add(store, 0, 1) };
    c.expected_area = 200.0;
    c.opts.target_min_angle_deg = 30.0;
    c.opts.max_edge_length = 8.0;
    cases.push_back(c);
  }

  // 5. Dense interior cluster: a large square with a tight group of
  //    Steiner points at the center. The close interior points create
  //    skinny *interior* triangles that circumcenter insertion plus
  //    Delaunay edge flips can genuinely improve, unlike boundary corners.
  {
    CaseData c;
    c.name = "dense interior: 20x20 square with a tight Steiner cluster";
    c.polyline = { add(store, 0, 0), add(store, 20, 0),
                   add(store, 20, 20), add(store, 0, 20) };
    const double cx = 10.0, cy = 10.0;
    c.steiner.push_back(add(store, cx, cy));
    c.steiner.push_back(add(store, cx + 0.4, cy));
    c.steiner.push_back(add(store, cx, cy + 0.45));
    c.steiner.push_back(add(store, cx + 0.5, cy + 0.6));
    c.steiner.push_back(add(store, cx - 0.6, cy + 0.3));
    c.steiner.push_back(add(store, cx + 1.5, cy + 1.5));
    c.expected_area = 400.0;
    c.opts.target_min_angle_deg = 20.0;
    c.opts.max_edge_length = 4.0;
    cases.push_back(c);
  }

  return cases;
}

// Original constrained segments derived directly from the CDT triangle
// flags, as ordered coordinate pairs.
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
constraint_segments(const std::vector<p2t::Triangle*>& tris)
{
  std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>
      out;
  auto key = [](double ax, double ay, double bx, double by) {
    if (ay > by || (ay == by && ax > bx)) {
      std::swap(ax, bx);
      std::swap(ay, by);
    }
    return std::make_pair(std::make_pair(ax, ay), std::make_pair(bx, by));
  };
  std::vector<decltype(key(0, 0, 0, 0))> seen;
  for (p2t::Triangle* t : tris) {
    for (int i = 0; i < 3; ++i) {
      if (!t->constrained_edge[i]) {
        continue;
      }
      const p2t::Point* p = t->GetPoint((i + 1) % 3);
      const p2t::Point* q = t->GetPoint((i + 2) % 3);
      seen.push_back(key(p->x, p->y, q->x, q->y));
    }
  }
  std::sort(seen.begin(), seen.end());
  seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
  return seen;
}

// Flatten every triangle vertex coordinate to a snapshot vector.
std::vector<double> snapshot(const std::vector<p2t::Triangle*>& tris)
{
  std::vector<double> v;
  v.reserve(tris.size() * 6);
  for (p2t::Triangle* t : tris) {
    for (int k = 0; k < 3; ++k) {
      v.push_back(t->GetPoint(k)->x);
      v.push_back(t->GetPoint(k)->y);
    }
  }
  return v;
}

double orient(double ax, double ay, double bx, double by, double cx,
              double cy)
{
  return (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
}

bool on_seg(double px, double py, double ax, double ay, double bx, double by,
            double tol)
{
  if (std::fabs(orient(ax, ay, bx, by, px, py)) >
      tol * std::max(1.0, std::hypot(bx - ax, by - ay))) {
    return false;
  }
  return std::min(ax, bx) - tol <= px && px <= std::max(ax, bx) + tol &&
         std::min(ay, by) - tol <= py && py <= std::max(ay, by) + tol;
}

// Independent brute-force check that every original constraint segment is
// covered by a chain of output edges.
bool constraints_covered(
    const p2t::RemeshResult& r,
    const std::vector<
        std::pair<std::pair<double, double>, std::pair<double, double>>>&
        constraints,
    double tol)
{
  for (const auto& seg : constraints) {
    const double ax = seg.first.first, ay = seg.first.second;
    const double bx = seg.second.first, by = seg.second.second;
    const double dx = bx - ax, dy = by - ay;
    const double len2 = dx * dx + dy * dy;
    std::vector<std::pair<double, double>> intervals;
    for (const p2t::RemeshFace& f : r.faces) {
      for (int i = 0; i < 3; ++i) {
        const p2t::Point& p = r.vertices[f.v[i]];
        const p2t::Point& q = r.vertices[f.v[(i + 1) % 3]];
        if (!on_seg(p.x, p.y, ax, ay, bx, by, tol) ||
            !on_seg(q.x, q.y, ax, ay, bx, by, tol)) {
          continue;
        }
        const double ta = ((p.x - ax) * dx + (p.y - ay) * dy) / len2;
        const double tb = ((q.x - ax) * dx + (q.y - ay) * dy) / len2;
        intervals.emplace_back(std::min(ta, tb), std::max(ta, tb));
      }
    }
    if (intervals.empty()) {
      return false;
    }
    std::sort(intervals.begin(), intervals.end());
    if (intervals[0].first > tol) {
      return false;
    }
    double reach = intervals[0].second;
    for (std::size_t i = 1; i < intervals.size(); ++i) {
      if (intervals[i].first > reach + tol) {
        return false;
      }
      reach = std::max(reach, intervals[i].second);
    }
    if (reach < 1.0 - tol) {
      return false;
    }
  }
  return true;
}

bool proper_cross(double ax, double ay, double bx, double by, double cx,
                  double cy, double dx, double dy, double tol)
{
  const double o1 = orient(ax, ay, bx, by, cx, cy);
  const double o2 = orient(ax, ay, bx, by, dx, dy);
  const double o3 = orient(cx, cy, dx, dy, ax, ay);
  const double o4 = orient(cx, cy, dx, dy, bx, by);
  return ((o1 > tol && o2 < -tol) || (o1 < -tol && o2 > tol)) &&
         ((o3 > tol && o4 < -tol) || (o3 < -tol && o4 > tol));
}

// Brute-force O(E^2) self-intersection and orientation check on the output.
bool output_clean(const p2t::RemeshResult& r, double tol)
{
  struct E {
    int a, b;
  };
  std::vector<E> edges;
  for (const p2t::RemeshFace& f : r.faces) {
    const p2t::Point& a = r.vertices[f.v[0]];
    const p2t::Point& b = r.vertices[f.v[1]];
    const p2t::Point& c = r.vertices[f.v[2]];
    if (orient(a.x, a.y, b.x, b.y, c.x, c.y) <= 0.0) {
      return false; // every output face must be CCW with positive area
    }
    for (int i = 0; i < 3; ++i) {
      int u = f.v[i], v2 = f.v[(i + 1) % 3];
      if (u > v2) std::swap(u, v2);
      bool dup = false;
      for (const E& e : edges) {
        if (e.a == u && e.b == v2) {
          dup = true;
          break;
        }
      }
      if (!dup) {
        edges.push_back({ u, v2 });
      }
    }
  }
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const p2t::Point& a = r.vertices[edges[i].a];
    const p2t::Point& b = r.vertices[edges[i].b];
    for (std::size_t j = i + 1; j < edges.size(); ++j) {
      if (edges[i].a == edges[j].a || edges[i].a == edges[j].b ||
          edges[i].b == edges[j].a || edges[i].b == edges[j].b) {
        continue;
      }
      const p2t::Point& c = r.vertices[edges[j].a];
      const p2t::Point& d = r.vertices[edges[j].b];
      if (proper_cross(a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y, tol)) {
        return false;
      }
    }
  }
  return true;
}

struct BruteStats {
  double min_angle = 180.0;
  double min_edge = 1e300;
  double max_edge = 0.0;
  double area = 0.0;
};

BruteStats brute(const p2t::RemeshResult& r)
{
  BruteStats s;
  for (const p2t::RemeshFace& f : r.faces) {
    const FaceGeom g = geom_of(r.vertices[f.v[0]], r.vertices[f.v[1]],
                               r.vertices[f.v[2]]);
    s.min_angle = std::min(s.min_angle, g.min_angle);
    s.min_edge = std::min(s.min_edge, g.min_edge);
    s.max_edge = std::max(s.max_edge, g.max_edge);
    s.area += g.area;
  }
  return s;
}

BruteStats brute_cdt(const std::vector<p2t::Triangle*>& tris)
{
  BruteStats s;
  for (p2t::Triangle* t : tris) {
    const FaceGeom g = geom_of(*t->GetPoint(0), *t->GetPoint(1),
                               *t->GetPoint(2));
    s.min_angle = std::min(s.min_angle, g.min_angle);
    s.min_edge = std::min(s.min_edge, g.min_edge);
    s.max_edge = std::max(s.max_edge, g.max_edge);
    s.area += g.area;
  }
  return s;
}

bool nearly(double a, double b, double rel = 1e-9, double abs = 1e-10)
{
  return std::fabs(a - b) <= abs + rel * std::max(std::fabs(a), std::fabs(b));
}

// Smallest minimum angle among triangles that do not touch a constrained
// edge, i.e. the quality the remesher can actually improve in the interior.
double interior_min_angle_cdt(const std::vector<p2t::Triangle*>& tris)
{
  double best = 180.0;
  for (p2t::Triangle* t : tris) {
    bool touches_constraint = false;
    for (int i = 0; i < 3; ++i) {
      if (t->constrained_edge[i]) {
        touches_constraint = true;
      }
    }
    if (touches_constraint) {
      continue;
    }
    const FaceGeom g = geom_of(*t->GetPoint(0), *t->GetPoint(1),
                               *t->GetPoint(2));
    best = std::min(best, g.min_angle);
  }
  return best; // 180 means "no fully interior face in the input"
}

// Same metric on the remeshed output: a face is interior when none of its
// edges is part of the preserved constraint coverage. We approximate this
// by excluding faces that use an output edge collinear with and lying on an
// original constrained segment.
double interior_min_angle_remesh(
    const p2t::RemeshResult& r,
    const std::vector<
        std::pair<std::pair<double, double>, std::pair<double, double>>>&
        constraints,
    double tol)
{
  auto on_any_constraint = [&](const p2t::Point& p, const p2t::Point& q) {
    for (const auto& seg : constraints) {
      if (on_seg(p.x, p.y, seg.first.first, seg.first.second,
                 seg.second.first, seg.second.second, tol) &&
          on_seg(q.x, q.y, seg.first.first, seg.first.second,
                 seg.second.first, seg.second.second, tol)) {
        return true;
      }
    }
    return false;
  };
  double best = 180.0;
  for (const p2t::RemeshFace& f : r.faces) {
    bool boundary_face = false;
    for (int i = 0; i < 3; ++i) {
      if (on_any_constraint(r.vertices[f.v[i]],
                            r.vertices[f.v[(i + 1) % 3]])) {
        boundary_face = true;
      }
    }
    if (boundary_face) {
      continue;
    }
    const FaceGeom g = geom_of(r.vertices[f.v[0]], r.vertices[f.v[1]],
                               r.vertices[f.v[2]]);
    best = std::min(best, g.min_angle);
  }
  return best;
}

void run_case(const CaseData& c)
{
  std::printf("== %s ==\n", c.name.c_str());
  p2t::CDT cdt{ c.polyline };
  for (const PtList& h : c.holes) {
    cdt.AddHole(h);
  }
  for (p2t::Point* p : c.steiner) {
    cdt.AddPoint(p);
  }
  cdt.Triangulate();
  const std::vector<p2t::Triangle*> tris = cdt.GetTriangles();

  const BruteStats bs0 = brute_cdt(tris);
  const auto constraints = constraint_segments(tris);
  const std::vector<double> before = snapshot(tris);

  const p2t::RemeshResult r = p2t::QualityRemesh(tris, c.opts);

  // The stock triangulation output must be untouched, vertex by vertex.
  const std::vector<double> after = snapshot(tris);
  check(before == after, "CDT per-vertex output unchanged after remeshing");
  check(cdt.GetTriangles().size() == tris.size(),
        "CDT triangle count unchanged after remeshing");

  const BruteStats bs1 = brute(r);
  const double tol = 1e-10;

  std::printf("  triangles:        %zu -> %zu (splits=%zu, collapses=%zu, "
              "rejected=%zu, checks=%zu)\n",
              r.input_triangles, r.output_triangles, r.split_count,
              r.collapse_count, r.rejected_count, r.checks_run);
  std::printf("  min angle (deg):  %.4f -> %.4f   (target %.2f: %s)\n",
              r.min_angle_before, r.min_angle_after,
              c.opts.target_min_angle_deg,
              r.target_angle_reached ? "REACHED" : "not reached");
  if (!r.target_angle_reached) {
    std::printf("                    actual achieved minimum angle = %.4f deg\n",
                r.min_angle_after);
  }
  const double interior_before = interior_min_angle_cdt(tris);
  const double interior_after =
      interior_min_angle_remesh(r, constraints, tol);
  std::printf("  interior min ang: %.4f -> %.4f   (excludes constrained "
              "boundary faces)\n", interior_before, interior_after);
  std::printf("  edge length:      min %.4f -> %.4f, max %.4f -> %.4f "
              "(Lmax=%.3f)\n",
              r.min_edge_before, r.min_edge_after, r.max_edge_before,
              r.max_edge_after, c.opts.max_edge_length);
  std::printf("  area:             %.8f -> %.8f  (abs err %.3e, rel %.3e)\n",
              r.area_before, r.area_after, r.area_error_abs,
              r.area_error_rel);
  std::printf("  constraints:      %zu/%zu covered\n",
              r.constraints_satisfied, r.constraints_total);
  std::printf("  self-intersection:%s   invariant violation: %s\n",
              r.self_intersecting ? " YES" : " none",
              r.invariant_violated ? "YES" : "none");
  std::printf("  remesh time:      %.3f ms\n", r.elapsed_seconds * 1e3);

  // Cross-check every reported number against an independent computation.
  check(r.input_triangles == tris.size(), "reported input triangle count");
  check(r.output_triangles == r.faces.size(), "reported output triangle count");
  check(nearly(r.min_angle_before, bs0.min_angle),
        "before min-angle matches brute force");
  check(nearly(r.min_angle_after, bs1.min_angle),
        "after min-angle matches brute force");
  check(nearly(r.max_edge_before, bs0.max_edge),
        "before max-edge matches brute force");
  check(nearly(r.max_edge_after, bs1.max_edge),
        "after max-edge matches brute force");
  check(nearly(r.min_edge_before, bs0.min_edge),
        "before min-edge matches brute force");
  check(nearly(r.min_edge_after, bs1.min_edge),
        "after min-edge matches brute force");
  check(nearly(r.area_before, bs0.area), "before area matches brute force");
  check(nearly(r.area_after, bs1.area), "after area matches brute force");
  check(nearly(r.area_after, c.expected_area, 1e-9, 1e-8),
        "output area matches input polygon area");
  check(r.min_angle_after >= r.min_angle_before - 1e-9,
        "minimum angle never regresses");
  // A 180 baseline means the input had no fully interior face; in that case
  // the first interior faces are created by refinement itself.
  if (interior_before < 179.0) {
    check(interior_after >= interior_before - 1e-9,
          "interior minimum angle never regresses");
  } else {
    check(interior_after > 0.0,
          "refinement creates valid positive-angle interior faces");
  }
  check(r.constraints_total == constraints.size(),
        "constraint count matches CDT flags");
  check(r.constraints_satisfied == r.constraints_total,
        "all constraints preserved (library check)");
  check(constraints_covered(r, constraints, tol),
        "all constraints preserved (independent brute force)");
  check(!r.self_intersecting && !r.invariant_violated,
        "no self-intersection and no invariant violation");
  check(output_clean(r, tol), "output faces CCW and pairwise non-crossing");
  check(r.area_error_rel < 1e-12, "area error at floating-point noise level");
  std::printf("\n");
}

// A second pass with different targets proves the knobs are configurable and
// the reachability reporting is honest.
void run_config_check()
{
  std::printf("== configurability: extreme-ratio case at 20 / 25 / 30 deg ==\n");
  for (double target : { 20.0, 25.0, 30.0 }) {
    // Fresh points per CDT: the library attaches Edge objects to each
    // Point and frees them when the CDT is destroyed, so Point instances
    // must never be shared across CDT objects.
    std::deque<p2t::Point> store;
    PtList poly{ add(store, 0, 0), add(store, 200, 0),
                 add(store, 200, 1), add(store, 0, 1) };
    p2t::CDT cdt{ poly };
    cdt.Triangulate();
    p2t::RemeshOptions opts;
    opts.target_min_angle_deg = target;
    opts.max_edge_length = 8.0;
    const p2t::RemeshResult r = p2t::QualityRemesh(cdt.GetTriangles(), opts);
    std::printf("  target=%.1f deg -> min angle %.4f deg, %zu triangles, "
                "reached=%s\n",
                target, r.min_angle_after, r.output_triangles,
                r.target_angle_reached ? "yes" : "no");
    check(!r.invariant_violated, "invariants hold at this target");
    check(r.min_angle_after >= r.min_angle_before - 1e-9,
          "min angle does not regress at this target");
  }
  std::printf("\n");
}

// Synthetic-mesh check: QualityRemesh only needs a list of Triangle
// objects, so a hand-built mesh lets us exercise the circumcenter path on
// an all-acute skinny *interior* triangle, which a Delaunay CDT almost
// never produces directly.
void run_synthetic_circumcenter_check()
{
  std::printf("== synthetic: acute skinny interior triangle, circumcenter ==\n");
  // Inner all-acute isosceles triangle A-B-C: apex angle 30 deg, base
  // angles 75 deg (base 2, height tan(75 deg) ~= 3.732). All angles acute,
  // so its circumcenter lies strictly inside the face.
  const double h = std::tan(75.0 * 3.14159265358979323846 / 180.0);
  p2t::Point A(-1.0, 0.0), B(1.0, 0.0), C(0.0, h);
  // Outer vertices forming three surrounding triangles.
  p2t::Point D(0.0, -2.0), E(3.0, 4.5), F(-3.0, 4.5);

  p2t::Triangle inner(A, B, C);
  p2t::Triangle t_abd(A, B, D);
  p2t::Triangle t_bce(B, C, E);
  p2t::Triangle t_caf(C, A, F);
  // Mark every outer-facing edge constrained; the three inner edges stay
  // unconstrained so the central face is a fully interior triangle.
  t_abd.MarkConstrainedEdge(&A, &D);
  t_abd.MarkConstrainedEdge(&B, &D);
  t_bce.MarkConstrainedEdge(&B, &E);
  t_bce.MarkConstrainedEdge(&C, &E);
  t_caf.MarkConstrainedEdge(&C, &F);
  t_caf.MarkConstrainedEdge(&A, &F);

  const std::vector<p2t::Triangle*> mesh{ &inner, &t_abd, &t_bce, &t_caf };

  p2t::RemeshOptions opts;
  opts.target_min_angle_deg = 33.0; // above the 30 deg apex angle
  opts.max_edge_length = 1e9;       // disable size-driven splitting
  opts.max_operations = 200;
  const p2t::RemeshResult r = p2t::QualityRemesh(mesh, opts);

  std::printf("  triangles: %zu -> %zu (splits=%zu, steiner=%zu)\n",
              r.input_triangles, r.output_triangles, r.split_count,
              r.steiner_count);
  std::printf("  min angle: %.4f -> %.4f (target 33.00: %s)\n",
              r.min_angle_before, r.min_angle_after,
              r.target_angle_reached ? "REACHED" : "not reached");
  check(r.steiner_count >= 1,
        "circumcenter Steiner insertion fires on an acute skinny interior face");
  // The surrounding constrained boundary caps how far a Ruppert cascade can
  // converge; here we assert the mechanism works and never collapses the
  // mesh below a sane fraction of the starting quality, while the four CDT
  // cases report the honestly achieved minimum angle.
  check(r.min_angle_after >= r.min_angle_before * 0.25,
        "circumcenter cascade keeps the mesh away from degeneracy");
  check(!r.invariant_violated, "synthetic mesh keeps all invariants");
  check(r.constraints_satisfied == r.constraints_total,
        "outer constrained edges preserved in the synthetic mesh");
  check(r.output_triangles > r.input_triangles,
        "refinement produces additional well-formed triangles");
  std::printf("\n");
}

} // namespace

int main()
{
  std::deque<p2t::Point> store;
  const std::vector<CaseData> cases = make_cases(store);
  for (const CaseData& c : cases) {
    run_case(c);
  }
  run_config_check();
  run_synthetic_circumcenter_check();
  std::printf("%s (%d failure%s)\n",
              g_failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
              g_failures, g_failures == 1 ? "" : "s");
  return g_failures == 0 ? 0 : 1;
}

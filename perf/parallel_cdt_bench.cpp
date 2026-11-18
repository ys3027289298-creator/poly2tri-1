// Benchmark + geometric-consistency checker for p2t::ParallelCDT.
//
// For four input classes (large polyline, with hole, many contours,
// slender) this program runs the classic single-threaded p2t::CDT and
// the divide-and-conquer p2t::ParallelCDT, then reports:
//   - wall time and speedup
//   - triangle count difference
//   - covered-area relative error
//   - constrained-edge coverage (every original edge, split at the cut
//     lines, must appear as a triangle edge in the stitched result)
// It also exercises the sequential fallback path.
//
// Self-contained: inputs are generated deterministically, no files.

#include <poly2tri/poly2tri.h>
#include <poly2tri/sweep/parallel_cdt.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct Rng {
  explicit Rng(uint64_t seed) : s(seed) {}
  double next() {  // uniform in [0, 1)
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(s >> 11) * (1.0 / 9007199254740992.0);
  }
  uint64_t s;
};

struct Input {
  std::string name;
  std::vector<p2t::Point*> polyline;
  std::vector<std::vector<p2t::Point*>> holes;
  Input() = default;
  Input(Input&&) = default;
  Input& operator=(Input&&) = default;
  Input(const Input&) = delete;
  Input& operator=(const Input&) = delete;
  ~Input() {
    for (auto* p : polyline) delete p;
    for (auto& h : holes)
      for (auto* p : h) delete p;
  }
  size_t NumPoints() const {
    size_t n = polyline.size();
    for (auto& h : holes) n += h.size();
    return n;
  }
};

// Star-shaped (hence simple) large polygon.
Input MakeLarge(size_t n) {
  Input in;
  in.name = "large";
  Rng rng(12345);
  for (size_t i = 0; i < n; ++i) {
    const double th = 2.0 * M_PI * static_cast<double>(i) / n;
    const double r = 1000.0 * (1.0 + 0.25 * std::sin(7 * th) +
                               0.1 * std::sin(13 * th + 1.0) +
                               0.05 * (rng.next() - 0.5));
    in.polyline.push_back(new p2t::Point(r * std::cos(th), r * std::sin(th)));
  }
  return in;
}

// Large polygon with a single hole.
Input MakeWithHole(size_t n) {
  Input in = MakeLarge(n);
  in.name = "with_hole";
  std::vector<p2t::Point*> hole;
  const size_t m = 2000;
  for (size_t i = 0; i < m; ++i) {
    const double th = 2.0 * M_PI * static_cast<double>(i) / m;
    hole.push_back(new p2t::Point(300.0 * std::cos(th), 300.0 * std::sin(th)));
  }
  in.holes.push_back(hole);
  return in;
}

// Polygon with many inner contours (holes).
Input MakeManyContours(size_t n) {
  Input in = MakeLarge(n);
  in.name = "many_contours";
  const size_t m = 500;
  for (int gy = -1; gy <= 2; ++gy) {
    for (int gx = -1; gx <= 2; ++gx) {
      const double cx = 150.0 * gx;
      const double cy = 150.0 * gy;
      std::vector<p2t::Point*> hole;
      for (size_t i = 0; i < m; ++i) {
        const double th = 2.0 * M_PI * static_cast<double>(i) / m;
        hole.push_back(
            new p2t::Point(cx + 40.0 * std::cos(th), cy + 40.0 * std::sin(th)));
      }
      in.holes.push_back(hole);
    }
  }
  return in;
}

// Very slender wavy strip (aspect ratio ~1e5 : 1).
Input MakeSlender(size_t n) {
  Input in;
  in.name = "slender";
  const size_t m = n / 2;
  const double length = 100000.0;
  for (size_t i = 0; i < m; ++i) {
    const double x = length * static_cast<double>(i) / (m - 1);
    in.polyline.push_back(new p2t::Point(x, 0.2 * std::sin(i * 0.37)));
  }
  for (size_t i = m; i-- > 0;) {
    const double x = length * static_cast<double>(i) / (m - 1);
    in.polyline.push_back(new p2t::Point(x, 1.0 + 0.2 * std::cos(i * 0.21)));
  }
  return in;
}

long double TriArea(const p2t::Triangle* t) {
  p2t::Triangle* tc = const_cast<p2t::Triangle*>(t);
  const p2t::Point* a = tc->GetPoint(0);
  const p2t::Point* b = tc->GetPoint(1);
  const p2t::Point* c = tc->GetPoint(2);
  const long double v =
      ((long double)b->x - a->x) * ((long double)c->y - a->y) -
      ((long double)b->y - a->y) * ((long double)c->x - a->x);
  return v < 0 ? -v * 0.5L : v * 0.5L;
}

long double TotalArea(const std::vector<p2t::Triangle*>& tris) {
  long double a = 0.0L;
  for (auto* t : tris) a += TriArea(t);
  return a;
}

double NowMs() {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// Bit-exact geometric key for a point (crossing points are recomputed
// with the exact same formula as inside ParallelCDT, so keys match).
struct PKey {
  uint64_t x, y;
  bool operator==(const PKey& o) const { return x == o.x && y == o.y; }
};
struct PKeyHash {
  size_t operator()(const PKey& k) const {
    return std::hash<uint64_t>()(k.x) ^ (std::hash<uint64_t>()(k.y) << 1);
  }
};
uint64_t Bits(double v) {
  uint64_t u;
  static_assert(sizeof u == sizeof v, "size");
  __builtin_memcpy(&u, &v, sizeof u);
  return u;
}
PKey KeyOf(double x, double y) { return PKey{Bits(x), Bits(y)}; }

struct EKey {
  PKey a, b;
  bool operator==(const EKey& o) const { return a == o.a && b == o.b; }
};
struct EKeyHash {
  size_t operator()(const EKey& k) const {
    PKeyHash h;
    return h(k.a) ^ (h(k.b) << 1);
  }
};
EKey EdgeKey(PKey p, PKey q) {  // order independent
  if (q.x < p.x || (q.x == p.x && q.y < p.y)) std::swap(p, q);
  return EKey{p, q};
}

// Verify that every original constrained edge (outer + holes), split at
// the parallel cut lines, is present as a triangle edge in the result.
bool CheckConstraintEdges(const Input& in,
                          const std::vector<p2t::Triangle*>& tris, int axis,
                          const std::vector<double>& splits,
                          size_t* missing_out) {
  std::unordered_multiset<EKey, EKeyHash> edges;
  edges.reserve(tris.size() * 3);
  for (auto* t : tris) {
    for (int i = 0; i < 3; ++i) {
      const p2t::Point* p = t->GetPoint(i);
      const p2t::Point* q = t->GetPoint((i + 1) % 3);
      edges.insert(EdgeKey(KeyOf(p->x, p->y), KeyOf(q->x, q->y)));
    }
  }

  size_t missing = 0;
  std::vector<const std::vector<p2t::Point*>*> loops{&in.polyline};
  for (auto& h : in.holes) loops.push_back(&h);
  for (auto* loop : loops) {
    const size_t m = loop->size();
    for (size_t i = 0; i < m; ++i) {
      const p2t::Point* a = (*loop)[i];
      const p2t::Point* b = (*loop)[(i + 1) % m];
      const double ua = axis == 0 ? a->x : a->y;
      const double ub = axis == 0 ? b->x : b->y;
      std::vector<double> ts{0.0, 1.0};
      if (axis >= 0 && ua != ub) {
        for (double B : splits) {
          const double t = (B - ua) / (ub - ua);
          if (t > 0.0 && t < 1.0) ts.push_back(t);
        }
        std::sort(ts.begin(), ts.end());
      }
      for (size_t k = 0; k + 1 < ts.size(); ++k) {
        const double t0 = ts[k], t1 = ts[k + 1];
        // t==0/1 must use the exact endpoint coordinates: a + 1.0*(b-a)
        // is not guaranteed to round back to b.
        const PKey p0 =
            t0 == 0.0 ? KeyOf(a->x, a->y)
                      : KeyOf(a->x + t0 * (b->x - a->x),
                              a->y + t0 * (b->y - a->y));
        const PKey p1 =
            t1 == 1.0 ? KeyOf(b->x, b->y)
                      : KeyOf(a->x + t1 * (b->x - a->x),
                              a->y + t1 * (b->y - a->y));
        if (edges.find(EdgeKey(p0, p1)) == edges.end()) ++missing;
      }
    }
  }
  *missing_out = missing;
  return missing == 0;
}

struct RunResult {
  double ms = 0;
  size_t tris = 0;
  long double area = 0;
};

RunResult RunSequential(const Input& in) {
  p2t::CDT cdt(in.polyline);
  for (auto& h : in.holes) cdt.AddHole(h);
  const double t0 = NowMs();
  cdt.Triangulate();
  const double t1 = NowMs();
  RunResult r;
  r.ms = t1 - t0;
  r.tris = cdt.GetTriangles().size();
  r.area = TotalArea(cdt.GetTriangles());
  return r;
}

void PrintRow(const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  std::vprintf(fmt, ap);
  va_end(ap);
  std::printf("\n");
}

}  // namespace

int main() {
  const unsigned hw = std::thread::hardware_concurrency();
  std::printf("hardware_concurrency = %u\n\n", hw);

  bool all_ok = true;

  struct Config {
    size_t segments;
    size_t threads;
  };
  const std::vector<Config> configs{
      {hw, hw},
      {4 * hw, hw},
  };

  std::vector<Input> inputs;
  inputs.push_back(MakeLarge(300000));
  inputs.push_back(MakeWithHole(100000));
  inputs.push_back(MakeManyContours(50000));
  inputs.push_back(MakeSlender(200000));

  for (auto& in : inputs) {
    const RunResult seq = RunSequential(in);
    std::printf(
        "== %-14s points=%-7zu | sequential: %8.1f ms  tris=%-7zu ==\n",
        in.name.c_str(), in.NumPoints(), seq.ms, seq.tris);
    for (const auto& cfg : configs) {
      p2t::ParallelOptions opt;
      opt.num_segments = cfg.segments;
      opt.num_threads = cfg.threads;
      p2t::ParallelCDT par(in.polyline, opt);
      for (auto& h : in.holes) par.AddHole(h);
      const double t0 = NowMs();
      par.Triangulate();
      const double t1 = NowMs();
      const auto tris = par.GetTriangles();
      const long double area = TotalArea(tris);
      const long double area_err =
          seq.area != 0 ? (area - seq.area) / seq.area : 0.0L;
      size_t missing = 0;
      const bool edges_ok = CheckConstraintEdges(
          in, tris, par.SplitAxis(), par.SplitPositions(), &missing);
      const bool area_ok = area_err > -1e-9L && area_err < 1e-9L;
      const long tri_diff =
          static_cast<long>(tris.size()) - static_cast<long>(seq.tris);
      all_ok = all_ok && edges_ok && area_ok;
      PrintRow(
          "  parallel S=%-3zu T=%-2zu used=%-2zu fb=%d | %7.1f ms  "
          "speedup=%5.2fx  tris=%-7zu (d=%+ld)  area_err=%+.2e  "
          "edges=%s",
          cfg.segments, cfg.threads, par.NumSegments(),
          par.UsedFallback() ? 1 : 0, t1 - t0, seq.ms / (t1 - t0),
          tris.size(), tri_diff, static_cast<double>(area_err),
          edges_ok ? "OK" : ("MISSING=" + std::to_string(missing)).c_str());
    }
    std::printf("\n");
  }

  // --- Fallback behaviour -------------------------------------------------
  std::printf("== fallback checks ==\n");
  {
    // (a) tiny input: below min_points_for_parallel -> sequential path.
    std::vector<p2t::Point*> tri{new p2t::Point(0, 0), new p2t::Point(1, 0),
                                 new p2t::Point(0, 1)};
    p2t::ParallelCDT par(tri);
    par.Triangulate();
    const bool ok = par.UsedFallback() && par.GetTriangles().size() == 1;
    std::printf("  tiny input        -> sequential path, tris=%zu : %s\n",
                par.GetTriangles().size(), ok ? "OK" : "FAIL");
    all_ok = all_ok && ok;
    for (auto* p : tri) delete p;
  }
  {
    // (b) impossible segmentation request -> recover via single thread.
    Input in = MakeLarge(50000);
    p2t::ParallelOptions opt;
    opt.num_segments = 1000000;  // more cuts than distinct coordinates
    opt.num_threads = hw;
    p2t::ParallelCDT par(in.polyline, opt);
    par.Triangulate();
    const RunResult seq = RunSequential(in);
    const long double err =
        (TotalArea(par.GetTriangles()) - seq.area) / seq.area;
    const bool ok = par.UsedFallback() && err > -1e-12L && err < 1e-12L;
    std::printf(
        "  S=1000000 request -> fallback=%d, area_err=%.2e : %s\n",
        par.UsedFallback() ? 1 : 0, static_cast<double>(err),
        ok ? "OK" : "FAIL");
    all_ok = all_ok && ok;
  }

  std::printf("\n%s\n", all_ok ? "ALL CHECKS PASSED" : "CHECKS FAILED");
  return all_ok ? 0 : 1;
}

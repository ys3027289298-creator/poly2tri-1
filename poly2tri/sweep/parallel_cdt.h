/*
 * Poly2Tri Copyright (c) 2009-2022, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include "cdt.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * Divide-and-conquer parallel triangulation on top of the classic
 * single-threaded p2t::CDT sweep.
 *
 * The input polyline (+ holes) is intersected with S axis-aligned slabs
 * whose cut lines are placed between vertex coordinates.  Every slab is
 * decomposed into a set of simple component polygons (outer loops plus
 * the holes they contain) which are triangulated independently - and in
 * parallel - by the unmodified sequential CDT.  Because the slabs
 * partition the plane, the union of all slab triangulations covers
 * exactly the same region as a single-threaded run and keeps every
 * original constrained edge (split, where a cut line crosses it, into
 * collinear pieces).  Triangles on both sides of a cut line share the
 * generated divider edge, so "stitching" is simply the concatenation of
 * the per-slab triangle sets.
 *
 * Each slab works on private copies of the points it uses: p2t::Point
 * carries mutable per-triangulation state (Point::edge_list) so Point
 * objects must never be shared between concurrent CDT instances.  The
 * sequential fallback below is the only path that touches the caller's
 * original Point objects and it is vertex-for-vertex identical to using
 * p2t::CDT directly.
 *
 * Any internal failure (clipping inconsistency, loop extraction,
 * per-component or global area validation, exceptions) discards all
 * partial work and re-runs the plain single-threaded CDT, so the result
 * is always geometrically consistent with the sequential library.
 *
 * Usage mirrors p2t::CDT:
 *
 *   p2t::ParallelCDT cdt(polyline);   // or with ParallelOptions
 *   cdt.AddHole(hole);
 *   cdt.Triangulate();
 *   auto triangles = cdt.GetTriangles();
 *
 * Triangles and generated Steiner points (on the cut lines) are owned by
 * this object, exactly like p2t::CDT owns its triangles.
 *
 * @author (parallel divide-and-conquer extension)
 */

namespace p2t {

struct ParallelOptions {
  /// Number of slabs the input is divided into. 0 = auto (= num_threads).
  size_t num_segments;
  /// Worker threads used to triangulate the slabs. 0 = auto
  /// (= std::thread::hardware_concurrency).
  size_t num_threads;
  /// Inputs with fewer total points always take the sequential path.
  size_t min_points_for_parallel;

  ParallelOptions()
      : num_segments(0), num_threads(0), min_points_for_parallel(4096) {}
};

class ParallelCDT {
public:
  explicit ParallelCDT(const std::vector<Point*>& polyline,
                       const ParallelOptions& options = ParallelOptions())
      : polyline_(polyline),
        options_(options),
        used_fallback_(false),
        num_segments_(1),
        split_axis_(-1) {}

  /// Add a hole (same contract as p2t::CDT::AddHole).
  void AddHole(const std::vector<Point*>& polyline) {
    holes_.push_back(polyline);
  }

  /// Triangulate - do this AFTER you've added the polyline and holes.
  void Triangulate();

  /// Get triangles (valid until this object is destroyed).
  std::vector<Triangle*> GetTriangles() { return triangles_; }

  /// True when the result was produced by the single-threaded CDT path
  /// (small input, num_segments < 2, or recovery from a failed divide
  /// and conquer pass).
  bool UsedFallback() const { return used_fallback_; }

  /// Number of slabs actually used (1 when the sequential path was taken).
  size_t NumSegments() const { return num_segments_; }

  /// Verification helpers: split axis (0 = x, 1 = y, -1 = sequential) and
  /// the split coordinates along that axis.
  int SplitAxis() const { return split_axis_; }
  const std::vector<double>& SplitPositions() const { return splits_; }

private:
  struct Segment {
    Point* a;
    Point* b;
  };

  struct SlabResult {
    SlabResult() : ok(false), area(0.0L) {}
    bool ok;
    std::vector<std::unique_ptr<CDT>> cdts;
    std::vector<std::unique_ptr<Point>> points;  // slab-local point copies
    long double area;                            // sum of triangle areas
  };

  void RunSequential();
  bool TryParallel(size_t num_segments, size_t num_threads);
  SlabResult ProcessSlab(int axis, double lo, double hi,
                         double min_dist) const;

  static double Coord(const Point* p, int axis) {
    return axis == 0 ? p->x : p->y;
  }
  static long double SignedArea(const std::vector<Point*>& loop);
  static long double TriangleArea(const Triangle* t);
  static bool PointInLoop(const Point* p, const std::vector<Point*>& loop);

  std::vector<Point*> polyline_;
  std::vector<std::vector<Point*>> holes_;
  ParallelOptions options_;

  bool used_fallback_;
  size_t num_segments_;
  int split_axis_;
  std::vector<double> splits_;

  // Keeps every slab CDT (and thereby its triangles) alive, plus the
  // generated point copies.  Never touches the caller's Point objects.
  std::vector<std::unique_ptr<CDT>> contexts_;
  std::vector<std::unique_ptr<Point>> owned_points_;
  std::vector<Triangle*> triangles_;
};

inline void ParallelCDT::RunSequential() {
  contexts_.clear();
  owned_points_.clear();
  triangles_.clear();
  splits_.clear();
  num_segments_ = 1;
  split_axis_ = -1;
  used_fallback_ = true;

  std::unique_ptr<CDT> cdt(new CDT(polyline_));
  for (size_t i = 0; i < holes_.size(); ++i) {
    cdt->AddHole(holes_[i]);
  }
  cdt->Triangulate();
  triangles_ = cdt->GetTriangles();
  contexts_.push_back(std::move(cdt));
}

inline void ParallelCDT::Triangulate() {
  size_t num_points = polyline_.size();
  for (size_t i = 0; i < holes_.size(); ++i) {
    num_points += holes_[i].size();
  }

  size_t num_threads = options_.num_threads;
  if (num_threads == 0) {
    num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) {
      num_threads = 1;
    }
  }
  size_t num_segments =
      options_.num_segments != 0 ? options_.num_segments : num_threads;

  if (num_points < options_.min_points_for_parallel || num_segments < 2 ||
      num_threads < 2 || polyline_.size() < 3) {
    RunSequential();
    return;
  }

  bool ok = false;
  try {
    ok = TryParallel(num_segments, num_threads);
  } catch (...) {
    ok = false;
  }
  if (!ok) {
    // Stitching (or any other step) failed: discard partial work and
    // re-run the unmodified single-threaded path.
    RunSequential();
  }
}

inline bool ParallelCDT::TryParallel(size_t num_segments,
                                     size_t num_threads) {
  // Bounding box over every loop; split along the longer axis so slender
  // inputs are cut across their short dimension.
  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();
  std::vector<double> all;  // every vertex coordinate along the split axis
  {
    std::vector<const std::vector<Point*>*> loops;
    loops.push_back(&polyline_);
    for (size_t i = 0; i < holes_.size(); ++i) loops.push_back(&holes_[i]);
    for (size_t l = 0; l < loops.size(); ++l) {
      for (size_t i = 0; i < loops[l]->size(); ++i) {
        const Point* p = (*loops[l])[i];
        xmin = std::min(xmin, p->x);
        xmax = std::max(xmax, p->x);
        ymin = std::min(ymin, p->y);
        ymax = std::max(ymax, p->y);
      }
    }
  }
  const int axis = (xmax - xmin >= ymax - ymin) ? 0 : 1;
  {
    std::vector<const std::vector<Point*>*> loops;
    loops.push_back(&polyline_);
    for (size_t i = 0; i < holes_.size(); ++i) loops.push_back(&holes_[i]);
    for (size_t l = 0; l < loops.size(); ++l) {
      for (size_t i = 0; i < loops[l]->size(); ++i) {
        all.push_back(Coord((*loops[l])[i], axis));
      }
    }
  }
  std::sort(all.begin(), all.end());
  std::vector<double> uniq(all);
  uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

  if (uniq.size() < num_segments + 1) {
    return false;  // not enough distinct coordinates to cut between
  }

  // Split positions at vertex-count quantiles so every slab receives a
  // comparable number of points.  Each split sits at the midpoint of a
  // gap between distinct vertex coordinates, therefore no vertex and no
  // cut-parallel edge can ever lie exactly on a cut line.  Only gaps
  // wider than a minimum width are eligible: a cut line passing within a
  // few ulps of a vertex would produce degenerate sliver components
  // (clipped parameters rounding to exactly 0/1, crossing points
  // coinciding with existing points) that the sequential sweep cannot
  // handle.  The search walks outward from the ideal quantile gap, so
  // load balance is preserved whenever the input density allows it.
  std::vector<double> splits;
  const size_t n = all.size();
  const double range = uniq.back() - uniq.front();
  double min_gap = range * 1e-6;
  for (int attempt = 0; attempt < 3 && splits.empty(); ++attempt) {
    for (size_t k = 1; k < num_segments; ++k) {
      size_t rank = (k * n) / num_segments;
      if (rank >= n) rank = n - 1;
      const double target = all[rank];
      std::vector<double>::const_iterator it =
          std::upper_bound(uniq.begin(), uniq.end(), target);
      if (it == uniq.begin() || it == uniq.end()) continue;
      const size_t g0 = static_cast<size_t>(it - uniq.begin()) - 1;
      const size_t no_gap = uniq.size();
      size_t best = no_gap;
      for (size_t d = 0; d < uniq.size(); ++d) {
        bool any = false;
        if (d <= g0) {
          any = true;
          if (uniq[g0 - d + 1] - uniq[g0 - d] >= min_gap) {
            best = g0 - d;
            break;
          }
        }
        if (g0 + d < uniq.size() - 1) {
          any = true;
          if (uniq[g0 + d + 1] - uniq[g0 + d] >= min_gap) {
            best = g0 + d;
            break;
          }
        }
        if (!any) break;
      }
      if (best == no_gap) continue;
      splits.push_back(uniq[best] + (uniq[best + 1] - uniq[best]) * 0.5);
    }
    std::sort(splits.begin(), splits.end());
    splits.erase(std::unique(splits.begin(), splits.end()), splits.end());
    if (splits.empty()) {
      min_gap *= 1e-3;  // relax and retry; eventually the plain CDT runs
    }
  }
  if (splits.empty()) {
    return false;
  }

  const size_t num_slabs = splits.size() + 1;
  // Slab points closer together than this distance would risk degenerate
  // configurations inside the sequential sweep; such slabs are rejected
  // (and the whole run falls back to the plain CDT) instead.  The
  // threshold sits ~3 orders of magnitude above the floating-point
  // danger zone (~1e-13 * range) so legitimate close features - e.g.
  // cut-line crossings from a boundary running nearly parallel to a cut
  // - do not trigger it.
  const double min_dist =
      1e-10 * std::max(xmax - xmin, ymax - ymin);
  std::vector<SlabResult> results(num_slabs);
  std::atomic<size_t> next(0);
  const size_t workers = std::min(num_threads, num_slabs);

  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (size_t w = 0; w < workers; ++w) {
    pool.push_back(std::thread([this, &results, &next, axis, &splits,
                                min_dist]() {
      const size_t num_slabs = splits.size() + 1;
      for (;;) {
        const size_t i = next.fetch_add(1);
        if (i >= num_slabs) break;
        const double lo = (i == 0) ? -std::numeric_limits<double>::infinity()
                                   : splits[i - 1];
        const double hi = (i + 1 == num_slabs)
                              ? std::numeric_limits<double>::infinity()
                              : splits[i];
        SlabResult r;
        try {
          r = this->ProcessSlab(axis, lo, hi, min_dist);
        } catch (...) {
          r.ok = false;
        }
        results[i] = std::move(r);
      }
    }));
  }
  for (size_t w = 0; w < pool.size(); ++w) {
    pool[w].join();
  }

  long double total_area = 0.0L;
  for (size_t i = 0; i < results.size(); ++i) {
    if (!results[i].ok) {
      return false;
    }
    total_area += results[i].area;
  }

  // Global validation: the stitched triangle set must cover exactly the
  // same area as the input region (outer loop minus holes).
  long double expected = SignedArea(polyline_);
  if (expected < 0) expected = -expected;
  for (size_t i = 0; i < holes_.size(); ++i) {
    long double h = SignedArea(holes_[i]);
    expected -= (h < 0) ? -h : h;
  }
  if (expected <= 0.0L) {
    return false;
  }
  const long double rel =
      (total_area - expected) / expected;
  if (rel > 1e-9L || rel < -1e-9L) {
    return false;
  }

  // Commit.
  contexts_.clear();
  owned_points_.clear();
  triangles_.clear();
  for (size_t i = 0; i < results.size(); ++i) {
    for (size_t j = 0; j < results[i].cdts.size(); ++j) {
      std::vector<Triangle*> tris = results[i].cdts[j]->GetTriangles();
      triangles_.insert(triangles_.end(), tris.begin(), tris.end());
      contexts_.push_back(std::move(results[i].cdts[j]));
    }
    for (size_t j = 0; j < results[i].points.size(); ++j) {
      owned_points_.push_back(std::move(results[i].points[j]));
    }
  }
  splits_ = splits;
  split_axis_ = axis;
  num_segments_ = num_slabs;
  used_fallback_ = false;
  return true;
}

inline ParallelCDT::SlabResult ParallelCDT::ProcessSlab(int axis, double lo,
                                                        double hi,
                                                        double min_dist) const {
  SlabResult result;

  std::vector<Segment> segments;
  std::vector<Point*> cross_lo;  // cut-line crossings on the lo boundary
  std::vector<Point*> cross_hi;  // cut-line crossings on the hi boundary
  std::unordered_map<const Point*, Point*> copies;

  // Every point a slab uses is a private copy: p2t::Point carries mutable
  // sweep state (edge_list) and must not be shared between CDT instances.
  struct PointBag {
    std::vector<std::unique_ptr<Point>>* store;
    Point* Make(double x, double y) {
      store->push_back(std::unique_ptr<Point>(new Point(x, y)));
      return store->back().get();
    }
  } bag{&result.points};

  const bool lo_finite = std::isfinite(lo);
  const bool hi_finite = std::isfinite(hi);

  // Clip every loop edge to the slab, emitting in-slab pieces and
  // registering cut-line crossings.
  std::vector<const std::vector<Point*>*> loops;
  loops.push_back(&polyline_);
  for (size_t i = 0; i < holes_.size(); ++i) loops.push_back(&holes_[i]);

  for (size_t l = 0; l < loops.size(); ++l) {
    const std::vector<Point*>& loop = *loops[l];
    const size_t m = loop.size();
    for (size_t i = 0; i < m; ++i) {
      const Point* a = loop[i];
      const Point* b = loop[(i + 1) % m];
      const double ua = Coord(a, axis);
      const double ub = Coord(b, axis);

      Point* pa = NULL;
      Point* pb = NULL;

      if (ua == ub) {
        // Edge parallel to the cut lines: fully inside or outside.
        if (!(ua > lo && ua < hi)) continue;
        std::unordered_map<const Point*, Point*>::iterator ia = copies.find(a);
        if (ia == copies.end()) {
          pa = bag.Make(a->x, a->y);
          copies[a] = pa;
        } else {
          pa = ia->second;
        }
        std::unordered_map<const Point*, Point*>::iterator ib = copies.find(b);
        if (ib == copies.end()) {
          pb = bag.Make(b->x, b->y);
          copies[b] = pb;
        } else {
          pb = ib->second;
        }
      } else {
        // Parametric clip of a->b against lo < u < hi (infinite bounds
        // naturally produce infinite parameters).
        const double d = ub - ua;
        const double tlo = (lo - ua) / d;
        const double thi = (hi - ua) / d;
        const double te = std::max(0.0, std::min(tlo, thi));
        const double tl = std::min(1.0, std::max(tlo, thi));
        if (!(te < tl)) continue;

        for (int end = 0; end < 2; ++end) {
          const double t = (end == 0) ? te : tl;
          const Point* orig = (end == 0) ? a : b;
          Point** out = (end == 0) ? &pa : &pb;
          if (t <= 0.0 || t >= 1.0) {
            std::unordered_map<const Point*, Point*>::iterator it =
                copies.find(orig);
            if (it == copies.end()) {
              *out = bag.Make(orig->x, orig->y);
              copies[orig] = *out;
            } else {
              *out = it->second;
            }
          } else {
            Point* p = bag.Make(a->x + t * (b->x - a->x),
                                a->y + t * (b->y - a->y));
            // t equals exactly one of tlo/thi (same bits after clamping).
            if (t == tlo) {
              cross_lo.push_back(p);
            } else {
              cross_hi.push_back(p);
            }
            *out = p;
          }
        }
      }

      if (pa == pb) continue;  // degenerate piece
      Segment seg;
      seg.a = pa;
      seg.b = pb;
      segments.push_back(seg);
    }
  }

  // Divider edges on each cut line: sort crossings along the line and
  // pair them up (even-odd fill rule over outer loop + holes).
  for (int side = 0; side < 2; ++side) {
    if (side == 0 ? !lo_finite : !hi_finite) continue;
    std::vector<Point*>& cross = (side == 0) ? cross_lo : cross_hi;
    if (cross.size() % 2 != 0) return result;  // inconsistent clipping
    std::sort(cross.begin(), cross.end(), [axis](const Point* p,
                                                 const Point* q) {
      return axis == 0 ? (p->y < q->y) : (p->x < q->x);
    });
    for (size_t i = 0; i < cross.size(); i += 2) {
      if (Coord(cross[i], 1 - axis) == Coord(cross[i + 1], 1 - axis)) {
        return result;  // zero-length divider
      }
      Segment seg;
      seg.a = cross[i];
      seg.b = cross[i + 1];
      segments.push_back(seg);
    }
  }

  if (segments.empty()) return result;

  // Sanity: no two points handed to one CDT instance may be (nearly)
  // coincident - the sequential sweep cannot handle those, so reject the
  // slab and let the whole run fall back to the plain single-threaded
  // CDT instead.
  if (min_dist > 0.0) {
    std::vector<Point*> sorted(result.points.size());
    for (size_t i = 0; i < result.points.size(); ++i) {
      sorted[i] = result.points[i].get();
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const Point* p, const Point* q) { return p->x < q->x; });
    const double md2 = min_dist * min_dist;
    for (size_t i = 0; i < sorted.size(); ++i) {
      for (size_t j = i + 1;
           j < sorted.size() && sorted[j]->x - sorted[i]->x < min_dist;
           ++j) {
        const double dx = sorted[i]->x - sorted[j]->x;
        const double dy = sorted[i]->y - sorted[j]->y;
        if (dx * dx + dy * dy < md2) return result;
      }
    }
  }

  // The clipped pieces plus dividers form a graph where every vertex has
  // degree exactly 2, so loops can be extracted unambiguously.
  typedef std::pair<Point*, size_t> AdjEntry;
  std::unordered_map<Point*, std::vector<AdjEntry>> adj;
  adj.reserve(segments.size() * 2);
  for (size_t i = 0; i < segments.size(); ++i) {
    adj[segments[i].a].push_back(AdjEntry(segments[i].b, i));
    adj[segments[i].b].push_back(AdjEntry(segments[i].a, i));
  }
  for (std::unordered_map<Point*, std::vector<AdjEntry>>::const_iterator it =
           adj.begin();
       it != adj.end(); ++it) {
    if (it->second.size() != 2) return result;  // not a clean arrangement
  }

  std::vector<char> used(segments.size(), 0);
  std::vector<std::vector<Point*>> extracted;
  for (size_t start = 0; start < segments.size(); ++start) {
    if (used[start]) continue;
    std::vector<Point*> loop;
    Point* const first = segments[start].a;
    Point* cur = first;
    size_t prev_seg = start;
    used[start] = 1;
    loop.push_back(cur);
    cur = (segments[start].a == cur) ? segments[start].b : segments[start].a;
    size_t guard = 0;
    while (cur != first) {
      if (++guard > segments.size()) return result;  // corrupt graph
      loop.push_back(cur);
      std::vector<AdjEntry>& lst = adj[cur];
      const size_t s =
          (lst[0].second == prev_seg) ? lst[1].second : lst[0].second;
      used[s] = 1;
      prev_seg = s;
      cur = (segments[s].a == cur) ? segments[s].b : segments[s].a;
    }
    if (loop.size() < 3) return result;
    extracted.push_back(std::move(loop));
  }

  // Classify loops by nesting parity: even depth = component boundary,
  // odd depth = hole.
  std::vector<long double> areas(extracted.size());
  std::vector<int> depth(extracted.size(), 0);
  for (size_t i = 0; i < extracted.size(); ++i) {
    areas[i] = SignedArea(extracted[i]);
    if (areas[i] < 0) areas[i] = -areas[i];
    if (areas[i] == 0.0L) return result;  // degenerate loop
  }
  for (size_t i = 0; i < extracted.size(); ++i) {
    for (size_t j = 0; j < extracted.size(); ++j) {
      if (i != j && PointInLoop(extracted[i][0], extracted[j])) {
        ++depth[i];
      }
    }
  }

  // Group holes with their innermost containing boundary loop.
  std::vector<std::vector<size_t>> comp_holes(extracted.size());
  std::vector<size_t> outers;
  for (size_t i = 0; i < extracted.size(); ++i) {
    if (depth[i] % 2 == 0) {
      outers.push_back(i);
    }
  }
  if (outers.empty()) return result;
  for (size_t i = 0; i < extracted.size(); ++i) {
    if (depth[i] % 2 == 0) continue;
    size_t parent = extracted.size();
    for (size_t j = 0; j < outers.size(); ++j) {
      const size_t o = outers[j];
      if (PointInLoop(extracted[i][0], extracted[o]) &&
          (parent == extracted.size() || depth[o] > depth[parent])) {
        parent = o;
      }
    }
    if (parent == extracted.size()) return result;  // orphan hole
    comp_holes[parent].push_back(i);
  }

  // Triangulate every component with the unmodified sequential CDT.
  for (size_t c = 0; c < outers.size(); ++c) {
    const size_t o = outers[c];
    std::unique_ptr<CDT> cdt(new CDT(extracted[o]));
    long double expect = areas[o];
    for (size_t h = 0; h < comp_holes[o].size(); ++h) {
      cdt->AddHole(extracted[comp_holes[o][h]]);
      expect -= areas[comp_holes[o][h]];
    }
    if (expect <= 0.0L) return result;
    cdt->Triangulate();

    std::vector<Triangle*> tris = cdt->GetTriangles();
    long double sum = 0.0L;
    for (size_t t = 0; t < tris.size(); ++t) {
      sum += TriangleArea(tris[t]);
    }
    // Per-component validation: the triangulation must cover exactly the
    // component region, otherwise the stitch is rejected.
    const long double rel = (sum - expect) / expect;
    if (rel > 1e-9L || rel < -1e-9L) return result;

    result.area += sum;
    result.cdts.push_back(std::move(cdt));
  }

  result.ok = true;
  return result;
}

inline long double ParallelCDT::SignedArea(const std::vector<Point*>& loop) {
  // Fan triangulation around a local origin (the first vertex): keeps
  // every cross product at local scale so large coordinates do not
  // introduce catastrophic cancellation (matters for tiny sliver
  // components clipped off near a cut line).
  const size_t n = loop.size();
  const long double ox = loop[0]->x;
  const long double oy = loop[0]->y;
  long double area = 0.0L;
  for (size_t i = 1; i + 1 < n; ++i) {
    const long double ax = (long double)loop[i]->x - ox;
    const long double ay = (long double)loop[i]->y - oy;
    const long double bx = (long double)loop[i + 1]->x - ox;
    const long double by = (long double)loop[i + 1]->y - oy;
    area += ax * by - bx * ay;
  }
  return area * 0.5L;
}

inline long double ParallelCDT::TriangleArea(const Triangle* t) {
  const Point* a = const_cast<Triangle*>(t)->GetPoint(0);
  const Point* b = const_cast<Triangle*>(t)->GetPoint(1);
  const Point* c = const_cast<Triangle*>(t)->GetPoint(2);
  const long double v =
      ((long double)b->x - a->x) * ((long double)c->y - a->y) -
      ((long double)b->y - a->y) * ((long double)c->x - a->x);
  return v < 0 ? -v * 0.5L : v * 0.5L;
}

inline bool ParallelCDT::PointInLoop(const Point* p,
                                     const std::vector<Point*>& loop) {
  bool inside = false;
  const size_t n = loop.size();
  for (size_t i = 0, j = n - 1; i < n; j = i++) {
    const Point* a = loop[i];
    const Point* b = loop[j];
    if ((a->y > p->y) != (b->y > p->y) &&
        p->x < (b->x - a->x) * (p->y - a->y) / (b->y - a->y) + a->x) {
      inside = !inside;
    }
  }
  return inside;
}

} // namespace p2t

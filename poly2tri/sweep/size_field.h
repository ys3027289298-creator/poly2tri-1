/*
 * Poly2Tri Copyright (c) 2009-2026, Poly2Tri Contributors
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

#include "../common/dll_symbol.h"
#include "cdt.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace p2t {

/**
 * Size field: target edge length as a function of a point of the domain.
 * Every evaluated value must be finite and strictly positive; a zero,
 * negative or non-finite value aborts the run with std::runtime_error.
 */
typedef std::function<double(const Point&)> SizeFieldFunction;

/**
 * Options of SizeFieldRefiner::Triangulate().
 *
 * An edge is measured against the size field evaluated at its midpoint
 * (multiplied by `scale`):
 *  - an edge longer than split_ratio * target is split at its midpoint,
 *  - an edge shorter than merge_ratio * target is collapsed when one of its
 *    endpoints was created by the refiner (or is a dense Steiner point and
 *    merge_dense_steiner_points is set).
 *
 * merge_ratio must not exceed split_ratio / 2, so a freshly split edge (now
 * about split_ratio / 2 * target long) is never immediately mergeable: split
 * and merge cannot oscillate between rounds.
 */
struct P2T_DLL_SYMBOL SizeFieldOptions {
  /// Multiplier applied to every value returned by the size field.
  double scale = 1.0;
  /// Maximum number of split/merge + re-triangulation rounds.
  int max_iterations = 8;
  /// Edges longer than split_ratio * target length are split.
  double split_ratio = 1.25;
  /// Edges shorter than merge_ratio * target length are merge candidates.
  double merge_ratio = 0.55;
  /**
   * When true, caller supplied Steiner points with at least one incident edge
   * shorter than merge_ratio * target are dropped from the mesh (removed from
   * the triangulation only; ownership stays with the caller). This thins an
   * over-dense input point cloud. When false (default) caller points are
   * never removed.
   */
  bool merge_dense_steiner_points = false;
};

/**
 * Outcome of one SizeFieldRefiner::Triangulate() run.
 */
struct P2T_DLL_SYMBOL SizeFieldReport {
  /// True when no split/merge was pending after the last round.
  bool converged = true;
  /// Number of split/merge + re-triangulation rounds that ran.
  int iterations = 0;
  /// Iteration cap the run was given.
  int max_iterations = 0;
  /// Number of triangles in the final mesh.
  std::size_t triangle_count = 0;
  /// Number of distinct mesh vertices.
  std::size_t vertex_count = 0;
  /// Points inserted into the contours (boundary/hole segments) by the refiner.
  std::size_t boundary_points_added = 0;
  /// Interior points inserted by the refiner.
  std::size_t interior_points_added = 0;
  /// Points removed because all their edges were too short.
  std::size_t points_merged = 0;
  /// Of points_merged, the caller supplied Steiner points that were dropped
  /// (only possible with merge_dense_steiner_points).
  std::size_t steiner_points_merged = 0;
  /**
   * Deviation between the produced edge length distribution and the target
   * distribution. For every unique edge e with length L and target h(mid(e))
   * the relative deviation is |L / h - 1|; the fields below summarize it.
   */
  double mean_relative_deviation = 0.0;
  double median_relative_deviation = 0.0;
  double p90_relative_deviation = 0.0;
  /// Fraction of edges with merge_ratio <= L/h <= split_ratio.
  double fraction_within_band = 0.0;
  /// Smallest and largest observed L / h ratio.
  double min_length_ratio = 0.0;
  double max_length_ratio = 0.0;
  /// Mean target length over the unique edges.
  double target_mean_edge_length = 0.0;
  /// Mean produced edge length over the unique edges.
  double actual_mean_edge_length = 0.0;
  /// Area of the input polygon (outer minus holes) and of the final mesh.
  double polygon_area = 0.0;
  double mesh_area = 0.0;
  /// |mesh_area - polygon_area| / polygon_area.
  double relative_area_error = 0.0;
  /// Wall clock time of the initial triangulation and of the whole run.
  double initial_triangulation_seconds = 0.0;
  double total_seconds = 0.0;
  /// Wall clock time of every refinement round (split + rebuild + checks).
  std::vector<double> iteration_seconds;
};

/**
 * Constrained Delaunay triangulator driven by a size field.
 *
 * Usage mirrors p2t::CDT:
 *
 *   p2t::SizeFieldRefiner refiner(polyline);
 *   refiner.AddHole(hole);
 *   refiner.AddPoint(steiner);
 *   auto field = [](const p2t::Point& p) { return 0.2 + 0.1 * p.x; };
 *   const auto& report = refiner.Triangulate(field, options);
 *   auto triangles = refiner.GetTriangles();
 *
 * Each round starts from a fresh constrained Delaunay triangulation of the
 * current points, so the constrained boundary is always re-imposed and a
 * split point lying on a boundary/hole segment is projected back onto that
 * segment before it is inserted into the contour. Interior split points are
 * kept away from the contours to avoid degenerate triangles. The boundary
 * is verified not to be crossed after every round; a crossing aborts the run
 * with std::runtime_error.
 *
 * Points passed to the constructor/AddHole/AddPoint stay owned by the caller
 * (same convention as p2t::CDT). Points inserted by the refiner are owned and
 * deleted by it.
 *
 * When the field and all edges already agree (no split or merge is pending),
 * the first round is returned untouched, i.e. the output is vertex-by-vertex
 * identical to a plain p2t::CDT triangulation of the same input.
 */
class P2T_DLL_SYMBOL SizeFieldRefiner {
public:
  /// Constructor - add polyline with non repeating points
  explicit SizeFieldRefiner(const std::vector<Point*>& polyline)
  {
    polylines_.push_back(polyline);
  }

  ~SizeFieldRefiner()
  {
    delete cdt_;
    for (Point* point : owned_) {
      delete point;
    }
  }

  SizeFieldRefiner(const SizeFieldRefiner&) = delete;
  SizeFieldRefiner& operator=(const SizeFieldRefiner&) = delete;

  /// Add a hole
  void AddHole(const std::vector<Point*>& polyline)
  {
    polylines_.push_back(polyline);
  }

  /// Add a steiner point
  void AddPoint(Point* point)
  {
    steiner_.push_back(point);
  }

  /**
   * Triangulate and refine according to the size field.
   * Throws std::runtime_error when the field or the options are invalid.
   */
  const SizeFieldReport& Triangulate(const SizeFieldFunction& field,
                                     const SizeFieldOptions& options = SizeFieldOptions())
  {
    Validate(field, options);
    report_ = SizeFieldReport();
    report_.max_iterations = options.max_iterations;

    const auto start = Clock::now();
    RebuildCdt();
    report_.initial_triangulation_seconds = SecondsSince(start);

    while (true) {
      const auto round_start = Clock::now();
      const Changes changes = FindChanges(field, options);
      if (changes.Empty()) {
        report_.converged = true;
        break;
      }
      if (report_.iterations >= options.max_iterations) {
        report_.converged = false;
        break;
      }
      const std::size_t applied = ApplyChanges(changes);
      if (applied == 0) {
        // Every candidate was a duplicate/too close to be inserted.
        report_.converged = false;
        break;
      }
      ++report_.iterations;
      RebuildCdt();
      report_.iteration_seconds.push_back(SecondsSince(round_start));
    }

    ComputeReport(field, options);
    report_.total_seconds = SecondsSince(start);
    return report_;
  }

  /// Get triangles of the final mesh
  std::vector<Triangle*> GetTriangles()
  {
    return cdt_ ? cdt_->GetTriangles() : std::vector<Triangle*>();
  }

  /// Report of the last Triangulate call
  const SizeFieldReport& report() const
  {
    return report_;
  }

private:
  typedef std::chrono::steady_clock Clock;

  /// Unordered pair of point pointers; the canonical smaller pointer is first.
  struct EdgeKey {
    Point* first;
    Point* second;

    EdgeKey(Point* a, Point* b) : first(a), second(b)
    {
      if (std::less<Point*>()(second, first)) {
        std::swap(first, second);
      }
    }

    bool operator<(const EdgeKey& other) const
    {
      if (std::less<Point*>()(first, other.first)) {
        return true;
      }
      if (std::less<Point*>()(other.first, first)) {
        return false;
      }
      return std::less<Point*>()(second, other.second);
    }
  };

  struct Split {
    Point* a;
    Point* b;
    bool constrained;

    Split(Point* first, Point* second, bool is_constrained)
        : a(first), b(second), constrained(is_constrained)
    {
    }
  };

  struct Changes {
    std::vector<Split> splits;
    std::vector<Point*> merged;

    bool Empty() const
    {
      return splits.empty() && merged.empty();
    }
  };

  static double Distance(const Point& a, const Point& b)
  {
    return std::hypot(a.x - b.x, a.y - b.y);
  }

  /**
   * Project `point` onto the segment a-b. The projection parameter is clamped
   * away from the endpoints so a split never reproduces an existing vertex.
   */
  static Point ProjectOntoSegment(const Point& point, const Point& a, const Point& b)
  {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double length2 = dx * dx + dy * dy;
    if (!(length2 > 0.0)) {
      return Point(a.x, a.y);
    }
    double t = ((point.x - a.x) * dx + (point.y - a.y) * dy) / length2;
    const double t_min = 1e-6;
    t = std::max(t_min, std::min(1.0 - t_min, t));
    return Point(a.x + t * dx, a.y + t * dy);
  }

  static double DistanceToSegment(const Point& point, const Point& a, const Point& b)
  {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double length2 = dx * dx + dy * dy;
    if (!(length2 > 0.0)) {
      return Distance(point, a);
    }
    double t = ((point.x - a.x) * dx + (point.y - a.y) * dy) / length2;
    t = std::max(0.0, std::min(1.0, t));
    const Point closest(a.x + t * dx, a.y + t * dy);
    return Distance(point, closest);
  }

  static double SignedArea(const std::vector<Point*>& polyline)
  {
    double area = 0.0;
    const std::size_t n = polyline.size();
    for (std::size_t i = 0; i < n; ++i) {
      const Point& p = *polyline[i];
      const Point& q = *polyline[(i + 1) % n];
      area += p.x * q.y - q.x * p.y;
    }
    return 0.5 * area;
  }

  static double TriangleArea(Triangle& triangle)
  {
    const Point& a = *triangle.GetPoint(0);
    const Point& b = *triangle.GetPoint(1);
    const Point& c = *triangle.GetPoint(2);
    return 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) -
                          (c.x - a.x) * (b.y - a.y));
  }

  static double SecondsSince(const Clock::time_point& start)
  {
    return std::chrono::duration<double>(Clock::now() - start).count();
  }

  std::vector<Point*> AllPoints() const
  {
    std::vector<Point*> all;
    for (const std::vector<Point*>& polyline : polylines_) {
      all.insert(all.end(), polyline.begin(), polyline.end());
    }
    all.insert(all.end(), steiner_.begin(), steiner_.end());
    all.insert(all.end(), interior_.begin(), interior_.end());
    return all;
  }

  double BoundingDiagonal() const
  {
    const std::vector<Point*> points = AllPoints();
    if (points.empty()) {
      return 0.0;
    }
    double xmin = points[0]->x;
    double xmax = points[0]->x;
    double ymin = points[0]->y;
    double ymax = points[0]->y;
    for (const Point* point : points) {
      xmin = std::min(xmin, point->x);
      xmax = std::max(xmax, point->x);
      ymin = std::min(ymin, point->y);
      ymax = std::max(ymax, point->y);
    }
    return std::hypot(xmax - xmin, ymax - ymin);
  }

  /// Evaluate and validate the field at (x, y); returns scale * value.
  double Evaluate(const SizeFieldFunction& field, double x, double y, double scale) const
  {
    const Point sample(x, y);
    const double value = field(sample);
    if (!std::isfinite(value) || value <= 0.0) {
      std::ostringstream message;
      message << "SizeFieldRefiner: size field returned " << value << " at ("
              << x << ", " << y
              << "); the target edge length must be positive and finite";
      throw std::runtime_error(message.str());
    }
    return value * scale;
  }

  void Validate(const SizeFieldFunction& field, const SizeFieldOptions& options) const
  {
    if (!field) {
      throw std::runtime_error("SizeFieldRefiner: size field function is empty");
    }
    if (!std::isfinite(options.scale) || !(options.scale > 0.0)) {
      throw std::runtime_error("SizeFieldRefiner: scale must be positive and finite");
    }
    if (options.max_iterations < 0) {
      throw std::runtime_error("SizeFieldRefiner: max_iterations must not be negative");
    }
    if (!std::isfinite(options.split_ratio) || !(options.split_ratio > 1.0)) {
      throw std::runtime_error("SizeFieldRefiner: split_ratio must be greater than 1");
    }
    // An edge is split when it is longer than split_ratio * target; the two
    // halves are then split_ratio / 2 * target long. Requiring
    // merge_ratio <= split_ratio / 2 guarantees a freshly split edge is not
    // immediately mergeable again, i.e. split and merge never oscillate.
    if (!std::isfinite(options.merge_ratio) || options.merge_ratio < 0.0 ||
        options.merge_ratio > 0.5 * options.split_ratio) {
      throw std::runtime_error(
          "SizeFieldRefiner: merge_ratio must be >= 0 and at most split_ratio / 2");
    }
    // Contours need at least three points (something to close) and non-empty
    // point sets are needed for the sampling below.
    if (polylines_.empty() || polylines_[0].size() < 3) {
      throw std::runtime_error(
          "SizeFieldRefiner: the outer polyline needs at least 3 points");
    }
    for (std::size_t i = 1; i < polylines_.size(); ++i) {
      if (polylines_[i].size() < 3) {
        std::ostringstream message;
        message << "SizeFieldRefiner: hole " << (i - 1)
                << " needs at least 3 points";
        throw std::runtime_error(message.str());
      }
    }

    // Reject fields which are illegal anywhere on a coarse sample of the
    // domain, not only where the mesh happens to ask for a value.
    const std::vector<Point*> points = AllPoints();
    double xmin = points[0]->x;
    double xmax = points[0]->x;
    double ymin = points[0]->y;
    double ymax = points[0]->y;
    for (const Point* point : points) {
      xmin = std::min(xmin, point->x);
      xmax = std::max(xmax, point->x);
      ymin = std::min(ymin, point->y);
      ymax = std::max(ymax, point->y);
    }
    const int samples = 8;
    for (int i = 0; i <= samples; ++i) {
      for (int j = 0; j <= samples; ++j) {
        const double x = xmin + (xmax - xmin) * static_cast<double>(i) / samples;
        const double y = ymin + (ymax - ymin) * static_cast<double>(j) / samples;
        Evaluate(field, x, y, options.scale);
      }
    }
    for (const Point* point : points) {
      Evaluate(field, point->x, point->y, options.scale);
    }
  }

  void RebuildCdt()
  {
    // The previous SweepContext owned the Edge objects stored in
    // Point::edge_list and freed them with itself; clear the stale pointers
    // before the points are inserted into a new context.
    delete cdt_;
    cdt_ = nullptr;
    for (Point* point : AllPoints()) {
      point->edge_list.clear();
    }

    cdt_ = new CDT(polylines_[0]);
    for (std::size_t i = 1; i < polylines_.size(); ++i) {
      cdt_->AddHole(polylines_[i]);
    }
    for (Point* point : steiner_) {
      cdt_->AddPoint(point);
    }
    for (Point* point : interior_) {
      cdt_->AddPoint(point);
    }
    cdt_->Triangulate();
    VerifyConstraints();
  }

  /// Collect the constrained segments of all contours.
  std::set<EdgeKey> ContourSegments() const
  {
    std::set<EdgeKey> segments;
    for (const std::vector<Point*>& polyline : polylines_) {
      const std::size_t n = polyline.size();
      for (std::size_t i = 0; i < n; ++i) {
        segments.insert(EdgeKey(polyline[i], polyline[(i + 1) % n]));
      }
    }
    return segments;
  }

  Changes FindChanges(const SizeFieldFunction& field,
                      const SizeFieldOptions& options) const
  {
    Changes changes;
    const std::set<EdgeKey> contour = ContourSegments();
    std::set<EdgeKey> seen;
    // Per point: (number of incident unique edges, number of short ones).
    std::map<Point*, std::pair<std::size_t, std::size_t>, std::less<Point*>> incident;

    const std::vector<Triangle*>& triangles = cdt_->GetTriangles();
    for (Triangle* triangle : triangles) {
      for (int i = 0; i < 3; ++i) {
        Point* a = triangle->GetPoint(i);
        Point* b = triangle->GetPoint((i + 1) % 3);
        const EdgeKey key(a, b);
        if (!seen.insert(key).second) {
          continue;
        }
        const double length = Distance(*a, *b);
        const double target =
            Evaluate(field, 0.5 * (a->x + b->x), 0.5 * (a->y + b->y), options.scale);
        const bool constrained = contour.find(key) != contour.end();
        if (length > options.split_ratio * target) {
          changes.splits.push_back(Split(key.first, key.second, constrained));
        }
        const bool is_short = length < options.merge_ratio * target;
        std::pair<std::size_t, std::size_t>& stats_a = incident[a];
        std::pair<std::size_t, std::size_t>& stats_b = incident[b];
        ++stats_a.first;
        ++stats_b.first;
        if (is_short) {
          ++stats_a.second;
          ++stats_b.second;
        }
      }
    }

    // Interior points created by an earlier round: merge when every
    // incident edge is short. Points supplied by the caller as Steiner points
    // are handled by the option below.
    for (Point* point : interior_) {
      const auto found = incident.find(point);
      if (found != incident.end() && found->second.first > 0 &&
          found->second.second == found->second.first) {
        changes.merged.push_back(point);
      }
    }
    if (options.merge_dense_steiner_points) {
      // Explicitly requested thinning: a caller supplied Steiner point is
      // dropped as soon as one of its incident edges is short, so an
      // over-dense input cloud relaxes towards the requested spacing instead
      // of keeping pairs of short edges that no point can absorb.
      for (Point* point : steiner_) {
        const auto found = incident.find(point);
        if (found != incident.end() && found->second.second > 0) {
          changes.merged.push_back(point);
        }
      }
    }

    // Contour points created by an earlier round: merge when both adjacent
    // segments are short and the merged segment would not need splitting.
    // At most every other point is merged per round so a chain of short
    // segments coarsens by one level per round instead of all at once.
    for (const std::vector<Point*>& polyline : polylines_) {
      const std::size_t n = polyline.size();
      if (n < 4) {
        continue;
      }
      std::size_t index = 0;
      while (index < n) {
        Point* point = polyline[index];
        if (owned_.find(point) == owned_.end()) {
          ++index;
          continue;
        }
        Point* prev = polyline[(index + n - 1) % n];
        Point* next = polyline[(index + 1) % n];
        const double length_prev = Distance(*prev, *point);
        const double length_next = Distance(*point, *next);
        const double target_prev = Evaluate(
            field, 0.5 * (prev->x + point->x), 0.5 * (prev->y + point->y), options.scale);
        const double target_next = Evaluate(
            field, 0.5 * (point->x + next->x), 0.5 * (point->y + next->y), options.scale);
        const bool mergeable =
            length_prev < options.merge_ratio * target_prev &&
            length_next < options.merge_ratio * target_next;
        if (mergeable) {
          const double merged_length = Distance(*prev, *next);
          const double merged_target = Evaluate(
              field, 0.5 * (prev->x + next->x), 0.5 * (prev->y + next->y),
              options.scale);
          if (merged_length <= options.split_ratio * merged_target) {
            changes.merged.push_back(point);
            index += 2; // never remove two adjacent contour points at once
            continue;
          }
        }
        ++index;
      }
    }
    return changes;
  }

  /// Remove a refiner-created contour/interior point, or a dense Steiner point.
  /// `was_steiner` reports whether the removed point was caller supplied.
  bool RemovePoint(Point* point, bool& was_steiner)
  {
    bool found = false;
    const auto interior_it = std::find(interior_.begin(), interior_.end(), point);
    if (interior_it != interior_.end()) {
      interior_.erase(interior_it);
      found = true;
    } else {
      const auto steiner_it = std::find(steiner_.begin(), steiner_.end(), point);
      if (steiner_it != steiner_.end()) {
        steiner_.erase(steiner_it);
        found = true;
        was_steiner = true;
      } else {
        for (std::vector<Point*>& polyline : polylines_) {
          const auto point_it = std::find(polyline.begin(), polyline.end(), point);
          if (point_it != polyline.end()) {
            polyline.erase(point_it);
            found = true;
            break;
          }
        }
      }
    }
    if (owned_.erase(point) > 0) {
      delete point;
    }
    return found;
  }

  /// Locate the contour segment (a, b) in either direction.
  bool FindSegment(Point* a, Point* b, std::size_t& polyline_index,
                   std::size_t& position) const
  {
    for (std::size_t k = 0; k < polylines_.size(); ++k) {
      const std::vector<Point*>& polyline = polylines_[k];
      const std::size_t n = polyline.size();
      for (std::size_t i = 0; i < n; ++i) {
        Point* first = polyline[i];
        Point* second = polyline[(i + 1) % n];
        if ((first == a && second == b) || (first == b && second == a)) {
          polyline_index = k;
          position = i;
          return true;
        }
      }
    }
    return false;
  }

  bool TooCloseToPoint(const Point& candidate, const std::vector<Point*>& points,
                       double epsilon) const
  {
    for (const Point* point : points) {
      if (Distance(candidate, *point) <= epsilon) {
        return true;
      }
    }
    return false;
  }

  bool TooCloseToContour(const Point& candidate, double epsilon) const
  {
    for (const std::vector<Point*>& polyline : polylines_) {
      const std::size_t n = polyline.size();
      for (std::size_t i = 0; i < n; ++i) {
        if (DistanceToSegment(candidate, *polyline[i], *polyline[(i + 1) % n]) <= epsilon) {
          return true;
        }
      }
    }
    return false;
  }

  std::size_t ApplyChanges(const Changes& changes)
  {
    std::size_t applied = 0;

    // Removals first: the split candidates were measured on the old mesh.
    std::set<Point*> removed;
    for (Point* point : changes.merged) {
      bool removable = true;
      for (const std::vector<Point*>& polyline : polylines_) {
        const std::size_t n = polyline.size();
        const auto position = std::find(polyline.begin(), polyline.end(), point);
        if (position == polyline.end()) {
          continue;
        }
        // Never shrink a contour to fewer than three points (a closed ring
        // needs three) and never remove two adjacent contour points in the
        // same round.
        const std::size_t index =
            static_cast<std::size_t>(position - polyline.begin());
        if (n <= 3 ||
            removed.count(polyline[(index + n - 1) % n]) > 0 ||
            removed.count(polyline[(index + 1) % n]) > 0) {
          removable = false;
        }
        break;
      }
      if (!removable) {
        continue;
      }
      bool was_steiner = false;
      if (RemovePoint(point, was_steiner)) {
        removed.insert(point);
        ++applied;
        ++report_.points_merged;
        if (was_steiner) {
          ++report_.steiner_points_merged;
        }
      }
    }

    const double diagonal = BoundingDiagonal();
    const double epsilon = std::max(1e-12, 1e-8 * diagonal);
    const double contour_gap = std::max(1e-12, 1e-6 * diagonal);
    std::vector<Point*> known = AllPoints();

    for (const Split& split : changes.splits) {
      // A split of an edge touching a point merged above would use a freed
      // point; skip it, the next round measures the new mesh anyway.
      if (removed.count(split.a) > 0 || removed.count(split.b) > 0) {
        continue;
      }
      if (split.constrained) {
        std::size_t polyline_index = 0;
        std::size_t position = 0;
        if (!FindSegment(split.a, split.b, polyline_index, position)) {
          continue;
        }
        std::vector<Point*>& polyline = polylines_[polyline_index];
        Point* a = polyline[position];
        Point* b = polyline[(position + 1) % polyline.size()];
        const Point midpoint(0.5 * (a->x + b->x), 0.5 * (a->y + b->y));
        // A split point of a constrained edge belongs to the contour: put it
        // exactly back onto the segment before it is inserted.
        const Point projected = ProjectOntoSegment(midpoint, *a, *b);
        if (TooCloseToPoint(projected, known, epsilon)) {
          continue;
        }
        Point* point = new Point(projected.x, projected.y);
        polyline.insert(polyline.begin() + static_cast<std::ptrdiff_t>(position) + 1,
                        point);
        owned_.insert(point);
        known.push_back(point);
        ++report_.boundary_points_added;
        ++applied;
      } else {
        const Point midpoint(0.5 * (split.a->x + split.b->x),
                             0.5 * (split.a->y + split.b->y));
        if (TooCloseToPoint(midpoint, known, epsilon) ||
            TooCloseToContour(midpoint, contour_gap)) {
          continue;
        }
        Point* point = new Point(midpoint.x, midpoint.y);
        interior_.push_back(point);
        owned_.insert(point);
        known.push_back(point);
        ++report_.interior_points_added;
        ++applied;
      }
    }
    return applied;
  }

  /**
   * Abort if any constrained contour segment is crossed by a triangle edge.
   * Edges that share a point with the segment are not crossings.
   */
  void VerifyConstraints() const
  {
    const double diagonal = std::max(1.0, BoundingDiagonal());
    const double orientation_epsilon = 1e-12 * diagonal * diagonal;

    struct Segment {
      const Point* a;
      const Point* b;
    };
    std::vector<Segment> segments;
    for (const std::vector<Point*>& polyline : polylines_) {
      const std::size_t n = polyline.size();
      for (std::size_t i = 0; i < n; ++i) {
        segments.push_back(Segment{polyline[i], polyline[(i + 1) % n]});
      }
    }

    for (Triangle* triangle : cdt_->GetTriangles()) {
      for (int i = 0; i < 3; ++i) {
        const Point* a = triangle->GetPoint(i);
        const Point* b = triangle->GetPoint((i + 1) % 3);
        for (const Segment& segment : segments) {
          if (a == segment.a || a == segment.b || b == segment.a || b == segment.b) {
            continue;
          }
          const double o1 = Orient(*a, *b, *segment.a);
          const double o2 = Orient(*a, *b, *segment.b);
          const double o3 = Orient(*segment.a, *segment.b, *a);
          const double o4 = Orient(*segment.a, *segment.b, *b);
          const bool crosses =
              ((o1 > orientation_epsilon && o2 < -orientation_epsilon) ||
               (o1 < -orientation_epsilon && o2 > orientation_epsilon)) &&
              ((o3 > orientation_epsilon && o4 < -orientation_epsilon) ||
               (o3 < -orientation_epsilon && o4 > orientation_epsilon));
          if (crosses) {
            std::ostringstream message;
            message << "SizeFieldRefiner: constrained segment (" << segment.a->x
                    << ", " << segment.a->y << ")-(" << segment.b->x << ", "
                    << segment.b->y << ") is crossed by triangle edge ("
                    << a->x << ", " << a->y << ")-(" << b->x << ", " << b->y
                    << ")";
            throw std::runtime_error(message.str());
          }
        }
      }
    }
  }

  static double Orient(const Point& a, const Point& b, const Point& c)
  {
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
  }

  void ComputeReport(const SizeFieldFunction& field, const SizeFieldOptions& options)
  {
    const std::vector<Triangle*>& triangles = cdt_->GetTriangles();
    report_.triangle_count = triangles.size();

    std::set<Point*> vertices;
    std::set<EdgeKey> edges;
    double mesh_area = 0.0;
    for (Triangle* triangle : triangles) {
      for (int i = 0; i < 3; ++i) {
        vertices.insert(triangle->GetPoint(i));
      }
      edges.insert(EdgeKey(triangle->GetPoint(0), triangle->GetPoint(1)));
      edges.insert(EdgeKey(triangle->GetPoint(1), triangle->GetPoint(2)));
      edges.insert(EdgeKey(triangle->GetPoint(2), triangle->GetPoint(0)));
      mesh_area += TriangleArea(*triangle);
    }
    report_.vertex_count = vertices.size();
    report_.mesh_area = mesh_area;

    const double outer_area = std::fabs(SignedArea(polylines_[0]));
    double holes_area = 0.0;
    for (std::size_t i = 1; i < polylines_.size(); ++i) {
      holes_area += std::fabs(SignedArea(polylines_[i]));
    }
    report_.polygon_area = std::max(0.0, outer_area - holes_area);
    if (report_.polygon_area > 0.0) {
      report_.relative_area_error =
          std::fabs(mesh_area - report_.polygon_area) / report_.polygon_area;
    }

    std::vector<double> deviations;
    deviations.reserve(edges.size());
    double target_sum = 0.0;
    double actual_sum = 0.0;
    double min_ratio = 0.0;
    double max_ratio = 0.0;
    std::size_t within_band = 0;
    bool first = true;
    for (const EdgeKey& edge : edges) {
      const double length = Distance(*edge.first, *edge.second);
      const double target = Evaluate(field, 0.5 * (edge.first->x + edge.second->x),
                                     0.5 * (edge.first->y + edge.second->y),
                                     options.scale);
      const double ratio = length / target;
      deviations.push_back(std::fabs(ratio - 1.0));
      target_sum += target;
      actual_sum += length;
      if (first) {
        min_ratio = ratio;
        max_ratio = ratio;
        first = false;
      } else {
        min_ratio = std::min(min_ratio, ratio);
        max_ratio = std::max(max_ratio, ratio);
      }
      if (ratio >= options.merge_ratio && ratio <= options.split_ratio) {
        ++within_band;
      }
    }

    const std::size_t count = deviations.size();
    if (count > 0) {
      double sum = 0.0;
      for (double deviation : deviations) {
        sum += deviation;
      }
      report_.mean_relative_deviation = sum / count;
      std::vector<double> sorted = deviations;
      std::sort(sorted.begin(), sorted.end());
      report_.median_relative_deviation = sorted[(count - 1) / 2];
      report_.p90_relative_deviation =
          sorted[std::min(count - 1, static_cast<std::size_t>(std::floor(0.9 * (count - 1))))];
      report_.fraction_within_band =
          static_cast<double>(within_band) / static_cast<double>(count);
      report_.target_mean_edge_length = target_sum / count;
      report_.actual_mean_edge_length = actual_sum / count;
    }
    report_.min_length_ratio = min_ratio;
    report_.max_length_ratio = max_ratio;
  }

  std::vector<std::vector<Point*>> polylines_;
  std::vector<Point*> steiner_;
  std::vector<Point*> interior_;
  std::set<Point*> owned_;
  CDT* cdt_ = nullptr;
  SizeFieldReport report_;
};

} // namespace p2t

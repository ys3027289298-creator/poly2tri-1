/*
 * Poly2Tri Copyright (c) 2009-2018, Poly2Tri Contributors
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

#include "advancing_front.h"
#include "sweep_context.h"
#include "sweep.h"

#include "../common/dll_symbol.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <utility>

/**
 *
 * @author Mason Green <mason.green@gmail.com>
 *
 */

namespace p2t {

class P2T_DLL_SYMBOL CDT
{
public:

  /**
   * Constructor - add polyline with non repeating points
   *
   * @param polyline
   */
  CDT(const std::vector<Point*>& polyline);

   /**
   * Destructor - clean up memory
   */
  ~CDT();

  /**
   * Add a hole
   *
   * @param polyline
   */
  void AddHole(const std::vector<Point*>& polyline);

  /**
   * Add a steiner point
   *
   * @param point
   */
  void AddPoint(Point* point);

  /**
   * Triangulate - do this AFTER you've added the polyline, holes, and Steiner points
   */
  void Triangulate();

  /**
   * Get CDT triangles
   */
  std::vector<Triangle*> GetTriangles();

  /**
   * Get triangle map
   */
  std::list<Triangle*> GetMap();

  private:

  /**
   * Internals
   */

  SweepContext* sweep_context_;
  Sweep* sweep_;

};

/**
 * One independent triangulation job for TriangulateBatch().
 *
 * The points are owned by the caller and must stay alive (and unmodified)
 * until the returned BatchResult is destroyed. As with CDT, a Point object
 * must not be triangulated twice: Point::edge_list keeps raw Edge pointers
 * from the previous run. Create fresh Point objects per run.
 */
struct BatchGroup {
  /// Outer boundary, at least 3 distinct points, no repeated points.
  std::vector<Point*> polyline;
  /// Optional holes, each at least 3 distinct points.
  std::vector<std::vector<Point*>> holes;
  /// Optional Steiner points.
  std::vector<Point*> steiner_points;
};

/**
 * Outcome of one BatchGroup. A failure in one group never affects another.
 */
struct BatchGroupResult {
  bool success = false;
  /// Empty on success, otherwise a description of the failure.
  std::string error;
  /// Triangles of a successful group. Valid while the owning BatchResult lives.
  std::vector<Triangle*> triangles;

  size_t triangle_count() const { return triangles.size(); }
};

/**
 * Options for TriangulateBatch().
 */
struct BatchOptions {
  /// Number of worker threads; 0 means std::thread::hardware_concurrency().
  /// The value is clamped to the number of groups. 1 runs serially.
  size_t thread_count = 0;
};

namespace batch_detail {

inline bool ValidateContour(const std::vector<Point*>& contour, const std::string& label,
                            std::string& error)
{
  if (contour.size() < 3) {
    error = label + " has fewer than 3 points";
    return false;
  }
  std::set<std::pair<double, double>> seen;
  for (const Point* point : contour) {
    if (!point) {
      error = label + " contains a null point";
      return false;
    }
    if (!std::isfinite(point->x) || !std::isfinite(point->y)) {
      error = label + " contains a non-finite point";
      return false;
    }
    if (!seen.insert(std::make_pair(point->x, point->y)).second) {
      error = label + " contains repeated points";
      return false;
    }
  }
  return true;
}

// CDT asserts (rather than throws) when a Steiner point duplicates a point that
// is already in the context, which would abort the whole batch. Points are
// deduplicated by identity in SweepContext, so this mirrors that check.
inline bool ContourContains(const std::vector<Point*>& contour, const Point* point)
{
  for (const Point* candidate : contour) {
    if (candidate == point) {
      return true;
    }
  }
  return false;
}

inline bool ValidateGroup(const BatchGroup& group, std::string& error)
{
  if (!ValidateContour(group.polyline, "polyline", error)) {
    return false;
  }
  for (size_t i = 0; i < group.holes.size(); ++i) {
    if (!ValidateContour(group.holes[i], "hole " + std::to_string(i), error)) {
      return false;
    }
  }
  std::set<const Point*> steiner_seen;
  for (size_t i = 0; i < group.steiner_points.size(); ++i) {
    const Point* point = group.steiner_points[i];
    const std::string label = "steiner point " + std::to_string(i);
    if (!point) {
      error = label + " is null";
      return false;
    }
    if (!std::isfinite(point->x) || !std::isfinite(point->y)) {
      error = label + " is non-finite";
      return false;
    }
    if (!steiner_seen.insert(point).second || ContourContains(group.polyline, point)) {
      error = label + " duplicates an existing point";
      return false;
    }
    for (const auto& hole : group.holes) {
      if (ContourContains(hole, point)) {
        error = label + " duplicates an existing point";
        return false;
      }
    }
  }
  return true;
}

inline void TriangulateOne(const BatchGroup& group, BatchGroupResult& result,
                           std::unique_ptr<CDT>& cdt_out)
{
  std::string error;
  if (!ValidateGroup(group, error)) {
    result.success = false;
    result.error = error;
    return;
  }

  try {
    std::unique_ptr<CDT> cdt(new CDT(group.polyline));
    for (const auto& hole : group.holes) {
      cdt->AddHole(hole);
    }
    for (Point* point : group.steiner_points) {
      cdt->AddPoint(point);
    }
    cdt->Triangulate();
    result.triangles = cdt->GetTriangles();
    result.success = true;
    result.error.clear();
    cdt_out = std::move(cdt);
  } catch (const std::exception& e) {
    result.success = false;
    result.error = e.what();
  } catch (...) {
    result.success = false;
    result.error = "unknown triangulation error";
  }
}

inline void RunGroupIsolated(const BatchGroup& group, BatchGroupResult& result,
                             std::unique_ptr<CDT>& cdt)
{
  try {
    TriangulateOne(group, result, cdt);
  } catch (...) {
    // Never let one group's failure escape: it would terminate a worker thread
    // or hide the results of the groups that already finished.
    cdt.reset();
    result.triangles.clear();
    result.success = false;
    result.error = "unknown triangulation error";
  }
}

inline size_t ResolveThreadCount(size_t requested, size_t tasks)
{
  if (tasks == 0) {
    return 0;
  }
  size_t count = requested;
  if (count == 0) {
    count = std::thread::hardware_concurrency();
  }
  if (count == 0) {
    count = 1;
  }
  return count < tasks ? count : tasks;
}

} // namespace batch_detail

class BatchResult;

/**
 * Triangulate many independent point sets in one call.
 *
 * Every group is processed on its own CDT, so a group that throws or fails
 * validation is reported as a failed BatchGroupResult while the remaining
 * groups are still triangulated. Results are returned in input order. The
 * per-group triangle output is identical to calling CDT directly on the same
 * points.
 *
 * Groups with fewer than 3 points, non-finite points, repeated points within
 * one contour, or a Steiner point that duplicates any other point of the group
 * are rejected before the sweep starts (the plain CDT interface is unchanged
 * and still does not validate input). Errors thrown by the sweep itself are
 * caught per group. Note that poly2tri also uses assert() for some malformed
 * geometry; an assertion fires inside the triangulating thread and cannot be
 * caught in-process, so wholly degenerate input should still be validated by
 * the caller if assertions are enabled.
 *
 * Triangulation is run on up to options.thread_count threads. When
 * thread_count is 1 everything runs on the calling thread.
 */
inline BatchResult TriangulateBatch(const std::vector<BatchGroup>& groups,
                                    const BatchOptions& options = BatchOptions());

class BatchResult
{
public:

  BatchResult() : elapsed_seconds_(0.0) {}

  /// Number of submitted groups.
  size_t size() const { return groups_.size(); }

  bool empty() const { return groups_.empty(); }

  /// Result for one group; throws std::out_of_range if index >= size().
  const BatchGroupResult& at(size_t index) const { return groups_.at(index); }

  const BatchGroupResult& operator[](size_t index) const { return groups_[index]; }

  const std::vector<BatchGroupResult>& groups() const { return groups_; }

  /// Number of groups that triangulated successfully.
  size_t succeeded() const
  {
    size_t count = 0;
    for (const auto& group : groups_) {
      if (group.success) {
        ++count;
      }
    }
    return count;
  }

  /// Number of groups that failed validation or triangulation.
  size_t failed() const { return groups_.size() - succeeded(); }

  /// Sum of triangle_count() over all groups.
  size_t total_triangles() const
  {
    size_t count = 0;
    for (const auto& group : groups_) {
      count += group.triangle_count();
    }
    return count;
  }

  /// Wall-clock duration of the whole TriangulateBatch() call, in seconds.
  double elapsed_seconds() const { return elapsed_seconds_; }

private:

  friend BatchResult TriangulateBatch(const std::vector<BatchGroup>& groups,
                                      const BatchOptions& options);

  std::vector<BatchGroupResult> groups_;
  // Keeps the Triangle objects referenced by groups_ alive.
  std::vector<std::unique_ptr<CDT>> cdts_;
  double elapsed_seconds_;
};

inline BatchResult TriangulateBatch(const std::vector<BatchGroup>& groups,
                                    const BatchOptions& options)
{
  BatchResult result;
  result.groups_.resize(groups.size());
  result.cdts_.resize(groups.size());

  const auto start = std::chrono::steady_clock::now();
  const size_t thread_count = batch_detail::ResolveThreadCount(options.thread_count, groups.size());

  if (thread_count <= 1) {
    for (size_t i = 0; i < groups.size(); ++i) {
      batch_detail::RunGroupIsolated(groups[i], result.groups_[i], result.cdts_[i]);
    }
  } else {
    std::atomic<size_t> next(0);
    auto work = [&groups, &result, &next]() {
      for (;;) {
        const size_t index = next.fetch_add(1, std::memory_order_relaxed);
        if (index >= groups.size()) {
          break;
        }
        batch_detail::RunGroupIsolated(groups[index], result.groups_[index], result.cdts_[index]);
      }
    };
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    try {
      for (size_t t = 0; t < thread_count; ++t) {
        workers.emplace_back(std::ref(work));
      }
    } catch (...) {
      // Thread creation failed (e.g. resource limits). Join what started, then
      // finish the remaining groups on this thread so every group still gets a
      // result instead of the whole batch failing.
    }
    for (auto& worker : workers) {
      if (worker.joinable()) {
        worker.join();
      }
    }
    work(); // no-op when the workers already consumed every group
  }
  const auto end = std::chrono::steady_clock::now();
  result.elapsed_seconds_ = std::chrono::duration<double>(end - start).count();
  return result;
}

}

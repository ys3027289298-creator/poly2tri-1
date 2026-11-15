#pragma once

#include "../common/dll_symbol.h"
#include "cdt.h"

#include <cstddef>
#include <vector>

namespace p2t {

/// Dual refinement constraint: maximum edge length and minimum triangle area.
struct P2T_DLL_SYMBOL RefinementCriteria {
  /// Triangles whose longest edge exceeds this length are refined.
  /// A value <= 0 disables the max-edge-length constraint.
  double max_edge_length;
  /// Triangles whose area falls below this value are refined.
  /// A value <= 0 disables the min-area constraint.
  double min_area;
  /// Maximum number of point-insertion/re-triangulation rounds.
  int max_iterations;
  /// When false, no refinement is performed: the output is vertex-by-vertex
  /// identical to a plain p2t::CDT triangulation.
  bool enabled;

  RefinementCriteria();
};

/// Statistics of a Refiner::Triangulate run.
struct P2T_DLL_SYMBOL RefinementReport {
  /// Number of point-insertion/re-triangulation rounds actually performed.
  int iterations = 0;
  /// Number of triangles in the final mesh.
  std::size_t triangle_count = 0;
  /// Longest triangle edge in the final mesh.
  double max_edge_length = 0.0;
  /// True when both constraints are satisfied; false when the iteration
  /// cap was reached (or no further points could be inserted) first.
  bool converged = true;
};

/**
 * Iteratively refining constrained Delaunay triangulator.
 *
 * Mirrors the p2t::CDT usage pattern (polyline, holes, Steiner points),
 * then Triangulate() checks every triangle against the dual
 * max-edge-length/min-area constraint. Any violated constraint inserts a
 * point at the offending triangle (midpoint of the longest edge, splitting
 * the polyline segment when that edge is constrained, otherwise the
 * triangle centroid) and the domain is re-triangulated so all constrained
 * edges are re-imposed. The loop runs until convergence or the iteration
 * cap. After every round the mesh is verified: a constrained edge crossed
 * by any triangle edge is reported via std::runtime_error instead of
 * being silently accepted.
 *
 * Points passed to the constructor/AddHole/AddPoint stay owned by the
 * caller (same convention as p2t::CDT); points inserted during refinement
 * are owned by the Refiner and deleted with it.
 */
class P2T_DLL_SYMBOL Refiner {
public:
  /// Constructor - add polyline with non repeating points
  explicit Refiner(const std::vector<Point*>& polyline);
  /// Destructor - clean up memory
  ~Refiner();

  Refiner(const Refiner&) = delete;
  Refiner& operator=(const Refiner&) = delete;

  /// Add a hole
  void AddHole(const std::vector<Point*>& polyline);
  /// Add a steiner point
  void AddPoint(Point* point);

  /// Triangulate - do this AFTER you've added the polyline, holes, and
  /// Steiner points. Iteratively refines according to the criteria.
  const RefinementReport& Triangulate(const RefinementCriteria& criteria);

  /// Get triangles of the final mesh
  std::vector<Triangle*> GetTriangles();

  /// Report of the last Triangulate call
  const RefinementReport& report() const;

private:
  /// Midpoint insertion splitting a constrained polyline segment a-b.
  struct EdgeSplit {
    Point point;
    const Point* a;
    const Point* b;
  };

  void RebuildCdt();
  bool FindViolations(const std::vector<Triangle*>& triangles,
                      const RefinementCriteria& criteria,
                      std::vector<EdgeSplit>& edge_splits,
                      std::vector<Point>& interior_points) const;
  std::size_t InsertPoints(const std::vector<EdgeSplit>& edge_splits,
                           const std::vector<Point>& interior_points);
  bool SplitSegment(const Point* a, const Point* b, Point* point);
  std::vector<Point*> AllPoints() const;
  std::vector<std::vector<Point*>> AllPolylines() const;

  std::vector<Point*> boundary_;
  std::vector<std::vector<Point*>> holes_;
  std::vector<Point*> steiner_points_;
  /// Interior points created by the refinement itself (owned).
  std::vector<Point*> inserted_points_;
  /// Points inserted into polylines by segment splits (owned).
  std::vector<Point*> split_points_;
  CDT* cdt_ = nullptr;
  RefinementReport report_;
};

/**
 * Verifies that no triangle edge properly crosses any of the constrained
 * edges defined by the given closed polylines (boundary and holes).
 * Throws std::runtime_error describing the crossing when one is found.
 */
void P2T_DLL_SYMBOL
CheckConstraintsIntact(const std::vector<Triangle*>& triangles,
                       const std::vector<std::vector<Point*>>& polylines);

} // namespace p2t

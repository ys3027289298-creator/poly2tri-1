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

#include "dll_symbol.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <stdexcept>
#include <utility>
#include <vector>

namespace p2t {

struct Edge;

struct P2T_DLL_SYMBOL Point {

  double x, y;

  /// Default constructor does nothing (for performance).
  Point()
  {
    x = 0.0;
    y = 0.0;
  }

  /// The edges this point constitutes an upper ending point
  std::vector<Edge*> edge_list;

  /// Construct using coordinates.
  Point(double x, double y);

  /// Set this point to all zeros.
  void set_zero()
  {
    x = 0.0;
    y = 0.0;
  }

  /// Set this point to some specified coordinates.
  void set(double x_, double y_)
  {
    x = x_;
    y = y_;
  }

  /// Negate this point.
  Point operator -() const
  {
    Point v;
    v.set(-x, -y);
    return v;
  }

  /// Add a point to this point.
  void operator +=(const Point& v)
  {
    x += v.x;
    y += v.y;
  }

  /// Subtract a point from this point.
  void operator -=(const Point& v)
  {
    x -= v.x;
    y -= v.y;
  }

  /// Multiply this point by a scalar.
  void operator *=(double a)
  {
    x *= a;
    y *= a;
  }

  /// Get the length of this point (the norm).
  double Length() const
  {
    return sqrt(x * x + y * y);
  }

  /// Convert this point into a unit point. Returns the Length.
  double Normalize()
  {
    const double len = Length();
    x /= len;
    y /= len;
    return len;
  }

};

P2T_DLL_SYMBOL std::ostream& operator<<(std::ostream&, const Point&);

// Represents a simple polygon's edge
struct P2T_DLL_SYMBOL Edge {

  Point* p, *q;

  /// Constructor
  Edge(Point& p1, Point& p2) : p(&p1), q(&p2)
  {
    if (p1.y > p2.y) {
      q = &p1;
      p = &p2;
    } else if (p1.y == p2.y) {
      if (p1.x > p2.x) {
        q = &p1;
        p = &p2;
      } else if (p1.x == p2.x) {
        // Repeat points
        throw std::runtime_error("Edge::Edge: p1 == p2");
      }
    }

    q->edge_list.push_back(this);
  }
};

// Triangle-based data structures are know to have better performance than quad-edge structures
// See: J. Shewchuk, "Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator"
//      "Triangulations in CGAL"
class P2T_DLL_SYMBOL Triangle {
public:

/// Constructor
Triangle(Point& a, Point& b, Point& c);

/// Flags to determine if an edge is a Constrained edge
bool constrained_edge[3];
/// Flags to determine if an edge is a Delauney edge
bool delaunay_edge[3];

Point* GetPoint(int index);
Point* PointCW(const Point& point);
Point* PointCCW(const Point& point);
Point* OppositePoint(Triangle& t, const Point& p);

Triangle* GetNeighbor(int index);
void MarkNeighbor(Point* p1, Point* p2, Triangle* t);
void MarkNeighbor(Triangle& t);

void MarkConstrainedEdge(int index);
void MarkConstrainedEdge(Edge& edge);
void MarkConstrainedEdge(Point* p, Point* q);

int Index(const Point* p);
int EdgeIndex(const Point* p1, const Point* p2);

Triangle* NeighborAcross(const Point& point);
Triangle* NeighborCW(const Point& point);
Triangle* NeighborCCW(const Point& point);
bool GetConstrainedEdgeCCW(const Point& p);
bool GetConstrainedEdgeCW(const Point& p);
void SetConstrainedEdgeCCW(const Point& p, bool ce);
void SetConstrainedEdgeCW(const Point& p, bool ce);
bool GetDelunayEdgeCCW(const Point& p);
bool GetDelunayEdgeCW(const Point& p);
void SetDelunayEdgeCCW(const Point& p, bool e);
void SetDelunayEdgeCW(const Point& p, bool e);

bool Contains(const Point* p);
bool Contains(const Edge& e);
bool Contains(const Point* p, const Point* q);
void Legalize(Point& point);
void Legalize(Point& opoint, Point& npoint);
/**
 * Clears all references to all other triangles and points
 */
void Clear();
void ClearNeighbor(const Triangle *triangle);
void ClearNeighbors();
void ClearDelunayEdges();

inline bool IsInterior();
inline void IsInterior(bool b);

void DebugPrint();

bool CircumcicleContains(const Point&) const;

private:

bool IsCounterClockwise() const;

/// Triangle points
Point* points_[3];
/// Neighbor list
Triangle* neighbors_[3];

/// Has this triangle been marked as an interior triangle?
bool interior_;
};

inline bool cmp(const Point* a, const Point* b)
{
  if (a->y < b->y) {
    return true;
  } else if (a->y == b->y) {
    // Make sure q is point with greater x value
    if (a->x < b->x) {
      return true;
    }
  }
  return false;
}

/// Add two points_ component-wise.
inline Point operator +(const Point& a, const Point& b)
{
  return Point(a.x + b.x, a.y + b.y);
}

/// Subtract two points_ component-wise.
inline Point operator -(const Point& a, const Point& b)
{
  return Point(a.x - b.x, a.y - b.y);
}

/// Multiply point by scalar
inline Point operator *(double s, const Point& a)
{
  return Point(s * a.x, s * a.y);
}

inline bool operator ==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}

inline bool operator !=(const Point& a, const Point& b)
{
  return !(a.x == b.x) || !(a.y == b.y);
}

/// Peform the dot product on two vectors.
inline double Dot(const Point& a, const Point& b)
{
  return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline double Cross(const Point& a, const Point& b)
{
  return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a point and a scalar. In 2D this produces
/// a point.
inline Point Cross(const Point& a, double s)
{
  return Point(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a point. In 2D this produces
/// a point.
inline Point Cross(double s, const Point& a)
{
  return Point(-s * a.y, s * a.x);
}

inline Point* Triangle::GetPoint(int index)
{
  return points_[index];
}

inline Triangle* Triangle::GetNeighbor(int index)
{
  return neighbors_[index];
}

inline bool Triangle::Contains(const Point* p)
{
  return p == points_[0] || p == points_[1] || p == points_[2];
}

inline bool Triangle::Contains(const Edge& e)
{
  return Contains(e.p) && Contains(e.q);
}

inline bool Triangle::Contains(const Point* p, const Point* q)
{
  return Contains(p) && Contains(q);
}

inline bool Triangle::IsInterior()
{
  return interior_;
}

inline void Triangle::IsInterior(bool b)
{
  interior_ = b;
}

/// Is this set a valid delaunay triangulation?
P2T_DLL_SYMBOL bool IsDelaunay(const std::vector<p2t::Triangle*>&);

/// Options controlling polyline topology repair.
struct RepairOptions {
  /// Points closer together than this distance are merged into one
  /// (the first point of a run is kept). A value of 0 disables
  /// distance-based merging; exact adjacent duplicates and the
  /// resulting zero-length edges are always removed.
  double merge_threshold = 0.0;
  /// When true, RepairReport::original_indices is filled with the
  /// input index of every surviving point.
  bool keep_original_indices = false;
};

/// Statistics of a RepairPolyline call.
struct RepairReport {
  std::size_t points_before = 0; ///< Points in the input polyline.
  std::size_t points_after = 0;  ///< Points in the repaired polyline.
  double area_before = 0.0;      ///< Signed area before repair.
  double area_after = 0.0;       ///< Signed area after repair (same sign as before).
  /// Maps each repaired point to its index in the input polyline.
  /// Only filled when RepairOptions::keep_original_indices is true.
  std::vector<std::size_t> original_indices;
};

/// Signed (shoelace) area of a closed polyline. Positive means
/// counter-clockwise.
inline double SignedArea(const std::vector<Point*>& polyline)
{
  double area = 0.0;
  const std::size_t n = polyline.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Point& a = *polyline[i];
    const Point& b = *polyline[(i + 1) % n];
    area += a.x * b.y - b.x * a.y;
  }
  return 0.5 * area;
}

/// Repair the topology of an imported polygon polyline in place:
///  * adjacent duplicate points are removed (including wrap-around),
///  * dangling edges (spikes of the form p, q, p) are removed,
///  * points closer than options.merge_threshold are merged,
///  * zero-length edges are thereby eliminated,
///  * the orientation (signed-area sign) of the contour is preserved.
/// At least three points are always kept. The Point objects themselves
/// are not modified; surviving points keep their identity, so the
/// triangulation output of an already-clean polyline is unchanged.
inline RepairReport RepairPolyline(std::vector<Point*>& polyline,
                                   const RepairOptions& options = RepairOptions())
{
  RepairReport report;
  report.points_before = polyline.size();
  report.area_before = SignedArea(polyline);

  const std::size_t n = polyline.size();
  if (n == 0) {
    return report;
  }

  typedef std::pair<Point*, std::size_t> IndexedPoint;

  // Pair every point with its original index.
  std::vector<IndexedPoint> pts;
  pts.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    pts.push_back(IndexedPoint(polyline[i], i));
  }

  // 1. Drop adjacent duplicate points (zero-length edges), wrap-around included.
  std::vector<IndexedPoint> dedup;
  dedup.reserve(n);
  for (const auto& e : pts) {
    if (dedup.empty() || *e.first != *dedup.back().first) {
      dedup.push_back(e);
    }
  }
  while (dedup.size() > 1 && *dedup.front().first == *dedup.back().first) {
    dedup.pop_back();
  }

  // 2. Remove dangling edges: spikes p, q, p (edge traversed out and back).
  std::vector<IndexedPoint> despiked;
  despiked.reserve(dedup.size());
  for (const auto& e : dedup) {
    if (despiked.size() >= 2 && *e.first == *despiked[despiked.size() - 2].first) {
      despiked.pop_back(); // drop the spike tip q and the returning p
    } else {
      despiked.push_back(e);
    }
  }

  // 3. Merge points closer than the threshold (greedy), always keeping
  //    at least three points.
  std::vector<IndexedPoint> merged;
  merged.reserve(despiked.size());
  const double tol = options.merge_threshold;
  const double tol2 = tol * tol;
  const std::size_t m = despiked.size();
  for (std::size_t i = 0; i < m; ++i) {
    const bool must_keep = merged.size() + (m - i) <= 3;
    if (!must_keep && !merged.empty()) {
      const Point& a = *despiked[i].first;
      const Point& b = *merged.back().first;
      const double dx = a.x - b.x, dy = a.y - b.y;
      if (dx * dx + dy * dy < tol2) {
        continue; // merge into the previously kept point
      }
    }
    merged.push_back(despiked[i]);
  }
  // Wrap-around: first and last kept points too close -> drop the last.
  while (merged.size() > 3) {
    const Point& a = *merged.front().first;
    const Point& b = *merged.back().first;
    const double dx = a.x - b.x, dy = a.y - b.y;
    if (dx * dx + dy * dy >= tol2) {
      break;
    }
    merged.pop_back();
  }

  // 4. Preserve orientation: if the repair flipped the signed-area sign,
  //    reverse the repaired polyline.
  double area = 0.0;
  for (std::size_t i = 0; i < merged.size(); ++i) {
    const Point& a = *merged[i].first;
    const Point& b = *merged[(i + 1) % merged.size()].first;
    area += a.x * b.y - b.x * a.y;
  }
  area *= 0.5;
  if ((report.area_before > 0.0 && area < 0.0) ||
      (report.area_before < 0.0 && area > 0.0)) {
    std::reverse(merged.begin(), merged.end());
  }

  // Write the repaired polyline back and fill the report.
  polyline.clear();
  polyline.reserve(merged.size());
  for (const auto& e : merged) {
    polyline.push_back(e.first);
  }
  report.points_after = polyline.size();
  report.area_after = SignedArea(polyline);
  if (options.keep_original_indices) {
    report.original_indices.reserve(merged.size());
    for (const auto& e : merged) {
      report.original_indices.push_back(e.second);
    }
  }
  return report;
}

}

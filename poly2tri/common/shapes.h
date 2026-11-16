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
#include <map>
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

// ---------------------------------------------------------------------------
// Polygon offsetting (contour expansion / shrinking before triangulation)
//
// OffsetPolygon() shifts every edge of a simple closed polygon along its
// outward normal by a configurable distance, joins adjacent offset edges
// (miter or bevel, configurable), clips the result at self-intersections and
// returns the cleaned loops, ready to be triangulated with p2t::CDT.
// ---------------------------------------------------------------------------

/// Join style used where two adjacent offset edges meet.
enum class OffsetJoin {
  Miter,  ///< Extend both offset edges to their intersection (clamped by the miter limit).
  Bevel   ///< Connect the two offset edge endpoints with a straight segment.
};

/// Options controlling OffsetPolygon().
struct OffsetOptions {
  /// Offset distance along the edge normals. Positive expands a CCW polygon
  /// outward, negative shrinks it inward, zero returns the input unchanged.
  double distance = 0.0;
  /// How adjacent offset edges are connected at vertices.
  OffsetJoin join = OffsetJoin::Miter;
  /// Maximum miter length as a multiple of |distance| (must be >= 1).
  /// Spikes longer than this fall back to a bevel join.
  double miter_limit = 2.0;
};

/// Signed area of a closed loop (positive for counter-clockwise winding).
inline double SignedArea(const std::vector<Point>& loop)
{
  double area = 0.0;
  for (std::size_t i = 0; i < loop.size(); ++i) {
    const Point& a = loop[i];
    const Point& b = loop[(i + 1) % loop.size()];
    area += a.x * b.y - b.x * a.y;
  }
  return 0.5 * area;
}

/// True if the closed loop intersects or touches itself anywhere
/// (non-adjacent edges crossing, touching or overlapping).
inline bool HasSelfIntersections(const std::vector<Point>& loop);

namespace detail {

inline bool OffsetFinite(double v)
{
  return std::isfinite(v);
}

/// Counter-clockwise test of the three points.
inline double OffsetOrient(const Point& a, const Point& b, const Point& c)
{
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool OffsetOnSegment(const Point& a, const Point& b, const Point& p, double eps)
{
  return std::fabs(OffsetOrient(a, b, p)) <= eps &&
         p.x >= std::min(a.x, b.x) - eps && p.x <= std::max(a.x, b.x) + eps &&
         p.y >= std::min(a.y, b.y) - eps && p.y <= std::max(a.y, b.y) + eps;
}

/// Intersection of segments ab and cd. Returns true if the segments cross or
/// touch (including collinear overlap). For a proper crossing the point is
/// written to `out` and `proper` is set to true.
inline bool OffsetSegmentIntersect(const Point& a, const Point& b,
                                   const Point& c, const Point& d,
                                   Point* out, bool* proper)
{
  const double d1 = OffsetOrient(c, d, a);
  const double d2 = OffsetOrient(c, d, b);
  const double d3 = OffsetOrient(a, b, c);
  const double d4 = OffsetOrient(a, b, d);
  const double eps = 1e-12;
  if (((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
      ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps))) {
    // Proper crossing: solve a + t*(b-a) on line cd.
    const double denom = d1 - d2;
    const double t = d1 / denom;
    if (out) {
      out->set(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
    }
    if (proper) {
      *proper = true;
    }
    return true;
  }
  if (proper) {
    *proper = false;
  }
  // Touching or collinear cases.
  if (std::fabs(d1) <= eps && OffsetOnSegment(c, d, a, eps)) {
    if (out) *out = a;
    return true;
  }
  if (std::fabs(d2) <= eps && OffsetOnSegment(c, d, b, eps)) {
    if (out) *out = b;
    return true;
  }
  if (std::fabs(d3) <= eps && OffsetOnSegment(a, b, c, eps)) {
    if (out) *out = c;
    return true;
  }
  if (std::fabs(d4) <= eps && OffsetOnSegment(a, b, d, eps)) {
    if (out) *out = d;
    return true;
  }
  return false;
}

/// Winding number of `p` with respect to the closed `loop`.
inline int OffsetWindingNumber(const std::vector<Point>& loop, const Point& p)
{
  int winding = 0;
  for (std::size_t i = 0; i < loop.size(); ++i) {
    const Point& a = loop[i];
    const Point& b = loop[(i + 1) % loop.size()];
    if (a.y <= p.y) {
      if (b.y > p.y && OffsetOrient(a, b, p) > 0) {
        ++winding;
      }
    } else {
      if (b.y <= p.y && OffsetOrient(a, b, p) < 0) {
        --winding;
      }
    }
  }
  return winding;
}

struct OffsetSubSegment {
  Point u, v;  ///< Directed from u to v, following the raw polygon direction.
};

/// Split the raw offset loop at every self-intersection into directed
/// sub-segments whose endpoints all coincide exactly at crossings.
inline std::vector<OffsetSubSegment> OffsetSplitAtIntersections(
    const std::vector<Point>& loop)
{
  const std::size_t n = loop.size();
  // Split points (parameter in [0,1] plus the exact point) per edge; the
  // endpoints are always present. Crossing points are stored as the exact
  // same value on both edges so that node merging can compare exactly.
  std::vector<std::vector<std::pair<double, Point>>> splits(n);
  for (std::size_t i = 0; i < n; ++i) {
    splits[i].push_back(std::make_pair(0.0, loop[i]));
    splits[i].push_back(std::make_pair(1.0, loop[(i + 1) % n]));
  }
  auto parameter = [&](std::size_t i, const Point& p) {
    const Point& a = loop[i];
    const Point& b = loop[(i + 1) % n];
    const double len = (b - a).Length();
    if (len == 0.0) {
      return 0.0;
    }
    return (p - a).Length() / len * (Dot(p - a, b - a) >= 0 ? 1.0 : -1.0);
  };
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      // Skip adjacent edges (they only share a vertex).
      if (j == i + 1 || (i == 0 && j == n - 1)) {
        continue;
      }
      const Point& a = loop[i];
      const Point& b = loop[(i + 1) % n];
      const Point& c = loop[j];
      const Point& d = loop[(j + 1) % n];
      Point x;
      bool proper = false;
      if (!OffsetSegmentIntersect(a, b, c, d, &x, &proper)) {
        continue;
      }
      if (proper) {
        splits[i].push_back(std::make_pair(parameter(i, x), x));
        splits[j].push_back(std::make_pair(parameter(j, x), x));
      } else {
        // Touching / collinear: add the touching endpoint(s) to the other edge.
        const double eps = 1e-12;
        auto add_if_on = [&](std::size_t edge, const Point& p) {
          const Point& p0 = loop[edge];
          const Point& p1 = loop[(edge + 1) % n];
          if (OffsetOnSegment(p0, p1, p, eps)) {
            splits[edge].push_back(std::make_pair(parameter(edge, p), p));
          }
        };
        add_if_on(i, c);
        add_if_on(i, d);
        add_if_on(j, a);
        add_if_on(j, b);
      }
    }
  }
  std::vector<OffsetSubSegment> result;
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<std::pair<double, Point>>& ts = splits[i];
    std::stable_sort(ts.begin(), ts.end(),
                     [](const std::pair<double, Point>& a,
                        const std::pair<double, Point>& b) {
                       return a.first < b.first;
                     });
    ts.erase(std::unique(ts.begin(), ts.end(),
                         [](const std::pair<double, Point>& a,
                            const std::pair<double, Point>& b) {
                           return std::fabs(a.first - b.first) < 1e-9;
                         }),
             ts.end());
    for (std::size_t k = 0; k + 1 < ts.size(); ++k) {
      const double t0 = ts[k].first;
      const double t1 = ts[k + 1].first;
      if (t1 - t0 < 1e-12) {
        continue;
      }
      OffsetSubSegment seg;
      seg.u = ts[k].second;
      seg.v = ts[k + 1].second;
      result.push_back(seg);
    }
  }
  return result;
}

/// Keep only the sub-segments that border the offset region (positive winding
/// number on the left, non-positive on the right), oriented so the region
/// stays on the left, and stitch them into closed loops. Loops with positive
/// signed area are outer boundaries, loops with negative area are holes.
inline std::vector<std::vector<Point>> OffsetExtractBoundary(
    const std::vector<Point>& raw)
{
  const std::vector<OffsetSubSegment> subs = OffsetSplitAtIntersections(raw);
  if (subs.empty()) {
    return {};
  }
  // Characteristic length used for the winding sample offset.
  double min_len = 0.0;
  for (const auto& s : subs) {
    const double len = (s.v - s.u).Length();
    if (len > 0 && (min_len == 0.0 || len < min_len)) {
      min_len = len;
    }
  }
  const double eps = min_len * 1e-6;

  // Kept directed edges, region on the left.
  std::vector<OffsetSubSegment> kept;
  for (const auto& s : subs) {
    const Point mid((s.u.x + s.v.x) * 0.5, (s.u.y + s.v.y) * 0.5);
    const Point dir = s.v - s.u;
    const double len = dir.Length();
    if (len == 0.0) {
      continue;
    }
    const Point left_normal(-dir.y / len, dir.x / len);
    const Point pl(mid.x + eps * left_normal.x, mid.y + eps * left_normal.y);
    const Point pr(mid.x - eps * left_normal.x, mid.y - eps * left_normal.y);
    const int wl = OffsetWindingNumber(raw, pl);
    const int wr = OffsetWindingNumber(raw, pr);
    if (wl > 0 && wr <= 0) {
      kept.push_back(s);
    } else if (wr > 0 && wl <= 0) {
      kept.push_back(OffsetSubSegment{s.v, s.u});
    }
  }
  if (kept.empty()) {
    return {};
  }

  // Merge identical endpoints into nodes (exact arithmetic is shared between
  // the two sides of every split, so exact comparison is sufficient).
  std::vector<Point> nodes;
  auto node_index = [&](const Point& p) {
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      if (nodes[i] == p) {
        return i;
      }
    }
    nodes.push_back(p);
    return nodes.size() - 1;
  };
  std::vector<std::pair<int, int>> edges;  // (from, to) node indices
  std::vector<std::vector<int>> outgoing;
  for (const auto& s : kept) {
    const int u = node_index(s.u);
    const int v = node_index(s.v);
    if (u == v) {
      continue;
    }
    if (static_cast<std::size_t>(std::max(u, v) + 1) > outgoing.size()) {
      outgoing.resize(std::max(u, v) + 1);
    }
    // Skip duplicate directed edges (collinear overlaps produce them).
    bool duplicate = false;
    for (const int e : outgoing[u]) {
      if (edges[e].second == v) {
        duplicate = true;
        break;
      }
    }
    if (duplicate) {
      continue;
    }
    outgoing[u].push_back(static_cast<int>(edges.size()));
    edges.push_back(std::make_pair(u, v));
  }

  // Stitch edges into loops: at every node continue with the unused outgoing
  // edge that makes the sharpest right turn (smallest clockwise angle)
  // relative to the incoming direction.
  const double kTwoPi = 6.28318530717958647692;
  std::vector<char> used(edges.size(), 0);
  std::vector<std::vector<Point>> loops;
  for (std::size_t start = 0; start < edges.size(); ++start) {
    if (used[start]) {
      continue;
    }
    std::vector<Point> loop;
    std::size_t current = start;
    while (!used[current]) {
      used[current] = 1;
      const int u = edges[current].first;
      const int v = edges[current].second;
      loop.push_back(nodes[u]);
      const Point in_dir = nodes[v] - nodes[u];
      std::size_t best = edges.size();
      double best_angle = 0.0;
      for (const int e : outgoing[v]) {
        if (used[e]) {
          continue;
        }
        const Point out_dir = nodes[edges[e].second] - nodes[v];
        // Clockwise angle from the incoming direction, in [0, 2*pi).
        double angle = std::atan2(-Cross(in_dir, out_dir), Dot(in_dir, out_dir));
        if (angle < 0) {
          angle += kTwoPi;
        }
        if (best == edges.size() || angle < best_angle) {
          best = e;
          best_angle = angle;
        }
      }
      if (best == edges.size()) {
        break;  // Dead end (should not happen for a valid arrangement).
      }
      current = best;
    }
    if (loop.size() >= 3) {
      loops.push_back(loop);
    }
  }
  return loops;
}

}  // namespace detail

inline bool HasSelfIntersections(const std::vector<Point>& loop)
{
  const std::size_t n = loop.size();
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (j == i + 1 || (i == 0 && j == n - 1)) {
        continue;
      }
      if (detail::OffsetSegmentIntersect(loop[i], loop[(i + 1) % n],
                                         loop[j], loop[(j + 1) % n],
                                         nullptr, nullptr)) {
        return true;
      }
    }
  }
  return false;
}

/// Offset a simple closed polygon along its edge normals.
///
/// Every edge is translated along its outward normal by `options.distance`,
/// adjacent offset edges are joined (miter or bevel, see OffsetOptions) and
/// self-intersections of the resulting raw curve are clipped away. The
/// returned loops can be triangulated with p2t::CDT: loops with positive
/// SignedArea() are outer boundaries, loops with negative area are holes.
///
/// With `options.distance == 0` the input vertices are returned unchanged.
///
/// @param polyline  Closed polygon, non-repeating points, any winding.
/// @param options   Offset distance, join style and miter limit.
/// @throws std::invalid_argument on an illegal polygon (fewer than 3 points,
///         zero-length edges) or an illegal distance / miter limit (NaN,
///         infinity, miter limit < 1).
inline std::vector<std::vector<Point>> OffsetPolygon(
    const std::vector<Point*>& polyline,
    const OffsetOptions& options = OffsetOptions())
{
  if (polyline.size() < 3) {
    throw std::invalid_argument("OffsetPolygon: need at least 3 points");
  }
  if (!detail::OffsetFinite(options.distance)) {
    throw std::invalid_argument("OffsetPolygon: distance must be finite");
  }
  if (!detail::OffsetFinite(options.miter_limit) || options.miter_limit < 1.0) {
    throw std::invalid_argument("OffsetPolygon: miter_limit must be finite and >= 1");
  }

  std::vector<Point> pts;
  pts.reserve(polyline.size());
  for (const Point* p : polyline) {
    if (p == nullptr) {
      throw std::invalid_argument("OffsetPolygon: null point");
    }
    pts.push_back(*p);
  }
  const std::size_t n = pts.size();
  for (std::size_t i = 0; i < n; ++i) {
    if (pts[i] == pts[(i + 1) % n]) {
      throw std::invalid_argument("OffsetPolygon: zero-length edge");
    }
  }

  // Zero distance: return the input vertex-for-vertex unchanged.
  if (options.distance == 0.0) {
    return std::vector<std::vector<Point>>{pts};
  }

  // Normalize to counter-clockwise so "outward" is well defined.
  if (SignedArea(pts) < 0) {
    std::reverse(pts.begin(), pts.end());
  }

  const double d = options.distance;
  // Outward (right-hand) unit normal per edge.
  std::vector<Point> normals(n);
  for (std::size_t i = 0; i < n; ++i) {
    Point e = pts[(i + 1) % n] - pts[i];
    const double len = e.Normalize();
    (void)len;
    normals[i] = Point(e.y, -e.x);
  }

  // Raw offset loop: join the offset edges at every vertex.
  std::vector<Point> raw;
  raw.reserve(2 * n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t prev = (i + n - 1) % n;
    const Point& v = pts[i];
    // Offset endpoints of the incoming and outgoing edges at this vertex.
    const Point a2 = v + d * normals[prev];  // end of previous offset edge
    const Point b1 = v + d * normals[i];     // start of next offset edge
    bool emit_miter = false;
    Point miter;
    if (options.join == OffsetJoin::Miter) {
      const Point dir_prev = pts[i] - pts[prev];
      const Point dir_next = pts[(i + 1) % n] - pts[i];
      const double denom = Cross(dir_prev, dir_next);
      if (std::fabs(denom) > 1e-30) {
        // Intersection of the two offset edge lines.
        const Point delta = b1 - a2;
        const double t = Cross(delta, dir_next) / denom;
        miter = a2 + t * dir_prev;
        const double miter_len = (miter - v).Length();
        emit_miter = miter_len <= options.miter_limit * std::fabs(d) * (1.0 + 1e-9);
      }
    }
    if (emit_miter) {
      raw.push_back(miter);
    } else {
      raw.push_back(a2);
      if (b1 != a2) {
        raw.push_back(b1);
      }
    }
  }

  // Clip self-intersections and stitch the region boundary.
  std::vector<std::vector<Point>> loops = detail::OffsetExtractBoundary(raw);
  if (loops.empty()) {
    return loops;  // Polygon collapsed (e.g. inward offset too large).
  }
  // Drop degenerate loops and order: outer boundaries first, then holes.
  std::vector<std::vector<Point>> outers, holes;
  for (auto& loop : loops) {
    const double area = SignedArea(loop);
    if (area > 0) {
      outers.push_back(std::move(loop));
    } else if (area < 0) {
      holes.push_back(std::move(loop));
    }
  }
  outers.insert(outers.end(), holes.begin(), holes.end());
  return outers;
}

}

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
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iosfwd>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
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

/*
 * Convex decomposition.
 *
 * Collision detection that treats a polygon as a single convex shape is wrong
 * as soon as its boundary has a reflex (concave) vertex.  DecomposePolygon()
 * splits a simple polygon, optionally carrying holes, into convex pieces in
 * three stages:
 *
 *   1. Ear clipping walks the boundary and clips ears.  A reflex vertex can
 *      never be an ear apex, so reflex vertices are identified here and stay
 *      in the ring until the ears around them have been trimmed away.
 *   2. Bridging re-joins the clipped ears across their shared diagonals.  A
 *      diagonal is only kept when the union stays convex, so the result is a
 *      small set of grown convex blocks rather than a full triangulation.
 *   3. The piece set is validated: every piece must be convex and the piece
 *      areas must add up to the input polygon's area.
 *
 * Anything that prevents a stage from finishing (self-intersecting or
 * degenerate input, no legal hole bridge, a piece that cannot honour
 * min_piece_area) is reported by throwing std::runtime_error.  This is an
 * independent analysis pass: the CDT triangulation path is untouched.
 */

struct ConvexPiece {
  /// Points of the piece in CCW order.  They alias the caller's points.
  std::vector<Point*> points;

  double Area() const
  {
    double a = 0.0;
    const size_t n = points.size();
    for (size_t i = 0; i < n; ++i) {
      const Point& p = *points[i];
      const Point& q = *points[(i + 1) % n];
      a += p.x * q.y - q.x * p.y;
    }
    return a * 0.5;
  }
};

/// How bridge diagonals are chosen (hole bridges and convex re-joinings).
enum class BridgeStrategy {
  /// Shortest valid bridge (Euclidean length).
  Shortest,
  /// Bridge whose endpoints are closest in the input vertex ordering.
  NearestIndex,
};

struct ConvexDecompositionOptions {
  BridgeStrategy bridge = BridgeStrategy::Shortest;

  /// Convex pieces below this area (same units as the input) may not remain:
  /// DecomposePolygon() merges them into a convex neighbour, and throws if no
  /// neighbour can absorb them.  Zero disables the check.
  double min_piece_area = 0.0;

  /// Relative tolerance of the stage-3 area check:
  /// |sum(pieces) - polygon| <= area_epsilon * max(1, |polygon area|).
  double area_epsilon = 1e-9;

  /// Pieces whose |area| is at most this are degenerate and are never emitted.
  double degenerate_area = 1e-12;
};

struct ConvexDecompositionResult {
  std::vector<ConvexPiece> pieces;

  /// |area of the outer boundary| minus the |area of each hole|.
  double polygon_area = 0.0;
  /// Sum of the (positive) piece areas.
  double pieces_area = 0.0;
  /// |polygon_area - pieces_area|.
  double area_error = 0.0;

  /// Boundary vertices, holes included.
  size_t input_vertices = 0;
  /// Vertices over all emitted pieces (a vertex shared by k pieces counts k).
  size_t output_vertices = 0;

  /// Reflex vertices on the outer boundary (holes are not counted).
  size_t reflex_vertices = 0;
  /// Triangles produced by stage 1 (ear clipping).
  size_t ears_clipped = 0;
  /// Successful convex merges performed by stage 2 (bridging).
  size_t bridges = 0;
};

namespace detail {

struct ConvexDecompImpl {
  typedef std::array<int, 3> Tri;

  struct EdgeMerge {
    int a;
    int b;
    int u; // piece a contains u -> v, piece b contains v -> u
    int v;
    double score;
  };

  const ConvexDecompositionOptions& opt;
  ConvexDecompositionResult& result;

  /// Node id -> point.  The two sides of a hole bridge are separate nodes at
  /// the same coordinates, so the slit is not read as an internal edge.
  std::vector<Point*> node_point;
  /// Node id -> position in the caller's input (used by NearestIndex).
  std::vector<size_t> node_orig;

  ConvexDecompImpl(const ConvexDecompositionOptions& o, ConvexDecompositionResult& r)
      : opt(o), result(r)
  {
  }

  // ---- geometry helpers --------------------------------------------------

  static double Cross3(const Point& a, const Point& b, const Point& c)
  {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  }

  static double SignedArea(const std::vector<Point*>& p)
  {
    double a = 0.0;
    const size_t n = p.size();
    for (size_t i = 0; i < n; ++i)
      a += p[i]->x * p[(i + 1) % n]->y - p[(i + 1) % n]->x * p[i]->y;
    return a * 0.5;
  }

  static bool SamePoint(const Point* a, const Point* b)
  {
    return a->x == b->x && a->y == b->y;
  }

  static bool PointInTriangle(const Point& p, const Point& a, const Point& b, const Point& c)
  {
    const double d1 = Cross3(a, b, p);
    const double d2 = Cross3(b, c, p);
    const double d3 = Cross3(c, a, p);
    const bool neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    const bool pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    return !(neg && pos);
  }

  static bool ProperlyIntersect(const Point& a, const Point& b, const Point& c, const Point& d)
  {
    const double d1 = Cross3(c, d, a);
    const double d2 = Cross3(c, d, b);
    const double d3 = Cross3(a, b, c);
    const double d4 = Cross3(a, b, d);
    return ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) &&
           ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0));
  }

  static bool OnSegment(const Point& p, const Point& a, const Point& b)
  {
    if (Cross3(a, b, p) != 0.0)
      return false;
    return p.x >= std::min(a.x, b.x) && p.x <= std::max(a.x, b.x) &&
           p.y >= std::min(a.y, b.y) && p.y <= std::max(a.y, b.y);
  }

  // ---- node-level helpers ------------------------------------------------

  double NodeArea(const std::vector<int>& ring) const
  {
    double a = 0.0;
    const size_t n = ring.size();
    for (size_t i = 0; i < n; ++i) {
      const Point& p = *node_point[ring[i]];
      const Point& q = *node_point[ring[(i + 1) % n]];
      a += p.x * q.y - q.x * p.y;
    }
    return a * 0.5;
  }

  bool IsConvexNodes(const std::vector<int>& ring, double tolerance) const
  {
    const size_t n = ring.size();
    if (n < 3)
      return false;
    for (size_t i = 0; i < n; ++i)
      if (Cross3(*node_point[ring[i]], *node_point[ring[(i + 1) % n]],
                 *node_point[ring[(i + 2) % n]]) < -tolerance)
        return false;
    return true;
  }

  // Removes consecutive duplicates (bridge slits) and the wrap-around repeat.
  static void CollapseDuplicates(std::vector<int>& ring)
  {
    std::vector<int> out;
    out.reserve(ring.size());
    for (int id : ring)
      if (out.empty() || out.back() != id)
        out.push_back(id);
    if (out.size() > 1 && out.front() == out.back())
      out.pop_back();
    ring.swap(out);
  }

  // ---- hole bridging -----------------------------------------------------

  // Turns the outer boundary plus holes into one CCW node ring.  Each hole is
  // attached by a zero-width bridge (slit), which is what lets stage 1 treat a
  // multiply connected region as a single simple ring.
  std::vector<int> PrepareRing(const std::vector<Point*>& outer,
                               const std::vector<std::vector<Point*>>& holes,
                               const std::map<const Point*, size_t>& orig_of)
  {
    std::vector<Point*> ring = outer;
    if (SignedArea(ring) < 0.0)
      std::reverse(ring.begin(), ring.end());

    std::vector<std::vector<Point*>> prepared = holes;
    for (auto& h : prepared) {
      if (h.size() < 3)
        throw std::runtime_error("ConvexDecomposition: hole needs at least 3 vertices");
      if (SignedArea(h) > 0.0) // holes wind CW inside the CCW outer ring
        std::reverse(h.begin(), h.end());
    }

    // Right to left so a later bridge cannot cross an earlier one.
    std::sort(prepared.begin(), prepared.end(),
              [](const std::vector<Point*>& a, const std::vector<Point*>& b) {
                double ax = -std::numeric_limits<double>::infinity();
                double bx = -std::numeric_limits<double>::infinity();
                for (auto* p : a)
                  ax = std::max(ax, p->x);
                for (auto* p : b)
                  bx = std::max(bx, p->x);
                return ax > bx;
              });

    for (const auto& hole : prepared)
      if (!MergeHole(ring, hole))
        throw std::runtime_error(
            "ConvexDecomposition: no legal bridge for a hole (self-intersecting or touching input?)");

    // Materialise nodes.  Identical coordinates on the two sides of a bridge
    // get distinct ids because they are distinct occurrences in the ring.
    std::vector<int> nodes;
    nodes.reserve(ring.size());
    node_point.clear();
    node_orig.clear();
    for (auto* p : ring) {
      nodes.push_back(static_cast<int>(node_point.size()));
      node_point.push_back(p);
      auto it = orig_of.find(p);
      node_orig.push_back(it == orig_of.end() ? 0 : it->second);
    }
    return nodes;
  }

  bool MergeHole(std::vector<Point*>& ring, const std::vector<Point*>& hole)
  {
    struct Candidate {
      size_t ri;
      size_t hi;
      double score;
    };
    std::vector<Candidate> candidates;

    for (size_t ri = 0; ri < ring.size(); ++ri) {
      for (size_t hi = 0; hi < hole.size(); ++hi) {
        const Point& a = *hole[hi];
        const Point& b = *ring[ri];
        if (a == b)
          continue;
        if (!BridgeIsValid(ring, hole, hi, ri))
          continue;
        candidates.push_back({ ri, hi, BridgeScore(*ring[ri], *hole[hi], ring, ri, hi) });
      }
    }
    if (candidates.empty())
      return false;

    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.score < b.score; });
    const Candidate& best = candidates.front();

    std::vector<Point*> merged;
    merged.reserve(ring.size() + hole.size() + 2);
    for (size_t i = 0; i <= best.ri; ++i)
      merged.push_back(ring[i]);
    for (size_t k = 0; k < hole.size(); ++k)
      merged.push_back(hole[(best.hi + k) % hole.size()]);
    merged.push_back(hole[best.hi]); // back along the bridge
    merged.push_back(ring[best.ri]); // other slit side
    for (size_t i = best.ri + 1; i < ring.size(); ++i)
      merged.push_back(ring[i]);

    ring.swap(merged);
    return true;
  }

  double BridgeScore(const Point& a, const Point& b, const std::vector<Point*>& ring, size_t ia,
                     size_t ib) const
  {
    if (opt.bridge == BridgeStrategy::Shortest) {
      const double dx = a.x - b.x;
      const double dy = a.y - b.y;
      return dx * dx + dy * dy;
    }
    // Index preference: hole vertex first, then the boundary vertex.
    return static_cast<double>(ia) + static_cast<double>(ib) / (ring.size() + 1.0);
  }

  static bool BridgeIsValid(const std::vector<Point*>& ring, const std::vector<Point*>& hole,
                            size_t hi, size_t ri)
  {
    const Point& a = *hole[hi];
    const Point& b = *ring[ri];

    // No other vertex may sit on the bridge.
    for (size_t i = 0; i < ring.size(); ++i)
      if (i != ri && OnSegment(*ring[i], a, b))
        return false;
    for (size_t i = 0; i < hole.size(); ++i)
      if (i != hi && OnSegment(*hole[i], a, b))
        return false;

    // The bridge must not cross an edge.
    for (size_t i = 0; i < ring.size(); ++i) {
      const size_t k = (i + 1) % ring.size();
      if (i == ri || k == ri)
        continue;
      if (ProperlyIntersect(a, b, *ring[i], *ring[k]))
        return false;
    }
    for (size_t i = 0; i < hole.size(); ++i) {
      const size_t k = (i + 1) % hole.size();
      if (i == hi || k == hi)
        continue;
      if (ProperlyIntersect(a, b, *hole[i], *hole[k]))
        return false;
    }
    return true;
  }

  // ---- stage 1: ear clipping ---------------------------------------------

  std::vector<Tri> EarClip(const std::vector<int>& ring)
  {
    const size_t n = ring.size();
    if (n < 3)
      throw std::runtime_error("ConvexDecomposition: polygon needs at least 3 vertices");

    std::vector<int> next(n), prev(n);
    for (size_t i = 0; i < n; ++i) {
      next[i] = static_cast<int>((i + 1) % n);
      prev[i] = static_cast<int>((i + n - 1) % n);
    }

    std::vector<Tri> triangles;
    triangles.reserve(n - 2);
    size_t remaining = n;
    int start = 0;

    while (remaining > 3) {
      bool clipped = false;
      int i = start;
      do {
        const int p = prev[i];
        const int q = next[i];
        if (IsEar(next, prev, i)) {
          triangles.push_back({{ p, i, q }});
          next[p] = q;
          prev[q] = p;
          if (i == start)
            start = q;
          remaining--;
          clipped = true;
          break;
        }
        i = next[i];
      } while (i != start);

      if (!clipped) {
        // No strict ear left (collinear or coincident bridge nodes).  Clip a
        // convex corner of positive area; stage 3 catches any area mistake.
        i = start;
        do {
          const int p = prev[i];
          const int q = next[i];
          if (Cross3(*node_point[p], *node_point[i], *node_point[q]) > 0.0) {
            triangles.push_back({{ p, i, q }});
            next[p] = q;
            prev[q] = p;
            if (i == start)
              start = q;
            remaining--;
            clipped = true;
            break;
          }
          i = next[i];
        } while (i != start);
        if (!clipped)
          throw std::runtime_error(
              "ConvexDecomposition: ear clipping failed (degenerate or self-intersecting polygon)");
      }
    }

    triangles.push_back({{ prev[start], start, next[start] }});
    return triangles;
  }

  bool IsEar(const std::vector<int>& next, const std::vector<int>& prev, int i) const
  {
    const int p = prev[i];
    const int q = next[i];
    const Point& a = *node_point[p];
    const Point& b = *node_point[i];
    const Point& c = *node_point[q];
    if (Cross3(a, b, c) <= 0.0)
      return false;

    for (int j = next[q]; j != p; j = next[j]) {
      if (j == p || j == i || j == q)
        continue;
      const Point& v = *node_point[j];
      if (SamePoint(&v, &a) || SamePoint(&v, &b) || SamePoint(&v, &c))
        continue;
      if (PointInTriangle(v, a, b, c))
        return false;
    }
    return true;
  }

  // ---- stage 2: bridging / convex merging --------------------------------

  // Re-joins ear triangles across their shared diagonals.  A bridge is only
  // taken when the union stays convex, so the pieces grow into maximal convex
  // blocks.  The minimum piece size is enforced afterwards, in stage 3.
  void BridgeIntoConvexBlocks(const std::vector<Tri>& triangles)
  {
    std::vector<std::vector<int>> pieces;
    pieces.reserve(triangles.size());
    for (const auto& t : triangles) {
      std::vector<int> piece(t.begin(), t.end());
      CollapseDuplicates(piece);
      if (piece.size() < 3)
        continue;
      if (NodeArea(piece) < 0.0)
        std::reverse(piece.begin(), piece.end());
      pieces.push_back(piece);
    }

    bool merged_any = true;
    while (merged_any) {
      merged_any = false;

      std::vector<EdgeMerge> candidates;
      CollectMergeCandidates(pieces, candidates);
      std::sort(candidates.begin(), candidates.end(),
                [](const EdgeMerge& x, const EdgeMerge& y) { return x.score < y.score; });

      for (const auto& cand : candidates) {
        if (pieces[cand.a].empty() || pieces[cand.b].empty())
          continue;
        std::vector<int> merged = MergeAcross(pieces[cand.a], pieces[cand.b], cand.u, cand.v);
        if (merged.empty() || !IsConvexNodes(merged, 1e-12))
          continue;
        pieces[cand.a] = merged;
        pieces[cand.b].clear();
        result.bridges++;
        merged_any = true;
        break; // adjacency changed: rebuild
      }
    }

    for (const auto& piece : pieces)
      EmitPiece(piece);
  }

  void CollectMergeCandidates(const std::vector<std::vector<int>>& pieces,
                              std::vector<EdgeMerge>& out) const
  {
    struct Dir {
      int piece;
      int from;
      int to;
    };
    std::map<std::pair<int, int>, std::vector<Dir>> edges;

    for (int pi = 0; pi < static_cast<int>(pieces.size()); ++pi) {
      const auto& piece = pieces[pi];
      if (piece.empty())
        continue;
      for (size_t i = 0; i < piece.size(); ++i) {
        const int u = piece[i];
        const int v = piece[(i + 1) % piece.size()];
        edges[std::make_pair(std::min(u, v), std::max(u, v))].push_back({ pi, u, v });
      }
    }

    for (const auto& entry : edges) {
      const auto& dirs = entry.second;
      if (dirs.size() != 2)
        continue;
      const Dir& d0 = dirs[0];
      const Dir& d1 = dirs[1];
      if (d0.from != d1.to || d0.to != d1.from)
        continue; // same direction => not a manifold internal diagonal
      out.push_back({ d0.piece, d1.piece, d0.from, d0.to,
                      DiagonalScore(d0.from, d0.to) });
    }
  }

  double DiagonalScore(int u, int v) const
  {
    if (opt.bridge == BridgeStrategy::Shortest) {
      const Point& a = *node_point[u];
      const Point& b = *node_point[v];
      return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
    }
    const long long d = std::llabs(static_cast<long long>(node_orig[u]) -
                                   static_cast<long long>(node_orig[v]));
    return static_cast<double>(d);
  }

  // Union of two CCW pieces that share a->u->v and b->v->u.
  std::vector<int> MergeAcross(const std::vector<int>& a, const std::vector<int>& b, int u,
                               int v) const
  {
    auto pos_of = [](const std::vector<int>& p, int id) -> long long {
      for (size_t i = 0; i < p.size(); ++i)
        if (p[i] == id)
          return static_cast<long long>(i);
      return -1;
    };
    const long long av = pos_of(a, v);
    const long long bu = pos_of(b, u);
    if (av < 0 || bu < 0 || pos_of(a, u) < 0 || pos_of(b, v) < 0)
      return {};

    // Both pieces are CCW and share the seam: a contains u -> v and b contains
    // v -> u.  Starting at v, follow a forward to u, then follow b forward
    // from u to v (both seam endpoints dropped on the b side).
    std::vector<int> merged;
    for (size_t k = 0; k < a.size(); ++k) {
      merged.push_back(a[(static_cast<size_t>(av) + k) % a.size()]);
      if (merged.back() == u)
        break;
    }
    for (size_t k = 1; k < b.size(); ++k) {
      const int id = b[(static_cast<size_t>(bu) + k) % b.size()];
      if (id == v)
        break;
      merged.push_back(id);
    }
    return merged;
  }

  void EmitPiece(const std::vector<int>& ring)
  {
    std::vector<int> ids = ring;
    CollapseDuplicates(ids);
    if (ids.size() < 3)
      return;
    if (NodeArea(ids) < 0.0)
      std::reverse(ids.begin(), ids.end());
    if (!IsConvexNodes(ids, 1e-12))
      throw std::runtime_error("ConvexDecomposition: produced a non-convex piece");
    const double area = NodeArea(ids);
    if (area <= opt.degenerate_area)
      return;
    ConvexPiece piece;
    piece.points.reserve(ids.size());
    for (int id : ids)
      piece.points.push_back(node_point[id]);
    result.pieces.push_back(piece);
    result.pieces_area += area;
  }

  // ---- stage 3: validate and enforce the piece area floor ----------------

  void EnforceMinimumPieceArea()
  {
    if (!(opt.min_piece_area > 0.0) || result.pieces.size() < 2)
      return;

    bool changed = true;
    while (changed) {
      changed = false;
      for (size_t i = 0; i < result.pieces.size(); ++i) {
        if (result.pieces[i].Area() >= opt.min_piece_area)
          continue;
        bool absorbed = false;
        for (size_t j = 0; j < result.pieces.size() && !absorbed; ++j) {
          if (i == j)
            continue;
          std::vector<Point*> merged = UnionConvex(result.pieces[i], result.pieces[j]);
          if (merged.empty())
            continue;
          ConvexPiece piece;
          piece.points = merged;
          if (!IsConvexPoints(piece.points))
            continue;
          result.pieces[j] = piece;
          result.pieces.erase(result.pieces.begin() + static_cast<std::ptrdiff_t>(i));
          absorbed = true;
          changed = true;
        }
        if (!absorbed)
          throw std::runtime_error(
              "ConvexDecomposition: piece smaller than min_piece_area could not be merged");
        break;
      }
    }

    result.pieces_area = 0.0;
    result.output_vertices = 0;
    for (const auto& piece : result.pieces) {
      result.pieces_area += piece.Area();
      result.output_vertices += piece.points.size();
    }
    result.area_error = std::fabs(result.polygon_area - result.pieces_area);
  }

  static bool IsConvexPoints(const std::vector<Point*>& p)
  {
    const size_t n = p.size();
    if (n < 3)
      return false;
    for (size_t i = 0; i < n; ++i)
      if (Cross3(*p[i], *p[(i + 1) % n], *p[(i + 2) % n]) < -1e-12)
        return false;
    return true;
  }

  // Concatenates two CCW pieces that share a full edge.
  static std::vector<Point*> UnionConvex(const ConvexPiece& a, const ConvexPiece& b)
  {
    const size_t na = a.points.size();
    const size_t nb = b.points.size();
    for (size_t i = 0; i < na; ++i) {
      const Point* u = a.points[i];
      const Point* v = a.points[(i + 1) % na];
      for (size_t j = 0; j < nb; ++j) {
        const Point* w = b.points[j];
        const Point* z = b.points[(j + 1) % nb];
        if (!SamePoint(u, z) || !SamePoint(v, w))
          continue;
        std::vector<Point*> out;
        out.reserve(na + nb - 2);
        // a: v -> ... -> u
        for (size_t k = 0; k < na; ++k) {
          const Point* p = a.points[(i + 1 + k) % na];
          out.push_back(const_cast<Point*>(p));
          if (p == u)
            break;
        }
        // b: u -> ... -> v walking forwards, dropping both bridge endpoints
        size_t start = (j + 1) % nb; // == u
        for (size_t k = 1; k + 1 < nb; ++k) {
          const Point* p = b.points[(start + k) % nb];
          out.push_back(const_cast<Point*>(p));
        }
        return out;
      }
    }
    return {};
  }

  // ---- driver -------------------------------------------------------------

  void Run(const std::vector<Point*>& outer, const std::vector<std::vector<Point*>>& holes)
  {
    std::map<const Point*, size_t> orig_of;
    size_t order = 0;
    for (auto* p : outer)
      orig_of[p] = order++;
    for (const auto& h : holes)
      for (auto* p : h)
        orig_of[p] = order++;

    std::vector<int> ring = PrepareRing(outer, holes, orig_of);

    // Reflex vertices of the outer boundary (holes are not counted).
    {
      std::vector<Point*> ccw = outer;
      if (SignedArea(ccw) < 0.0)
        std::reverse(ccw.begin(), ccw.end());
      const size_t n = ccw.size();
      for (size_t i = 0; i < n; ++i)
        if (Cross3(*ccw[(i + n - 1) % n], *ccw[i], *ccw[(i + 1) % n]) < 0.0)
          result.reflex_vertices++;
    }

    std::vector<Tri> triangles = EarClip(ring);
    result.ears_clipped = triangles.size();
    BridgeIntoConvexBlocks(triangles);
    EnforceMinimumPieceArea();

    // Stage 3: validate.
    result.output_vertices = 0;
    for (const auto& piece : result.pieces)
      result.output_vertices += piece.points.size();

    result.area_error = std::fabs(result.polygon_area - result.pieces_area);
    const double scale = std::max(1.0, std::fabs(result.polygon_area));
    if (result.pieces.empty())
      throw std::runtime_error("ConvexDecomposition: no convex pieces produced");
    if (result.area_error > opt.area_epsilon * scale)
      throw std::runtime_error("ConvexDecomposition: piece areas do not match polygon area");
    for (const auto& piece : result.pieces)
      if (!IsConvexPoints(piece.points))
        throw std::runtime_error("ConvexDecomposition: produced a non-convex piece");
  }
};

} // namespace detail

/**
 * Decompose a simple polygon, optionally with holes, into convex pieces.
 *
 * The outer boundary may wind either way; holes are re-oriented as needed.
 * Returned pieces alias the caller's points, which must outlive the result.
 * Throws std::runtime_error when the input cannot be decomposed.
 */
inline ConvexDecompositionResult DecomposePolygon(
    const std::vector<Point*>& outer,
    const std::vector<std::vector<Point*>>& holes = std::vector<std::vector<Point*>>(),
    const ConvexDecompositionOptions& options = ConvexDecompositionOptions())
{
  if (outer.size() < 3)
    throw std::runtime_error("ConvexDecomposition: polygon needs at least 3 vertices");
  for (const auto& h : holes)
    if (h.size() < 3)
      throw std::runtime_error("ConvexDecomposition: hole needs at least 3 vertices");

  ConvexDecompositionResult result;
  result.input_vertices = outer.size();
  for (const auto& h : holes)
    result.input_vertices += h.size();

  result.polygon_area = std::fabs(detail::ConvexDecompImpl::SignedArea(outer));
  for (const auto& h : holes)
    result.polygon_area -= std::fabs(detail::ConvexDecompImpl::SignedArea(h));
  if (!(result.polygon_area > 0.0))
    throw std::runtime_error(
        "ConvexDecomposition: polygon has zero or negative area (degenerate, self-intersecting, "
        "or holes cover it)");

  detail::ConvexDecompImpl impl(options, result);
  impl.Run(outer, holes);
  return result;
}

}

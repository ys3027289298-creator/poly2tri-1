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
#include "shapes.h"

#include "../sweep/cdt.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <utility>

namespace p2t {

Point::Point(double x, double y) : x(x), y(y)
{
}

std::ostream& operator<<(std::ostream& out, const Point& point) {
  return out << point.x << "," << point.y;
}

Triangle::Triangle(Point& a, Point& b, Point& c)
{
  points_[0] = &a; points_[1] = &b; points_[2] = &c;
  neighbors_[0] = nullptr; neighbors_[1] = nullptr; neighbors_[2] = nullptr;
  constrained_edge[0] = constrained_edge[1] = constrained_edge[2] = false;
  delaunay_edge[0] = delaunay_edge[1] = delaunay_edge[2] = false;
  interior_ = false;
}

// Update neighbor pointers
void Triangle::MarkNeighbor(Point* p1, Point* p2, Triangle* t)
{
  if ((p1 == points_[2] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[2]))
    neighbors_[0] = t;
  else if ((p1 == points_[0] && p2 == points_[2]) || (p1 == points_[2] && p2 == points_[0]))
    neighbors_[1] = t;
  else if ((p1 == points_[0] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[0]))
    neighbors_[2] = t;
  else
    assert(0);
}

// Exhaustive search to update neighbor pointers
void Triangle::MarkNeighbor(Triangle& t)
{
  if (t.Contains(points_[1], points_[2])) {
    neighbors_[0] = &t;
    t.MarkNeighbor(points_[1], points_[2], this);
  } else if (t.Contains(points_[0], points_[2])) {
    neighbors_[1] = &t;
    t.MarkNeighbor(points_[0], points_[2], this);
  } else if (t.Contains(points_[0], points_[1])) {
    neighbors_[2] = &t;
    t.MarkNeighbor(points_[0], points_[1], this);
  }
}

/**
 * Clears all references to all other triangles and points
 */
void Triangle::Clear()
{
    Triangle *t;
    for (auto& neighbor : neighbors_) {
      t = neighbor;
      if (t != nullptr) {
        t->ClearNeighbor(this);
      }
    }
    ClearNeighbors();
    points_[0]=points_[1]=points_[2] = nullptr;
}

void Triangle::ClearNeighbor(const Triangle *triangle )
{
    if( neighbors_[0] == triangle )
    {
        neighbors_[0] = nullptr;
    }
    else if( neighbors_[1] == triangle )
    {
        neighbors_[1] = nullptr;
    }
    else
    {
        neighbors_[2] = nullptr;
    }
}

void Triangle::ClearNeighbors()
{
  neighbors_[0] = nullptr;
  neighbors_[1] = nullptr;
  neighbors_[2] = nullptr;
}

void Triangle::ClearDelunayEdges()
{
  delaunay_edge[0] = delaunay_edge[1] = delaunay_edge[2] = false;
}

Point* Triangle::OppositePoint(Triangle& t, const Point& p)
{
  Point *cw = t.PointCW(p);
  return PointCW(*cw);
}

// Legalized triangle by rotating clockwise around point(0)
void Triangle::Legalize(Point& point)
{
  points_[1] = points_[0];
  points_[0] = points_[2];
  points_[2] = &point;
}

// Legalize triagnle by rotating clockwise around oPoint
void Triangle::Legalize(Point& opoint, Point& npoint)
{
  if (&opoint == points_[0]) {
    points_[1] = points_[0];
    points_[0] = points_[2];
    points_[2] = &npoint;
  } else if (&opoint == points_[1]) {
    points_[2] = points_[1];
    points_[1] = points_[0];
    points_[0] = &npoint;
  } else if (&opoint == points_[2]) {
    points_[0] = points_[2];
    points_[2] = points_[1];
    points_[1] = &npoint;
  } else {
    assert(0);
  }
}

int Triangle::Index(const Point* p)
{
  if (p == points_[0]) {
    return 0;
  } else if (p == points_[1]) {
    return 1;
  } else if (p == points_[2]) {
    return 2;
  }
  assert(0);
  return -1;
}

int Triangle::EdgeIndex(const Point* p1, const Point* p2)
{
  if (points_[0] == p1) {
    if (points_[1] == p2) {
      return 2;
    } else if (points_[2] == p2) {
      return 1;
    }
  } else if (points_[1] == p1) {
    if (points_[2] == p2) {
      return 0;
    } else if (points_[0] == p2) {
      return 2;
    }
  } else if (points_[2] == p1) {
    if (points_[0] == p2) {
      return 1;
    } else if (points_[1] == p2) {
      return 0;
    }
  }
  return -1;
}

void Triangle::MarkConstrainedEdge(int index)
{
  constrained_edge[index] = true;
}

void Triangle::MarkConstrainedEdge(Edge& edge)
{
  MarkConstrainedEdge(edge.p, edge.q);
}

// Mark edge as constrained
void Triangle::MarkConstrainedEdge(Point* p, Point* q)
{
  if ((q == points_[0] && p == points_[1]) || (q == points_[1] && p == points_[0])) {
    constrained_edge[2] = true;
  } else if ((q == points_[0] && p == points_[2]) || (q == points_[2] && p == points_[0])) {
    constrained_edge[1] = true;
  } else if ((q == points_[1] && p == points_[2]) || (q == points_[2] && p == points_[1])) {
    constrained_edge[0] = true;
  }
}

// The point counter-clockwise to given point
Point* Triangle::PointCW(const Point& point)
{
  if (&point == points_[0]) {
    return points_[2];
  } else if (&point == points_[1]) {
    return points_[0];
  } else if (&point == points_[2]) {
    return points_[1];
  }
  assert(0);
  return nullptr;
}

// The point counter-clockwise to given point
Point* Triangle::PointCCW(const Point& point)
{
  if (&point == points_[0]) {
    return points_[1];
  } else if (&point == points_[1]) {
    return points_[2];
  } else if (&point == points_[2]) {
    return points_[0];
  }
  assert(0);
  return nullptr;
}

// The neighbor across to given point
Triangle* Triangle::NeighborAcross(const Point& point)
{
  if (&point == points_[0]) {
    return neighbors_[0];
  } else if (&point == points_[1]) {
    return neighbors_[1];
  }
  return neighbors_[2];
}

// The neighbor clockwise to given point
Triangle* Triangle::NeighborCW(const Point& point)
{
  if (&point == points_[0]) {
    return neighbors_[1];
  } else if (&point == points_[1]) {
    return neighbors_[2];
  }
  return neighbors_[0];
}

// The neighbor counter-clockwise to given point
Triangle* Triangle::NeighborCCW(const Point& point)
{
  if (&point == points_[0]) {
    return neighbors_[2];
  } else if (&point == points_[1]) {
    return neighbors_[0];
  }
  return neighbors_[1];
}

bool Triangle::GetConstrainedEdgeCCW(const Point& p)
{
  if (&p == points_[0]) {
    return constrained_edge[2];
  } else if (&p == points_[1]) {
    return constrained_edge[0];
  }
  return constrained_edge[1];
}

bool Triangle::GetConstrainedEdgeCW(const Point& p)
{
  if (&p == points_[0]) {
    return constrained_edge[1];
  } else if (&p == points_[1]) {
    return constrained_edge[2];
  }
  return constrained_edge[0];
}

void Triangle::SetConstrainedEdgeCCW(const Point& p, bool ce)
{
  if (&p == points_[0]) {
    constrained_edge[2] = ce;
  } else if (&p == points_[1]) {
    constrained_edge[0] = ce;
  } else {
    constrained_edge[1] = ce;
  }
}

void Triangle::SetConstrainedEdgeCW(const Point& p, bool ce)
{
  if (&p == points_[0]) {
    constrained_edge[1] = ce;
  } else if (&p == points_[1]) {
    constrained_edge[2] = ce;
  } else {
    constrained_edge[0] = ce;
  }
}

bool Triangle::GetDelunayEdgeCCW(const Point& p)
{
  if (&p == points_[0]) {
    return delaunay_edge[2];
  } else if (&p == points_[1]) {
    return delaunay_edge[0];
  }
  return delaunay_edge[1];
}

bool Triangle::GetDelunayEdgeCW(const Point& p)
{
  if (&p == points_[0]) {
    return delaunay_edge[1];
  } else if (&p == points_[1]) {
    return delaunay_edge[2];
  }
  return delaunay_edge[0];
}

void Triangle::SetDelunayEdgeCCW(const Point& p, bool e)
{
  if (&p == points_[0]) {
    delaunay_edge[2] = e;
  } else if (&p == points_[1]) {
    delaunay_edge[0] = e;
  } else {
    delaunay_edge[1] = e;
  }
}

void Triangle::SetDelunayEdgeCW(const Point& p, bool e)
{
  if (&p == points_[0]) {
    delaunay_edge[1] = e;
  } else if (&p == points_[1]) {
    delaunay_edge[2] = e;
  } else {
    delaunay_edge[0] = e;
  }
}

void Triangle::DebugPrint()
{
  std::cout << *points_[0] << " " << *points_[1] << " " << *points_[2] << std::endl;
}

bool Triangle::CircumcicleContains(const Point& point) const
{
  assert(IsCounterClockwise());
  const double dx = points_[0]->x - point.x;
  const double dy = points_[0]->y - point.y;
  const double ex = points_[1]->x - point.x;
  const double ey = points_[1]->y - point.y;
  const double fx = points_[2]->x - point.x;
  const double fy = points_[2]->y - point.y;

  const double ap = dx * dx + dy * dy;
  const double bp = ex * ex + ey * ey;
  const double cp = fx * fx + fy * fy;

  return (dx * (fy * bp - cp * ey) - dy * (fx * bp - cp * ex) + ap * (fx * ey - fy * ex)) < 0;
}

bool Triangle::IsCounterClockwise() const
{
  return (points_[1]->x - points_[0]->x) * (points_[2]->y - points_[0]->y) -
             (points_[2]->x - points_[0]->x) * (points_[1]->y - points_[0]->y) >
         0;
}

bool IsDelaunay(const std::vector<p2t::Triangle*>& triangles)
{
  for (const auto triangle : triangles) {
    for (const auto other : triangles) {
      if (triangle == other) {
        continue;
      }
      for (int i = 0; i < 3; ++i) {
        if (triangle->CircumcicleContains(*other->GetPoint(i))) {
          return false;
        }
      }
    }
  }
  return true;
}

namespace boolean_approx {

namespace {

double Cross2(const Point& o, const Point& a, const Point& b)
{
  return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Proper-or-boundary segment intersection, used by the simple-polygon check.
bool SegmentsIntersect(const Point& p1, const Point& p2,
                       const Point& q1, const Point& q2)
{
  const double d1 = Cross2(q1, q2, p1);
  const double d2 = Cross2(q1, q2, p2);
  const double d3 = Cross2(p1, p2, q1);
  const double d4 = Cross2(p1, p2, q2);
  if (((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) &&
      ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))) {
    return true;
  }
  const auto on_segment = [](const Point& a, const Point& b, const Point& c) {
    return std::fabs(Cross2(a, b, c)) <= 1e-12 &&
           std::min(a.x, b.x) - 1e-12 <= c.x && c.x <= std::max(a.x, b.x) + 1e-12 &&
           std::min(a.y, b.y) - 1e-12 <= c.y && c.y <= std::max(a.y, b.y) + 1e-12;
  };
  return on_segment(q1, q2, p1) || on_segment(q1, q2, p2) ||
         on_segment(p1, p2, q1) || on_segment(p1, p2, q2);
}

// Ray-casting containment in one ordered ring. The boundary itself is not
// treated as inside: StrictlyInsideRegion() filters boundary points before
// points are handed to the CDT, and an operand centroid on the other
// polygon's edge is an ambiguous case that ray casting resolves one way.
bool PointInRing(const Point& p, const std::vector<Point>& ring)
{
  bool inside = false;
  const size_t n = ring.size();
  for (size_t i = 0, j = n - 1; i < n; j = i++) {
    const Point& a = ring[j];
    const Point& b = ring[i];
    const bool crosses = ((a.y > p.y) != (b.y > p.y));
    if (crosses) {
      const double x_hit = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
      if (x_hit > p.x) {
        inside = !inside;
      }
    }
  }
  return inside;
}

double RingArea(const std::vector<Point>& ring)
{
  double twice = 0.0;
  for (size_t i = 0; i < ring.size(); ++i) {
    const Point& a = ring[i];
    const Point& b = ring[(i + 1) % ring.size()];
    twice += a.x * b.y - b.x * a.y;
  }
  return 0.5 * std::fabs(twice);
}

bool OnRing(const Point& p, const std::vector<Point>& ring)
{
  for (size_t i = 0; i < ring.size(); ++i) {
    const Point& a = ring[i];
    const Point& b = ring[(i + 1) % ring.size()];
    if (std::fabs(Cross2(a, b, p)) <= 1e-12 &&
        std::min(a.x, b.x) - 1e-12 <= p.x && p.x <= std::max(a.x, b.x) + 1e-12 &&
        std::min(a.y, b.y) - 1e-12 <= p.y && p.y <= std::max(a.y, b.y) + 1e-12) {
      return true;
    }
  }
  return false;
}

// Strictly inside the outer ring and outside every ring boundary/hole.
bool StrictlyInsideRegion(const Point& p, const Region& region)
{
  if (!PointInRing(p, region.outer) || OnRing(p, region.outer)) {
    return false;
  }
  for (const auto& hole : region.holes) {
    if (PointInRing(p, hole)) {
      return false;
    }
  }
  return true;
}

std::vector<Tri> TriangulateRegion(const Region& region, double spacing)
{
  RequireSimpleRing(region.outer);
  std::vector<Point*> outline;
  for (const Point& p : region.outer) {
    outline.push_back(new Point(p.x, p.y));
  }
  std::vector<std::vector<Point*>> hole_ptrs;
  for (const auto& hole : region.holes) {
    RequireSimpleRing(hole);
    hole_ptrs.emplace_back();
    for (const Point& p : hole) {
      hole_ptrs.back().push_back(new Point(p.x, p.y));
    }
  }
  // Interior Steiner grid, inserted through the existing CDT::AddPoint path,
  // refines the triangulation so narrow overlaps get their own triangles.
  std::vector<Point*> steiners;
  if (spacing > 0.0) {
    double xmin = region.outer[0].x, xmax = xmin;
    double ymin = region.outer[0].y, ymax = ymin;
    for (const Point& p : region.outer) {
      xmin = std::min(xmin, p.x);
      xmax = std::max(xmax, p.x);
      ymin = std::min(ymin, p.y);
      ymax = std::max(ymax, p.y);
    }
    for (double x = xmin + spacing * 0.5; x < xmax; x += spacing) {
      for (double y = ymin + spacing * 0.5; y < ymax; y += spacing) {
        Point candidate(x, y);
        if (StrictlyInsideRegion(candidate, region)) {
          steiners.push_back(new Point(x, y));
        }
      }
    }
  }
  std::vector<Tri> out;
  try {
    CDT cdt(outline);
    for (const auto& hole : hole_ptrs) {
      cdt.AddHole(hole);
    }
    for (Point* s : steiners) {
      cdt.AddPoint(s);
    }
    cdt.Triangulate();
    for (Triangle* t : cdt.GetTriangles()) {
      out.push_back(Tri{*t->GetPoint(0), *t->GetPoint(1), *t->GetPoint(2)});
    }
  } catch (...) {
    for (Point* p : outline) {
      delete p;
    }
    for (auto& hole : hole_ptrs) {
      for (Point* p : hole) {
        delete p;
      }
    }
    for (Point* p : steiners) {
      delete p;
    }
    throw;
  }
  for (Point* p : outline) {
    delete p;
  }
  for (auto& hole : hole_ptrs) {
    for (Point* p : hole) {
      delete p;
    }
  }
  for (Point* p : steiners) {
    delete p;
  }
  return out;
}

// Key on exact coordinates so coincident edges from the two triangulations
// collapse onto the same undirected edge.
std::pair<std::pair<double, double>, std::pair<double, double>>
EdgeKey(const Point& a, const Point& b)
{
  auto ka = std::make_pair(a.x, a.y);
  auto kb = std::make_pair(b.x, b.y);
  return ka < kb ? std::make_pair(ka, kb) : std::make_pair(kb, ka);
}

// Walk edges that occur exactly once among the kept triangles; each such
// edge borders an unkept triangle and belongs to the merged outline. The
// triangles always come from a single conforming triangulation, so boundary
// edges form closed walks; a multimap with edge consumption also handles
// several walks that touch at one vertex.
std::vector<std::vector<Point>> ExtractContours(const std::vector<Tri>& kept)
{
  struct EdgeInfo {
    Point from, to; // Directed as first seen.
    int count = 0;
  };
  std::map<std::pair<std::pair<double, double>, std::pair<double, double>>,
           EdgeInfo>
      edges;
  const auto record = [&](const Point& from, const Point& to) {
    auto it = edges.find(EdgeKey(from, to));
    if (it == edges.end()) {
      edges.emplace(EdgeKey(from, to), EdgeInfo{from, to, 1});
    } else {
      ++it->second.count;
    }
  };
  for (const Tri& t : kept) {
    record(t.a, t.b);
    record(t.b, t.c);
    record(t.c, t.a);
  }

  using Key = std::pair<double, double>;
  std::multimap<Key, Key> outgoing;
  for (const auto& entry : edges) {
    if (entry.second.count == 1) {
      const EdgeInfo& e = entry.second;
      outgoing.emplace(std::make_pair(e.from.x, e.from.y),
                       std::make_pair(e.to.x, e.to.y));
    }
  }

  std::vector<std::vector<Point>> contours;
  while (!outgoing.empty()) {
    std::vector<Point> loop;
    const Key start = outgoing.begin()->first;
    Key cur = start;
    for (;;) {
      loop.push_back(Point(cur.first, cur.second));
      auto it = outgoing.find(cur);
      if (it == outgoing.end()) {
        break; // Non-closed chain: keep what was walked.
      }
      const Key nxt = it->second;
      outgoing.erase(it);
      cur = nxt;
      if (cur == start) {
        break;
      }
    }
    contours.push_back(loop);
  }
  return contours;
}

} // namespace

Point Tri::Centroid() const
{
  return Point((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
}

double Tri::Area() const
{
  return 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) -
                         (c.x - a.x) * (b.y - a.y));
}

void RequireSimpleRing(const std::vector<Point>& ring)
{
  const size_t n = ring.size();
  if (n < 3) {
    throw std::runtime_error("boolean_approx: ring needs at least 3 points");
  }
  for (size_t i = 0; i < n; ++i) {
    if (ring[i].x == ring[(i + 1) % n].x && ring[i].y == ring[(i + 1) % n].y) {
      throw std::runtime_error("boolean_approx: repeated point in ring");
    }
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      // Adjacent edges share exactly their common endpoint; all other pairs
      // must be disjoint.
      if ((j == i + 1) || (i == 0 && j == n - 1)) {
        continue;
      }
      if (SegmentsIntersect(ring[i], ring[(i + 1) % n], ring[j],
                            ring[(j + 1) % n])) {
        throw std::runtime_error(
            "boolean_approx: self-intersecting ring is not supported");
      }
    }
  }
}

bool PointInRegion(const Point& p, const Region& region)
{
  if (!PointInRing(p, region.outer)) {
    return false;
  }
  for (const auto& hole : region.holes) {
    if (PointInRing(p, hole)) {
      return false;
    }
  }
  return true;
}

double RegionArea(const Region& region)
{
  double area = RingArea(region.outer);
  for (const auto& hole : region.holes) {
    area -= RingArea(hole);
  }
  return area;
}

std::vector<Tri> Triangulate(const Region& region, double spacing)
{
  return TriangulateRegion(region, spacing);
}

Result Apply(const Region& lhs, const Region& rhs, Op op, double spacing)
{
  const std::vector<Tri> tris_a = TriangulateRegion(lhs, spacing);
  const std::vector<Tri> tris_b = TriangulateRegion(rhs, spacing);

  // Classify every triangle of both triangulations by its centroid.
  std::vector<Tri> in_a, out_a, in_b, out_b;
  for (const Tri& t : tris_a) {
    (PointInRegion(t.Centroid(), rhs) ? in_a : out_a).push_back(t);
  }
  for (const Tri& t : tris_b) {
    (PointInRegion(t.Centroid(), lhs) ? in_b : out_b).push_back(t);
  }

  Result result;
  // The retained triangles form disjoint-interior sets, so result.area is
  // exactly the sum of their areas:
  //  - intersection: A's triangles whose centroid is inside B;
  //  - difference:   A's triangles whose centroid is outside B;
  //  - union:        A's outside triangles + B's outside triangles + the
  //                  intersection piece (A's inside triangles). The overlap
  //                  is added once, which makes the area identities hold:
  //                  |A u B| = |A| + |B| - |A n B|,
  //                  |A \\ B| = |A| - |A n B|.
  std::vector<Tri> outline_a, outline_b;
  switch (op) {
  case Op::Intersection:
    outline_a = in_a;
    break;
  case Op::Difference:
    outline_a = out_a;
    break;
  case Op::Union:
    outline_a.reserve(out_a.size() + in_a.size());
    for (const Tri& t : out_a) {
      outline_a.push_back(t);
    }
    for (const Tri& t : in_a) {
      outline_a.push_back(t);
    }
    outline_b = out_b;
    break;
  }
  result.triangles.reserve(outline_a.size() + outline_b.size());
  for (const Tri& t : outline_a) {
    result.triangles.push_back(t);
    result.area += t.Area();
  }
  for (const Tri& t : outline_b) {
    result.triangles.push_back(t);
    result.area += t.Area();
  }

  // Walk the boundary of each side's own conforming triangle set; the two
  // outlines coincide geometrically where the regions meet.
  std::vector<std::vector<Point>> contours_a = ExtractContours(outline_a);
  std::vector<std::vector<Point>> contours_b = ExtractContours(outline_b);
  result.contours.reserve(contours_a.size() + contours_b.size());
  for (auto& loop : contours_a) {
    result.contours.push_back(std::move(loop));
  }
  for (auto& loop : contours_b) {
    result.contours.push_back(std::move(loop));
  }
  return result;
}

} // namespace boolean_approx

} // namespace p2t

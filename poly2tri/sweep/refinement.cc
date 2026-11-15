#include "refinement.h"

#include "../common/shapes.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace p2t {

namespace {

double EdgeLength(const Point& a, const Point& b)
{
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  return std::sqrt(dx * dx + dy * dy);
}

double TriangleArea(const Point& a, const Point& b, const Point& c)
{
  return 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) -
                         (c.x - a.x) * (b.y - a.y));
}

double MaxEdgeLength(const std::vector<Triangle*>& triangles)
{
  double max_length = 0.0;
  for (Triangle* t : triangles) {
    for (int i = 0; i < 3; i++) {
      max_length =
          std::max(max_length, EdgeLength(*t->GetPoint(i), *t->GetPoint((i + 1) % 3)));
    }
  }
  return max_length;
}

/// Sign of the orientation of (a, b, c) with tolerance eps.
int Orient(const Point& a, const Point& b, const Point& c, double eps)
{
  const double v =
      (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
  if (v > eps) {
    return 1;
  }
  if (v < -eps) {
    return -1;
  }
  return 0;
}

bool SamePoint(const Point& a, const Point& b, double eps)
{
  return std::fabs(a.x - b.x) <= eps && std::fabs(a.y - b.y) <= eps;
}

} // namespace

RefinementCriteria::RefinementCriteria()
    : max_edge_length(0.0), min_area(0.0), max_iterations(16), enabled(true)
{
}

Refiner::Refiner(const std::vector<Point*>& polyline) : boundary_(polyline)
{
}

Refiner::~Refiner()
{
  delete cdt_;
  for (Point* p : inserted_points_) {
    delete p;
  }
  for (Point* p : split_points_) {
    delete p;
  }
}

void Refiner::AddHole(const std::vector<Point*>& polyline)
{
  holes_.push_back(polyline);
}

void Refiner::AddPoint(Point* point)
{
  steiner_points_.push_back(point);
}

std::vector<Point*> Refiner::AllPoints() const
{
  std::vector<Point*> all = boundary_;
  for (const std::vector<Point*>& hole : holes_) {
    all.insert(all.end(), hole.begin(), hole.end());
  }
  all.insert(all.end(), steiner_points_.begin(), steiner_points_.end());
  all.insert(all.end(), inserted_points_.begin(), inserted_points_.end());
  return all;
}

std::vector<std::vector<Point*>> Refiner::AllPolylines() const
{
  std::vector<std::vector<Point*>> polylines{boundary_};
  polylines.insert(polylines.end(), holes_.begin(), holes_.end());
  return polylines;
}

void Refiner::RebuildCdt()
{
  // Each CDT creates Edge objects which register themselves in
  // Point::edge_list and are destroyed together with their SweepContext.
  // Clear the stale pointers before the points are reused.
  for (Point* p : AllPoints()) {
    p->edge_list.clear();
  }
  delete cdt_;
  cdt_ = new CDT(boundary_);
  for (const std::vector<Point*>& hole : holes_) {
    cdt_->AddHole(hole);
  }
  for (Point* p : steiner_points_) {
    cdt_->AddPoint(p);
  }
  for (Point* p : inserted_points_) {
    cdt_->AddPoint(p);
  }
  cdt_->Triangulate();
}

bool Refiner::FindViolations(const std::vector<Triangle*>& triangles,
                             const RefinementCriteria& criteria,
                             std::vector<EdgeSplit>& edge_splits,
                             std::vector<Point>& interior_points) const
{
  bool any = false;
  for (Triangle* t : triangles) {
    const Point* p[3]{t->GetPoint(0), t->GetPoint(1), t->GetPoint(2)};
    // Edge i is opposite vertex i.
    double length[3];
    int longest = 0;
    for (int i = 0; i < 3; i++) {
      length[i] = EdgeLength(*p[(i + 1) % 3], *p[(i + 2) % 3]);
      if (length[i] > length[longest]) {
        longest = i;
      }
    }
    const bool edge_violation = criteria.max_edge_length > 0.0 &&
                                length[longest] > criteria.max_edge_length;
    const bool area_violation =
        criteria.min_area > 0.0 &&
        TriangleArea(*p[0], *p[1], *p[2]) < criteria.min_area;
    if (!edge_violation && !area_violation) {
      continue;
    }
    any = true;
    if (edge_violation) {
      const Point* a = p[(longest + 1) % 3];
      const Point* b = p[(longest + 2) % 3];
      const Point midpoint(0.5 * (a->x + b->x), 0.5 * (a->y + b->y));
      if (t->constrained_edge[longest]) {
        // The longest edge lies on the boundary or a hole: split the
        // polyline segment itself so the constraint edge gets shorter too.
        edge_splits.push_back(EdgeSplit{midpoint, a, b});
      } else {
        interior_points.push_back(midpoint);
      }
    } else {
      interior_points.emplace_back((p[0]->x + p[1]->x + p[2]->x) / 3.0,
                                   (p[0]->y + p[1]->y + p[2]->y) / 3.0);
    }
  }
  return any;
}

bool Refiner::SplitSegment(const Point* a, const Point* b, Point* point)
{
  std::vector<std::vector<Point*>*> polylines{&boundary_};
  for (std::vector<Point*>& hole : holes_) {
    polylines.push_back(&hole);
  }
  for (std::vector<Point*>* polyline : polylines) {
    const std::size_t n = polyline->size();
    for (std::size_t i = 0; i < n; i++) {
      const std::size_t j = (i + 1) % n;
      const Point* pi = (*polyline)[i];
      const Point* pj = (*polyline)[j];
      if ((pi == a && pj == b) || (pi == b && pj == a)) {
        polyline->insert(polyline->begin() + static_cast<std::ptrdiff_t>(i) + 1,
                         point);
        return true;
      }
    }
  }
  return false;
}

std::size_t Refiner::InsertPoints(const std::vector<EdgeSplit>& edge_splits,
                                  const std::vector<Point>& interior_points)
{
  std::vector<Point*> known_points = AllPoints();
  double max_coord = 1.0;
  for (const Point* p : known_points) {
    max_coord = std::max(max_coord,
                         std::max(std::fabs(p->x), std::fabs(p->y)));
  }
  const double eps = 1e-9 * max_coord;
  const double eps2 = eps * eps;
  const auto is_known = [&](const Point& candidate) {
    for (const Point* p : known_points) {
      const double dx = p->x - candidate.x;
      const double dy = p->y - candidate.y;
      if (dx * dx + dy * dy <= eps2) {
        return true;
      }
    }
    return false;
  };

  std::size_t added = 0;
  for (const EdgeSplit& split : edge_splits) {
    if (is_known(split.point)) {
      continue;
    }
    Point* point = new Point(split.point.x, split.point.y);
    if (!SplitSegment(split.a, split.b, point)) {
      delete point;
      continue;
    }
    split_points_.push_back(point);
    known_points.push_back(point);
    ++added;
  }
  for (const Point& candidate : interior_points) {
    if (is_known(candidate)) {
      continue;
    }
    Point* point = new Point(candidate.x, candidate.y);
    inserted_points_.push_back(point);
    known_points.push_back(point);
    ++added;
  }
  return added;
}

const RefinementReport& Refiner::Triangulate(const RefinementCriteria& criteria)
{
  report_ = RefinementReport();
  const int max_iterations = std::max(0, criteria.max_iterations);
  RebuildCdt();
  while (true) {
    const std::vector<Triangle*> triangles = cdt_->GetTriangles();
    // A constrained edge crossed by the triangulation is an error and
    // must never be silently accepted.
    CheckConstraintsIntact(triangles, AllPolylines());
    if (!criteria.enabled) {
      break;
    }
    std::vector<EdgeSplit> edge_splits;
    std::vector<Point> interior_points;
    FindViolations(triangles, criteria, edge_splits, interior_points);
    if (edge_splits.empty() && interior_points.empty()) {
      break; // both constraints satisfied: converged
    }
    if (report_.iterations >= max_iterations ||
        InsertPoints(edge_splits, interior_points) == 0) {
      report_.converged = false;
      break;
    }
    ++report_.iterations;
    RebuildCdt();
  }
  const std::vector<Triangle*> triangles = cdt_->GetTriangles();
  report_.triangle_count = triangles.size();
  report_.max_edge_length = MaxEdgeLength(triangles);
  return report_;
}

std::vector<Triangle*> Refiner::GetTriangles()
{
  return cdt_->GetTriangles();
}

const RefinementReport& Refiner::report() const
{
  return report_;
}

void CheckConstraintsIntact(
    const std::vector<Triangle*>& triangles,
    const std::vector<std::vector<Point*>>& polylines)
{
  struct Segment {
    const Point* p;
    const Point* q;
  };
  std::vector<Segment> segments;
  double max_coord = 1.0;
  for (const std::vector<Point*>& polyline : polylines) {
    const std::size_t n = polyline.size();
    for (std::size_t i = 0; i < n; i++) {
      segments.push_back(Segment{polyline[i], polyline[(i + 1) % n]});
    }
    for (const Point* p : polyline) {
      max_coord = std::max(max_coord,
                           std::max(std::fabs(p->x), std::fabs(p->y)));
    }
  }
  const double orient_eps = 1e-9 * max_coord * max_coord;
  const double point_eps = 1e-9 * max_coord;
  for (Triangle* t : triangles) {
    for (int i = 0; i < 3; i++) {
      const Point* c = t->GetPoint(i);
      const Point* d = t->GetPoint((i + 1) % 3);
      for (const Segment& s : segments) {
        // Sharing an endpoint with the constrained edge is fine.
        if (SamePoint(*c, *s.p, point_eps) || SamePoint(*c, *s.q, point_eps) ||
            SamePoint(*d, *s.p, point_eps) || SamePoint(*d, *s.q, point_eps)) {
          continue;
        }
        const int o1 = Orient(*c, *d, *s.p, orient_eps);
        const int o2 = Orient(*c, *d, *s.q, orient_eps);
        const int o3 = Orient(*s.p, *s.q, *c, orient_eps);
        const int o4 = Orient(*s.p, *s.q, *d, orient_eps);
        if (o1 * o2 < 0 && o3 * o4 < 0) {
          std::ostringstream oss;
          oss << "CheckConstraintsIntact: constrained edge (" << s.p->x << ", "
              << s.p->y << ")-(" << s.q->x << ", " << s.q->y
              << ") is crossed by triangle edge (" << c->x << ", " << c->y
              << ")-(" << d->x << ", " << d->y << ")";
          throw std::runtime_error(oss.str());
        }
      }
    }
  }
}

} // namespace p2t

/*
 * Poly2Tri Copyright (c) 2009-2024, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "self_intersection.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

namespace p2t {

namespace {

double Cross2(const Point& a, const Point& b, const Point& c)
{
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool PointsEqual(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}

bool OnSegment(const Point& p, const Point& a, const Point& b)
{
  return std::fpclassify(Cross2(a, b, p)) == FP_ZERO
    && std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x)
    && std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
}

Point IntersectionPoint(const Point& a1, const Point& a2,
                        const Point& b1, const Point& b2)
{
  const double dx = a2.x - a1.x;
  const double dy = a2.y - a1.y;
  const double sx = b2.x - b1.x;
  const double sy = b2.y - b1.y;
  const double denom = dx * sy - dy * sx;
  if (std::fpclassify(denom) != FP_ZERO) {
    const double t = ((b1.x - a1.x) * sy - (b1.y - a1.y) * sx) / denom;
    return Point(a1.x + t * dx, a1.y + t * dy);
  }
  // Collinear overlap: split at both endpoints of the overlap box.
  return b1;
}

struct ParamPoint {
  double parameter;
  Point point;
};

bool ParamLess(const ParamPoint& a, const ParamPoint& b)
{
  if (a.parameter != b.parameter) {
    return a.parameter < b.parameter;
  }
  if (a.point.x != b.point.x) {
    return a.point.x < b.point.x;
  }
  return a.point.y < b.point.y;
}

/// Classifies how edge (a1,a2) meets edge (b1,b2). Shared endpoints between
/// adjacent edges are reported as SI_NONE by the caller through the index
/// filtering done in DetectSelfIntersection.
SelfIntersectionCode ClassifyPair(const Point& a1, const Point& a2,
                                  const Point& b1, const Point& b2)
{
  const double d1 = Cross2(b1, b2, a1);
  const double d2 = Cross2(b1, b2, a2);
  const double d3 = Cross2(a1, a2, b1);
  const double d4 = Cross2(a1, a2, b2);

  const bool a_straddles = (d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0);
  const bool b_straddles = (d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0);
  if (a_straddles && b_straddles) {
    return SI_PROPER_CROSSING;
  }

  // Endpoint of one edge landing in the other edge's interior: T-junction.
  if (OnSegment(a1, b1, b2) && !PointsEqual(a1, b1) && !PointsEqual(a1, b2)) {
    return SI_VERTEX_ON_EDGE;
  }
  if (OnSegment(a2, b1, b2) && !PointsEqual(a2, b1) && !PointsEqual(a2, b2)) {
    return SI_VERTEX_ON_EDGE;
  }
  if (OnSegment(b1, a1, a2) && !PointsEqual(b1, a1) && !PointsEqual(b1, a2)) {
    return SI_VERTEX_ON_EDGE;
  }
  if (OnSegment(b2, a1, a2) && !PointsEqual(b2, a1) && !PointsEqual(b2, a2)) {
    return SI_VERTEX_ON_EDGE;
  }

  const bool share = PointsEqual(a1, b1) || PointsEqual(a1, b2)
    || PointsEqual(a2, b1) || PointsEqual(a2, b2);
  if (share) {
    return SI_NONE; // Common endpoint of two edges; legal polygon join.
  }

  // No shared endpoint and collinear: the segments must overlap to matter.
  if (std::fpclassify(d1) == FP_ZERO && std::fpclassify(d2) == FP_ZERO) {
    const bool separated =
      std::max(a1.x, a2.x) < std::min(b1.x, b2.x)
      || std::max(b1.x, b2.x) < std::min(a1.x, a2.x)
      || std::max(a1.y, a2.y) < std::min(b1.y, b2.y)
      || std::max(b1.y, b2.y) < std::min(a1.y, a2.y);
    if (!separated) {
      return SI_COLLINEAR_OVERLAP;
    }
  }
  return SI_NONE;
}

} // namespace

SelfIntersectionInfo
DetectSelfIntersection(const std::vector<Point*>& contour, int contour_index)
{
  SelfIntersectionInfo result;
  result.contour_index = contour_index;

  const size_t n = contour.size();
  if (n < 3) {
    result.code = SI_DEGENERATE_CONTOUR;
    result.edge_index_a = static_cast<int>(n);
    return result;
  }

  // Repeated coordinates make the sweep edge construction throw or crash.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      if (PointsEqual(*contour[i], *contour[j])) {
        result.code = SI_DUPLICATE_POINT;
        result.edge_index_a = static_cast<int>(i);
        result.edge_index_b = static_cast<int>(j);
        return result;
      }
    }
  }

  for (size_t i = 0; i < n; ++i) {
    const size_t i_next = (i + 1) % n;
    for (size_t j = i + 1; j < n; ++j) {
      const size_t j_next = (j + 1) % n;
      const bool adjacent = (i == 0 && j == n - 1) || (j_next == i)
        || (i_next == j);
      if (adjacent) {
        continue;
      }
      const SelfIntersectionCode code = ClassifyPair(
        *contour[i], *contour[i_next], *contour[j], *contour[j_next]);
      if (code != SI_NONE) {
        result.code = code;
        result.edge_index_a = static_cast<int>(i);
        result.edge_index_b = static_cast<int>(j);
        return result;
      }
    }
  }

  if (std::fpclassify(SignedArea(contour)) == FP_ZERO) {
    result.code = SI_DEGENERATE_CONTOUR;
    result.edge_index_a = 0;
  }
  return result;
}

SelfIntersectionInfo
DetectSelfIntersections(const std::vector<std::vector<Point*> >& contours)
{
  for (size_t c = 0; c < contours.size(); ++c) {
    const SelfIntersectionInfo info = DetectSelfIntersection(contours[c],
      static_cast<int>(c));
    if (info.code != SI_NONE) {
      return info;
    }
  }
  return SelfIntersectionInfo();
}

double SignedArea(const std::vector<Point*>& contour)
{
  double twice_area = 0.0;
  const size_t n = contour.size();
  for (size_t i = 0; i < n; ++i) {
    const Point& a = *contour[i];
    const Point& b = *contour[(i + 1) % n];
    twice_area += a.x * b.y - b.x * a.y;
  }
  return 0.5 * twice_area;
}

std::string SelfIntersectionError::FormatMessage(const SelfIntersectionInfo& info)
{
  static const char* const messages[] = {
    "no self intersection",
    "proper edge crossing",
    "vertex on non-adjacent edge (T-junction)",
    "collinear overlapping edges",
    "duplicate point",
    "degenerate contour"
  };
  const int index = static_cast<int>(info.code);
  std::ostringstream stream;
  stream << "poly2tri: " << messages[index] << " in contour "
         << info.contour_index << " at edge " << info.edge_index_a;
  if (info.edge_index_b >= 0) {
    stream << " vs edge " << info.edge_index_b;
  }
  return stream.str();
}

EvenOddResolver::EvenOddResolver(const std::vector<Point*>& contour)
  : contour_(contour)
{
}

EvenOddResolver::~EvenOddResolver()
{
  for (Point* point : owned_points_) {
    delete point;
  }
}

size_t EvenOddResolver::FindOrAddPoint(const Point& point)
{
  for (size_t i = 0; i < owned_points_.size(); ++i) {
    if (PointsEqual(*owned_points_[i], point)) {
      return i;
    }
  }
  owned_points_.push_back(new Point(point.x, point.y));
  return owned_points_.size() - 1;
}

void EvenOddResolver::BuildSplitEdges()
{
  const size_t n = contour_.size();

  // For every original edge collect its endpoints plus every intersection
  // with any non-adjacent edge, expressed in that edge's parameter space.
  std::vector<std::vector<ParamPoint> > splits(n);
  for (size_t i = 0; i < n; ++i) {
    splits[i].push_back(ParamPoint{0.0, *contour_[i]});
    splits[i].push_back(ParamPoint{1.0, *contour_[(i + 1) % n]});
  }

  for (size_t i = 0; i < n; ++i) {
    const size_t i_next = (i + 1) % n;
    for (size_t j = i + 1; j < n; ++j) {
      const size_t j_next = (j + 1) % n;
      const bool adjacent = (i == 0 && j == n - 1) || (j_next == i)
        || (i_next == j);
      if (adjacent) {
        continue;
      }
      const Point& a1 = *contour_[i];
      const Point& a2 = *contour_[i_next];
      const Point& b1 = *contour_[j];
      const Point& b2 = *contour_[j_next];
      const SelfIntersectionCode code = ClassifyPair(a1, a2, b1, b2);
      if (code == SI_NONE) {
        continue;
      }
      if (code == SI_COLLINEAR_OVERLAP || code == SI_DUPLICATE_POINT) {
        SelfIntersectionInfo info;
        info.code = code;
        info.edge_index_a = static_cast<int>(i);
        info.edge_index_b = static_cast<int>(j);
        info.contour_index = 0;
        throw SelfIntersectionError(info);
      }

      const Point cross = IntersectionPoint(a1, a2, b1, b2);
      const double ti = std::hypot(cross.x - a1.x, cross.y - a1.y)
        / std::hypot(a2.x - a1.x, a2.y - a1.y);
      const double tj = std::hypot(cross.x - b1.x, cross.y - b1.y)
        / std::hypot(b2.x - b1.x, b2.y - b1.y);
      splits[i].push_back(ParamPoint{ti, cross});
      splits[j].push_back(ParamPoint{tj, cross});

      if (code == SI_VERTEX_ON_EDGE) {
        // Record the vertex endpoint too so the parameter list keeps both
        // representations; FindOrAddPoint collapses identical coordinates.
        if (OnSegment(a1, b1, b2)) {
          splits[j].push_back(ParamPoint{tj, a1});
        } else if (OnSegment(a2, b1, b2)) {
          splits[j].push_back(ParamPoint{tj, a2});
        } else if (OnSegment(b1, a1, a2)) {
          splits[i].push_back(ParamPoint{ti, b1});
        } else if (OnSegment(b2, a1, a2)) {
          splits[i].push_back(ParamPoint{ti, b2});
        }
      }
    }
  }

  for (size_t i = 0; i < n; ++i) {
    std::sort(splits[i].begin(), splits[i].end(), ParamLess);
    for (size_t k = 1; k < splits[i].size(); ++k) {
      if (PointsEqual(splits[i][k - 1].point, splits[i][k].point)) {
        continue;
      }
      const size_t first = FindOrAddPoint(splits[i][k - 1].point);
      const size_t second = FindOrAddPoint(splits[i][k].point);
      const auto directed = std::make_pair(first, second);
      if (std::find(split_edges_.begin(), split_edges_.end(), directed)
          == split_edges_.end()) {
        split_edges_.push_back(directed);
      }
    }
  }
}

void EvenOddResolver::ExtractLoops()
{
  // Undirected split-edge graph expressed as directed half-edges. The
  // circular-order face walk below enumerates every planar face, one per
  // side of every split edge; even-odd membership is tested afterwards.
  std::map<size_t, std::vector<size_t> > outgoing;
  for (const auto& edge : split_edges_) {
    outgoing[edge.first].push_back(edge.second);
    outgoing[edge.second].push_back(edge.first);
  }
  // Circular angular order of outgoing half-edges at every node. Face
  // traversal then always takes the edge immediately before the reverse of
  // the incoming edge, which enumerates every planar face exactly once.
  for (auto& entry : outgoing) {
    const Point& origin = *owned_points_[entry.first];
    std::sort(entry.second.begin(), entry.second.end(),
      [&](size_t lhs, size_t rhs) {
        const double al = std::atan2(owned_points_[lhs]->y - origin.y,
                                     owned_points_[lhs]->x - origin.x);
        const double ar = std::atan2(owned_points_[rhs]->y - origin.y,
                                     owned_points_[rhs]->x - origin.x);
        return al < ar;
      });
  }
  std::map<std::pair<size_t, size_t>, bool> used;
  for (const auto& edge : split_edges_) {
    used[edge] = false;
    used[std::make_pair(edge.second, edge.first)] = false;
  }

  for (auto& usage : used) {
    if (usage.second) {
      continue;
    }
    const size_t start = usage.first.first;
    size_t from = start;
    size_t to = usage.first.second;
    std::vector<Point*> loop;
    usage.second = true;
    loop.push_back(owned_points_[from]);

    while (to != start) {
      loop.push_back(owned_points_[to]);
      const std::vector<size_t>& candidates = outgoing[to];
      // Locate the reverse half-edge and take its circular predecessor.
      const auto reverse = std::find(candidates.begin(), candidates.end(),
        from);
      size_t reverse_index = static_cast<size_t>(reverse - candidates.begin());
      size_t best_index = reverse_index == 0
        ? candidates.size() - 1 : reverse_index - 1;
      const size_t best = candidates[best_index];
      from = to;
      to = best;
      auto& mark = used[std::make_pair(from, to)];
      if (mark) {
        break; // Defensive: malformed topology should not loop forever.
      }
      mark = true;
    }

    while (loop.size() > 1 && loop.front() == loop.back()) {
      loop.pop_back();
    }
    // A genuine planar face never revisits a node. A repeated node marks a
    // self-touching exterior walk produced at crossing vertices; those are
    // unions of real faces and must not be triangulated directly.
    bool simple = loop.size() >= 3;
    if (simple) {
      for (size_t i = 0; i < loop.size() && simple; ++i) {
        for (size_t j = i + 1; j < loop.size(); ++j) {
          if (loop[i] == loop[j]) {
            simple = false;
            break;
          }
        }
      }
    }
    if (simple && std::fpclassify(SignedArea(loop)) != FP_ZERO
        && InsideEvenOddFill(loop)) {
      loops_.push_back(loop);
    }
  }
}

bool EvenOddResolver::InsideEvenOddFill(const std::vector<Point*>& loop) const
{
  // Try several strictly interior samples; any odd-parity point proves the
  // face is filled under the even-odd rule. Multiple samples avoid ray
  // degeneracies where a sample happens to sit on an original contour edge.
  const size_t n = loop.size();
  const double third = 1.0 / 3.0;
  for (size_t i = 0; i < n; ++i) {
    const Point& a = *loop[i];
    const Point& b = *loop[(i + 1) % n];
    const Point& c = *loop[(i + 2) % n];
    Point sample(third * (a.x + b.x + c.x), third * (a.y + b.y + c.y));

    bool on_edge = false;
    for (size_t j = 0; j < contour_.size(); ++j) {
      if (OnSegment(sample, *contour_[j],
                    *contour_[(j + 1) % contour_.size()])) {
        on_edge = true;
        break;
      }
    }
    if (on_edge) {
      continue;
    }

    int crossings = 0;
    for (size_t j = 0; j < contour_.size(); ++j) {
      const Point& ea = *contour_[j];
      const Point& eb = *contour_[(j + 1) % contour_.size()];
      const bool straddles = (ea.y > sample.y) != (eb.y > sample.y);
      if (!straddles) {
        continue;
      }
      const double x = ea.x
        + (sample.y - ea.y) * (eb.x - ea.x) / (eb.y - ea.y);
      if (x > sample.x) {
        ++crossings;
      }
    }
    return (crossings % 2) == 1;
  }
  return false;
}

void EvenOddResolver::Resolve()
{
  SelfIntersectionInfo info = DetectSelfIntersection(contour_, 0);
  if (info.code == SI_DUPLICATE_POINT || info.code == SI_DEGENERATE_CONTOUR
      || info.code == SI_COLLINEAR_OVERLAP) {
    throw SelfIntersectionError(info);
  }
  if (info.code == SI_NONE) {
    loops_.push_back(contour_);
    return;
  }
  BuildSplitEdges();
  ExtractLoops();
  if (loops_.empty()) {
    SelfIntersectionInfo empty;
    empty.code = SI_DEGENERATE_CONTOUR;
    empty.edge_index_a = 0;
    empty.contour_index = 0;
    throw SelfIntersectionError(empty);
  }
}

double EvenOddResolver::FilledArea() const
{
  double total = 0.0;
  for (const auto& loop : loops_) {
    total += std::fabs(SignedArea(loop));
  }
  return total;
}

} // namespace p2t

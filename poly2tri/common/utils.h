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

// Otherwise #defines like M_PI are undeclared under Visual Studio
#define _USE_MATH_DEFINES

#include "shapes.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <utility>
#include <vector>

// C99 removes M_PI from math.h
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace p2t {

const double PI_3div4 = 3 * M_PI / 4;
const double PI_div2 = 1.57079632679489661923;
const double EPSILON = 1e-12;

enum Orientation { CW, CCW, COLLINEAR };

/**
 * Forumla to calculate signed area<br>
 * Positive if CCW<br>
 * Negative if CW<br>
 * 0 if collinear<br>
 * <pre>
 * A[P1,P2,P3]  =  (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + (x3*y1 - y3*x1)
 *              =  (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3)
 * </pre>
 */
inline Orientation Orient2d(const Point& pa, const Point& pb, const Point& pc)
{
  double detleft = (pa.x - pc.x) * (pb.y - pc.y);
  double detright = (pa.y - pc.y) * (pb.x - pc.x);
  double val = detleft - detright;

// Using a tolerance here fails on concave-by-subepsilon boundaries
//   if (val > -EPSILON && val < EPSILON) {
// Using == on double makes -Wfloat-equal warnings yell at us
  if (std::fpclassify(val) == FP_ZERO) {
    return COLLINEAR;
  } else if (val > 0) {
    return CCW;
  }
  return CW;
}

/*
bool InScanArea(Point& pa, Point& pb, Point& pc, Point& pd)
{
  double pdx = pd.x;
  double pdy = pd.y;
  double adx = pa.x - pdx;
  double ady = pa.y - pdy;
  double bdx = pb.x - pdx;
  double bdy = pb.y - pdy;

  double adxbdy = adx * bdy;
  double bdxady = bdx * ady;
  double oabd = adxbdy - bdxady;

  if (oabd <= EPSILON) {
    return false;
  }

  double cdx = pc.x - pdx;
  double cdy = pc.y - pdy;

  double cdxady = cdx * ady;
  double adxcdy = adx * cdy;
  double ocad = cdxady - adxcdy;

  if (ocad <= EPSILON) {
    return false;
  }

  return true;
}

*/

inline bool InScanArea(const Point& pa, const Point& pb, const Point& pc, const Point& pd)
{
  double oadb = (pa.x - pb.x)*(pd.y - pb.y) - (pd.x - pb.x)*(pa.y - pb.y);
  if (oadb >= -EPSILON) {
    return false;
  }

  double oadc = (pa.x - pc.x)*(pd.y - pc.y) - (pd.x - pc.x)*(pa.y - pc.y);
  if (oadc <= EPSILON) {
    return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// Triangulation quality report
//
// Pure post-processing pass over the triangles returned by CDT::GetTriangles.
// It never modifies the triangulation: with every metric disabled it is a
// no-op and the per-vertex output of the triangulation is left untouched.
// ---------------------------------------------------------------------------

/// Which statistics to gather, and how (the "statistics口径" is configurable).
struct QualityReportOptions {
  bool min_angle = true;      ///< per-triangle minimum interior angle (degrees)
  bool aspect_ratio = true;   ///< longest edge / (2*sqrt(3)*inradius); equilateral == 1
  bool area = true;           ///< triangle area distribution
  bool constraints = true;    ///< verify constrained edges appear in the mesh
  std::size_t histogram_bins = 10; ///< bins for the min-angle and area histograms
};

struct QualityReport {
  std::size_t triangle_count = 0;

  // Minimum interior angle per triangle, in degrees (always <= 60).
  double min_angle_deg = 0.0;      ///< smallest min-angle over all triangles
  double max_min_angle_deg = 0.0;  ///< largest min-angle over all triangles
  double mean_min_angle_deg = 0.0;
  std::vector<std::size_t> min_angle_hist; ///< histogram over [0, 60] degrees

  // Aspect ratio: longest_edge / (2*sqrt(3) * inradius). 1.0 == equilateral.
  double min_aspect_ratio = 0.0;
  double max_aspect_ratio = 0.0;
  double mean_aspect_ratio = 0.0;

  // Area distribution.
  double min_area = 0.0;
  double max_area = 0.0;
  double mean_area = 0.0;
  double total_area = 0.0;
  std::vector<std::size_t> area_hist; ///< histogram over [min_area, max_area]

  // Constrained-edge satisfaction.
  std::size_t constraints_total = 0;
  std::size_t constraints_satisfied = 0;

  // Worst triangle (the one with the smallest min-angle).
  std::size_t worst_triangle_index = 0;
  double worst_centroid_x = 0.0; ///< centroid of the worst triangle
  double worst_centroid_y = 0.0;
  double worst_min_angle_deg = 0.0;
};

namespace detail {

inline double TriangleArea(const Point& a, const Point& b, const Point& c)
{
  return 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

inline double EdgeLength(const Point& a, const Point& b)
{
  const double dx = a.x - b.x, dy = a.y - b.y;
  return std::sqrt(dx * dx + dy * dy);
}

/// Smallest interior angle of the triangle, in degrees.
inline double TriangleMinAngleDeg(const Point& a, const Point& b, const Point& c)
{
  const double la = EdgeLength(b, c); // opposite a
  const double lb = EdgeLength(a, c); // opposite b
  const double lc = EdgeLength(a, b); // opposite c
  // The smallest angle is opposite the shortest edge; use the law of cosines
  // on the two shorter edges for numerical robustness.
  double edges[3] = { la, lb, lc };
  std::sort(edges, edges + 3);
  const double s0 = edges[0], s1 = edges[1], s2 = edges[2];
  if (s0 <= 0.0) {
    return 0.0; // degenerate triangle
  }
  double cos_min = (s1 * s1 + s2 * s2 - s0 * s0) / (2.0 * s1 * s2);
  cos_min = std::max(-1.0, std::min(1.0, cos_min));
  return std::acos(cos_min) * 180.0 / M_PI;
}

} // namespace detail

/// Build a quality report for a triangulation.
///
/// @param triangles         triangles as returned by CDT::GetTriangles()
/// @param constrained_edges input constrained edges (boundary and hole
///                          polylines, as consecutive point pairs); only
///                          checked when options.constraints is true
/// @param options           which statistics to compute and their resolution
inline QualityReport MakeQualityReport(
    const std::vector<Triangle*>& triangles,
    const std::vector<std::pair<Point*, Point*>>& constrained_edges,
    const QualityReportOptions& options)
{
  QualityReport report;
  report.triangle_count = triangles.size();
  const std::size_t bins = std::max<std::size_t>(options.histogram_bins, 1);
  if (options.min_angle) {
    report.min_angle_hist.assign(bins, 0);
  }
  if (options.area) {
    report.area_hist.assign(bins, 0);
  }
  if (triangles.empty()) {
    report.constraints_total = constrained_edges.size();
    return report;
  }

  const double inf = std::numeric_limits<double>::infinity();
  double min_angle = inf, max_min_angle = -inf, sum_min_angle = 0.0;
  double min_aspect = inf, max_aspect = -inf, sum_aspect = 0.0;
  double min_area = inf, max_area = -inf, sum_area = 0.0;
  std::vector<double> areas(options.area ? triangles.size() : 0);

  for (std::size_t i = 0; i < triangles.size(); ++i) {
    const Point& a = *triangles[i]->GetPoint(0);
    const Point& b = *triangles[i]->GetPoint(1);
    const Point& c = *triangles[i]->GetPoint(2);

    double min_angle_deg = 0.0;
    if (options.min_angle || options.aspect_ratio) {
      min_angle_deg = detail::TriangleMinAngleDeg(a, b, c);
    }
    if (options.min_angle) {
      min_angle = std::min(min_angle, min_angle_deg);
      max_min_angle = std::max(max_min_angle, min_angle_deg);
      sum_min_angle += min_angle_deg;
      // Min-angle histogram over [0, 60] degrees.
      std::size_t bin = static_cast<std::size_t>(min_angle_deg / 60.0 * bins);
      bin = std::min(bin, bins - 1);
      ++report.min_angle_hist[bin];
    }
    if (i == 0 || min_angle_deg < report.worst_min_angle_deg) {
      report.worst_min_angle_deg = min_angle_deg;
      report.worst_triangle_index = i;
      report.worst_centroid_x = (a.x + b.x + c.x) / 3.0;
      report.worst_centroid_y = (a.y + b.y + c.y) / 3.0;
    }

    const double area = detail::TriangleArea(a, b, c);
    if (options.area) {
      areas[i] = area;
      min_area = std::min(min_area, area);
      max_area = std::max(max_area, area);
      sum_area += area;
    }

    if (options.aspect_ratio) {
      const double la = detail::EdgeLength(b, c);
      const double lb = detail::EdgeLength(a, c);
      const double lc = detail::EdgeLength(a, b);
      const double longest = std::max(la, std::max(lb, lc));
      const double semiperimeter = 0.5 * (la + lb + lc);
      const double inradius = semiperimeter > 0.0 ? area / semiperimeter : 0.0;
      const double aspect = inradius > 0.0
          ? longest / (2.0 * std::sqrt(3.0) * inradius)
          : inf;
      min_aspect = std::min(min_aspect, aspect);
      max_aspect = std::max(max_aspect, aspect);
      sum_aspect += aspect;
    }
  }

  const double n = static_cast<double>(triangles.size());
  if (options.min_angle) {
    report.min_angle_deg = min_angle;
    report.max_min_angle_deg = max_min_angle;
    report.mean_min_angle_deg = sum_min_angle / n;
  }
  if (options.aspect_ratio) {
    report.min_aspect_ratio = min_aspect;
    report.max_aspect_ratio = max_aspect;
    report.mean_aspect_ratio = sum_aspect / n;
  }
  if (options.area) {
    report.min_area = min_area;
    report.max_area = max_area;
    report.mean_area = sum_area / n;
    report.total_area = sum_area;
    const double range = max_area - min_area;
    for (std::size_t i = 0; i < areas.size(); ++i) {
      std::size_t bin = range > 0.0
          ? static_cast<std::size_t>((areas[i] - min_area) / range * bins)
          : 0;
      bin = std::min(bin, bins - 1);
      ++report.area_hist[bin];
    }
  }

  if (options.constraints) {
    report.constraints_total = constrained_edges.size();
    for (const auto& edge : constrained_edges) {
      for (auto* t : triangles) {
        if (t->Contains(edge.first) && t->Contains(edge.second)) {
          ++report.constraints_satisfied;
          break;
        }
      }
    }
  }
  return report;
}

}

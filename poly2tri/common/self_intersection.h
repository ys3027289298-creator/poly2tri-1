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

#pragma once

#include "shapes.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace p2t {

/// What should happen when a contour proves not to be a simple polygon.
enum SelfIntersectionPolicy {
  /// Reject the input with a SelfIntersectionError (default behaviour).
  REJECT_SELF_INTERSECTION = 0,
  /// Rewrite the contour into simple even-odd fill loops, then triangulate
  /// those. The resulting filled area equals the winding-independent
  /// (even-odd) area of the input within a floating point tolerance.
  RESOLVE_EVEN_ODD = 1
};

/// Machine readable reason for rejecting a contour before the sweep runs.
enum SelfIntersectionCode {
  SI_NONE = 0,
  /// Two non-adjacent edges cross properly in their interiors.
  SI_PROPER_CROSSING = 1,
  /// A contour vertex lies in the interior of a non-adjacent edge
  /// (a "T-junction"); poly2tri cannot express this constraint directly.
  SI_VERTEX_ON_EDGE = 2,
  /// Two edges overlap along a collinear segment.
  SI_COLLINEAR_OVERLAP = 3,
  /// The same coordinate occurs more than once inside one contour.
  SI_DUPLICATE_POINT = 4,
  /// A contour has fewer than three points or collapses to a line.
  SI_DEGENERATE_CONTOUR = 5
};

struct SelfIntersectionInfo {
  SelfIntersectionCode code;
  /// Indices of the two offending edges, local to their contour. For
  /// SI_DUPLICATE_POINT and SI_DEGENERATE_CONTOUR the second value is -1.
  int edge_index_a;
  int edge_index_b;
  /// Zero-based index of the offending contour (0 is the outer polyline).
  int contour_index;

  SelfIntersectionInfo()
  : code(SI_NONE), edge_index_a(-1), edge_index_b(-1), contour_index(-1)
  {
  }
};

/// Thrown by CDT::Triangulate when a contour is not a simple polygon and
/// the active policy is REJECT_SELF_INTERSECTION (or the defect cannot be
/// resolved automatically).
class SelfIntersectionError : public std::runtime_error {
public:
  explicit SelfIntersectionError(const SelfIntersectionInfo& info)
  : std::runtime_error(FormatMessage(info)), info_(info)
  {
  }

  const SelfIntersectionInfo& info() const
  {
    return info_;
  }

  SelfIntersectionCode code() const
  {
    return info_.code;
  }

private:
  static std::string FormatMessage(const SelfIntersectionInfo& info);
  SelfIntersectionInfo info_;
};

/// Checks one closed contour (outer polyline or hole) for self intersection.
/// Returns SI_NONE for a simple polygon, otherwise a description of the
/// first defect found. Adjacent edges are allowed to share their common
/// endpoint, which is the normal polygon representation.
P2T_DLL_SYMBOL SelfIntersectionInfo
DetectSelfIntersection(const std::vector<Point*>& contour, int contour_index = 0);

/// Convenience wrapper for the outer polyline plus every hole.
P2T_DLL_SYMBOL SelfIntersectionInfo
DetectSelfIntersections(const std::vector<std::vector<Point*> >& contours);

/// Signed (shoelace) area of a closed contour; positive for CCW input.
P2T_DLL_SYMBOL double SignedArea(const std::vector<Point*>& contour);

/// Rewrites a self-intersecting contour into simple (non-self-intersecting)
/// loops following the even-odd fill rule. Every intersection vertex is
/// inserted into the affected edges, then directed half-edges are paired
/// into closed loops. All returned points and their ownership stay inside
/// the resolver instance; the loops reference those internal points.
class P2T_DLL_SYMBOL EvenOddResolver {
public:
  explicit EvenOddResolver(const std::vector<Point*>& contour);
  ~EvenOddResolver();

  /// Runs the split and loop extraction. Throws SelfIntersectionError when
  /// the defect cannot be expressed under the even-odd rule, for example
  /// collinear overlaps or duplicate points.
  void Resolve();

  /// Simple loops after Resolve(). Outer loops are CCW, hole loops CW.
  const std::vector<std::vector<Point*> >& loops() const
  {
    return loops_;
  }

  /// Sum of absolute signed loop areas, i.e. the even-odd filled area.
  double FilledArea() const;

private:
  EvenOddResolver(const EvenOddResolver&);
  EvenOddResolver& operator=(const EvenOddResolver&);

  size_t FindOrAddPoint(const Point& point);
  void BuildSplitEdges();
  void ExtractLoops();
  bool InsideEvenOddFill(const std::vector<Point*>& loop) const;

  const std::vector<Point*>& contour_;
  std::vector<Point*> owned_points_;
  std::vector<std::pair<size_t, size_t> > split_edges_;
  std::vector<std::vector<Point*> > loops_;
};

} // namespace p2t

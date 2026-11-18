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
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include "dll_symbol.h"
#include "shapes.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace p2t {

/// Thrown by ExtractBoundaryLoops when an edge is shared by more than two
/// triangles, or by two triangles that wind it in the same direction.
struct P2T_DLL_SYMBOL NonManifoldEdgeError : public std::runtime_error {

  /// Index of the first half-edge (3 * triangle + corner) belonging to the
  /// non-manifold edge.
  std::size_t half_edge;
  /// Index of the triangle in the input array: half_edge / 3.
  std::size_t triangle;
  /// Corner of the triangle the half-edge starts at: half_edge % 3.
  int corner;

  explicit NonManifoldEdgeError(std::size_t half_edge);
};

/// One closed boundary loop of a triangle mesh.
struct P2T_DLL_SYMBOL BoundaryLoop {

  /// Ordered loop vertices (open: the last point connects back to the first).
  /// Outer loops wind counter-clockwise, holes clockwise.
  std::vector<Point*> points;
  /// Shoelace area: positive for outer loops, negative for holes.
  double signed_area;
  /// True when signed_area < 0 (clockwise loop).
  bool is_hole;
  /// Index of the smallest containing loop in BoundaryLoops::loops, -1 if none.
  int parent;
  /// Nesting depth: 0 = outer boundary, 1 = hole, 2 = island inside a hole, ...
  int depth;
};

/// Result of ExtractBoundaryLoops.
struct P2T_DLL_SYMBOL BoundaryLoops {

  std::vector<BoundaryLoop> loops;

  /// Indices of the depth-0 (outer boundary) loops.
  std::vector<std::size_t> OuterLoops() const;
};

/**
 * Extracts the boundary loops of a consistently oriented triangle mesh.
 *
 * Builds a half-edge structure over the input triangles, chains boundary
 * half-edges (those without a twin) into closed loops, classifies each loop
 * as outer boundary or hole by its signed area, and finally derives the
 * nesting hierarchy (parent/depth) from point-in-polygon containment tests.
 *
 * Throws NonManifoldEdgeError if an edge is used by more than two triangles
 * or by two triangles with the same winding.
 */
P2T_DLL_SYMBOL BoundaryLoops
ExtractBoundaryLoops(const std::vector<Triangle*>& triangles);

}

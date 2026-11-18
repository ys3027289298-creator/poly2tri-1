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

#include "shapes.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

// Winding order normalization for triangulated meshes.
//
// Downstream renderers derive face normals from vertex order, so a mesh whose
// triangles mix clockwise and counter-clockwise winding renders with inverted
// (black) faces. The helpers here follow a three step pipeline:
//
//   1. AnalyzeWinding: estimate the overall orientation from the sum of the
//      per-face signed double areas (a positive sum means the mesh as a whole
//      winds CCW) and count the CW/CCW/degenerate faces.
//   2. NormalizeWinding: flip each individual triangle so every face is CCW
//      (or CW when requested). Flipping only swaps indices, so face areas are
//      preserved bit for bit.
//   3. VerifyWinding: walk every directed edge and check that neighboring
//      faces traverse their shared edge in opposite directions, which is
//      exactly the condition that the winding is consistent across the whole
//      mesh. Boundary edges are traversed once and are not checked.
//
// Every routine is linear in the number of faces (the edge verification uses
// a hash table with one entry per unique edge).

namespace p2t {

// Orientation of a face or of the mesh as a whole.
enum class Winding {
  CCW,
  CW,
  Degenerate,
};

// Outcome of AnalyzeWinding / VerifyWinding. All area quantities are *double*
// areas: the signed area of a triangle is half of the signed double area.
struct WindingStats {
  std::size_t face_count = 0;
  std::size_t ccw_faces = 0;
  std::size_t cw_faces = 0;
  // Zero-area faces, which carry no orientation and are never flipped.
  std::size_t degenerate_faces = 0;
  // Sum of the per-face signed double areas; twice the signed polygon area.
  double signed_double_area = 0.0;
  // Sum of the absolute per-face double areas; twice the polygon area.
  double unsigned_double_area = 0.0;
  // Edges shared by exactly two faces traversed in the same direction.
  std::size_t inconsistent_edges = 0;
  // Edges shared by more than two faces (not a 2-manifold mesh).
  std::size_t non_manifold_edges = 0;
  // Edges belonging to a single face.
  std::size_t boundary_edges = 0;

  // Orientation of the mesh as a whole, estimated from the signed area sum.
  Winding OverallOrientation() const
  {
    if (signed_double_area > 0.0) {
      return Winding::CCW;
    }
    if (signed_double_area < 0.0) {
      return Winding::CW;
    }
    return Winding::Degenerate;
  }

  // True when every interior edge is traversed in opposite directions and no
  // edge is shared by more than two faces.
  bool IsManifoldConsistent() const
  {
    return inconsistent_edges == 0 && non_manifold_edges == 0;
  }
};

// Signed double area of a single triangle: positive for CCW winding, negative
// for CW, zero for a degenerate (collinear or repeated) triangle.
inline double SignedDoubleArea(const Point& a, const Point& b, const Point& c)
{
  return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

template<class Index>
inline double SignedDoubleArea(const std::vector<Point>& points, const Index* tri)
{
  return SignedDoubleArea(points[tri[0]], points[tri[1]], points[tri[2]]);
}

// Reverse a triangle's vertex order in place.
template<class Index>
inline void FlipTriangle(Index* tri)
{
  Index tmp = tri[1];
  tri[1] = tri[2];
  tri[2] = tmp;
}

// Step 1: estimate the overall orientation and tally face windings. The mesh
// is not modified. The edge counters stay zero here; use VerifyWinding for the
// shared-edge check.
template<class Index>
WindingStats AnalyzeWinding(const std::vector<Index>& triangles,
                            const std::vector<Point>& points)
{
  WindingStats stats;
  stats.face_count = triangles.size() / 3;

  for (std::size_t f = 0; f < stats.face_count; ++f) {
    const Index* tri = &triangles[3 * f];
    const double area2 = SignedDoubleArea(points, tri);
    stats.signed_double_area += area2;
    stats.unsigned_double_area += area2 >= 0.0 ? area2 : -area2;
    if (area2 > 0.0) {
      ++stats.ccw_faces;
    } else if (area2 < 0.0) {
      ++stats.cw_faces;
    } else {
      ++stats.degenerate_faces;
    }
  }
  return stats;
}

// Step 2: rewrite every triangle so its winding matches target. Degenerate
// triangles are left alone. Returns the number of triangles actually flipped.
//
// Vertex indices must be non-negative and fit into 32 bits, which is also the
// limit imposed by the edge hash used in VerifyWinding.
template<class Index>
std::size_t NormalizeWinding(std::vector<Index>& triangles,
                             const std::vector<Point>& points,
                             Winding target = Winding::CCW)
{
  std::size_t flipped = 0;
  const std::size_t face_count = triangles.size() / 3;

  for (std::size_t f = 0; f < face_count; ++f) {
    Index* tri = &triangles[3 * f];
    const double area2 = SignedDoubleArea(points, tri);
    if (area2 == 0.0) {
      continue;
    }
    const bool face_is_ccw = area2 > 0.0;
    if ((target == Winding::CCW && !face_is_ccw) ||
        (target == Winding::CW && face_is_ccw)) {
      FlipTriangle(tri);
      ++flipped;
    }
  }
  return flipped;
}

namespace winding_detail {

// Deterministic 64 bit finalizer (SplitMix64). The standard library's default
// hash for integers is the identity on major implementations, which leaves the
// packed edge keys clustered in the low buckets and makes the edge table's
// runtime depend on the platform allocator. Mixing the bits gives stable,
// cache-friendly behavior across standard libraries and build configurations.
struct SplitMix64Hash {
  std::size_t operator()(std::uint64_t x) const
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return static_cast<std::size_t>(x ^ (x >> 31));
  }
};

// Packs an undirected edge into a 64 bit key with the smaller index first, so
// that both (a, b) and (b, a) map to the same entry.
inline std::uint64_t EdgeKey(std::uint32_t a, std::uint32_t b)
{
  if (a > b) {
    const std::uint32_t tmp = a;
    a = b;
    b = tmp;
  }
  return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
}

// Occurrences of one undirected edge across the triangles.
struct EdgeInfo {
  // How many triangles contain this edge.
  std::uint32_t count = 0;
  // Direction of the first occurrence along the canonical (low, high) key:
  // 1 for low -> high, 0 for high -> low.
  std::uint32_t first_direction = 0;
  // Once two occurrences have been seen: 1 when both run the same way,
  // 0 when they oppose each other (the manifold-consistent case).
  std::uint32_t same_direction = 0;
};

// Records one directed edge (a -> b) into the edge table.
inline void RecordEdge(std::unordered_map<std::uint64_t, EdgeInfo, SplitMix64Hash>& edges,
                       std::uint32_t a,
                       std::uint32_t b)
{
  const std::uint64_t key = EdgeKey(a, b);
  const std::uint32_t direction = (a < b) ? 1u : 0u;
  EdgeInfo& info = edges[key];

  if (info.count == 0) {
    info.first_direction = direction;
    info.count = 1;
  } else if (info.count == 1) {
    info.same_direction = (direction == info.first_direction) ? 1u : 0u;
    info.count = 2;
  } else {
    // Third or later occurrence: the edge is non-manifold. The direction
    // fields stay as recorded for the first two faces.
    ++info.count;
  }
}

} // namespace winding_detail

// Step 3: walk every directed edge and check that interior edges are
// traversed in opposite directions by their two incident faces. Also counts
// boundary and non-manifold edges. Does not modify the mesh; run it after
// NormalizeWinding, or on any mesh to diagnose its winding.
template<class Index>
WindingStats VerifyWinding(const std::vector<Index>& triangles,
                           const std::vector<Point>& points)
{
  using namespace winding_detail;

  WindingStats stats;
  stats.face_count = triangles.size() / 3;
  std::unordered_map<std::uint64_t, EdgeInfo, SplitMix64Hash> edges;
  // Keep the load factor low enough that inserts stay O(1) without rehashing
  // on every growth step; the table holds at most 3 entries per face.
  edges.max_load_factor(0.7f);
  edges.reserve(stats.face_count * 2);

  for (std::size_t f = 0; f < stats.face_count; ++f) {
    const Index* tri = &triangles[3 * f];
    const double area2 = SignedDoubleArea(points, tri);
    stats.signed_double_area += area2;
    stats.unsigned_double_area += area2 >= 0.0 ? area2 : -area2;
    if (area2 > 0.0) {
      ++stats.ccw_faces;
    } else if (area2 < 0.0) {
      ++stats.cw_faces;
    } else {
      ++stats.degenerate_faces;
    }
    // Degenerate triangles can contain self-loops (a == b); they are not real
    // edges and would poison the shared-edge statistics, so skip them.
    for (int e = 0; e < 3; ++e) {
      const std::uint32_t a = static_cast<std::uint32_t>(tri[e]);
      const std::uint32_t b = static_cast<std::uint32_t>(tri[(e + 1) % 3]);
      if (a != b) {
        RecordEdge(edges, a, b);
      }
    }
  }

  for (const auto& entry : edges) {
    const EdgeInfo& info = entry.second;
    if (info.count == 1) {
      ++stats.boundary_edges;
    } else if (info.count == 2) {
      // Opposite directions are consistent; equal directions mean the two
      // faces wind the same way and one of them must be flipped.
      if (info.same_direction == 1) {
        ++stats.inconsistent_edges;
      }
    } else {
      ++stats.non_manifold_edges;
    }
  }
  return stats;
}

// Convenience: run the full analyze / normalize / verify pipeline and return
// the post-normalization statistics. A correctly oriented manifold mesh
// reports IsManifoldConsistent() == true.
template<class Index>
WindingStats NormalizeAndVerifyWinding(std::vector<Index>& triangles,
                                       const std::vector<Point>& points,
                                       Winding target = Winding::CCW)
{
  NormalizeWinding(triangles, points, target);
  return VerifyWinding(triangles, points);
}

// Build an indexed mesh from a CDT triangle list. The Point objects keep their
// identity, so shared vertices are merged into one index via pointer equality.
// The output can be fed straight into NormalizeWinding before export.
inline void ExtractTriangles(const std::vector<Triangle*>& tris,
                             std::vector<Point>& points,
                             std::vector<int>& triangles)
{
  points.clear();
  triangles.clear();
  triangles.reserve(tris.size() * 3);

  std::unordered_map<const Point*, int> index;
  index.reserve(tris.size() * 2);

  for (Triangle* tri : tris) {
    for (int i = 0; i < 3; ++i) {
      const Point* p = tri->GetPoint(i);
      auto found = index.find(p);
      if (found == index.end()) {
        const int id = static_cast<int>(points.size());
        index.emplace(p, id);
        points.push_back(*p);
        triangles.push_back(id);
      } else {
        triangles.push_back(found->second);
      }
    }
  }
}

} // namespace p2t

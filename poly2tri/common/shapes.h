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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iosfwd>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
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

// ===========================================================================
// Triangle mesh post-processing pipeline
//
// Stage 1 (weldVertices): merge corners whose coordinates are within
//   `tolerance` of an already accepted vertex (uniform grid hashing).
// Stage 2 (reorder): greedily grow triangle strips, then triangle fans, out
//   of the shared-edge adjacency; everything left over is emitted as plain
//   3-index segments, so no triangle is ever lost or duplicated.
// Stage 3 (emitIndices): concatenate every segment into one index buffer and
//   record one entry per segment (kind, offset, length) in the segment table.
//
// Every stage can be switched off independently through TrianglePipelineOptions.
// All temporary storage is linear in the number of input corners / unique
// vertices; a thread-local allocation tracker (MemoryTracker) reports the
// measured peak byte usage without any third-party dependency.
// ===========================================================================

namespace mesh_detail {

// Thread-local allocation accounting used by TrackedAllocator.
class MemoryTracker
{
public:
  static MemoryTracker& Instance()
  {
    static thread_local MemoryTracker tracker;
    return tracker;
  }

  void OnAllocate(std::size_t bytes)
  {
    live_bytes_ += bytes;
    if (live_bytes_ > peak_bytes_) {
      peak_bytes_ = live_bytes_;
    }
    allocation_count_ += 1;
  }

  void OnDeallocate(std::size_t bytes) { live_bytes_ -= bytes; }

  std::size_t PeakBytes() const { return peak_bytes_; }
  std::size_t LiveBytes() const { return live_bytes_; }
  std::size_t AllocationCount() const { return allocation_count_; }

  void Reset()
  {
    live_bytes_ = 0;
    peak_bytes_ = 0;
    allocation_count_ = 0;
  }

private:
  MemoryTracker() = default;
  std::size_t live_bytes_ = 0;
  std::size_t peak_bytes_ = 0;
  std::size_t allocation_count_ = 0;
};

// Allocator adapter that forwards every (de)allocation to MemoryTracker.
template <typename T>
struct TrackedAllocator
{
  using value_type = T;

  TrackedAllocator() noexcept = default;
  template <typename U>
  TrackedAllocator(const TrackedAllocator<U>&) noexcept
  {
  }

  T* allocate(std::size_t n)
  {
    // std::allocator keeps the standard overflow / bad_alloc semantics.
    T* p = std::allocator<T>().allocate(n);
    MemoryTracker::Instance().OnAllocate(n * sizeof(T));
    return p;
  }

  void deallocate(T* p, std::size_t n)
  {
    MemoryTracker::Instance().OnDeallocate(n * sizeof(T));
    std::allocator<T>().deallocate(p, n);
  }

  template <typename U>
  bool operator==(const TrackedAllocator<U>&) const noexcept
  {
    return true;
  }
  template <typename U>
  bool operator!=(const TrackedAllocator<U>&) const noexcept
  {
    return false;
  }
};

template <typename T>
using TrackedVector = std::vector<T, TrackedAllocator<T>>;

struct CellKey
{
  std::int64_t cx;
  std::int64_t cy;

  bool operator==(const CellKey& other) const
  {
    return cx == other.cx && cy == other.cy;
  }
};

struct CellKeyHash
{
  std::size_t operator()(const CellKey& key) const
  {
    // SplitMix-style mixing of the two 64-bit cell coordinates.
    std::uint64_t h = static_cast<std::uint64_t>(key.cx) * 0x9e3779b97f4a7c15ULL;
    h ^= static_cast<std::uint64_t>(key.cy) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
    h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
    return static_cast<std::size_t>(h ^ (h >> 31));
  }
};

struct PairHash
{
  std::size_t operator()(const std::pair<double, double>& p) const
  {
    std::uint64_t hx = 0;
    std::uint64_t hy = 0;
    std::memcpy(&hx, &p.first, sizeof(double));
    std::memcpy(&hy, &p.second, sizeof(double));
    std::uint64_t h = hx * 0x9e3779b97f4a7c15ULL;
    h ^= hy + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
    return static_cast<std::size_t>(h ^ (h >> 27));
  }
};

struct PairEqual
{
  bool operator()(const std::pair<double, double>& a,
                  const std::pair<double, double>& b) const
  {
    return a.first == b.first && a.second == b.second;
  }
};

using CellMap = std::unordered_map<CellKey, TrackedVector<std::uint32_t>, CellKeyHash,
                                   std::equal_to<CellKey>,
                                   TrackedAllocator<std::pair<const CellKey,
                                                             TrackedVector<std::uint32_t>>>>;

using ExactMap =
  std::unordered_map<std::pair<double, double>, std::uint32_t, PairHash, PairEqual,
                     TrackedAllocator<std::pair<const std::pair<double, double>,
                                                 std::uint32_t>>>;

using EdgeMap = std::unordered_map<std::uint64_t, TrackedVector<std::uint32_t>,
                                   std::hash<std::uint64_t>, std::equal_to<std::uint64_t>,
                                   TrackedAllocator<std::pair<const std::uint64_t,
                                                               TrackedVector<std::uint32_t>>>>;

inline std::uint64_t EdgeKey(std::uint32_t a, std::uint32_t b)
{
  if (a > b) {
    std::swap(a, b);
  }
  return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
}

inline std::int64_t QuantizeCell(double v, double inv_tol)
{
  return static_cast<std::int64_t>(std::floor(v * inv_tol));
}

} // namespace mesh_detail

/// Kind of an emitted index segment.
enum class SegmentKind
{
  Triangle, ///< Plain 3-index triangle (also used for degenerate input).
  Strip,    ///< Triangle strip: n+2 indices decode n triangles, alternating
            ///< winding.
  Fan       ///< Triangle fan: hub index first, then n+1 rim indices for n
            ///< triangles.
};

/// One entry of the per-segment length table.
struct Segment
{
  SegmentKind kind;
  std::uint32_t offset; ///< First index of this segment in the index buffer.
  std::uint32_t length; ///< Number of indices belonging to this segment.
};

/// Independent on/off switches for the three pipeline stages.
struct TrianglePipelineOptions
{
  /// Stage 1: merge vertices closer than `tolerance` (L-infinity distance).
  bool weldVertices = true;
  /// Stage 2: reorder triangles into strips and fans.
  bool reorder = true;
  /// Stage 3: materialize the concatenated index buffer and segment table.
  bool emitIndices = true;
  /// Stage 2 ordering policy. When false (default) strips are grown first,
  /// which compresses grids and long bands best. When true, edge-continuous
  /// triangle groups around a common hub are emitted as fans first, which
  /// compresses radial meshes best. Both policies always keep both encodings
  /// available and reconstruct the identical triangle set.
  bool preferFans = false;
  /// Coordinate tolerance of stage 1. Must be >= 0; 0 means exact welding.
  double tolerance = 0.0;
};

/// Measured statistics of a pipeline run.
struct PipelineMetrics
{
  std::size_t input_vertex_count = 0;   ///< Corners handed in (3 per triangle).
  std::size_t output_vertex_count = 0;  ///< Unique vertices after stage 1.
  std::size_t merged_vertex_count = 0;  ///< input_vertex_count - output.
  std::size_t triangle_count = 0;
  std::size_t strip_count = 0;
  std::size_t fan_count = 0;
  std::size_t plain_segment_count = 0;
  std::size_t degenerate_count = 0;
  std::size_t index_count = 0;          ///< Total emitted indices.
  /// Saved corners relative to independent triangles: 1 - index_count/(3T).
  double vertex_reuse_rate = 0.0;
  /// Every input triangle is decoded from the output exactly once.
  bool reconstruction_consistent = true;
  std::size_t reconstruction_mismatches = 0;
  /// Peak working-memory bytes measured by the tracked allocator. This covers
  /// every internal container of the three stages (weld hash grid, edge map,
  /// segment/index working buffers); it excludes the returned PipelineResult,
  /// whose size is itself O(V) and is observable directly by the caller.
  std::size_t peak_memory_bytes = 0;
};

/// Output of TriangleMeshPipeline::Build.
struct PipelineResult
{
  std::vector<Point> vertices; ///< Welded vertex table (stage 1).
  std::vector<std::uint32_t> indices; ///< Concatenated index buffer (stage 3).
  std::vector<Segment> segments;      ///< Per-segment length/kind table.
  /// Convenience copy of segment lengths, one entry per segment.
  std::vector<std::uint32_t> segment_lengths;
  PipelineMetrics metrics;
};

///
/// Three-stage triangle-soup post-processing pipeline. The class is stateless;
/// each Build() call is independent and re-entrant.
///
class TriangleMeshPipeline
{
public:
  /// Build from a flat corner array: [x0,y0,x1,y1,...], 3 corners per
  /// triangle. Throws std::invalid_argument on a negative tolerance or a
  /// corner array whose size is not a multiple of 6, and
  /// std::overflow_error if the mesh needs more than 2^32-1 vertices.
  static PipelineResult Build(const std::vector<double>& flat_corners,
                              const TrianglePipelineOptions& options = {})
  {
    if (flat_corners.size() % 6 != 0) {
      throw std::invalid_argument(
        "TriangleMeshPipeline: flat corner array size must be a multiple of 6");
    }
    const std::size_t triangle_count = flat_corners.size() / 6;
    mesh_detail::TrackedVector<Point> corners;
    corners.reserve(triangle_count * 3);
    for (std::size_t i = 0; i + 1 < flat_corners.size(); i += 2) {
      corners.emplace_back(flat_corners[i], flat_corners[i + 1]);
    }
    return Run(corners, options);
  }

  /// Convenience overload taking the triangulation output of a p2t::CDT.
  /// Reads only coordinates; CDT ownership is untouched.
  static PipelineResult Build(const std::vector<Triangle*>& triangles,
                              const TrianglePipelineOptions& options = {})
  {
    mesh_detail::TrackedVector<Point> corners;
    corners.reserve(triangles.size() * 3);
    for (Triangle* t : triangles) {
      for (int i = 0; i < 3; ++i) {
        const Point* p = t->GetPoint(i);
        corners.emplace_back(p->x, p->y);
      }
    }
    return Run(corners, options);
  }

private:
  using U32 = std::uint32_t;

  struct State
  {
    explicit State(std::size_t corner_count)
    {
      triangles.reserve(corner_count); // 3 ids per triangle
      used.reserve(corner_count / 3);
    }
    mesh_detail::TrackedVector<Point> vertices;
    mesh_detail::TrackedVector<U32> triangles; // welded ids, 3 per triangle
    mesh_detail::TrackedVector<unsigned char> used;
    mesh_detail::TrackedVector<unsigned char> degenerate;
    mesh_detail::TrackedVector<Segment> segments;
    mesh_detail::TrackedVector<U32> indices;
  };

  static PipelineResult Run(const mesh_detail::TrackedVector<Point>& corners,
                            const TrianglePipelineOptions& options)
  {
    if (options.tolerance < 0.0) {
      throw std::invalid_argument(
        "TriangleMeshPipeline: tolerance must be non-negative");
    }

    PipelineResult result;
    auto& tracker = mesh_detail::MemoryTracker::Instance();
    tracker.Reset();

    const std::size_t corner_count = corners.size();
    const std::size_t triangle_count = corner_count / 3;
    if (corner_count == 0) {
      // Explicit empty-input handling: valid, empty output, zero memory.
      result.metrics.reconstruction_consistent = true;
      result.metrics.peak_memory_bytes = 0;
      return result;
    }
    if (corner_count > (std::numeric_limits<U32>::max)()) {
      throw std::overflow_error(
        "TriangleMeshPipeline: mesh exceeds the 32-bit index range");
    }

    State state(corner_count);
    WeldVertices(corners, options, state);
    state.used.assign(triangle_count, 0);
    state.degenerate.assign(triangle_count, 0);
    for (std::size_t t = 0; t < triangle_count; ++t) {
      const U32* ids = &state.triangles[t * 3];
      if (ids[0] == ids[1] || ids[1] == ids[2] || ids[0] == ids[2]) {
        state.degenerate[t] = 1; // repeated vertex after welding
        continue;
      }
      // Collinear (zero-area) triangle: cross product of two edge vectors.
      const Point& p0 = state.vertices[ids[0]];
      const Point& p1 = state.vertices[ids[1]];
      const Point& p2 = state.vertices[ids[2]];
      const double cross = (p1.x - p0.x) * (p2.y - p0.y) -
                           (p2.x - p0.x) * (p1.y - p0.y);
      if (cross == 0.0) {
        state.degenerate[t] = 1;
      }
    }

    if (options.reorder) {
      if (options.preferFans) {
        // Fan-first policy for radial meshes; strips take the leftovers.
        ExtractFans(state);
        ExtractStrips(state);
      } else {
        // Strip-first policy: best for regular grids, bands and holed meshes;
        // the fan stage groups any remaining edge-continuous hub clusters.
        ExtractStrips(state);
        ExtractFans(state);
      }
    }
    EmitPlain(state); // remaining triangles, plus degenerates if reordered

    if (options.emitIndices) {
      MaterializeIndices(state);
    }

    FillResult(corners, options, state, result);
    result.metrics.peak_memory_bytes = tracker.PeakBytes();
    return result;
  }

  // ---------------------------------------------------------------- Stage 1
  static void WeldVertices(const mesh_detail::TrackedVector<Point>& corners,
                           const TrianglePipelineOptions& options, State& state)
  {
    const std::size_t corner_count = corners.size();
    state.vertices.reserve(corner_count);
    state.triangles.reserve(corner_count);

    if (!options.weldVertices) {
      // Every corner keeps its own identity; indices are simply 0..N-1.
      for (std::size_t i = 0; i < corner_count; ++i) {
        state.vertices.push_back(corners[i]);
        state.triangles.push_back(static_cast<U32>(i));
      }
      return;
    }

    if (options.tolerance == 0.0) {
      // Exact welding: bit-identical coordinate pair -> same vertex.
      mesh_detail::ExactMap exact;
      exact.reserve(corner_count * 2);
      for (std::size_t i = 0; i < corner_count; ++i) {
        const Point& p = corners[i];
        const std::pair<double, double> key(p.x, p.y);
        auto found = exact.find(key);
        if (found == exact.end()) {
          const U32 id = static_cast<U32>(state.vertices.size());
          exact.emplace(key, id);
          state.vertices.push_back(p);
          state.triangles.push_back(id);
        } else {
          state.triangles.push_back(found->second);
        }
      }
      return;
    }

    // Tolerance welding on a uniform grid whose cell size equals the
    // tolerance. Any point within tolerance lives in one of the 3x3 cells
    // around the lookup cell, so candidates are found in expected O(1).
    const double inv_tol = 1.0 / options.tolerance;
    mesh_detail::CellMap cells;
    cells.reserve(corner_count * 2);
    for (std::size_t i = 0; i < corner_count; ++i) {
      const Point& p = corners[i];
      const mesh_detail::CellKey home{ mesh_detail::QuantizeCell(p.x, inv_tol),
                                       mesh_detail::QuantizeCell(p.y, inv_tol) };
      U32 matched = static_cast<U32>(state.vertices.size());
      for (std::int64_t dx = -1; dx <= 1 && matched == state.vertices.size(); ++dx) {
        for (std::int64_t dy = -1; dy <= 1; ++dy) {
          const mesh_detail::CellKey q{ home.cx + dx, home.cy + dy };
          auto bucket = cells.find(q);
          if (bucket == cells.end()) {
            continue;
          }
          for (U32 candidate : bucket->second) {
            const Point& qp = state.vertices[candidate];
            if (std::fabs(qp.x - p.x) <= options.tolerance &&
                std::fabs(qp.y - p.y) <= options.tolerance) {
              matched = candidate;
              break;
            }
          }
          if (matched != state.vertices.size()) {
            break;
          }
        }
      }
      if (matched == state.vertices.size()) {
        matched = static_cast<U32>(state.vertices.size());
        state.vertices.push_back(p);
        cells[home].push_back(matched);
      }
      state.triangles.push_back(matched);
    }
  }

  // ------------------------------------------------- Stage 2: triangle strips
  static mesh_detail::EdgeMap BuildEdgeMap(const State& state)
  {
    const std::size_t triangle_count = state.used.size();
    mesh_detail::EdgeMap edges;
    edges.reserve(triangle_count * 3);
    for (std::size_t t = 0; t < triangle_count; ++t) {
      if (state.used[t] || state.degenerate[t]) {
        continue;
      }
      const U32* ids = &state.triangles[t * 3];
      for (int e = 0; e < 3; ++e) {
        const U32 a = ids[e];
        const U32 b = ids[(e + 1) % 3];
        if (a != b) {
          edges[mesh_detail::EdgeKey(a, b)].push_back(static_cast<U32>(t));
        }
      }
    }
    return edges;
  }

  static U32 PickUnusedNeighbor(const mesh_detail::EdgeMap& edges, U32 a, U32 b,
                                const mesh_detail::TrackedVector<unsigned char>& used)
  {
    auto found = edges.find(mesh_detail::EdgeKey(a, b));
    if (found == edges.end()) {
      return static_cast<U32>(used.size()); // "none" sentinel
    }
    for (U32 t : found->second) {
      if (!used[t]) {
        return t;
      }
    }
    return static_cast<U32>(used.size());
  }

  static U32 ThirdVertex(const State& state, U32 triangle, U32 a, U32 b)
  {
    const U32* ids = &state.triangles[triangle * 3];
    for (int i = 0; i < 3; ++i) {
      if (ids[i] != a && ids[i] != b) {
        return ids[i];
      }
    }
    return static_cast<U32>(state.vertices.size()); // cannot happen for valid edge
  }

  static void ExtractStrips(State& state)
  {
    mesh_detail::EdgeMap edges = BuildEdgeMap(state);
    const std::size_t triangle_count = state.used.size();
    const U32 none = static_cast<U32>(triangle_count);

    for (std::size_t seed = 0; seed < triangle_count; ++seed) {
      if (state.used[seed] || state.degenerate[seed]) {
        continue;
      }
      const U32* seed_ids = &state.triangles[seed * 3];
      // Seed the strip on edge (seed_ids[0], seed_ids[1]); the third vertex
      // is the first decoded apex.
      mesh_detail::TrackedVector<U32> strip;
      strip.push_back(seed_ids[0]);
      strip.push_back(seed_ids[1]);
      strip.push_back(seed_ids[2]);
      state.used[seed] = 1;

      // Grow at the front: the shared edge is (strip[0], strip[1]).
      bool grew = true;
      while (grew) {
        grew = false;
        const U32 t = PickUnusedNeighbor(edges, strip[0], strip[1], state.used);
        if (t != none) {
          const U32 x = ThirdVertex(state, t, strip[0], strip[1]);
          if (x != strip[0] && x != strip[1] &&
              static_cast<std::size_t>(x) < state.vertices.size()) {
            strip.insert(strip.begin(), x);
            state.used[t] = 1;
            grew = true;
          }
        }
      }

      // Grow at the back: the shared edge is the last two strip vertices.
      grew = true;
      while (grew) {
        grew = false;
        const U32 a = strip[strip.size() - 2];
        const U32 b = strip[strip.size() - 1];
        const U32 t = PickUnusedNeighbor(edges, a, b, state.used);
        if (t != none) {
          const U32 x = ThirdVertex(state, t, a, b);
          if (x != a && x != b && static_cast<std::size_t>(x) < state.vertices.size()) {
            strip.push_back(x);
            state.used[t] = 1;
            grew = true;
          }
        }
      }

      // A strip needs at least 2 triangles (4 indices); a lone triangle is
      // left for the fan/plain stages.
      if (strip.size() >= 4) {
        state.segments.push_back({ SegmentKind::Strip, 0,
                                   static_cast<U32>(strip.size()) });
        state.indices.insert(state.indices.end(), strip.begin(), strip.end());
      } else {
        state.used[seed] = 0;
      }
    }
  }

  // --------------------------------------------------- Stage 2: triangle fans

  // Try to build one fan component starting from triangle `seed` around
  // `hub`. A candidate triangle joins the fan only if it (a) really contains
  // the hub and (b) shares a rim edge with the current open end of the fan,
  // so no triangle can ever be invented. Returns true on emission.
  static bool GrowFanFromSeed(State& state, U32 hub, U32 seed)
  {
    // The two non-hub vertices are the seed triangle's rim edge.
    mesh_detail::TrackedVector<U32> rim;
    const U32* seed_ids = &state.triangles[seed * 3];
    for (int i = 0; i < 3; ++i) {
      if (seed_ids[i] != hub) {
        rim.push_back(seed_ids[i]);
      }
    }

    mesh_detail::TrackedVector<U32> fan_triangles;
    mesh_detail::TrackedVector<unsigned char> in_fan(state.used.size(), 0);
    mesh_detail::TrackedVector<unsigned char> rim_used(state.vertices.size(), 0);
    fan_triangles.push_back(seed);
    in_fan[seed] = 1;
    for (U32 r : rim) {
      rim_used[r] = 1;
    }

    // A triangle contains both given vertices.
    auto has_both = [&](U32 t, U32 a, U32 b) {
      const U32* ids = &state.triangles[t * 3];
      bool has_a = false, has_b = false;
      for (int i = 0; i < 3; ++i) {
        has_a |= (ids[i] == a);
        has_b |= (ids[i] == b);
      }
      return has_a && has_b;
    };
    // Third vertex of a triangle known to contain hub and `end`.
    auto other_rim = [&](U32 t, U32 end) -> U32 {
      const U32* ids = &state.triangles[t * 3];
      for (int i = 0; i < 3; ++i) {
        if (ids[i] != hub && ids[i] != end) {
          return ids[i];
        }
      }
      return static_cast<U32>(state.vertices.size());
    };
    // Find an unused, non-degenerate triangle sharing the full edge
    // (hub, end) with the current boundary triangle, whose third rim vertex
    // is not already part of this fan. Smallest triangle id wins on ties.
    auto find_extension = [&](U32 boundary_tri, U32 end) -> U32 {
      U32 best = static_cast<U32>(state.used.size());
      for (std::size_t t = 0; t < state.used.size(); ++t) {
        if (state.used[t] || in_fan[t] || state.degenerate[t]) {
          continue;
        }
        if (!has_both(static_cast<U32>(t), hub, end)) {
          continue; // must share the full (hub, end) edge
        }
        const U32 x = other_rim(static_cast<U32>(t), end);
        if (static_cast<std::size_t>(x) >= state.vertices.size() || rim_used[x]) {
          continue; // would revisit an existing rim vertex -> branch, reject
        }
        (void)boundary_tri;
        if (best == state.used.size() || t < best) {
          best = static_cast<U32>(t);
        }
      }
      return best;
    };

    // Grow the back end, then the front end, until neither can extend.
    bool grew = true;
    while (grew) {
      grew = false;
      const U32 end = rim.back();
      const U32 best = find_extension(fan_triangles.back(), end);
      if (best != state.used.size()) {
        const U32 x = other_rim(best, end);
        in_fan[best] = 1;
        fan_triangles.push_back(best);
        rim_used[x] = 1;
        rim.push_back(x);
        grew = true;
      }
    }
    grew = true;
    while (grew) {
      grew = false;
      const U32 end = rim.front();
      const U32 best = find_extension(fan_triangles.front(), end);
      if (best != state.used.size()) {
        const U32 x = other_rim(best, end);
        in_fan[best] = 1;
        fan_triangles.insert(fan_triangles.begin(), best);
        rim_used[x] = 1;
        rim.insert(rim.begin(), x);
        grew = true;
      }
    }

    if (fan_triangles.size() < 2) {
      return false; // single triangle is not a fan; leave it for strips/plain
    }
    // Emit canonical fan sequence [hub, rim0, rim1, ...]. Every consecutive
    // rim pair together with the hub is one real, verified triangle.
    mesh_detail::TrackedVector<U32> sequence;
    sequence.push_back(hub);
    for (U32 r : rim) {
      sequence.push_back(r);
    }
    if (sequence.size() > 3 * fan_triangles.size()) {
      return false; // never worse than plain triangles
    }
    for (U32 t : fan_triangles) {
      state.used[t] = 1;
    }
    state.segments.push_back({ SegmentKind::Fan, 0,
                               static_cast<U32>(sequence.size()) });
    state.indices.insert(state.indices.end(), sequence.begin(), sequence.end());
    return true;
  }

  static void ExtractFans(State& state)
  {
    const std::size_t triangle_count = state.used.size();
    const std::size_t vertex_count = state.vertices.size();

    // Per-vertex list of unused, non-degenerate triangles incident to it.
    mesh_detail::TrackedVector<mesh_detail::TrackedVector<U32>> incident(
      vertex_count, mesh_detail::TrackedVector<U32>());
    for (std::size_t t = 0; t < triangle_count; ++t) {
      if (state.used[t] || state.degenerate[t]) {
        continue;
      }
      const U32* ids = &state.triangles[t * 3];
      for (int i = 0; i < 3; ++i) {
        incident[ids[i]].push_back(static_cast<U32>(t));
      }
    }

    // Hubs with the most surrounding triangles first, so the largest fans
    // are taken while the most triangles are still available.
    mesh_detail::TrackedVector<U32> hubs;
    for (U32 hub = 0; hub < vertex_count; ++hub) {
      if (incident[hub].size() >= 2) {
        hubs.push_back(hub);
      }
    }
    std::sort(hubs.begin(), hubs.end(), [&](U32 a, U32 b) {
      if (incident[a].size() != incident[b].size()) {
        return incident[a].size() > incident[b].size();
      }
      return a < b;
    });

    // Grow a fan from every still-free incident triangle of each hub.
    // GrowFanFromSeed only accepts triangles sharing a full (hub, rim) edge,
    // so it can never invent a triangle; failures leave the input untouched
    // for the plain stage.
    for (U32 hub : hubs) {
      for (U32 seed : incident[hub]) {
        if (!state.used[seed] && !state.degenerate[seed]) {
          GrowFanFromSeed(state, hub, seed);
        }
      }
    }
  }

  // --------------------------- Leftovers: plain 3-index triangle segments
  static void EmitPlain(State& state)
  {
    const std::size_t triangle_count = state.used.size();
    for (std::size_t t = 0; t < triangle_count; ++t) {
      if (state.used[t]) {
        continue;
      }
      const U32* ids = &state.triangles[t * 3];
      state.segments.push_back({ SegmentKind::Triangle, 0, 3 });
      state.indices.push_back(ids[0]);
      state.indices.push_back(ids[1]);
      state.indices.push_back(ids[2]);
      state.used[t] = 1;
    }
  }

  // ---------------------------------------------------------------- Stage 3
  static void MaterializeIndices(State& state)
  {
    // state.indices already holds the concatenated per-segment data; only the
    // offset table still needs to be stamped.
    U32 offset = 0;
    for (Segment& segment : state.segments) {
      segment.offset = offset;
      offset += segment.length;
    }
  }

  // ------------------------------- Reconstruction proof + metrics/reporting
  static void FillResult(const mesh_detail::TrackedVector<Point>& corners,
                         const TrianglePipelineOptions& options, State& state,
                         PipelineResult& result)
  {
    const std::size_t corner_count = corners.size();
    const std::size_t triangle_count = corner_count / 3;

    result.vertices.assign(state.vertices.begin(), state.vertices.end());
    if (options.emitIndices) {
      result.indices.assign(state.indices.begin(), state.indices.end());
      result.segments.assign(state.segments.begin(), state.segments.end());
      result.segment_lengths.reserve(state.segments.size());
      for (const Segment& segment : state.segments) {
        result.segment_lengths.push_back(segment.length);
      }
    }

    PipelineMetrics& metrics = result.metrics;
    metrics.input_vertex_count = corner_count;
    metrics.output_vertex_count = state.vertices.size();
    metrics.merged_vertex_count = corner_count - state.vertices.size();
    metrics.triangle_count = triangle_count;
    for (const Segment& segment : state.segments) {
      switch (segment.kind) {
        case SegmentKind::Strip:
          metrics.strip_count += 1;
          break;
        case SegmentKind::Fan:
          metrics.fan_count += 1;
          break;
        case SegmentKind::Triangle:
          metrics.plain_segment_count += 1;
          break;
      }
    }
    for (unsigned char flag : state.degenerate) {
      metrics.degenerate_count += flag;
    }
    if (options.emitIndices) {
      metrics.index_count = state.indices.size();
      if (triangle_count > 0) {
        metrics.vertex_reuse_rate =
          1.0 - static_cast<double>(state.indices.size()) / (3.0 * triangle_count);
      }
    } else {
      // Stage 3 switched off: no index buffer is materialized, so there is no
      // externally observable reuse figure to report.
      metrics.index_count = 0;
      metrics.vertex_reuse_rate = 0.0;
    }

    VerifyReconstruction(corners, options, state, metrics);
  }

  // Decode every emitted segment back into triangles and compare the resulting
  // multiset with the welded input triangle multiset. Runs in O(T log T) and
  // uses only storage already bounded by the pipeline's linear memory budget.
  static void VerifyReconstruction(const mesh_detail::TrackedVector<Point>& corners,
                                   const TrianglePipelineOptions& options, State& state,
                                   PipelineMetrics& metrics)
  {
    using IdTri = std::pair<U32, std::pair<U32, U32>>;
    auto sorted_ids = [](U32 a, U32 b, U32 c) {
      U32 k[3] = { a, b, c };
      std::sort(k, k + 3);
      return std::make_pair(k[0], std::make_pair(k[1], k[2]));
    };

    // Canonical multiset of the triangles entering stage 2 (welded ids).
    mesh_detail::TrackedVector<IdTri> expected;
    expected.reserve(metrics.triangle_count);
    for (std::size_t t = 0; t < metrics.triangle_count; ++t) {
      const U32* ids = &state.triangles[t * 3];
      expected.push_back(sorted_ids(ids[0], ids[1], ids[2]));
    }

    mesh_detail::TrackedVector<IdTri> actual;
    actual.reserve(metrics.triangle_count);
    if (options.emitIndices) {
      std::size_t cursor = 0;
      for (const Segment& segment : state.segments) {
        const U32* ids = &state.indices[cursor];
        if (segment.kind == SegmentKind::Triangle) {
          actual.push_back(sorted_ids(ids[0], ids[1], ids[2]));
        } else if (segment.kind == SegmentKind::Strip) {
          // Strip: sliding window of three, alternating winding.
          const std::size_t decoded = segment.length - 2;
          for (std::size_t i = 0; i < decoded; ++i) {
            actual.push_back(sorted_ids(ids[i], ids[i + 1], ids[i + 2]));
          }
        } else {
          // Fan: ids[0] is the shared hub; triangle i is (hub, ids[i+1],
          // ids[i+2]). Length L encodes L-2 triangles.
          const U32 hub = ids[0];
          const std::size_t decoded = segment.length - 2;
          for (std::size_t i = 0; i < decoded; ++i) {
            actual.push_back(sorted_ids(hub, ids[i + 1], ids[i + 2]));
          }
        }
        cursor += segment.length;
      }
    } else {
      // Stage 3 switched off: prove stage 2 consumed every triangle once.
      for (std::size_t t = 0; t < metrics.triangle_count; ++t) {
        if (state.used[t]) {
          const U32* ids = &state.triangles[t * 3];
          actual.push_back(sorted_ids(ids[0], ids[1], ids[2]));
        }
      }
    }

    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
    std::size_t mismatches = 0;
    if (expected.size() != actual.size()) {
      mismatches += static_cast<std::size_t>(
        std::llabs(static_cast<long long>(expected.size()) -
                   static_cast<long long>(actual.size())));
    }
    const std::size_t common = (std::min)(expected.size(), actual.size());
    for (std::size_t i = 0; i < common; ++i) {
      if (expected[i] != actual[i]) {
        mismatches += 1;
      }
    }

    // Coordinate-level end-to-end check, only when welding is lossless:
    // exact welding (tolerance 0) or no welding preserve every input corner
    // coordinate, so decoded coordinate triples must equal the input ones.
    // Under a positive tolerance the merge is lossy by definition, so the
    // welded-id multiset check above is the authoritative consistency proof.
    if (mismatches == 0 && options.emitIndices &&
        (!options.weldVertices || options.tolerance == 0.0)) {
      using Key = std::pair<double, double>;
      using TriKey = std::pair<Key, std::pair<Key, Key>>;
      auto make_key = [](const Point& p) { return std::make_pair(p.x, p.y); };
      auto sorted_coords = [&](TriKey tri) {
        Key keys[3] = { tri.first, tri.second.first, tri.second.second };
        std::sort(keys, keys + 3);
        return std::make_pair(keys[0], std::make_pair(keys[1], keys[2]));
      };
      mesh_detail::TrackedVector<TriKey> expected_coords;
      mesh_detail::TrackedVector<TriKey> actual_coords;
      expected_coords.reserve(metrics.triangle_count);
      actual_coords.reserve(metrics.triangle_count);
      for (std::size_t t = 0; t < metrics.triangle_count; ++t) {
        expected_coords.push_back(sorted_coords(std::make_pair(
          make_key(corners[t * 3]),
          std::make_pair(make_key(corners[t * 3 + 1]), make_key(corners[t * 3 + 2])))));
      }
      std::size_t cursor = 0;
      for (const Segment& segment : state.segments) {
        const U32* ids = &state.indices[cursor];
        if (segment.kind == SegmentKind::Triangle) {
          actual_coords.push_back(sorted_coords(std::make_pair(
            make_key(state.vertices[ids[0]]),
            std::make_pair(make_key(state.vertices[ids[1]]),
                           make_key(state.vertices[ids[2]])))));
        } else if (segment.kind == SegmentKind::Strip) {
          const std::size_t decoded = segment.length - 2;
          for (std::size_t i = 0; i < decoded; ++i) {
            actual_coords.push_back(sorted_coords(std::make_pair(
              make_key(state.vertices[ids[i]]),
              std::make_pair(make_key(state.vertices[ids[i + 1]]),
                             make_key(state.vertices[ids[i + 2]])))));
          }
        } else {
          // Fan: fixed hub at ids[0].
          const Point& hub_point = state.vertices[ids[0]];
          const std::size_t decoded = segment.length - 2;
          for (std::size_t i = 0; i < decoded; ++i) {
            actual_coords.push_back(sorted_coords(std::make_pair(
              make_key(hub_point),
              std::make_pair(make_key(state.vertices[ids[i + 1]]),
                             make_key(state.vertices[ids[i + 2]])))));
          }
        }
        cursor += segment.length;
      }
      std::sort(expected_coords.begin(), expected_coords.end());
      std::sort(actual_coords.begin(), actual_coords.end());
      for (std::size_t i = 0; i < expected_coords.size(); ++i) {
        if (expected_coords[i] != actual_coords[i]) {
          mismatches += 1;
        }
      }
    }

    metrics.reconstruction_mismatches = mismatches;
    metrics.reconstruction_consistent = mismatches == 0;
  }
};

}

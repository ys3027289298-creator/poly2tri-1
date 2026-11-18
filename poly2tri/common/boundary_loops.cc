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
#include "boundary_loops.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace p2t {

NonManifoldEdgeError::NonManifoldEdgeError(std::size_t half_edge)
    : std::runtime_error("Non-manifold edge at half-edge " +
                         std::to_string(half_edge) + " (triangle " +
                         std::to_string(half_edge / 3) + ", corner " +
                         std::to_string(half_edge % 3) + ")"),
      half_edge(half_edge), triangle(half_edge / 3),
      corner(static_cast<int>(half_edge % 3))
{
}

std::vector<std::size_t> BoundaryLoops::OuterLoops() const
{
  std::vector<std::size_t> outer;
  for (std::size_t i = 0; i < loops.size(); ++i) {
    if (loops[i].parent == -1) {
      outer.push_back(i);
    }
  }
  return outer;
}

namespace {

constexpr std::size_t kInvalid = std::numeric_limits<std::size_t>::max();

// Ray-casting test against the +x axis; points exactly on the boundary are
// treated as outside.
bool PointInPolygon(const Point& p, const std::vector<Point*>& polygon)
{
  bool inside = false;
  const std::size_t n = polygon.size();
  for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
    const Point& a = *polygon[j];
    const Point& b = *polygon[i];
    if ((a.y > p.y) != (b.y > p.y)) {
      const double x_intersect = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
      if (x_intersect > p.x) {
        inside = !inside;
      }
    }
  }
  return inside;
}

// Stable LSD radix sort of 64-bit keys with 32-bit payloads, 16 bits per
// pass. Beats std::sort here because the comparison-based alternative
// dominates the runtime of ExtractBoundaryLoops on large meshes.
void RadixSort(std::vector<std::uint64_t>& keys,
               std::vector<std::uint32_t>& payload)
{
  const std::size_t n = keys.size();
  std::vector<std::uint64_t> tmp_keys(n);
  std::vector<std::uint32_t> tmp_payload(n);
  std::vector<std::uint32_t> count(1 << 16);
  for (int shift = 0; shift < 64; shift += 16) {
    std::fill(count.begin(), count.end(), 0);
    for (std::size_t k = 0; k < n; ++k) {
      ++count[(keys[k] >> shift) & 0xffff];
    }
    std::uint32_t sum = 0;
    for (std::uint32_t& c : count) {
      const std::uint32_t t = c;
      c = sum;
      sum += t;
    }
    for (std::size_t k = 0; k < n; ++k) {
      const std::uint32_t pos = count[(keys[k] >> shift) & 0xffff]++;
      tmp_keys[pos] = keys[k];
      tmp_payload[pos] = payload[k];
    }
    keys.swap(tmp_keys);
    payload.swap(tmp_payload);
  }
}

} // namespace

BoundaryLoops ExtractBoundaryLoops(const std::vector<Triangle*>& triangles)
{
  BoundaryLoops result;
  const std::size_t num_triangles = triangles.size();
  const std::size_t num_half_edges = 3 * num_triangles;
  if (num_triangles == 0) {
    return result;
  }

  // Collect the corner points and assign each vertex an index by radix
  // sorting the pointer keys (faster than a hash map or std::sort here).
  std::vector<Point*> corner(num_half_edges);
  for (std::size_t t = 0; t < num_triangles; ++t) {
    Triangle& triangle = *triangles[t];
    for (std::size_t c = 0; c < 3; ++c) {
      corner[3 * t + c] = triangle.GetPoint(static_cast<int>(c));
    }
  }
  std::vector<std::uint64_t> sort_keys(num_half_edges);
  std::vector<std::uint32_t> sort_payload(num_half_edges);
  for (std::size_t he = 0; he < num_half_edges; ++he) {
    sort_keys[he] = reinterpret_cast<std::uintptr_t>(corner[he]);
    sort_payload[he] = static_cast<std::uint32_t>(he);
  }
  RadixSort(sort_keys, sort_payload);
  std::vector<Point*> vertices;
  vertices.reserve(num_half_edges);
  std::vector<std::uint32_t> from(num_half_edges), to(num_half_edges);
  for (std::size_t k = 0; k < num_half_edges; ++k) {
    if (k == 0 || sort_keys[k] != sort_keys[k - 1]) {
      vertices.push_back(corner[sort_payload[k]]);
    }
    from[sort_payload[k]] = static_cast<std::uint32_t>(vertices.size() - 1);
  }
  corner.clear();
  for (std::size_t he = 0; he < num_half_edges; ++he) {
    to[he] = from[he - (he % 3) + (he % 3 + 1) % 3];
  }

  // Pair half-edges with their twins by radix sorting on the undirected
  // edge key. The sort is stable, so within one edge's group the half-edges
  // stay ordered by their index.
  for (std::size_t he = 0; he < num_half_edges; ++he) {
    std::uint32_t a = from[he], b = to[he];
    if (a > b) {
      std::swap(a, b);
    }
    sort_keys[he] = (static_cast<std::uint64_t>(a) << 32) | b;
    sort_payload[he] = static_cast<std::uint32_t>(he);
  }
  RadixSort(sort_keys, sort_payload);
  std::vector<std::size_t> twin(num_half_edges, kInvalid);
  std::size_t i = 0;
  while (i < num_half_edges) {
    std::size_t j = i + 1;
    while (j < num_half_edges && sort_keys[j] == sort_keys[i]) {
      ++j;
    }
    const std::size_t count = j - i;
    if (count == 2) {
      const std::size_t he0 = sort_payload[i], he1 = sort_payload[i + 1];
      if (from[he0] == from[he1]) {
        // Same undirected edge wound in the same direction twice.
        throw NonManifoldEdgeError(he0);
      }
      twin[he0] = he1;
      twin[he1] = he0;
    } else if (count > 2) {
      throw NonManifoldEdgeError(sort_payload[i]);
    }
    i = j;
  }

  // Chain boundary half-edges: the loop successor of h is found by rotating
  // around h's target vertex through adjacent triangles until the next
  // boundary half-edge is reached.
  const auto he_next = [](std::size_t he) {
    return he - (he % 3) + (he % 3 + 1) % 3;
  };
  std::vector<std::size_t> next_boundary(num_half_edges, kInvalid);
  for (std::size_t h = 0; h < num_half_edges; ++h) {
    if (twin[h] != kInvalid) {
      continue;
    }
    std::size_t g = he_next(h);
    for (std::size_t steps = 0; twin[g] != kInvalid; ++steps) {
      if (steps > num_half_edges) {
        // A boundary half-edge must terminate the rotation; this is only
        // reachable on degenerate input.
        throw NonManifoldEdgeError(h);
      }
      g = he_next(twin[g]);
    }
    next_boundary[h] = g;
  }

  // Trace the loops and classify them by signed area.
  std::vector<char> visited(num_half_edges, 0);
  for (std::size_t h = 0; h < num_half_edges; ++h) {
    if (twin[h] != kInvalid || visited[h]) {
      continue;
    }
    BoundaryLoop loop;
    double area2 = 0.0;
    std::size_t he = h;
    do {
      visited[he] = 1;
      const Point* p = vertices[from[he]];
      const Point* q = vertices[to[he]];
      loop.points.push_back(const_cast<Point*>(p));
      area2 += p->x * q->y - q->x * p->y;
      he = next_boundary[he];
    } while (he != h);
    loop.signed_area = 0.5 * area2;
    loop.is_hole = area2 < 0.0;
    loop.parent = -1;
    loop.depth = 0;
    result.loops.push_back(std::move(loop));
  }

  // Nesting hierarchy: the parent of a loop is the smallest-area loop that
  // contains its leftmost-lowest vertex.
  const std::size_t num_loops = result.loops.size();
  for (std::size_t a = 0; a < num_loops; ++a) {
    BoundaryLoop& loop = result.loops[a];
    const Point* representative = loop.points[0];
    for (const Point* p : loop.points) {
      if (p->x < representative->x ||
          (p->x == representative->x && p->y < representative->y)) {
        representative = p;
      }
    }
    const double abs_area = std::fabs(loop.signed_area);
    double best_area = std::numeric_limits<double>::max();
    for (std::size_t b = 0; b < num_loops; ++b) {
      if (b == a) {
        continue;
      }
      const double container_area = std::fabs(result.loops[b].signed_area);
      // A valid container always encloses a larger area.
      if (container_area <= abs_area || container_area >= best_area) {
        continue;
      }
      if (PointInPolygon(*representative, result.loops[b].points)) {
        loop.parent = static_cast<int>(b);
        best_area = container_area;
      }
    }
  }
  for (std::size_t a = 0; a < num_loops; ++a) {
    int depth = 0;
    for (int p = result.loops[a].parent; p != -1; p = result.loops[p].parent) {
      ++depth;
    }
    result.loops[a].depth = depth;
  }

  return result;
}

}

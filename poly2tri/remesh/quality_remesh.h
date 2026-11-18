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

// Quality-driven 2D remeshing on top of a poly2tri CDT result.
//
// The remesher never modifies the CDT output: the input Triangle*/Point*
// objects are read once, copied into a private mesh, and all edge splits and
// edge collapses operate on that copy. A full invariant check (constrained
// edge coverage, segment self-intersection, positive oriented area and
// manifold incidence) runs after every tentative operation; an operation
// whose check fails is rolled back and its edge is blacklisted.
//
// Boundary / constrained edges may only be split, never collapsed, so every
// original constrained segment stays present as an exact chain of sub-edges.
// Interior edge collapses are restricted to good-quality regions, so they
// can never undo a refinement performed to chase the target minimum angle.

#pragma once

#include "../common/dll_symbol.h"
#include "../common/shapes.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace p2t {

/// Configurable quality targets and safety limits.
struct P2T_DLL_SYMBOL RemeshOptions {
  /// Triangles whose minimum angle (degrees) is below this are "bad".
  double target_min_angle_deg = 25.0;
  /// Edges longer than this are split. A value <= 0 selects the median edge
  /// length of the input mesh automatically.
  double max_edge_length = 0.0;
  /// Interior edges shorter than this become collapse candidates (only in
  /// good-quality regions). A value <= 0 selects 0.5 * max_edge_length.
  double min_edge_length = 0.0;
  /// A triangle is also "bad" when longest_edge / shortest_edge exceeds this.
  double max_edge_ratio = 2.5;
  /// Hard stop on the number of attempted split/collapse operations.
  std::size_t max_operations = 3000;
  /// Run the full constraint/self-intersection validation after every
  /// tentative operation (required by the quality-driven procedure).
  bool check_each_operation = true;
};

/// One output triangle as compact vertex indices.
struct P2T_DLL_SYMBOL RemeshFace {
  int v[3];
};

/// Full remeshing output, including the measured quality report. All angle
/// values are in degrees.
struct P2T_DLL_SYMBOL RemeshResult {
  /// Owned output mesh (only alive vertices, indices match RemeshFace::v).
  std::vector<Point> vertices;
  std::vector<RemeshFace> faces;

  std::size_t input_triangles = 0;
  std::size_t output_triangles = 0;

  double min_angle_before = 0.0;
  double min_angle_after = 0.0;
  double max_edge_before = 0.0;
  double max_edge_after = 0.0;
  double min_edge_before = 0.0;
  double min_edge_after = 0.0;

  double area_before = 0.0;
  double area_after = 0.0;
  double area_error_abs = 0.0;
  double area_error_rel = 0.0;

  /// Number of original constrained segments and how many are still exactly
  /// covered by a chain of constrained sub-edges after remeshing.
  std::size_t constraints_total = 0;
  std::size_t constraints_satisfied = 0;

  bool self_intersecting = false;
  bool target_angle_reached = false;
  /// Stays false while every post-operation validation passed.
  bool invariant_violated = false;

  std::size_t split_count = 0;
  std::size_t steiner_count = 0;
  std::size_t collapse_count = 0;
  std::size_t rejected_count = 0;
  std::size_t checks_run = 0;

  double elapsed_seconds = 0.0;
};

namespace remesh_detail {

inline uint64_t EdgeKey(int a, int b)
{
  if (a > b) {
    const int t = a;
    a = b;
    b = t;
  }
  return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
         static_cast<uint32_t>(b);
}

// Decode helpers: go through uint32_t first so vertex ids never hit
// implementation-defined signed overflow when the id space grows large.
inline int KeyHi(uint64_t key)
{
  return static_cast<int>(static_cast<uint32_t>(key >> 32));
}

inline int KeyLo(uint64_t key)
{
  return static_cast<int>(static_cast<uint32_t>(key & 0xffffffffu));
}

struct RVertex {
  double x, y;
  bool alive;
};

struct RFace {
  int v[3];
  bool alive;
  double min_angle_deg;
  double shortest, longest;
};

/// Original constrained segment, stored in absolute coordinates so that the
/// coverage check is independent of internal vertex indexing.
struct RSegment {
  double ax, ay, bx, by;
};

struct RMesh {
  std::vector<RVertex> verts;
  std::vector<RFace> faces;
  std::unordered_set<uint64_t> constrained;
  std::vector<RSegment> constraints;
  /// Edge key -> incident face indices (one for boundary, two otherwise).
  std::unordered_map<uint64_t, std::vector<int>> edges;
  std::vector<bool> locked;
  double scale; // bounding-box diagonal, used for tolerances

  // Tight tolerance for topology predicates (intersection, orientation).
  double eps() const { return 1e-11 * std::max(1.0, scale); }
  // Looser tolerance for collinearity/coverage bookkeeping.
  double tol() const { return 1e-9 * std::max(1.0, scale); }
};

inline double EdgLen(const RMesh& m, int a, int b)
{
  const double dx = m.verts[a].x - m.verts[b].x;
  const double dy = m.verts[a].y - m.verts[b].y;
  return std::sqrt(dx * dx + dy * dy);
}

inline double Orient(double ax, double ay, double bx, double by,
                     double cx, double cy)
{
  return (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
}

inline double FaceArea(const RMesh& m, const RFace& f)
{
  return 0.5 * std::fabs(Orient(m.verts[f.v[0]].x, m.verts[f.v[0]].y,
                                m.verts[f.v[1]].x, m.verts[f.v[1]].y,
                                m.verts[f.v[2]].x, m.verts[f.v[2]].y));
}

inline double MinAngleDeg(const RMesh& m, const RFace& f)
{
  static const double kDeg = 180.0 / 3.14159265358979323846;
  double best = 180.0;
  for (int i = 0; i < 3; ++i) {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    const double ux = m.verts[f.v[j]].x - m.verts[f.v[i]].x;
    const double uy = m.verts[f.v[j]].y - m.verts[f.v[i]].y;
    const double vx = m.verts[f.v[k]].x - m.verts[f.v[i]].x;
    const double vy = m.verts[f.v[k]].y - m.verts[f.v[i]].y;
    const double dot = ux * vx + uy * vy;
    const double cross = ux * vy - uy * vx;
    best = std::min(best, std::atan2(std::fabs(cross), dot) * kDeg);
  }
  return best;
}

inline void FaceMetrics(RMesh& m, RFace& f)
{
  const double l0 = EdgLen(m, f.v[0], f.v[1]);
  const double l1 = EdgLen(m, f.v[1], f.v[2]);
  const double l2 = EdgLen(m, f.v[2], f.v[0]);
  f.shortest = std::min(l0, std::min(l1, l2));
  f.longest = std::max(l0, std::max(l1, l2));
  f.min_angle_deg = MinAngleDeg(m, f);
}

/// Rebuild the edge-incidence map and the locked (constrained-endpoint)
/// vertex flags from the live face set.
inline void RebuildTopology(RMesh& m)
{
  m.edges.clear();
  m.locked.assign(m.verts.size(), false);
  for (int fi = 0; fi < static_cast<int>(m.faces.size()); ++fi) {
    const RFace& f = m.faces[fi];
    if (!f.alive) {
      continue;
    }
    for (int i = 0; i < 3; ++i) {
      const uint64_t key = EdgeKey(f.v[i], f.v[(i + 1) % 3]);
      m.edges[key].push_back(fi);
    }
  }
  for (const auto& kv : m.edges) {
    // Boundary edges are implicitly constrained even if a source flag was
    // missing; they must be preserved just like explicit constraints.
    if (kv.second.size() == 1 || m.constrained.count(kv.first) != 0) {
      const int a = KeyHi(kv.first);
      const int b = KeyLo(kv.first);
      m.locked[a] = true;
      m.locked[b] = true;
    }
  }
}

inline bool OnSegment(double px, double py, const RSegment& s, double tol)
{
  const double edge_len = std::hypot(s.bx - s.ax, s.by - s.ay);
  const double o = Orient(s.ax, s.ay, s.bx, s.by, px, py);
  if (std::fabs(o) > tol * std::max(1.0, edge_len)) {
    return false;
  }
  return std::min(s.ax, s.bx) - tol <= px && px <= std::max(s.ax, s.bx) + tol &&
         std::min(s.ay, s.by) - tol <= py && py <= std::max(s.ay, s.by) + tol;
}

/// Every original constrained segment must be exactly covered by a chain of
/// live constrained sub-edges (no gaps, no overshoot).
inline bool ConstraintsCovered(const RMesh& m)
{
  const double tol = m.tol();
  for (const RSegment& s : m.constraints) {
    const double dx = s.bx - s.ax, dy = s.by - s.ay;
    const double len2 = dx * dx + dy * dy;
    if (len2 <= 0.0) {
      continue;
    }
    std::vector<std::pair<double, double>> intervals;
    for (const auto& kv : m.edges) {
      if (m.constrained.count(kv.first) == 0) {
        continue;
      }
      const int a = KeyHi(kv.first);
      const int b = KeyLo(kv.first);
      const RVertex& pa = m.verts[a];
      const RVertex& pb = m.verts[b];
      if (!OnSegment(pa.x, pa.y, s, tol) ||
          !OnSegment(pb.x, pb.y, s, tol)) {
        continue;
      }
      const double ta = ((pa.x - s.ax) * dx + (pa.y - s.ay) * dy) / len2;
      const double tb = ((pb.x - s.ax) * dx + (pb.y - s.ay) * dy) / len2;
      intervals.emplace_back(std::min(ta, tb), std::max(ta, tb));
    }
    if (intervals.empty()) {
      return false;
    }
    std::sort(intervals.begin(), intervals.end());
    if (intervals[0].first > tol) {
      return false; // gap before the first sub-edge
    }
    double reach = intervals[0].second;
    for (std::size_t i = 1; i < intervals.size(); ++i) {
      if (intervals[i].first > reach + tol) {
        return false; // gap between sub-edges
      }
      reach = std::max(reach, intervals[i].second);
    }
    if (reach < 1.0 - tol) {
      return false; // gap before the segment end
    }
  }
  return true;
}

inline bool ProperIntersect(double ax, double ay, double bx, double by,
                            double cx, double cy, double dx, double dy,
                            double tol)
{
  const double o1 = Orient(ax, ay, bx, by, cx, cy);
  const double o2 = Orient(ax, ay, bx, by, dx, dy);
  const double o3 = Orient(cx, cy, dx, dy, ax, ay);
  const double o4 = Orient(cx, cy, dx, dy, bx, by);
  const bool s1 = ((o1 > tol && o2 < -tol) || (o1 < -tol && o2 > tol));
  const bool s2 = ((o3 > tol && o4 < -tol) || (o3 < -tol && o4 > tol));
  if (s1 && s2) {
    return true; // proper crossing through both edge interiors
  }
  // Collinear overlap where a non-endpoint point of one edge lies strictly
  // inside the other (shared endpoints alone are allowed).
  auto interior_overlap = [&](double px, double py) {
    if (std::fabs(Orient(ax, ay, bx, by, px, py)) > tol) {
      return false;
    }
    const double ex = std::max(ax, bx) - std::min(ax, bx);
    const double ey = std::max(ay, by) - std::min(ay, by);
    const bool on =
        std::min(ax, bx) - tol <= px && px <= std::max(ax, bx) + tol &&
        std::min(ay, by) - tol <= py && py <= std::max(ay, by) + tol;
    const bool endpoint =
        (std::fabs(px - ax) <= tol && std::fabs(py - ay) <= tol) ||
        (std::fabs(px - bx) <= tol && std::fabs(py - by) <= tol);
    // Ignore endpoint-sized noise: overlaps shorter than the tolerance are
    // not treated as real geometry.
    return on && !endpoint && (ex > 4.0 * tol || ey > 4.0 * tol);
  };
  return interior_overlap(cx, cy) || interior_overlap(dx, dy);
}

/// O(E^2) global self-intersection test. Demo meshes are small, and running
/// the exhaustive test after every operation is exactly the requested guard.
inline bool HasSelfIntersection(const RMesh& m)
{
  const double tol = m.eps();
  std::vector<std::pair<int, int>> edges;
  edges.reserve(m.edges.size());
  for (const auto& kv : m.edges) {
    const int a = KeyHi(kv.first);
    const int b = KeyLo(kv.first);
    edges.emplace_back(a, b);
    if (kv.second.size() > 2) {
      return true; // non-manifold incidence
    }
  }
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const RVertex& a = m.verts[edges[i].first];
    const RVertex& b = m.verts[edges[i].second];
    for (std::size_t j = i + 1; j < edges.size(); ++j) {
      if (edges[i].first == edges[j].first ||
          edges[i].first == edges[j].second ||
          edges[i].second == edges[j].first ||
          edges[i].second == edges[j].second) {
        continue; // shared endpoint is fine
      }
      const RVertex& c = m.verts[edges[j].first];
      const RVertex& d = m.verts[edges[j].second];
      if (ProperIntersect(a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y, tol)) {
        return true;
      }
    }
  }
  return false;
}

inline bool AllFacesPositive(const RMesh& m)
{
  const double min_area = 1e-15 * std::max(1.0, m.scale * m.scale);
  for (const RFace& f : m.faces) {
    if (f.alive && FaceArea(m, f) < min_area) {
      return false;
    }
  }
  return true;
}

inline bool Validate(const RMesh& m)
{
  return AllFacesPositive(m) && ConstraintsCovered(m) &&
         !HasSelfIntersection(m);
}

struct Snapshot {
  std::vector<RVertex> verts;
  std::vector<RFace> faces;
  std::unordered_set<uint64_t> constrained;
};

inline Snapshot TakeSnapshot(const RMesh& m)
{
  return Snapshot{ m.verts, m.faces, m.constrained };
}

inline void Restore(RMesh& m, const Snapshot& s)
{
  m.verts = s.verts;
  m.faces = s.faces;
  m.constrained = s.constrained;
  RebuildTopology(m);
}

/// Build a CCW face from three vertex indices.
inline RFace MakeFace(RMesh& m, int a, int b, int c)
{
  RFace f;
  f.alive = true;
  f.v[0] = a;
  f.v[1] = b;
  f.v[2] = c;
  if (Orient(m.verts[a].x, m.verts[a].y, m.verts[b].x, m.verts[b].y,
             m.verts[c].x, m.verts[c].y) < 0.0) {
    f.v[1] = c;
    f.v[2] = b;
  }
  FaceMetrics(m, f);
  return f;
}

/// Circumcenter of triangle (a, b, c).
inline RVertex Circumcenter(const RMesh& m, int a, int b, int c)
{
  const double ax = m.verts[a].x, ay = m.verts[a].y;
  const double bx = m.verts[b].x, by = m.verts[b].y;
  const double cx = m.verts[c].x, cy = m.verts[c].y;
  const double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  const double a2 = ax * ax + ay * ay;
  const double b2 = bx * bx + by * by;
  const double c2 = cx * cx + cy * cy;
  RVertex o;
  o.x = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
  o.y = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
  o.alive = true;
  return o;
}

/// Choose the best split point on edge (a, b). Boundary/constrained edges
/// must keep their midpoint so the original segment stays exactly covered.
/// For interior edges the candidate set is the midpoint plus the projection
/// of each incident triangle's circumcenter onto the edge; whichever gives
/// the largest worst minimum angle in the two child triangles wins. This is
/// still an edge split (the new vertex lies on the edge) but behaves like
/// circumcenter insertion where that is favorable.
inline RVertex ChooseSplitPoint(RMesh& m, uint64_t key)
{
  const int a = KeyHi(key);
  const int b = KeyLo(key);
  RVertex mid;
  mid.x = 0.5 * (m.verts[a].x + m.verts[b].x);
  mid.y = 0.5 * (m.verts[a].y + m.verts[b].y);
  mid.alive = true;

  if (m.constrained.count(key) != 0 ||
      m.edges.find(key)->second.size() == 1) {
    return mid;
  }

  std::vector<RVertex> candidates{ mid };
  for (const int fi : m.edges.find(key)->second) {
    const RFace& f = m.faces[fi];
    int ic = -1;
    for (int i = 0; i < 3; ++i) {
      if (f.v[i] != a && f.v[i] != b) {
        ic = f.v[i];
      }
    }
    if (ic < 0) {
      continue;
    }
    const RVertex o = Circumcenter(m, a, b, ic);
    // Project the circumcenter onto the edge line.
    const double ex = m.verts[b].x - m.verts[a].x;
    const double ey = m.verts[b].y - m.verts[a].y;
    const double el2 = ex * ex + ey * ey;
    double t = ((o.x - m.verts[a].x) * ex + (o.y - m.verts[a].y) * ey) / el2;
    t = std::max(0.15, std::min(0.85, t)); // keep away from endpoints
    RVertex p;
    p.x = m.verts[a].x + t * ex;
    p.y = m.verts[a].y + t * ey;
    p.alive = true;
    candidates.push_back(p);
  }

  auto worst_child_angle = [&](const RVertex& p) {
    double worst = 180.0;
    for (const int fi : m.edges.find(key)->second) {
      const RFace& f = m.faces[fi];
      int ic = -1;
      for (int i = 0; i < 3; ++i) {
        if (f.v[i] != a && f.v[i] != b) {
          ic = f.v[i];
        }
      }
      RMesh probe = m; // cheap local copy for scoring
      RFace f1, f2;
      f1.alive = f2.alive = true;
      f1.v[0] = a; f1.v[1] = static_cast<int>(probe.verts.size());
      f1.v[2] = ic;
      f2.v[0] = static_cast<int>(probe.verts.size());
      f2.v[1] = b; f2.v[2] = ic;
      probe.verts.push_back(p);
      for (RFace* ff : { &f1, &f2 }) {
        if (Orient(probe.verts[ff->v[0]].x, probe.verts[ff->v[0]].y,
                   probe.verts[ff->v[1]].x, probe.verts[ff->v[1]].y,
                   probe.verts[ff->v[2]].x, probe.verts[ff->v[2]].y) < 0.0) {
          std::swap(ff->v[1], ff->v[2]);
        }
        worst = std::min(worst, MinAngleDeg(probe, *ff));
      }
    }
    return worst;
  };

  RVertex best = mid;
  double best_angle = -1.0;
  for (const RVertex& p : candidates) {
    const double ang = worst_child_angle(p);
    if (ang > best_angle) {
      best_angle = ang;
      best = p;
    }
  }
  return best;
}

/// Split edge (a, b) at the supplied point (which lies on the edge). Works
/// for both interior edges (one or two incident faces) and constrained
/// boundary edges.
inline bool SplitEdge(RMesh& m, uint64_t key, const RVertex& split_point)
{
  const auto it = m.edges.find(key);
  if (it == m.edges.end()) {
    return false;
  }
  const int a = KeyHi(key);
  const int b = KeyLo(key);
  const bool was_constrained = m.constrained.count(key) != 0;

  const int mi = static_cast<int>(m.verts.size());
  m.verts.push_back(split_point);

  // Copy the incident face indices first: push_back below may reallocate the
  // face vector and invalidate any RFace reference taken before it.
  const std::vector<int> incident = it->second;
  for (const int fi : incident) {
    const RFace f = m.faces[fi]; // snapshot by value; the slot is killed
    int ia = -1, ib = -1;
    for (int i = 0; i < 3; ++i) {
      if (f.v[i] == a) ia = i;
      if (f.v[i] == b) ib = i;
    }
    if (ia < 0 || ib < 0) {
      return false;
    }
    const int ic = 3 - ia - ib;
    const int c = f.v[ic];
    // Split the CCW triangle (a, b, c) into (a, mid, c) and (mid, b, c).
    m.faces.push_back(MakeFace(m, a, mi, c));
    m.faces.push_back(MakeFace(m, mi, b, c));
    m.faces[fi].alive = false;
  }

  m.constrained.erase(key);
  if (was_constrained) {
    m.constrained.insert(EdgeKey(a, mi));
    m.constrained.insert(EdgeKey(mi, b));
  }
  RebuildTopology(m);
  for (RFace& f : m.faces) {
    if (f.alive) {
      FaceMetrics(m, f);
    }
  }
  return true;
}

/// Insert an interior point into whichever live face contains it, replacing
/// that face by a three-triangle fan. Three outcomes are handled:
///   * point strictly inside a face -> three-way fan;
///   * point on an unconstrained interior edge -> that edge is split;
///   * point on a constrained edge, or outside the domain -> rejected.
/// Returns the new vertex index, or -1 when the point cannot be inserted.
inline int InsertInteriorPoint(RMesh& m, const RVertex& point)
{
  const double tol = m.eps();
  int host = -1;
  int on_edge_key_hi = -1, on_edge_key_lo = -1;
  for (int fi = 0; fi < static_cast<int>(m.faces.size()); ++fi) {
    const RFace& f = m.faces[fi];
    if (!f.alive) {
      continue;
    }
    bool inside = true;
    for (int i = 0; i < 3; ++i) {
      const int p = f.v[i], q = f.v[(i + 1) % 3];
      const double o = Orient(m.verts[p].x, m.verts[p].y,
                              m.verts[q].x, m.verts[q].y,
                              point.x, point.y);
      if (o <= tol) {
        inside = false;
      }
      if (std::fabs(o) <= tol) {
        // Only remember when the point actually lies within the edge span.
        const RVertex& vp = m.verts[p];
        const RVertex& vq = m.verts[q];
        const bool on =
            std::min(vp.x, vq.x) - tol <= point.x &&
            point.x <= std::max(vp.x, vq.x) + tol &&
            std::min(vp.y, vq.y) - tol <= point.y &&
            point.y <= std::max(vp.y, vq.y) + tol;
        if (on) {
          on_edge_key_hi = p;
          on_edge_key_lo = q;
        }
      }
    }
    if (inside) {
      host = fi;
      break;
    }
  }

  if (host < 0) {
    if (on_edge_key_hi < 0) {
      return -1; // point outside the domain
    }
    const uint64_t ekey = EdgeKey(on_edge_key_hi, on_edge_key_lo);
    if (m.constrained.count(ekey) != 0) {
      return -1; // never split a constrained edge with a Steiner point
    }
    if (!SplitEdge(m, ekey, point)) {
      return -1;
    }
    // The split point is the newest vertex.
    return static_cast<int>(m.verts.size()) - 1;
  }

  const RFace f = m.faces[host];
  // Refuse if the host touches a constrained edge: a Steiner vertex must
  // never land on or cross a constraint.
  for (int i = 0; i < 3; ++i) {
    if (m.constrained.count(EdgeKey(f.v[i], f.v[(i + 1) % 3])) != 0) {
      return -1;
    }
  }
  const int ci = static_cast<int>(m.verts.size());
  m.verts.push_back(point);
  for (int i = 0; i < 3; ++i) {
    m.faces.push_back(MakeFace(m, f.v[i], f.v[(i + 1) % 3], ci));
  }
  m.faces[host].alive = false;
  RebuildTopology(m);
  for (RFace& ff : m.faces) {
    if (ff.alive) {
      FaceMetrics(m, ff);
    }
  }
  return ci;
}

/// Ruppert-style circumcenter insertion: add a Steiner vertex at the
/// circumcenter of one interior triangle. For an acute face the center is
/// inside the face; for an obtuse face it lies in a neighboring face, which
/// InsertInteriorPoint locates and splits. Only fully interior faces (no
/// constrained edge) are eligible, so no constraint is ever touched.
/// Returns false (and changes nothing) when not eligible.
inline int InsertCircumcenter(RMesh& m, int face_index)
{
  if (face_index < 0 || face_index >= static_cast<int>(m.faces.size())) {
    return -1;
  }
  const RFace f = m.faces[face_index];
  if (!f.alive || f.min_angle_deg <= 0.0) {
    return -1;
  }
  for (int i = 0; i < 3; ++i) {
    const uint64_t key = EdgeKey(f.v[i], f.v[(i + 1) % 3]);
    if (m.constrained.count(key) != 0) {
      return -1; // never insert across a constrained edge
    }
  }
  const RVertex cc = Circumcenter(m, f.v[0], f.v[1], f.v[2]);
  return InsertInteriorPoint(m, cc);
}

/// Delaunay edge flip: replace a non-constrained interior edge shared by
/// two faces (a, b, c) and (b, a, d) with the opposite diagonal (c, d).
/// Returns false if the edge cannot be flipped (boundary, constrained, or
/// the new diagonal would create a degenerate/non-convex quad).
inline bool FlipInteriorEdge(RMesh& m, uint64_t key)
{
  const auto it = m.edges.find(key);
  if (it == m.edges.end() || it->second.size() != 2 ||
      m.constrained.count(key) != 0) {
    return false;
  }
  const int a = KeyHi(key);
  const int b = KeyLo(key);
  const RFace f0 = m.faces[it->second[0]];
  const RFace f1 = m.faces[it->second[1]];
  int c = -1, d = -1;
  for (int i = 0; i < 3; ++i) {
    if (f0.v[i] != a && f0.v[i] != b) c = f0.v[i];
    if (f1.v[i] != a && f1.v[i] != b) d = f1.v[i];
  }
  if (c < 0 || d < 0 || c == d) {
    return false;
  }
  // The quad (a, c, b, d) must be strictly convex for a valid diagonal.
  if (Orient(m.verts[a].x, m.verts[a].y, m.verts[c].x, m.verts[c].y,
             m.verts[b].x, m.verts[b].y) <= 0.0) {
    return false;
  }
  if (Orient(m.verts[b].x, m.verts[b].y, m.verts[d].x, m.verts[d].y,
             m.verts[a].x, m.verts[a].y) <= 0.0) {
    return false;
  }
  // Reject if the new diagonal already exists (would duplicate an edge).
  if (m.edges.count(EdgeKey(c, d)) != 0) {
    return false;
  }
  m.faces[it->second[0]].alive = false;
  m.faces[it->second[1]].alive = false;
  m.faces.push_back(MakeFace(m, a, c, d));
  m.faces.push_back(MakeFace(m, c, b, d));
  RebuildTopology(m);
  for (RFace& ff : m.faces) {
    if (ff.alive) {
      FaceMetrics(m, ff);
    }
  }
  return true;
}

/// In-circle predicate: positive when point d is strictly inside the
/// circumcircle of the CCW triangle (a, b, c).
inline double InCircle(const RMesh& m, int a, int b, int c, int d)
{
  const double ax = m.verts[a].x, ay = m.verts[a].y;
  const double bx = m.verts[b].x, by = m.verts[b].y;
  const double cx = m.verts[c].x, cy = m.verts[c].y;
  const double dx = m.verts[d].x, dy = m.verts[d].y;
  const double a1 = ax - dx, a2 = ay - dy;
  const double b1 = bx - dx, b2 = by - dy;
  const double c1 = cx - dx, c2 = cy - dy;
  return (a1 * a1 + a2 * a2) * (b1 * c2 - b2 * c1) -
         (b1 * b1 + b2 * b2) * (a1 * c2 - a2 * c1) +
         (c1 * c1 + c2 * c2) * (a1 * b2 - a2 * b1);
}

/// Re-establish the Delaunay property around a freshly inserted vertex by
/// walking outward and flipping every non-constrained interior edge whose
/// opposite vertex violates the empty-circumcircle test. A bounded work
/// count guards the walk.
inline void LegalizeAround(RMesh& m, int new_vertex, std::size_t max_flips)
{
  std::size_t flips = 0;
  bool changed = true;
  while (changed && flips < max_flips) {
    changed = false;
    RebuildTopology(m);
    for (const auto& kv : m.edges) {
      const uint64_t key = kv.first;
      if (kv.second.size() != 2 || m.constrained.count(key) != 0) {
        continue;
      }
      const int a = KeyHi(key);
      const int b = KeyLo(key);
      const RFace& f0 = m.faces[kv.second[0]];
      const RFace& f1 = m.faces[kv.second[1]];
      const bool n0 = f0.v[0] == new_vertex || f0.v[1] == new_vertex ||
                      f0.v[2] == new_vertex;
      const bool n1 = f1.v[0] == new_vertex || f1.v[1] == new_vertex ||
                      f1.v[2] == new_vertex;
      if (n0 == n1) {
        continue; // edge must border exactly one triangle with the new vertex
      }
      const RFace& owner = n0 ? f0 : f1;
      const RFace& other = n0 ? f1 : f0;
      int opp = -1;
      for (int i = 0; i < 3; ++i) {
        if (other.v[i] != a && other.v[i] != b) {
          opp = other.v[i];
        }
      }
      if (opp < 0) {
        continue;
      }
      // The owner triangle is (a, b, new_vertex); check whether opp enters
      // its circumcircle (positive determinant for a CCW face).
      const double det = InCircle(m, owner.v[0], owner.v[1], owner.v[2], opp);
      if (det > m.eps() * m.eps() * m.eps()) {
        if (FlipInteriorEdge(m, key)) {
          ++flips;
          changed = true;
          break; // restart the scan after mutating the mesh
        }
      }
    }
  }
}

/// Collapse an interior edge to its midpoint and rebuild the one-ring.
/// Boundary/constrained edges and edges touching a locked vertex are
/// rejected; callers additionally rely on post-operation validation.
inline bool CollapseEdge(RMesh& m, uint64_t key)
{
  const auto it = m.edges.find(key);
  if (it == m.edges.end() || it->second.size() != 2) {
    return false;
  }
  const int a = KeyHi(key);
  const int bb = KeyLo(key);
  if (m.constrained.count(key) != 0 || m.locked[a] || m.locked[bb]) {
    return false;
  }

  RVertex mid;
  mid.x = 0.5 * (m.verts[a].x + m.verts[bb].x);
  mid.y = 0.5 * (m.verts[a].y + m.verts[bb].y);
  mid.alive = true;
  const int mi = static_cast<int>(m.verts.size());
  m.verts.push_back(mid);

  // The two incident triangles are removed. The link vertices (opposite the
  // collapsed edge in each incident triangle) bound the one-ring that has to
  // be re-triangulated.
  const std::vector<int> incident = it->second;
  std::vector<int> links;
  for (const int fi : incident) {
    const RFace f = m.faces[fi];
    for (int i = 0; i < 3; ++i) {
      const int vi = f.v[i];
      if (vi != a && vi != bb) {
        links.push_back(vi);
      }
    }
    m.faces[fi].alive = false;
  }
  if (links.size() != 2 || links[0] == links[1]) {
    return false;
  }

  // Every other live face touching a or b is rewired to the midpoint.
  for (RFace& f : m.faces) {
    if (!f.alive) {
      continue;
    }
    for (int i = 0; i < 3; ++i) {
      if (f.v[i] == a || f.v[i] == bb) {
        f.v[i] = mi;
      }
    }
  }
  // Reconnect the two link vertices across the collapsed edge.
  m.faces.push_back(MakeFace(m, mi, links[0], links[1]));

  m.verts[a].alive = false;
  m.verts[bb].alive = false;
  RebuildTopology(m);
  for (RFace& f : m.faces) {
    if (f.alive) {
      FaceMetrics(m, f);
    }
  }
  return true;
}

struct MeshStats {
  double min_angle;
  double max_edge;
  double min_edge;
  double total_area;
  std::size_t alive_faces;
};

inline MeshStats ComputeStats(const RMesh& m)
{
  MeshStats s;
  s.min_angle = 180.0;
  s.max_edge = 0.0;
  s.min_edge = std::numeric_limits<double>::infinity();
  s.total_area = 0.0;
  s.alive_faces = 0;
  for (const RFace& f : m.faces) {
    if (!f.alive) {
      continue;
    }
    ++s.alive_faces;
    s.min_angle = std::min(s.min_angle, f.min_angle_deg);
    s.max_edge = std::max(s.max_edge, f.longest);
    s.min_edge = std::min(s.min_edge, f.shortest);
    s.total_area += FaceArea(m, f);
  }
  return s;
}

inline double MedianEdgeLength(const RMesh& m)
{
  std::vector<double> lengths;
  lengths.reserve(m.edges.size());
  for (const auto& kv : m.edges) {
    const int a = KeyHi(kv.first);
    const int b = KeyLo(kv.first);
    lengths.push_back(EdgLen(m, a, b));
  }
  if (lengths.empty()) {
    return 1.0;
  }
  std::sort(lengths.begin(), lengths.end());
  return lengths[lengths.size() / 2];
}

} // namespace remesh_detail

P2T_DLL_SYMBOL inline RemeshResult
QualityRemesh(const std::vector<Triangle*>& triangles,
              const RemeshOptions& options = RemeshOptions())
{
  using namespace remesh_detail;

  const auto wall_start = std::chrono::steady_clock::now();

  RemeshResult result;
  result.input_triangles = triangles.size();

  RMesh m;

  // --- Copy the CDT output read-only into a private mesh. -----------------
  std::unordered_map<const Point*, int> index;
  auto intern = [&](const Point* p) {
    const auto found = index.find(p);
    if (found != index.end()) {
      return found->second;
    }
    const int id = static_cast<int>(m.verts.size());
    m.verts.push_back(RVertex{ p->x, p->y, true });
    index.emplace(p, id);
    return id;
  };

  double minx = std::numeric_limits<double>::infinity();
  double miny = std::numeric_limits<double>::infinity();
  double maxx = -std::numeric_limits<double>::infinity();
  double maxy = -std::numeric_limits<double>::infinity();

  for (Triangle* t : triangles) {
    int ids[3];
    for (int i = 0; i < 3; ++i) {
      Point* p = t->GetPoint(i);
      ids[i] = intern(p);
      minx = std::min(minx, p->x);
      miny = std::min(miny, p->y);
      maxx = std::max(maxx, p->x);
      maxy = std::max(maxy, p->y);
    }
    RFace f;
    f.alive = true;
    for (int i = 0; i < 3; ++i) {
      f.v[i] = ids[i];
    }
    m.faces.push_back(f);
    // constrained_edge[i] is the edge opposite vertex i.
    for (int i = 0; i < 3; ++i) {
      if (t->constrained_edge[i]) {
        const Point* p = t->GetPoint((i + 1) % 3);
        const Point* q = t->GetPoint((i + 2) % 3);
        m.constrained.insert(EdgeKey(index[p], index[q]));
      }
    }
  }
  m.scale = std::hypot(maxx - minx, maxy - miny);

  // Orient every face counter-clockwise and compute initial metrics.
  for (RFace& f : m.faces) {
    if (Orient(m.verts[f.v[0]].x, m.verts[f.v[0]].y,
               m.verts[f.v[1]].x, m.verts[f.v[1]].y,
               m.verts[f.v[2]].x, m.verts[f.v[2]].y) < 0.0) {
      std::swap(f.v[1], f.v[2]);
    }
    FaceMetrics(m, f);
  }
  RebuildTopology(m);
  // Edges with a single incident face are domain boundary edges and are
  // constraints even if the source flags did not mark them.
  for (const auto& kv : m.edges) {
    if (kv.second.size() == 1) {
      m.constrained.insert(kv.first);
    }
  }
  for (uint64_t key : m.constrained) {
    const int a = KeyHi(key);
    const int b = KeyLo(key);
    m.constraints.push_back(
        RSegment{ m.verts[a].x, m.verts[a].y, m.verts[b].x, m.verts[b].y });
  }
  result.constraints_total = m.constraints.size();
  RebuildTopology(m);

  const double l_max = options.max_edge_length > 0.0 ?
      options.max_edge_length : MedianEdgeLength(m);
  const double l_min = options.min_edge_length > 0.0 ?
      options.min_edge_length : 0.5 * l_max;
  const double l_floor = 1e-3 * std::max(1.0, m.scale);

  const MeshStats before = ComputeStats(m);
  result.min_angle_before = before.min_angle;
  result.max_edge_before = before.max_edge;
  result.min_edge_before = before.min_edge;
  result.area_before = before.total_area;

  // --- Quality-driven split/collapse loop. --------------------------------
  // Topologically impossible operations are blacklisted permanently. Edges
  // that merely fail the quality guard go into a per-pass retry set: they
  // become eligible again once another operation changes the neighborhood,
  // which keeps the refinement from stalling.
  std::unordered_set<uint64_t> blacklist;
  std::unordered_set<uint64_t> retry;
  const double target = options.target_min_angle_deg;
  const double ratio_limit = options.max_edge_ratio;
  const double angle_guard = 1e-9; // degrees

  // Full post-operation acceptance: topology invariants plus the *local*
  // quality contract. A global non-regression rule would freeze refinement
  // whenever a different, untouched triangle owns the current global worst
  // angle; instead edge splits/collapses must not create a child whose
  // minimum angle is below the worst angle of the faces they replace.
  // Circumcenter (Ruppert) insertions are exempt from the angle guard:
  // their first children can legitimately be more skinny and only converge
  // to the target through the following insertion cascade.
  auto acceptable = [&](const std::vector<int>& replaced_faces,
                        double prev_min_angle, bool is_split,
                        bool is_steiner) {
    if (!options.check_each_operation) {
      return true;
    }
    ++result.checks_run;
    if (!Validate(m)) {
      return false;
    }
    if (!is_steiner) {
      double prev_worst = 180.0;
      for (const int fi : replaced_faces) {
        prev_worst = std::min(prev_worst, m.faces[fi].min_angle_deg);
      }
      // Find the new faces occupying the replaced region by checking every
      // face whose vertices include the newest vertex.
      double new_worst = 180.0;
      const int newest = static_cast<int>(m.verts.size()) - 1;
      for (const RFace& f : m.faces) {
        if (!f.alive || (f.v[0] != newest && f.v[1] != newest &&
                         f.v[2] != newest)) {
          continue;
        }
        new_worst = std::min(new_worst, f.min_angle_deg);
      }
      const double floor_angle = is_split ? prev_worst - angle_guard
                                          : prev_min_angle - angle_guard;
      if (new_worst + angle_guard < floor_angle) {
        return false;
      }
    }
    if (ComputeStats(m).min_edge < 0.5 * l_floor) {
      return false;
    }
    return true;
  };

  // A run of rejected operations with no success means the mesh has
  // converged: every remaining candidate is blocked by a constraint or the
  // quality guard, so further iterations cannot change it.
  std::size_t stall = 0;
  const std::size_t kMaxStall = 64;
  for (std::size_t op = 0; op < options.max_operations; ++op) {
    if (stall >= kMaxStall) {
      break;
    }
    RebuildTopology(m);
    for (RFace& f : m.faces) {
      if (f.alive) {
        FaceMetrics(m, f);
      }
    }

    // Locate the worst bad triangle (smallest minimum angle, then largest
    // edge ratio). Bad = small angle OR edge ratio.
    int worst_face = -1;
    double worst_angle = target;
    double worst_ratio = -1.0;
    for (int fi = 0; fi < static_cast<int>(m.faces.size()); ++fi) {
      const RFace& f = m.faces[fi];
      if (!f.alive) {
        continue;
      }
      const double ratio = f.longest / std::max(f.shortest, 1e-300);
      const bool bad = f.min_angle_deg < target || ratio > ratio_limit;
      if (bad && (f.min_angle_deg < worst_angle ||
                  (f.min_angle_deg <= worst_angle && ratio > worst_ratio))) {
        worst_angle = f.min_angle_deg;
        worst_ratio = ratio;
        worst_face = fi;
      }
    }

    // Priority 1: a fully interior bad face gets a Ruppert circumcenter
    // insertion. This is the only operation that can raise the minimum
    // angle of an interior skinny triangle; constrained edges are never
    // touched by it.
    if (worst_face >= 0) {
      const RFace& f = m.faces[worst_face];
      bool interior_face = true;
      for (int i = 0; i < 3; ++i) {
        if (m.constrained.count(EdgeKey(f.v[i], f.v[(i + 1) % 3])) != 0) {
          interior_face = false;
        }
      }
      if (interior_face && f.min_angle_deg < target) {
        const Snapshot snap = TakeSnapshot(m);
        const double prev_min_angle = ComputeStats(m).min_angle;
        const std::vector<int> replaced{ worst_face };
        const int new_vertex = InsertCircumcenter(m, worst_face);
        bool ok = new_vertex >= 0;
        if (ok) {
          // Restore the Delaunay property locally before validating; this is
          // what lets circumcenter insertion actually raise the min angle.
          LegalizeAround(m, new_vertex, 64);
          ok = acceptable(replaced, prev_min_angle, true, true);
        }
        if (ok) {
          ++result.steiner_count;
          retry.clear();
          stall = 0;
          continue;
        }
        Restore(m, snap);
        ++result.rejected_count;
        ++stall;
        // Fall through to edge splitting for this face below.
      }

      // Priority 2: split the longest edge of the worst face. Works on
      // boundary/constrained faces (midpoint) and as a fallback when the
      // circumcenter insertion was not applicable.
      uint64_t face_split_key = 0;
      double face_split_len = -1.0;
      for (int i = 0; i < 3; ++i) {
        const uint64_t key = EdgeKey(f.v[i], f.v[(i + 1) % 3]);
        if (blacklist.count(key) != 0 || retry.count(key) != 0) {
          continue;
        }
        const int ea = KeyHi(key);
        const int eb = KeyLo(key);
        const double len = EdgLen(m, ea, eb);
        if (len > face_split_len) {
          face_split_len = len;
          face_split_key = key;
        }
      }
      if (face_split_len >= 0.0 && 0.5 * face_split_len >= l_floor) {
        const Snapshot snap = TakeSnapshot(m);
        const double prev_min_angle = ComputeStats(m).min_angle;
        const std::vector<int> replaced =
            m.edges.find(face_split_key)->second;
        const RVertex split_point = ChooseSplitPoint(m, face_split_key);
        const bool applied = SplitEdge(m, face_split_key, split_point);
        const bool ok = applied &&
                        acceptable(replaced, prev_min_angle, true, false);
        if (ok) {
          ++result.split_count;
          retry.clear();
          stall = 0;
          continue;
        }
        Restore(m, snap);
        if (applied) {
          retry.insert(face_split_key);
        } else {
          blacklist.insert(face_split_key);
        }
        ++result.rejected_count;
      }
    }

    // Priority 3: longest edge over the global maximum (size control).
    uint64_t split_key = 0;
    double split_len = -1.0;
    for (const auto& kv : m.edges) {
      if (blacklist.count(kv.first) != 0 || retry.count(kv.first) != 0) {
        continue;
      }
      const int a = KeyHi(kv.first);
      const int b = KeyLo(kv.first);
      const double len = EdgLen(m, a, b);
      if (len > l_max && len > split_len) {
        split_len = len;
        split_key = kv.first;
      }
    }

    if (split_len >= 0.0 && 0.5 * split_len >= l_floor) {
      const Snapshot snap = TakeSnapshot(m);
      const double prev_min_angle = ComputeStats(m).min_angle;
      const std::vector<int> replaced = m.edges.find(split_key)->second;
      const RVertex split_point = ChooseSplitPoint(m, split_key);
      const bool applied = SplitEdge(m, split_key, split_point);
      const bool ok = applied &&
                      acceptable(replaced, prev_min_angle, true, false);
      if (!ok) {
        Restore(m, snap);
        // A topology failure can never become valid; a quality failure is
        // only deferred to the next pass.
        if (applied) {
          retry.insert(split_key);
        } else {
          blacklist.insert(split_key);
        }
        ++result.rejected_count;
        ++stall;
      } else {
        ++result.split_count;
        retry.clear(); // neighborhood changed; deferred edges re-eligible
        stall = 0;
      }
      continue;
    }

    // Candidate 3: collapse short interior edges in good-quality regions.
    uint64_t collapse_key = 0;
    double collapse_len = std::numeric_limits<double>::infinity();
    for (const auto& kv : m.edges) {
      const uint64_t key = kv.first;
      if (blacklist.count(key) != 0 || retry.count(key) != 0 ||
          kv.second.size() != 2 || m.constrained.count(key) != 0) {
        continue;
      }
      const int a = KeyHi(key);
      const int b = KeyLo(key);
      if (m.locked[a] || m.locked[b]) {
        continue;
      }
      bool good_region = true;
      for (const int fi : kv.second) {
        const RFace& f = m.faces[fi];
        const double ratio = f.longest / std::max(f.shortest, 1e-300);
        if (f.min_angle_deg < target || ratio > ratio_limit) {
          good_region = false;
          break;
        }
      }
      if (!good_region) {
        continue;
      }
      const double len = EdgLen(m, a, b);
      if (len < l_min && len < collapse_len) {
        collapse_len = len;
        collapse_key = key;
      }
    }

    if (collapse_len < std::numeric_limits<double>::infinity()) {
      const Snapshot snap = TakeSnapshot(m);
      const double prev_min_angle = ComputeStats(m).min_angle;
      const std::vector<int> replaced =
          m.edges.find(collapse_key)->second;
      const bool applied = CollapseEdge(m, collapse_key);
      const bool ok = applied &&
                      acceptable(replaced, prev_min_angle, false, false);
      if (!ok) {
        Restore(m, snap);
        if (applied) {
          retry.insert(collapse_key);
        } else {
          blacklist.insert(collapse_key);
        }
        ++result.rejected_count;
        ++stall;
      } else {
        ++result.collapse_count;
        retry.clear();
        stall = 0;
      }
      continue;
    }

    // No candidate edge remains: let deferred edges try once more, then
    // stop. The stall counter below also terminates a long run of rejects.
    if (!retry.empty()) {
      retry.clear();
      ++stall;
      if (stall < kMaxStall) {
        continue;
      }
    }
    break; // nothing left to do
  }

  // --- Final validation and report. ---------------------------------------
  RebuildTopology(m);
  for (RFace& f : m.faces) {
    if (f.alive) {
      FaceMetrics(m, f);
    }
  }
  ++result.checks_run;
  result.self_intersecting = HasSelfIntersection(m);
  result.invariant_violated =
      result.self_intersecting || !AllFacesPositive(m) ||
      !ConstraintsCovered(m);

  // Count satisfied constraints one segment at a time via a cheap view.
  result.constraints_satisfied = 0;
  for (const RSegment& s : m.constraints) {
    RMesh view;
    view.verts = m.verts;
    view.edges = m.edges;
    view.constrained = m.constrained;
    view.scale = m.scale;
    view.constraints.push_back(s);
    if (ConstraintsCovered(view)) {
      ++result.constraints_satisfied;
    }
  }

  const MeshStats after = ComputeStats(m);
  result.output_triangles = after.alive_faces;
  result.min_angle_after = after.min_angle;
  result.max_edge_after = after.max_edge;
  result.min_edge_after = after.min_edge;
  result.area_after = after.total_area;
  result.area_error_abs = std::fabs(after.total_area - before.total_area);
  result.area_error_rel = before.total_area > 0.0 ?
      result.area_error_abs / before.total_area : 0.0;
  result.target_angle_reached = after.min_angle >= target;

  // Compact the output to alive vertices.
  std::vector<int> remap(m.verts.size(), -1);
  for (const RFace& f : m.faces) {
    if (!f.alive) {
      continue;
    }
    for (int i = 0; i < 3; ++i) {
      const int old = f.v[i];
      if (remap[old] < 0) {
        remap[old] = static_cast<int>(result.vertices.size());
        result.vertices.push_back(Point(m.verts[old].x, m.verts[old].y));
      }
    }
  }
  for (const RFace& f : m.faces) {
    if (!f.alive) {
      continue;
    }
    RemeshFace out;
    for (int i = 0; i < 3; ++i) {
      out.v[i] = remap[f.v[i]];
    }
    result.faces.push_back(out);
  }

  const auto wall_end = std::chrono::steady_clock::now();
  result.elapsed_seconds =
      std::chrono::duration<double>(wall_end - wall_start).count();
  return result;
}

} // namespace p2t

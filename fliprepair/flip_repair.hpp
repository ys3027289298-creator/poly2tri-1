/*
 * Edge-flip repair for poly2tri triangulations.
 *
 * Use case: triangles adjacent to constrained edges got flipped over, so the
 * triangulation no longer respects its constraints ("形状切破"). This header
 * implements a repair pass that
 *
 *   1. finds edges violating the constraints (mesh edges crossing a required
 *      constrained segment, and interior edges failing the empty-circle
 *      criterion),
 *   2. tries to flip them, guided by the empty-circle (in-circle) criterion,
 *   3. never flips a constrained edge,
 *   4. iterates to a fixpoint or until the round budget is exhausted.
 *
 * Round count and criterion strictness are tunable via Options. The Report
 * always tells the truth: if the round budget ran out before convergence,
 * round_limit_hit is set and the mesh is returned as-is (no silent
 * half-repaired result).
 *
 * Self-contained: uses only the public poly2tri API, no third-party code.
 */
#pragma once

#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace fliprepair {

struct Options {
  /// Maximum number of flip rounds before giving up.
  int max_rounds = 64;
  /// Strictness of the empty-circle criterion: an edge is a Delaunay
  /// violation only if the normalized in-circle determinant exceeds this
  /// value. Larger = stricter (fewer edges flagged, fewer flips accepted).
  double incircle_epsilon = 1e-9;
};

struct Report {
  int rounds_used = 0;
  long flips_applied = 0;
  /// No violating edges remain.
  bool converged = false;
  /// A full round applied zero flips (no further progress possible).
  bool fixpoint_reached = false;
  /// Round budget exhausted without converging. The mesh is NOT fully
  /// repaired in this case and the caller must not treat it as such.
  bool round_limit_hit = false;
  /// Violations still present after the last round.
  int remaining_violations = 0;
  double min_angle_before_deg = 0.0;
  double min_angle_after_deg = 0.0;
  /// All required constrained segments exist as mesh edges afterwards.
  bool constraints_intact = false;
  int missing_constraints = 0;
};

using EdgeKey = std::pair<p2t::Point*, p2t::Point*>;

inline EdgeKey MakeEdgeKey(p2t::Point* a, p2t::Point* b)
{
  return a < b ? EdgeKey(a, b) : EdgeKey(b, a);
}

/// 2D cross product of (b - a) x (c - a).
inline double Cross2(const p2t::Point& a, const p2t::Point& b, const p2t::Point& c)
{
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

/// Do segments (a,b) and (c,d) properly intersect (cross in their interiors,
/// merely touching endpoints does not count)?
inline bool ProperlyIntersect(const p2t::Point& a, const p2t::Point& b,
                              const p2t::Point& c, const p2t::Point& d)
{
  const double abc = Cross2(a, b, c);
  const double abd = Cross2(a, b, d);
  const double cda = Cross2(c, d, a);
  const double cdb = Cross2(c, d, b);
  return ((abc > 0) != (abd > 0)) && abc != 0 && abd != 0 &&
         ((cda > 0) != (cdb > 0)) && cda != 0 && cdb != 0;
}

/// Normalized in-circle determinant of d against the CCW triangle (a,b,c).
/// Positive means d lies inside the circumcircle of (a,b,c). Dimensionless:
/// the raw determinant is divided by the product of the squared distances
/// so that Options::incircle_epsilon is scale independent.
inline double InCircleNormalized(const p2t::Point& a, const p2t::Point& b,
                                 const p2t::Point& c, const p2t::Point& d)
{
  const double adx = a.x - d.x, ady = a.y - d.y;
  const double bdx = b.x - d.x, bdy = b.y - d.y;
  const double cdx = c.x - d.x, cdy = c.y - d.y;

  const double alift = adx * adx + ady * ady;
  const double blift = bdx * bdx + bdy * bdy;
  const double clift = cdx * cdx + cdy * cdy;

  const double det = alift * (bdx * cdy - cdx * bdy) +
                     blift * (cdx * ady - adx * cdy) +
                     clift * (adx * bdy - bdx * ady);
  const double denom = alift * blift * clift;
  if (denom <= 0) {
    return 0; // degenerate configuration, treat as not violating
  }
  return det / denom;
}

/// Smallest interior angle over all triangles, in degrees.
inline double MinAngleDeg(const std::vector<p2t::Triangle*>& triangles)
{
  double best = 1e30;
  for (auto* t : triangles) {
    for (int i = 0; i < 3; ++i) {
      p2t::Point* p = t->GetPoint(i);
      p2t::Point* u = t->GetPoint((i + 1) % 3);
      p2t::Point* v = t->GetPoint((i + 2) % 3);
      const double ux = u->x - p->x, uy = u->y - p->y;
      const double vx = v->x - p->x, vy = v->y - p->y;
      const double lu = std::sqrt(ux * ux + uy * uy);
      const double lv = std::sqrt(vx * vx + vy * vy);
      if (lu == 0 || lv == 0) {
        return 0; // degenerate triangle
      }
      double cosang = (ux * vx + uy * vy) / (lu * lv);
      cosang = std::max(-1.0, std::min(1.0, cosang));
      best = std::min(best, std::acos(cosang) * 180.0 / M_PI);
    }
  }
  return best;
}

struct EdgeInfo {
  p2t::Triangle* tri[2] = {nullptr, nullptr};
  int count = 0;
};

inline std::map<EdgeKey, EdgeInfo> BuildEdgeMap(
    const std::vector<p2t::Triangle*>& triangles)
{
  std::map<EdgeKey, EdgeInfo> edges;
  for (auto* t : triangles) {
    for (int i = 0; i < 3; ++i) {
      p2t::Point* p1 = t->GetPoint((i + 1) % 3);
      p2t::Point* p2 = t->GetPoint((i + 2) % 3);
      EdgeInfo& info = edges[MakeEdgeKey(p1, p2)];
      if (info.count < 2) {
        info.tri[info.count] = t;
      }
      ++info.count;
    }
  }
  return edges;
}

/// Rotate the triangle pair sharing the edge opposite p in t and opposite op
/// in ot, replacing that edge with the (p, op) diagonal. Mirrors
/// Sweep::RotateTrianglePair, restricted to the public Triangle API.
inline void FlipEdge(p2t::Triangle& t, p2t::Point& p, p2t::Triangle& ot, p2t::Point& op)
{
  p2t::Triangle* n1 = t.NeighborCCW(p);
  p2t::Triangle* n2 = t.NeighborCW(p);
  p2t::Triangle* n3 = ot.NeighborCCW(op);
  p2t::Triangle* n4 = ot.NeighborCW(op);

  const bool ce1 = t.GetConstrainedEdgeCCW(p);
  const bool ce2 = t.GetConstrainedEdgeCW(p);
  const bool ce3 = ot.GetConstrainedEdgeCCW(op);
  const bool ce4 = ot.GetConstrainedEdgeCW(op);

  const bool de1 = t.GetDelunayEdgeCCW(p);
  const bool de2 = t.GetDelunayEdgeCW(p);
  const bool de3 = ot.GetDelunayEdgeCCW(op);
  const bool de4 = ot.GetDelunayEdgeCW(op);

  t.Legalize(p, op);
  ot.Legalize(op, p);

  ot.SetDelunayEdgeCCW(p, de1);
  t.SetDelunayEdgeCW(p, de2);
  t.SetDelunayEdgeCCW(op, de3);
  ot.SetDelunayEdgeCW(op, de4);

  ot.SetConstrainedEdgeCCW(p, ce1);
  t.SetConstrainedEdgeCW(p, ce2);
  t.SetConstrainedEdgeCCW(op, ce3);
  ot.SetConstrainedEdgeCW(op, ce4);

  t.ClearNeighbors();
  ot.ClearNeighbors();
  if (n1) {
    ot.MarkNeighbor(*n1);
  }
  if (n2) {
    t.MarkNeighbor(*n2);
  }
  if (n3) {
    t.MarkNeighbor(*n3);
  }
  if (n4) {
    ot.MarkNeighbor(*n4);
  }
  t.MarkNeighbor(ot);
}

namespace detail {

struct EdgeViolation {
  EdgeKey key;
  bool crosses_constraint = false;
  bool violates_empty_circle = false;
};

/// Count how many required constrained segments properly cross segment (a,b).
inline int CountConstraintCrossings(
    const p2t::Point& a, const p2t::Point& b,
    const std::set<EdgeKey>& required)
{
  int n = 0;
  for (const auto& seg : required) {
    if (ProperlyIntersect(a, b, *seg.first, *seg.second)) {
      ++n;
    }
  }
  return n;
}

/// In-circle margin of the edge shared by the quad with diagonal (a,b) and
/// opposite vertices c, d. Positive (beyond epsilon) means the edge should
/// be flipped per the empty-circle criterion.
inline double EdgeInCircleMargin(const p2t::Point& a, const p2t::Point& b,
                                 const p2t::Point& c, const p2t::Point& d)
{
  if (Cross2(a, b, c) > 0) {
    return InCircleNormalized(a, b, c, d);
  }
  return InCircleNormalized(b, a, c, d);
}

} // namespace detail

/**
 * Repair a triangulation whose constrained edges were flipped over.
 *
 * @param triangles           interior triangles of the mesh (mutated in place)
 * @param required_constraints segments that must exist as mesh edges and are
 *                            never allowed to be flipped
 * @param options             round budget and criterion strictness
 */
inline Report Repair(std::vector<p2t::Triangle*>& triangles,
                     const std::vector<EdgeKey>& required_constraints,
                     const Options& options)
{
  Report report;
  report.min_angle_before_deg = MinAngleDeg(triangles);

  std::set<EdgeKey> required;
  for (const auto& seg : required_constraints) {
    required.insert(MakeEdgeKey(seg.first, seg.second));
  }

  struct Scan {
    std::vector<detail::EdgeViolation> violations;
    int missing_constraints = 0;
  };

  // Detect all edges currently violating a constraint, plus the required
  // constrained segments missing from the mesh. A missing segment is a
  // violation in its own right, even when no flippable edge crosses it
  // (e.g. the only crossing edges are themselves constrained).
  const auto scan = [&]() {
    Scan result;
    const auto edges = BuildEdgeMap(triangles);
    std::set<EdgeKey> missing;
    for (const auto& seg : required) {
      if (edges.find(seg) == edges.end()) {
        missing.insert(seg);
      }
    }
    result.missing_constraints = static_cast<int>(missing.size());
    for (const auto& entry : edges) {
      const EdgeInfo& info = entry.second;
      if (info.count != 2 || required.count(entry.first)) {
        continue; // boundary edge or constrained edge: never a flip candidate
      }
      p2t::Point* a = entry.first.first;
      p2t::Point* b = entry.first.second;
      detail::EdgeViolation v;
      v.key = entry.first;
      for (const auto& seg : missing) {
        if (ProperlyIntersect(*a, *b, *seg.first, *seg.second)) {
          v.crosses_constraint = true;
          break;
        }
      }
      p2t::Point* c = info.tri[0]->GetPoint(info.tri[0]->EdgeIndex(a, b));
      p2t::Point* d = info.tri[1]->GetPoint(info.tri[1]->EdgeIndex(a, b));
      if (detail::EdgeInCircleMargin(*a, *b, *c, *d) > options.incircle_epsilon) {
        v.violates_empty_circle = true;
      }
      if (v.crosses_constraint || v.violates_empty_circle) {
        result.violations.push_back(v);
      }
    }
    return result;
  };

  while (report.rounds_used < options.max_rounds) {
    const Scan current = scan();
    if (current.violations.empty() && current.missing_constraints == 0) {
      report.converged = true;
      break;
    }
    if (current.violations.empty()) {
      // Constraints are missing but nothing is flippable: stuck.
      report.fixpoint_reached = true;
      break;
    }
    ++report.rounds_used;
    long flips_this_round = 0;
    for (const auto& violation : current.violations) {
      // Rebuild: earlier flips this round may have changed the neighborhood.
      const auto edges = BuildEdgeMap(triangles);
      const auto it = edges.find(violation.key);
      if (it == edges.end() || it->second.count != 2) {
        continue; // edge no longer interior
      }
      if (required.count(it->first)) {
        continue; // constrained edges are never flipped
      }
      p2t::Point* a = it->first.first;
      p2t::Point* b = it->first.second;
      p2t::Triangle* t0 = it->second.tri[0];
      p2t::Triangle* t1 = it->second.tri[1];
      p2t::Point* c = t0->GetPoint(t0->EdgeIndex(a, b));
      p2t::Point* d = t1->GetPoint(t1->EdgeIndex(a, b));

      // The quad must be convex, otherwise the flip is not a valid move.
      if (!ProperlyIntersect(*a, *b, *c, *d)) {
        continue;
      }

      const int old_crossings = detail::CountConstraintCrossings(*a, *b, required);
      const int new_crossings = detail::CountConstraintCrossings(*c, *d, required);
      const double margin_before = detail::EdgeInCircleMargin(*a, *b, *c, *d);
      const double margin_after = detail::EdgeInCircleMargin(*c, *d, *a, *b);

      // Accept the flip if it removes a constraint crossing, or if the
      // empty-circle criterion says the edge is illegal and flipping strictly
      // improves it (without introducing new constraint crossings).
      const bool recovers_constraint = new_crossings < old_crossings;
      const bool improves_empty_circle =
          new_crossings == old_crossings &&
          margin_before > options.incircle_epsilon && margin_after < margin_before;
      if (!recovers_constraint && !improves_empty_circle) {
        continue;
      }

      FlipEdge(*t0, *c, *t1, *d);
      ++flips_this_round;
      ++report.flips_applied;
    }
    if (flips_this_round == 0) {
      report.fixpoint_reached = true;
      break;
    }
  }

  // Final state: count what is still wrong, honestly.
  const Scan final_scan = scan();
  report.remaining_violations =
      static_cast<int>(final_scan.violations.size()) + final_scan.missing_constraints;
  if (!report.converged && !report.fixpoint_reached) {
    report.round_limit_hit = true;
  }

  // Reconcile constrained-edge flags with the required segment set (the
  // corruption may have scattered stale flags onto unrelated edges).
  const auto edges = BuildEdgeMap(triangles);
  for (auto* t : triangles) {
    for (int i = 0; i < 3; ++i) {
      p2t::Point* p1 = t->GetPoint((i + 1) % 3);
      p2t::Point* p2 = t->GetPoint((i + 2) % 3);
      t->constrained_edge[i] = required.count(MakeEdgeKey(p1, p2)) != 0;
    }
  }
  for (const auto& seg : required) {
    if (edges.find(seg) == edges.end()) {
      ++report.missing_constraints;
    }
  }
  report.constraints_intact = report.missing_constraints == 0;
  report.min_angle_after_deg = MinAngleDeg(triangles);
  return report;
}

} // namespace fliprepair

/*
 * Poly2Tri Copyright (c) 2009-2022, Poly2Tri Contributors
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

#include "../common/self_intersection.h"
#include <list>
#include <vector>
#include <cstddef>
#include <unordered_set>

namespace p2t {

// Inital triangle factor, seed triangle will extend 30% of
// PointSet width to both left and right.
const double kAlpha = 0.3;

struct Point;
class Triangle;
struct Node;
struct Edge;
class AdvancingFront;

class SweepContext {
public:

/// Constructor
explicit SweepContext(std::vector<Point*> polyline);
/// Destructor
~SweepContext();

/// Select how a non-simple contour is handled. Existing callers keep the
/// default reject behaviour and identical output for simple polygons.
void set_self_intersection_policy(SelfIntersectionPolicy policy)
{
  self_intersection_policy_ = policy;
}

SelfIntersectionPolicy self_intersection_policy() const
{
  return self_intersection_policy_;
}

void set_head(Point* p1);

Point* head() const;

void set_tail(Point* p1);

Point* tail() const;

size_t point_count() const;

Node* LocateNode(const Point& point);

void RemoveNode(Node* node);

void CreateAdvancingFront();

/// Try to map a node to all sides of this triangle that don't have a neighbor
void MapTriangleToNodes(Triangle& t);

void AddToMap(Triangle* triangle);

Point* GetPoint(size_t index);

Point* GetPoints();

void RemoveFromMap(Triangle* triangle);

void AddHole(const std::vector<Point*>& polyline);

void AddPoint(Point* point);

AdvancingFront* front() const;

void MeshClean(Triangle& triangle);

std::vector<Triangle*> &GetTriangles();
std::list<Triangle*> &GetMap();

/// Runs the pre-sweep contour validation and, when requested, resolves the
/// even-odd fill loops. Returns true when the regular sweep should run on
/// the original points; false when an already triangulated result is stored.
bool ValidateOrResolve();

/// Takes ownership of a triangle created while triangulating resolved loops.
void AdoptTriangle(Triangle* triangle);

/// Takes ownership of a point created while splitting intersection edges.
void AdoptPoint(Point* point);

std::vector<Edge*> edge_list;

struct Basin {
  Node* left_node;
  Node* bottom_node;
  Node* right_node;
  double width;
  bool left_highest;

  Basin()
  : left_node(nullptr), bottom_node(nullptr), right_node(nullptr), width(0.0), left_highest(false)
  {
  }

  void Clear()
  {
    left_node = nullptr;
    bottom_node = nullptr;
    right_node = nullptr;
    width = 0.0;
    left_highest = false;
  }
};

struct EdgeEvent {
  Edge* constrained_edge;
  bool right;

  EdgeEvent() : constrained_edge(NULL), right(false)
  {
  }
};

Basin basin;
EdgeEvent edge_event;

private:

friend class Sweep;

std::vector<Triangle*> triangles_;
std::list<Triangle*> map_;
std::vector<Point*> points_;
std::unordered_set<Point*> point_set_;

// Advancing front
AdvancingFront* front_;
// head point used with advancing front
Point* head_;
// tail point used with advancing front
Point* tail_;

Node *af_head_, *af_middle_, *af_tail_;

/// Original contours, kept for validation and even-odd resolution.
std::vector<std::vector<Point*> > contours_;
/// Triangles created for resolved loops and owned by this context.
std::vector<Triangle*> adopted_triangles_;
/// Points created for intersection splits and owned by this context.
std::vector<Point*> adopted_points_;
SelfIntersectionPolicy self_intersection_policy_ = REJECT_SELF_INTERSECTION;

void InitTriangulation();
void InitEdges(const std::vector<Point*>& polyline);

};

inline AdvancingFront* SweepContext::front() const
{
  return front_;
}

inline size_t SweepContext::point_count() const
{
  return points_.size();
}

inline void SweepContext::set_head(Point* p1)
{
  head_ = p1;
}

inline Point* SweepContext::head() const
{
  return head_;
}

inline void SweepContext::set_tail(Point* p1)
{
  tail_ = p1;
}

inline Point* SweepContext::tail() const
{
  return tail_;
}

}

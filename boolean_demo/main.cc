// Boolean approximation demo on top of p2t::CDT.
//
// Four small hand-made shape pairs; for each pair it prints the
// intersection, union and A-minus-B difference: kept-triangle count,
// merged-outline point count, resulting area, the area-conservation error
// against the identity built from the two triangulated input areas
// (union: A+B-I, difference: A-I), and the error against the exact
// analytic area of the operation. A final case shows a self-intersecting
// input being rejected through the library's std::runtime_error error path.
//
// Build (uses only the library's own sources, no extra dependencies):
//   g++ -std=c++14 -I.. main.cc ../poly2tri/common/shapes.cc
//       ../poly2tri/sweep/*.cc -o boolean_demo

#include <poly2tri/poly2tri.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using p2t::Point;
using p2t::boolean_approx::Op;
using p2t::boolean_approx::Region;
using p2t::boolean_approx::Result;

namespace ba = p2t::boolean_approx;

namespace {

// Steiner spacing for the interior refinement grid; small enough that the
// centroid test resolves the narrowest overlap in these pairs.
const double kSpacing = 0.25;

Region MakeRegion(std::vector<Point> outer,
                  std::vector<std::vector<Point>> holes = {})
{
  return Region{std::move(outer), std::move(holes)};
}

double TriangulatedArea(const Region& r)
{
  double area = 0.0;
  for (const auto& t : ba::Triangulate(r, kSpacing)) {
    area += t.Area();
  }
  return area;
}

size_t OutlinePointCount(const Result& r)
{
  size_t n = 0;
  for (const auto& loop : r.contours) {
    n += loop.size();
  }
  return n;
}

const char* OpName(Op op)
{
  switch (op) {
  case Op::Intersection:
    return "intersection";
  case Op::Union:
    return "union       ";
  case Op::Difference:
    return "difference  ";
  }
  return "?";
}

// exact: analytic area of the operation, used to show the approximation
// error of the centroid test.
void RunPair(const std::string& title, const Region& a, const Region& b,
             double exact_intersection)
{
  const double area_a = TriangulatedArea(a);
  const double area_b = TriangulatedArea(b);

  std::cout << "------------------------------------------------------------\n";
  std::cout << title << "\n";
  std::cout << "  input areas (sum of triangulated triangles): A="
            << area_a << "  B=" << area_b
            << "  exact overlap=" << exact_intersection << "\n";
  std::cout << "  op             triangles  outline pts   area"
            << "        conservation err   exact-area err\n";

  const double inter = ba::Apply(a, b, Op::Intersection, kSpacing).area;
  // Conservation identities built from the same approximate intersection
  // returned by the centroid test. They measure the triangle bookkeeping;
  // the exact-area column measures the approximation against analytic areas.
  const double identity[3] = {
      inter,                       // intersection: its own area
      area_a + area_b - inter,     // union:  |A u B| = |A| + |B| - |A n B|
      area_a - inter               // difference: |A \\ B| = |A| - |A n B|
  };
  const double exact[3] = {exact_intersection,
                           area_a + area_b - exact_intersection,
                           area_a - exact_intersection};
  const Op ops[3] = {Op::Intersection, Op::Union, Op::Difference};

  for (int i = 0; i < 3; ++i) {
    const Result r = ba::Apply(a, b, ops[i], kSpacing);
    const double cons_err = std::fabs(r.area - identity[i]);
    const double exact_err = std::fabs(r.area - exact[i]);
    std::cout << "  " << OpName(ops[i]) << "  " << std::setw(9)
              << r.triangles.size() << std::setw(13) << OutlinePointCount(r)
              << std::setw(11) << std::fixed << std::setprecision(6) << r.area
              << "   " << std::scientific << std::setprecision(3)
              << std::setw(10) << cons_err << "      " << std::setw(10)
              << exact_err << std::defaultfloat << "\n";
  }
}

} // namespace

int main()
{
  // Pair 1: two overlapping axis-aligned squares (integer-grid shapes).
  // A = [0,4]^2 (area 16), B = [3,7]x[1,5] (area 16), overlap = 1x3 = 3.
  const Region square_a = MakeRegion({
      Point(0, 0), Point(4, 0), Point(4, 4), Point(0, 4)});
  const Region square_b = MakeRegion({
      Point(3, 1), Point(7, 1), Point(7, 5), Point(3, 5)});
  RunPair("Pair 1: overlapping axis-aligned squares, 16 and 16", square_a,
          square_b, 3.0);

  // Pair 2: concave L-shape (area 12) overlapping a square (area 4);
  // B straddles the L's missing corner; exact overlap area is 3.
  const Region l_shape = MakeRegion({
      Point(0, 0), Point(4, 0), Point(4, 2), Point(2, 2),
      Point(2, 4), Point(0, 4)});
  const Region inner_square = MakeRegion({
      Point(1, 1), Point(3, 1), Point(3, 3), Point(1, 3)});
  RunPair("Pair 2: concave L-shape (12) and square (4) over its corner",
          l_shape, inner_square, 3.0);

  // Pair 3: square with a square hole vs. a small square sitting in the hole.
  // Outer 6x6 (36) minus 2x2 hole (4) = 32; B has area 1 and no overlap.
  const Region holed = MakeRegion(
      {Point(0, 0), Point(6, 0), Point(6, 6), Point(0, 6)},
      {{Point(2, 2), Point(4, 2), Point(4, 4), Point(2, 4)}});
  const Region in_hole = MakeRegion({
      Point(2.5, 2.5), Point(3.5, 2.5), Point(3.5, 3.5), Point(2.5, 3.5)});
  RunPair("Pair 3: holed square (32) and small square (1) inside its hole",
          holed, in_hole, 0.0);

  // Pair 4: two right triangles of area 8 whose overlap is the axis-aligned
  // triangle with vertices (0,0),(2,0),(0,2), area 4 (half of each).
  const Region tri_a = MakeRegion({
      Point(0, 0), Point(4, 0), Point(0, 4)});
  const Region tri_b = MakeRegion({
      Point(0, 0), Point(4, 4), Point(0, 4)});
  RunPair("Pair 4: two overlapping right triangles (8 and 8, overlap 4)",
          tri_a, tri_b, 4.0);

  // Self-intersecting input must take the existing error path.
  std::cout << "------------------------------------------------------------\n";
  std::cout << "Error path: bow-tie ring (0,0)->(4,4)->(0,4)->(4,0)\n";
  const Region bowtie = MakeRegion({
      Point(0, 0), Point(4, 4), Point(0, 4), Point(4, 0)});
  try {
    ba::Apply(bowtie, square_a, Op::Union, kSpacing);
    std::cout << "  UNEXPECTED: no exception was thrown\n";
    return 1;
  } catch (const std::runtime_error& e) {
    std::cout << "  rejected with std::runtime_error: \"" << e.what() << "\"\n";
  }
  std::cout << "------------------------------------------------------------\n";
  std::cout << "All operations completed.\n";
  return 0;
}

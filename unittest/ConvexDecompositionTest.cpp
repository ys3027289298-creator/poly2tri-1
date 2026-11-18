/*
 * Convex decomposition tests.
 *
 * These drive the whole three-stage pipeline on small hand-built polygons:
 * stage 1 (ear clipping) and stage 2 (bridging) are observed through the
 * result counters, stage 3 (area validation) through the reported area error.
 * The key invariant is that the pieces really tile the input: every piece is
 * convex, their areas add up to the polygon area, and sampled points land in a
 * piece exactly when they land in the polygon.
 *
 * The file follows the same Boost.Test layout as the other unit tests.  It is
 * built by tools/run_convex_decomposition_tests.sh so that the existing build
 * scripts stay untouched.
 */

#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif

#include <boost/test/unit_test.hpp>
#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

struct OwnedPolygon {
  std::vector<p2t::Point*> outer;
  std::vector<std::vector<p2t::Point*>> holes;
  std::vector<std::unique_ptr<p2t::Point>> storage;

  p2t::Point* add(double x, double y)
  {
    storage.emplace_back(new p2t::Point(x, y));
    return storage.back().get();
  }

  void ring(std::vector<p2t::Point*>& out, const std::vector<std::pair<double, double>>& pts)
  {
    for (const auto& p : pts)
      out.push_back(add(p.first, p.second));
  }
};

double signed_area(const std::vector<p2t::Point*>& ring)
{
  double a = 0.0;
  for (size_t i = 0; i < ring.size(); ++i) {
    const auto& p = *ring[i];
    const auto& q = *ring[(i + 1) % ring.size()];
    a += p.x * q.y - q.x * p.y;
  }
  return a * 0.5;
}

double piece_area(const p2t::ConvexPiece& piece)
{
  const auto& p = piece.points;
  double a = 0.0;
  for (size_t i = 0; i < p.size(); ++i) {
    const auto& q = *p[i];
    const auto& r = *p[(i + 1) % p.size()];
    a += q.x * r.y - r.x * q.y;
  }
  return a * 0.5;
}

bool is_convex(const p2t::ConvexPiece& piece)
{
  const auto& p = piece.points;
  if (p.size() < 3)
    return false;
  for (size_t i = 0; i < p.size(); ++i) {
    const auto& a = *p[i];
    const auto& b = *p[(i + 1) % p.size()];
    const auto& c = *p[(i + 2) % p.size()];
    if ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) < -1e-10)
      return false;
  }
  return true;
}

bool point_in_ring(const std::vector<p2t::Point*>& ring, double x, double y)
{
  bool inside = false;
  for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
    const double xi = ring[i]->x, yi = ring[i]->y;
    const double xj = ring[j]->x, yj = ring[j]->y;
    if (((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi))
      inside = !inside;
  }
  return inside;
}

// Checks the invariants that every successful decomposition must satisfy.
void check_result(const OwnedPolygon& poly, const p2t::ConvexDecompositionResult& result)
{
  BOOST_REQUIRE(!result.pieces.empty());

  double expected = std::fabs(signed_area(poly.outer));
  for (const auto& h : poly.holes)
    expected -= std::fabs(signed_area(h));

  double total = 0.0;
  for (const auto& piece : result.pieces) {
    BOOST_CHECK_MESSAGE(is_convex(piece), "piece is not convex");
    const double area = piece_area(piece);
    BOOST_CHECK_MESSAGE(area > 0.0, "piece does not wind CCW");
    total += area;
  }
  BOOST_CHECK_CLOSE(result.polygon_area, expected, 1e-9);
  BOOST_CHECK_LT(std::fabs(total - expected), 1e-9 * std::max(1.0, std::fabs(expected)));
  BOOST_CHECK_LT(result.area_error, 1e-9 * std::max(1.0, std::fabs(expected)));

  double minx = poly.outer[0]->x, maxx = minx;
  double miny = poly.outer[0]->y, maxy = miny;
  for (const auto* p : poly.outer) {
    minx = std::min(minx, p->x);
    maxx = std::max(maxx, p->x);
    miny = std::min(miny, p->y);
    maxy = std::max(maxy, p->y);
  }
  std::mt19937 rng(2024);
  std::uniform_real_distribution<double> dx(minx - 0.5, maxx + 0.5);
  std::uniform_real_distribution<double> dy(miny - 0.5, maxy + 0.5);

  for (int k = 0; k < 2000; ++k) {
    const double x = dx(rng);
    const double y = dy(rng);
    bool in_polygon = point_in_ring(poly.outer, x, y);
    for (const auto& h : poly.holes)
      if (point_in_ring(h, x, y))
        in_polygon = false;

    int hits = 0;
    for (const auto& piece : result.pieces) {
      const auto& p = piece.points;
      bool inside = true;
      for (size_t i = 0; i < p.size() && inside; ++i) {
        const auto& a = *p[i];
        const auto& b = *p[(i + 1) % p.size()];
        if ((b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x) < -1e-9)
          inside = false;
      }
      if (inside)
        hits++;
    }
    BOOST_CHECK_MESSAGE(hits <= 1, "point covered by more than one piece");
    BOOST_CHECK_MESSAGE(in_polygon == (hits == 1), "piece coverage does not match polygon");
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(ConvexPolygonIsOnePiece)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 0 }, { 4, 0 }, { 5, 2 }, { 4, 4 }, { 0, 4 }, { -1, 2 } });

  const auto result = p2t::DecomposePolygon(poly.outer);
  check_result(poly, result);
  BOOST_CHECK_EQUAL(result.pieces.size(), 1u);
  BOOST_CHECK_EQUAL(result.reflex_vertices, 0u);
}

BOOST_AUTO_TEST_CASE(SingleReflexIsTwoPieces)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 0 }, { 4, 0 }, { 4, 1 }, { 1, 1 }, { 1, 4 }, { 0, 4 } });

  const auto result = p2t::DecomposePolygon(poly.outer);
  check_result(poly, result);
  BOOST_CHECK_EQUAL(result.reflex_vertices, 1u);
  BOOST_CHECK_EQUAL(result.pieces.size(), 2u);
  BOOST_CHECK_GT(result.ears_clipped, 0u);
  BOOST_CHECK_GT(result.bridges, 0u);
  // 6 input vertices; the two pieces share the bridge vertices.
  BOOST_CHECK_GT(result.output_vertices, result.input_vertices);
}

BOOST_AUTO_TEST_CASE(MultipleReflexesUseMorePieces)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 0 },  { 4, 0 },  { 4, 3 },  { 3.2, 3 }, { 3.2, 1.4 },
                          { 2.4, 1.4 }, { 2.4, 3 }, { 1.6, 3 }, { 1.6, 1.4 }, { 0.8, 1.4 },
                          { 0.8, 3 }, { 0, 3 } });

  const auto result = p2t::DecomposePolygon(poly.outer);
  check_result(poly, result);
  BOOST_CHECK_EQUAL(result.reflex_vertices, 4u);
  BOOST_CHECK_GT(result.pieces.size(), 2u);
  BOOST_CHECK_LT(result.pieces.size(), poly.outer.size());
}

BOOST_AUTO_TEST_CASE(HoleIsBridged)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 0 }, { 4, 0 }, { 4, 4 }, { 0, 4 } });
  poly.holes.emplace_back();
  poly.ring(poly.holes.back(), { { 1.5, 1.5 }, { 2.5, 1.5 }, { 2.5, 2.5 }, { 1.5, 2.5 } });

  const auto result = p2t::DecomposePolygon(poly.outer, poly.holes);
  check_result(poly, result);
  BOOST_CHECK_GT(result.pieces.size(), 1u);
  BOOST_CHECK_GT(result.bridges, 0u);
}

BOOST_AUTO_TEST_CASE(CwInputIsAccepted)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 4 }, { 1, 4 }, { 1, 1 }, { 4, 1 }, { 4, 0 }, { 0, 0 } });

  const auto result = p2t::DecomposePolygon(poly.outer);
  check_result(poly, result);
  BOOST_CHECK_EQUAL(result.pieces.size(), 2u);
}

BOOST_AUTO_TEST_CASE(BridgeStrategiesBothWork)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 0 },  { 4, 0 },  { 4, 3 },  { 3.2, 3 }, { 3.2, 1.4 },
                          { 2.4, 1.4 }, { 2.4, 3 }, { 1.6, 3 }, { 1.6, 1.4 }, { 0.8, 1.4 },
                          { 0.8, 3 }, { 0, 3 } });

  p2t::ConvexDecompositionOptions shortest;
  shortest.bridge = p2t::BridgeStrategy::Shortest;
  p2t::ConvexDecompositionOptions nearest;
  nearest.bridge = p2t::BridgeStrategy::NearestIndex;

  const auto a = p2t::DecomposePolygon(poly.outer, {}, shortest);
  const auto b = p2t::DecomposePolygon(poly.outer, {}, nearest);
  check_result(poly, a);
  check_result(poly, b);
}

BOOST_AUTO_TEST_CASE(MinPieceAreaIsHonoured)
{
  OwnedPolygon poly;
  poly.ring(poly.outer, { { 0, 0 }, { 4, 0 }, { 4, 1 }, { 1, 1 }, { 1, 4 }, { 0, 4 } });

  p2t::ConvexDecompositionOptions options;
  options.min_piece_area = 0.5;
  const auto result = p2t::DecomposePolygon(poly.outer, {}, options);
  check_result(poly, result);
  for (const auto& piece : result.pieces)
    BOOST_CHECK_GE(piece_area(piece), 0.5);
}

BOOST_AUTO_TEST_CASE(UndecomposablePolygonsThrow)
{
  {
    // Self-intersecting bowtie: the two lobes cancel, so no valid piece set
    // exists and the area check must reject the input.
    OwnedPolygon poly;
    poly.ring(poly.outer, { { 0, 0 }, { 2, 2 }, { 2, 0 }, { 0, 2 } });
    BOOST_CHECK_THROW(p2t::DecomposePolygon(poly.outer), std::runtime_error);
  }
  {
    // An impossible minimum piece size cannot be satisfied by merging.
    OwnedPolygon poly;
    poly.ring(poly.outer, { { 0, 0 }, { 4, 0 }, { 4, 4 }, { 0, 4 } });
    poly.holes.emplace_back();
    poly.ring(poly.holes.back(), { { 0.5, 0.5 }, { 2.5, 0.5 }, { 2.5, 2.5 }, { 0.5, 2.5 } });
    p2t::ConvexDecompositionOptions options;
    options.min_piece_area = 100.0;
    BOOST_CHECK_THROW(p2t::DecomposePolygon(poly.outer, poly.holes, options),
                      std::runtime_error);
  }
  {
    // Too few vertices.
    OwnedPolygon poly;
    BOOST_CHECK_THROW(p2t::DecomposePolygon(poly.outer), std::runtime_error);
  }
}

#if !defined(_WIN32) && !defined(BOOST_TEST_DYN_LINK)
#define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>
#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <sstream>
#include <vector>

BOOST_AUTO_TEST_CASE(TriangleTest)
{
  p2t::Point a(0, 0);
  p2t::Point b(1, 0);
  p2t::Point c(0.5, .5);
  p2t::Triangle triangle(a, b, c);
  BOOST_CHECK(triangle.Contains(&a));
  BOOST_CHECK(triangle.Contains(&b));
  BOOST_CHECK(triangle.Contains(&c));
  BOOST_CHECK(triangle.CircumcicleContains(p2t::Point(0.5, 0.1)));
  BOOST_CHECK(!triangle.CircumcicleContains(p2t::Point(1, 0.4)));
}

// ===========================================================================
// TriangleMeshPipeline tests
//
// Helpers build small, hand-made triangle soups with deliberately duplicated
// corner vertices, then assert the four required properties:
//   1. vertex count change (input corners -> welded unique vertices),
//   2. strip count and fan count,
//   3. vertex reuse rate,
//   4. reconstruction identity (decoded triangle multiset == input multiset).
// ===========================================================================

namespace {

using Soup = std::vector<double>;

void PushTri(Soup& soup, double ax, double ay, double bx, double by, double cx,
             double cy)
{
  soup.insert(soup.end(), { ax, ay, bx, by, cx, cy });
}

// Regular nx*ny grid, two triangles per cell, every corner emitted per
// triangle (so shared vertices are duplicated on input).
Soup MakeGrid(int nx, int ny)
{
  Soup soup;
  auto vid = [&](int x, int y) { return y * (nx + 1) + x; };
  auto px = [&](int id) { return static_cast<double>(id % (nx + 1)); };
  auto py = [&](int id) { return static_cast<double>(id / (nx + 1)); };
  auto push = [&](int a, int b, int c) {
    for (int id : { a, b, c }) {
      soup.push_back(px(id));
      soup.push_back(py(id));
    }
  };
  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x) {
      const int a = vid(x, y), b = vid(x + 1, y), c = vid(x + 1, y + 1),
                d = vid(x, y + 1);
      push(a, b, c);
      push(a, c, d);
    }
  }
  return soup;
}

// Square ring (outer square minus inner square), triangulated into 8 triangles.
Soup MakeHoleRing()
{
  Soup soup;
  PushTri(soup, 0, 0, 4, 0, 3, 1);
  PushTri(soup, 0, 0, 3, 1, 1, 1);
  PushTri(soup, 4, 0, 4, 4, 3, 3);
  PushTri(soup, 4, 0, 3, 3, 3, 1);
  PushTri(soup, 4, 4, 0, 4, 1, 3);
  PushTri(soup, 4, 4, 1, 3, 3, 3);
  PushTri(soup, 0, 4, 0, 0, 1, 1);
  PushTri(soup, 0, 4, 1, 1, 1, 3);
  return soup;
}

// Long 1*n band of two triangles per cell.
Soup MakeLongStrip(int n)
{
  Soup soup;
  for (int i = 0; i < n; ++i) {
    PushTri(soup, i, 0, i + 1, 0, i + 1, 1);
    PushTri(soup, i, 0, i + 1, 1, i, 1);
  }
  return soup;
}

// Four triangles, two of which are zero-area / collinear degenerates.
Soup MakeCollinear()
{
  Soup soup;
  PushTri(soup, 0, 0, 1, 0, 2, 0); // degenerate: collinear, distinct points
  PushTri(soup, 0, 0, 2, 0, 2, 2); // valid
  PushTri(soup, 0, 0, 2, 2, 0, 2); // valid
  PushTri(soup, 0, 0, 0, 2, 0, 1); // degenerate: collinear, distinct points
  return soup;
}

// Pure radial fan: k triangles around a center hub, a closed rim ring.
Soup MakeRadialFan(int k)
{
  Soup soup;
  const double pi = 3.14159265358979323846;
  for (int i = 0; i < k; ++i) {
    const double a0 = 2 * pi * i / k;
    const double a1 = 2 * pi * (i + 1) / k;
    PushTri(soup, 0, 0, std::cos(a0), std::sin(a0), std::cos(a1), std::sin(a1));
  }
  return soup;
}

// Independent reconstruction check that does not trust the pipeline's own
// metrics: decode every emitted segment back into coordinate triples and
// compare with the input triangle multiset (orientation-insensitive).
bool Reconstructs(const Soup& input, const p2t::PipelineResult& r)
{
  std::map<std::pair<double, double>, std::uint32_t> id_of;
  for (std::size_t i = 0; i < r.vertices.size(); ++i) {
    id_of[{ r.vertices[i].x, r.vertices[i].y }] = static_cast<std::uint32_t>(i);
  }
  using Tri = std::array<std::uint32_t, 3>;
  std::set<Tri> expected, actual;
  const auto canon = [](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    Tri t{ a, b, c };
    std::sort(t.begin(), t.end());
    return t;
  };
  for (std::size_t t = 0; t < r.metrics.triangle_count; ++t) {
    auto ia = id_of.find({ input[6 * t], input[6 * t + 1] });
    auto ib = id_of.find({ input[6 * t + 2], input[6 * t + 3] });
    auto ic = id_of.find({ input[6 * t + 4], input[6 * t + 5] });
    if (ia == id_of.end() || ib == id_of.end() || ic == id_of.end()) {
      return false;
    }
    expected.insert(canon(ia->second, ib->second, ic->second));
  }
  std::size_t cursor = 0;
  for (const p2t::Segment& s : r.segments) {
    const auto* idx = &r.indices[cursor];
    if (s.kind == p2t::SegmentKind::Triangle) {
      actual.insert(canon(idx[0], idx[1], idx[2]));
    } else if (s.kind == p2t::SegmentKind::Strip) {
      for (std::uint32_t i = 0; i + 2 < s.length; ++i) {
        actual.insert(canon(idx[i], idx[i + 1], idx[i + 2]));
      }
    } else {
      const std::uint32_t hub = idx[0];
      for (std::uint32_t i = 0; i + 2 < s.length; ++i) {
        actual.insert(canon(hub, idx[i + 1], idx[i + 2]));
      }
    }
    cursor += s.length;
  }
  return expected == actual;
}

std::string Summarize(const char* name, const p2t::PipelineResult& r)
{
  std::ostringstream os;
  os << name << ": input_vertices=" << r.metrics.input_vertex_count
     << " output_vertices=" << r.metrics.output_vertex_count
     << " merged=" << r.metrics.merged_vertex_count
     << " triangles=" << r.metrics.triangle_count
     << " strips=" << r.metrics.strip_count << " fans=" << r.metrics.fan_count
     << " plain=" << r.metrics.plain_segment_count
     << " degenerate=" << r.metrics.degenerate_count
     << " indices=" << r.metrics.index_count
     << " reuse=" << r.metrics.vertex_reuse_rate
     << " consistent=" << r.metrics.reconstruction_consistent
     << " peak_bytes=" << r.metrics.peak_memory_bytes;
  return os.str();
}

} // namespace

BOOST_AUTO_TEST_CASE(PipelineRegularGridTest)
{
  const Soup soup = MakeGrid(4, 3);
  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(soup);
  BOOST_TEST_MESSAGE(Summarize("regular_grid_4x3", r));

  BOOST_CHECK_EQUAL(r.metrics.input_vertex_count, 72u); // 24 tris * 3 corners
  BOOST_CHECK_EQUAL(r.metrics.output_vertex_count, 20u); // 5*4 grid points
  BOOST_CHECK_EQUAL(r.metrics.merged_vertex_count, 52u);
  BOOST_CHECK_EQUAL(r.metrics.triangle_count, 24u);
  BOOST_CHECK_GE(r.metrics.strip_count, 1u); // grids are strip territory
  BOOST_CHECK_EQUAL(r.metrics.degenerate_count, 0u);
  BOOST_CHECK_GT(r.metrics.vertex_reuse_rate, 0.0);
  BOOST_CHECK(r.metrics.reconstruction_consistent);
  BOOST_CHECK_EQUAL(r.metrics.reconstruction_mismatches, 0u);
  BOOST_CHECK(Reconstructs(soup, r));
  // Segment table and index buffer agree on offsets/lengths.
  std::uint32_t total = 0;
  for (const p2t::Segment& s : r.segments) {
    BOOST_CHECK_EQUAL(s.offset, total);
    total += s.length;
  }
  BOOST_CHECK_EQUAL(total, r.indices.size());
  BOOST_CHECK_EQUAL(r.segment_lengths.size(), r.segments.size());
}

BOOST_AUTO_TEST_CASE(PipelineHoleTest)
{
  const Soup soup = MakeHoleRing();
  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(soup);
  BOOST_TEST_MESSAGE(Summarize("hole_ring", r));

  BOOST_CHECK_EQUAL(r.metrics.input_vertex_count, 24u);
  BOOST_CHECK_EQUAL(r.metrics.output_vertex_count, 8u);
  BOOST_CHECK_EQUAL(r.metrics.merged_vertex_count, 16u);
  BOOST_CHECK_EQUAL(r.metrics.triangle_count, 8u);
  BOOST_CHECK_GE(r.metrics.strip_count, 1u);
  BOOST_CHECK_GT(r.metrics.vertex_reuse_rate, 0.0);
  BOOST_CHECK(r.metrics.reconstruction_consistent);
  BOOST_CHECK(Reconstructs(soup, r));
}

BOOST_AUTO_TEST_CASE(PipelineLongStripTest)
{
  const Soup soup = MakeLongStrip(8);
  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(soup);
  BOOST_TEST_MESSAGE(Summarize("long_strip_1x8", r));

  BOOST_CHECK_EQUAL(r.metrics.input_vertex_count, 48u);
  BOOST_CHECK_EQUAL(r.metrics.output_vertex_count, 18u);
  BOOST_CHECK_EQUAL(r.metrics.triangle_count, 16u);
  BOOST_CHECK_GE(r.metrics.strip_count, 1u);
  BOOST_CHECK_GT(r.metrics.vertex_reuse_rate, 0.0);
  BOOST_CHECK(r.metrics.reconstruction_consistent);
  BOOST_CHECK(Reconstructs(soup, r));
}

BOOST_AUTO_TEST_CASE(PipelineCollinearDegenerateTest)
{
  const Soup soup = MakeCollinear();
  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(soup);
  BOOST_TEST_MESSAGE(Summarize("collinear_degenerate", r));

  BOOST_CHECK_EQUAL(r.metrics.triangle_count, 4u);
  BOOST_CHECK_EQUAL(r.metrics.degenerate_count, 2u);
  // Degenerate triangles are always preserved as plain 3-index segments.
  BOOST_CHECK_GE(r.metrics.plain_segment_count, 2u);
  BOOST_CHECK(r.metrics.reconstruction_consistent);
  BOOST_CHECK(Reconstructs(soup, r));
}

BOOST_AUTO_TEST_CASE(PipelineRadialFanTest)
{
  const Soup soup = MakeRadialFan(6);

  // Default strip-first policy: the ring is encoded as strips correctly.
  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(soup);
  BOOST_TEST_MESSAGE(Summarize("radial_fan_6_strip_first", r));
  BOOST_CHECK_EQUAL(r.metrics.triangle_count, 6u);
  BOOST_CHECK_EQUAL(r.metrics.output_vertex_count, 8u); // 1 hub + 6 rim
  BOOST_CHECK_GE(r.metrics.strip_count, 1u);
  BOOST_CHECK(r.metrics.reconstruction_consistent);
  BOOST_CHECK(Reconstructs(soup, r));

  // Fan-first policy: the same mesh becomes one fan segment and reaches a
  // strictly higher reuse rate, proving the fan output path end to end.
  p2t::TrianglePipelineOptions prefer_fans;
  prefer_fans.preferFans = true;
  const p2t::PipelineResult rf = p2t::TriangleMeshPipeline::Build(soup, prefer_fans);
  BOOST_TEST_MESSAGE(Summarize("radial_fan_6_prefer_fans", rf));
  BOOST_CHECK_EQUAL(rf.metrics.fan_count, 1u);
  BOOST_CHECK_EQUAL(rf.metrics.strip_count, 0u);
  BOOST_CHECK_EQUAL(rf.metrics.plain_segment_count, 0u);
  BOOST_CHECK_EQUAL(rf.segments.size(), 1u);
  BOOST_CHECK(rf.segments[0].kind == p2t::SegmentKind::Fan);
  BOOST_CHECK_GT(rf.metrics.vertex_reuse_rate, r.metrics.vertex_reuse_rate);
  BOOST_CHECK(rf.metrics.reconstruction_consistent);
  BOOST_CHECK(Reconstructs(soup, rf));

  // A 5-triangle ring under fan-first is likewise a single fan.
  const Soup penta = MakeRadialFan(5);
  const p2t::PipelineResult rp = p2t::TriangleMeshPipeline::Build(penta, prefer_fans);
  BOOST_TEST_MESSAGE(Summarize("radial_fan_5_prefer_fans", rp));
  BOOST_CHECK_EQUAL(rp.metrics.output_vertex_count, 7u);
  BOOST_CHECK_EQUAL(rp.metrics.fan_count, 1u);
  BOOST_CHECK(rp.metrics.reconstruction_consistent);
  BOOST_CHECK(Reconstructs(penta, rp));
}

BOOST_AUTO_TEST_CASE(PipelineStageSwitchesTest)
{
  const Soup soup = MakeGrid(4, 3);

  // Stage 1 off: no welding, every corner keeps its identity. The external
  // coordinate-based reconstruction helper needs unique vertices, so for this
  // case we rely on the pipeline's own welded-id consistency proof instead.
  p2t::TrianglePipelineOptions no_weld;
  no_weld.weldVertices = false;
  const p2t::PipelineResult rw = p2t::TriangleMeshPipeline::Build(soup, no_weld);
  BOOST_CHECK_EQUAL(rw.metrics.output_vertex_count, rw.metrics.input_vertex_count);
  BOOST_CHECK_EQUAL(rw.metrics.merged_vertex_count, 0u);
  BOOST_CHECK(rw.metrics.reconstruction_consistent);
  BOOST_CHECK_EQUAL(rw.metrics.reconstruction_mismatches, 0u);

  // Stage 2 off: every triangle is a plain 3-index segment, no strips/fans.
  p2t::TrianglePipelineOptions no_reorder;
  no_reorder.reorder = false;
  const p2t::PipelineResult rr = p2t::TriangleMeshPipeline::Build(soup, no_reorder);
  BOOST_CHECK_EQUAL(rr.metrics.strip_count, 0u);
  BOOST_CHECK_EQUAL(rr.metrics.fan_count, 0u);
  BOOST_CHECK_EQUAL(rr.metrics.plain_segment_count, 24u);
  BOOST_CHECK_EQUAL(rr.metrics.index_count, 72u);
  BOOST_CHECK_CLOSE(rr.metrics.vertex_reuse_rate, 0.0, 1e-9);
  BOOST_CHECK(rr.metrics.reconstruction_consistent);

  // Stage 3 off: index buffer/segment table stay empty, but stage 2 still ran.
  p2t::TrianglePipelineOptions no_emit;
  no_emit.emitIndices = false;
  const p2t::PipelineResult re = p2t::TriangleMeshPipeline::Build(soup, no_emit);
  BOOST_CHECK(re.indices.empty());
  BOOST_CHECK(re.segments.empty());
  BOOST_CHECK(re.segment_lengths.empty());
  BOOST_CHECK_EQUAL(re.metrics.index_count, 0u);
  BOOST_CHECK(re.metrics.reconstruction_consistent);

  // All three off: identity pass over the raw corner array.
  p2t::TrianglePipelineOptions all_off;
  all_off.weldVertices = false;
  all_off.reorder = false;
  all_off.emitIndices = false;
  const p2t::PipelineResult ro = p2t::TriangleMeshPipeline::Build(soup, all_off);
  BOOST_CHECK_EQUAL(ro.metrics.output_vertex_count, soup.size() / 2);
  BOOST_CHECK(ro.indices.empty());
  BOOST_CHECK(ro.segments.empty());
  BOOST_CHECK(ro.metrics.reconstruction_consistent);
}

BOOST_AUTO_TEST_CASE(PipelineToleranceWeldTest)
{
  // Jitter some duplicated corners by less than the tolerance; they must weld
  // back onto the grid vertices, while reconstruction stays exact on ids.
  Soup jittered = MakeGrid(4, 3);
  for (std::size_t i = 0; i < jittered.size(); i += 12) {
    jittered[i] += 0.01;
    jittered[i + 1] += 0.01;
  }
  p2t::TrianglePipelineOptions opts;
  opts.tolerance = 0.05;
  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(jittered, opts);
  BOOST_TEST_MESSAGE(Summarize("tolerance_weld_0.05", r));

  BOOST_CHECK_EQUAL(r.metrics.output_vertex_count, 20u);
  BOOST_CHECK_GT(r.metrics.merged_vertex_count, 0u);
  BOOST_CHECK(r.metrics.reconstruction_consistent);

  // Tolerance of zero still performs exact welding.
  p2t::TrianglePipelineOptions exact;
  exact.tolerance = 0.0;
  const p2t::PipelineResult re = p2t::TriangleMeshPipeline::Build(MakeGrid(2, 2), exact);
  BOOST_CHECK_EQUAL(re.metrics.output_vertex_count, 9u);
  BOOST_CHECK(re.metrics.reconstruction_consistent);
}

BOOST_AUTO_TEST_CASE(PipelineEmptyAndNegativeToleranceTest)
{
  // Empty input is a valid no-op with zero memory and empty outputs.
  const p2t::PipelineResult empty = p2t::TriangleMeshPipeline::Build(Soup{});
  BOOST_CHECK(empty.vertices.empty());
  BOOST_CHECK(empty.indices.empty());
  BOOST_CHECK(empty.segments.empty());
  BOOST_CHECK_EQUAL(empty.metrics.triangle_count, 0u);
  BOOST_CHECK_EQUAL(empty.metrics.peak_memory_bytes, 0u);
  BOOST_CHECK(empty.metrics.reconstruction_consistent);

  // Negative tolerance is rejected explicitly.
  p2t::TrianglePipelineOptions bad;
  bad.tolerance = -1.0;
  BOOST_CHECK_THROW(p2t::TriangleMeshPipeline::Build(MakeGrid(1, 1), bad),
                    std::invalid_argument);

  // Malformed corner array (not a multiple of 6) is rejected explicitly.
  Soup malformed{ 0, 0, 1, 0 };
  BOOST_CHECK_THROW(p2t::TriangleMeshPipeline::Build(malformed), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(PipelineMemoryIsLinearTest)
{
  // Peak bytes must scale linearly with the welded vertex count: bytes/V
  // stays bounded as the grid grows, which is the operational signature of
  // O(V) memory (no quadratic accumulation).
  double previous_ratio = 0.0;
  for (int n : { 4, 8, 16 }) {
    const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(MakeGrid(n, n));
    const double ratio =
      static_cast<double>(r.metrics.peak_memory_bytes) / r.metrics.output_vertex_count;
    BOOST_TEST_MESSAGE("grid " << n << "x" << n << " V="
                               << r.metrics.output_vertex_count
                               << " peak=" << r.metrics.peak_memory_bytes
                               << " bytes/V=" << ratio);
    BOOST_CHECK_GT(r.metrics.peak_memory_bytes, 0u);
    if (previous_ratio > 0.0) {
      BOOST_CHECK_LE(ratio, previous_ratio * 1.5); // no super-linear blow-up
    }
    previous_ratio = ratio;
  }
}

BOOST_AUTO_TEST_CASE(PipelineCdtOverloadUnchangedTest)
{
  // The CDT overload must read triangulation output without changing it:
  // run a CDT, snapshot its triangles, build a pipeline, and confirm the CDT
  // triangle count and point coordinates are byte-for-byte unchanged.
  std::vector<p2t::Point*> polyline{ new p2t::Point(0, 0), new p2t::Point(1, 0),
                                     new p2t::Point(1, 1), new p2t::Point(0, 1) };
  p2t::CDT cdt{ polyline };
  cdt.Triangulate();
  const auto triangles = cdt.GetTriangles();
  const std::size_t before = triangles.size();

  std::vector<std::pair<double, double>> snapshot;
  for (p2t::Triangle* t : triangles) {
    for (int i = 0; i < 3; ++i) {
      const p2t::Point* p = t->GetPoint(i);
      snapshot.emplace_back(p->x, p->y);
    }
  }

  const p2t::PipelineResult r = p2t::TriangleMeshPipeline::Build(triangles);
  BOOST_CHECK_EQUAL(r.metrics.triangle_count, before);
  BOOST_CHECK(r.metrics.reconstruction_consistent);

  std::size_t k = 0;
  for (p2t::Triangle* t : triangles) {
    for (int i = 0; i < 3; ++i) {
      const p2t::Point* p = t->GetPoint(i);
      BOOST_CHECK_EQUAL(p->x, snapshot[k].first);
      BOOST_CHECK_EQUAL(p->y, snapshot[k].second);
      ++k;
    }
  }
  for (p2t::Point* p : polyline) {
    delete p;
  }
}

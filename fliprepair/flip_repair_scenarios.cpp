/*
 * Flip-repair scenario driver.
 *
 * Builds four small hand-made triangulations, deliberately flips edges to
 * simulate "triangles flipped over constrained edges", then runs the repair
 * and reports per group: flip count, min-angle improvement, constraint
 * integrity, and convergence status. Non-convergence (fixpoint with
 * remaining violations or exhausted round budget) is reported, never
 * silently swallowed. Exit code is non-zero if any expectation is unmet.
 */
#include "flip_repair.hpp"

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

using fliprepair::EdgeKey;

struct Mesh {
  std::vector<std::unique_ptr<p2t::Point>> points;
  std::vector<std::unique_ptr<p2t::Triangle>> storage;
  std::vector<p2t::Triangle*> triangles;

  p2t::Point* AddPoint(double x, double y)
  {
    points.emplace_back(new p2t::Point(x, y));
    return points.back().get();
  }

  p2t::Triangle* AddTriangle(p2t::Point* a, p2t::Point* b, p2t::Point* c)
  {
    if (fliprepair::Cross2(*a, *b, *c) < 0) {
      std::swap(b, c); // keep triangles counter-clockwise
    }
    storage.emplace_back(new p2t::Triangle(*a, *b, *c));
    triangles.push_back(storage.back().get());
    return storage.back().get();
  }

  void LinkNeighbors()
  {
    for (auto* t : triangles) {
      for (auto* o : triangles) {
        if (t != o) {
          t->MarkNeighbor(*o); // no-op unless an edge is shared
        }
      }
    }
  }

  void MarkConstrained(p2t::Point* p, p2t::Point* q)
  {
    for (auto* t : triangles) {
      if (t->Contains(p, q)) {
        t->MarkConstrainedEdge(p, q);
      }
    }
  }
};

/// Structural sanity check: CCW orientation, symmetric neighbor links, no
/// edge shared by more than two triangles.
bool ValidateMesh(std::vector<p2t::Triangle*>& triangles, std::string& error)
{
  for (auto* t : triangles) {
    if (fliprepair::Cross2(*t->GetPoint(0), *t->GetPoint(1), *t->GetPoint(2)) <= 0) {
      error = "triangle not counter-clockwise";
      return false;
    }
    for (int i = 0; i < 3; ++i) {
      p2t::Point* p1 = t->GetPoint((i + 1) % 3);
      p2t::Point* p2 = t->GetPoint((i + 2) % 3);
      p2t::Triangle* n = t->GetNeighbor(i);
      if (n == nullptr) {
        continue;
      }
      if (!n->Contains(p1, p2)) {
        error = "neighbor does not share the edge";
        return false;
      }
      const int ni = n->EdgeIndex(p1, p2);
      if (ni < 0 || n->GetNeighbor(ni) != t) {
        error = "neighbor link not symmetric";
        return false;
      }
    }
  }
  for (const auto& entry : fliprepair::BuildEdgeMap(triangles)) {
    if (entry.second.count > 2) {
      error = "edge shared by more than two triangles";
      return false;
    }
  }
  return true;
}

/// Flip an arbitrary interior edge, ignoring constrained status. This is the
/// "corruption" used to fabricate broken inputs; the repair itself never
/// does this. Returns false if the edge is missing or the quad is concave.
bool CorruptFlipEdge(const std::vector<p2t::Triangle*>& triangles,
                     p2t::Point* a, p2t::Point* b)
{
  const auto edges = fliprepair::BuildEdgeMap(triangles);
  const auto it = edges.find(fliprepair::MakeEdgeKey(a, b));
  if (it == edges.end() || it->second.count != 2) {
    return false;
  }
  p2t::Triangle* t0 = it->second.tri[0];
  p2t::Triangle* t1 = it->second.tri[1];
  p2t::Point* c = t0->GetPoint(t0->EdgeIndex(a, b));
  p2t::Point* d = t1->GetPoint(t1->EdgeIndex(a, b));
  if (!fliprepair::ProperlyIntersect(*a, *b, *c, *d)) {
    return false;
  }
  fliprepair::FlipEdge(*t0, *c, *t1, *d);
  return true;
}

bool HasEdge(const std::vector<p2t::Triangle*>& triangles,
             p2t::Point* a, p2t::Point* b)
{
  const auto edges = fliprepair::BuildEdgeMap(triangles);
  return edges.find(fliprepair::MakeEdgeKey(a, b)) != edges.end();
}

void PrintReport(const fliprepair::Report& r, int required_count, int max_rounds)
{
  std::printf("  翻转次数: %ld (轮数 %d/%d)\n", r.flips_applied, r.rounds_used, max_rounds);
  std::printf("  最小角: %.4f° -> %.4f° (提升 %+.4f°)\n", r.min_angle_before_deg,
              r.min_angle_after_deg, r.min_angle_after_deg - r.min_angle_before_deg);
  std::printf("  约束完整性: %s (缺失 %d/%d)\n", r.constraints_intact ? "完好" : "已破坏",
              r.missing_constraints, required_count);
  std::printf("  收敛: %s, 不动点: %s, 轮数用完未收敛: %s, 剩余违反: %d\n",
              r.converged ? "是" : "否", r.fixpoint_reached ? "是" : "否",
              r.round_limit_hit ? "是" : "否", r.remaining_violations);
}

int g_failures = 0;

void Expect(bool condition, const std::string& what)
{
  if (!condition) {
    std::printf("  [失败] 预期未满足: %s\n", what.c_str());
    ++g_failures;
  }
}

// ---------------------------------------------------------------------------
// 组1: 带约束 — 单条内部约束边被翻走, 修复应将其翻回.
// ---------------------------------------------------------------------------
void ScenarioSingleConstraint()
{
  std::printf("=== 组1: 带约束 ===\n");
  Mesh mesh;
  p2t::Point* a = mesh.AddPoint(0.0, 0.0);
  p2t::Point* b = mesh.AddPoint(4.0, 0.0);
  p2t::Point* c = mesh.AddPoint(4.0, 4.0);
  p2t::Point* d = mesh.AddPoint(0.0, 4.0);
  p2t::Point* e = mesh.AddPoint(1.5, 1.8);
  mesh.AddTriangle(a, b, e);
  mesh.AddTriangle(b, c, e);
  mesh.AddTriangle(c, d, e);
  mesh.AddTriangle(d, a, e);
  mesh.LinkNeighbors();
  mesh.MarkConstrained(e, c); // required constrained segment E-C

  Expect(CorruptFlipEdge(mesh.triangles, e, c), "组1: 构造损坏(翻走约束边 E-C)");
  Expect(!HasEdge(mesh.triangles, e, c), "组1: 损坏后约束边确实缺失");

  fliprepair::Options options;
  const fliprepair::Report r =
      fliprepair::Repair(mesh.triangles, {EdgeKey(e, c)}, options);
  PrintReport(r, 1, options.max_rounds);

  std::string error;
  Expect(ValidateMesh(mesh.triangles, error), "组1: 网格结构有效 " + error);
  Expect(r.converged, "组1: 收敛");
  Expect(r.flips_applied >= 1, "组1: 至少翻转一次");
  Expect(r.constraints_intact && HasEdge(mesh.triangles, e, c),
         "组1: 约束边 E-C 被恢复");
  Expect(!r.round_limit_hit, "组1: 未触及轮数上限");
}

// ---------------------------------------------------------------------------
// 组2: 多约束 — 两条约束边先后被翻走, 修复应逐一恢复.
// ---------------------------------------------------------------------------
void ScenarioMultiConstraint()
{
  std::printf("=== 组2: 多约束 ===\n");
  Mesh mesh;
  p2t::Point* a = mesh.AddPoint(0.0, 0.0);
  p2t::Point* b = mesh.AddPoint(4.0, 0.0);
  p2t::Point* c = mesh.AddPoint(5.0, 3.0);
  p2t::Point* d = mesh.AddPoint(4.0, 6.0);
  p2t::Point* e = mesh.AddPoint(0.0, 6.0);
  p2t::Point* f = mesh.AddPoint(-1.0, 3.0);
  mesh.AddTriangle(a, b, c);
  mesh.AddTriangle(a, c, d);
  mesh.AddTriangle(a, d, e);
  mesh.AddTriangle(a, e, f);
  mesh.LinkNeighbors();
  mesh.MarkConstrained(a, c);
  mesh.MarkConstrained(a, d);

  Expect(CorruptFlipEdge(mesh.triangles, a, c), "组2: 构造损坏(翻走约束边 A-C)");
  Expect(CorruptFlipEdge(mesh.triangles, a, d), "组2: 构造损坏(翻走约束边 A-D)");
  Expect(!HasEdge(mesh.triangles, a, c) && !HasEdge(mesh.triangles, a, d),
         "组2: 损坏后两条约束边均缺失");

  fliprepair::Options options;
  const fliprepair::Report r =
      fliprepair::Repair(mesh.triangles, {EdgeKey(a, c), EdgeKey(a, d)}, options);
  PrintReport(r, 2, options.max_rounds);

  std::string error;
  Expect(ValidateMesh(mesh.triangles, error), "组2: 网格结构有效 " + error);
  Expect(r.converged, "组2: 收敛");
  Expect(r.flips_applied >= 2, "组2: 至少翻转两次");
  Expect(r.constraints_intact && HasEdge(mesh.triangles, a, c) &&
             HasEdge(mesh.triangles, a, d),
         "组2: 两条约束边均被恢复");
}

// ---------------------------------------------------------------------------
// 组3: 约束冲突 — 两条约束线段互相交叉, 不可能同时满足. 修复必须如实报告
// 无法收敛, 且另跑一组轮数预算不足的用例验证 "轮数用完" 也会被报出来.
// ---------------------------------------------------------------------------
void ScenarioConflictingConstraints()
{
  std::printf("=== 组3: 约束冲突 ===\n");
  Mesh mesh;
  p2t::Point* a = mesh.AddPoint(0.0, 0.0);
  p2t::Point* b = mesh.AddPoint(4.0, 0.0);
  p2t::Point* c = mesh.AddPoint(4.0, 4.0);
  p2t::Point* d = mesh.AddPoint(0.0, 4.0);
  mesh.AddTriangle(a, b, c);
  mesh.AddTriangle(a, c, d);
  mesh.LinkNeighbors();
  mesh.MarkConstrained(a, c);

  // A-C 与 B-D 互相交叉: 恢复 B-D 必须翻走 A-C, 而约束边一律不许翻.
  fliprepair::Options options;
  const fliprepair::Report r =
      fliprepair::Repair(mesh.triangles, {EdgeKey(a, c), EdgeKey(b, d)}, options);
  PrintReport(r, 2, options.max_rounds);

  std::string error;
  Expect(ValidateMesh(mesh.triangles, error), "组3: 网格结构有效 " + error);
  Expect(!r.converged, "组3: 冲突下不得报告收敛");
  Expect(r.fixpoint_reached, "组3: 应停在不动点(违反边不可翻)");
  Expect(!r.constraints_intact && r.missing_constraints == 1,
         "组3: 约束完整性被破坏且如实报告缺失 1 条");
  Expect(r.remaining_violations > 0, "组3: 剩余违反必须非零");
  Expect(HasEdge(mesh.triangles, a, c), "组3: 已有约束边 A-C 未被翻走");

  // 轮数预算不足: 用组2的损坏输入, 只给 1 轮, 必须报告 "轮数用完未收敛".
  std::printf("--- 组3b: 轮数预算不足 ---\n");
  Mesh mesh2;
  p2t::Point* a2 = mesh2.AddPoint(0.0, 0.0);
  p2t::Point* b2 = mesh2.AddPoint(4.0, 0.0);
  p2t::Point* c2 = mesh2.AddPoint(5.0, 3.0);
  p2t::Point* d2 = mesh2.AddPoint(4.0, 6.0);
  p2t::Point* e2 = mesh2.AddPoint(0.0, 6.0);
  p2t::Point* f2 = mesh2.AddPoint(-1.0, 3.0);
  mesh2.AddTriangle(a2, b2, c2);
  mesh2.AddTriangle(a2, c2, d2);
  mesh2.AddTriangle(a2, d2, e2);
  mesh2.AddTriangle(a2, e2, f2);
  mesh2.LinkNeighbors();
  mesh2.MarkConstrained(a2, c2);
  mesh2.MarkConstrained(a2, d2);
  Expect(CorruptFlipEdge(mesh2.triangles, a2, c2), "组3b: 构造损坏(翻走 A-C)");
  Expect(CorruptFlipEdge(mesh2.triangles, a2, d2), "组3b: 构造损坏(翻走 A-D)");

  fliprepair::Options tight;
  tight.max_rounds = 1; // 恢复两条约束边需要两轮, 一轮不够
  const fliprepair::Report r2 =
      fliprepair::Repair(mesh2.triangles, {EdgeKey(a2, c2), EdgeKey(a2, d2)}, tight);
  PrintReport(r2, 2, tight.max_rounds);

  Expect(ValidateMesh(mesh2.triangles, error), "组3b: 网格结构有效 " + error);
  Expect(!r2.converged, "组3b: 轮数不足时不得报告收敛");
  Expect(r2.round_limit_hit, "组3b: 必须报告轮数用完未收敛");
  Expect(!r2.constraints_intact, "组3b: 约束完整性必须报告为破坏");
}

// ---------------------------------------------------------------------------
// 组4: 无约束 — 纯空圆准则修复: 一条内部边违反空圆判据, 修复后最小角提升
// 且三角剖分恢复 Delaunay 性质.
// ---------------------------------------------------------------------------
void ScenarioNoConstraint()
{
  std::printf("=== 组4: 无约束 ===\n");
  Mesh mesh;
  p2t::Point* a = mesh.AddPoint(0.0, 0.0);
  p2t::Point* b = mesh.AddPoint(4.0, 0.0);
  p2t::Point* c = mesh.AddPoint(3.0, 2.0);
  p2t::Point* d = mesh.AddPoint(0.0, 3.0);
  // 正确的 Delaunay 对角线是 A-C; 先按 A-C 建网再翻成违反空圆判据的 B-D.
  mesh.AddTriangle(a, b, c);
  mesh.AddTriangle(a, c, d);
  mesh.LinkNeighbors();

  Expect(CorruptFlipEdge(mesh.triangles, a, c), "组4: 构造损坏(翻成非 Delaunay 对角线)");
  Expect(!p2t::IsDelaunay(mesh.triangles), "组4: 损坏后确实违反空圆判据");

  fliprepair::Options options;
  const fliprepair::Report r = fliprepair::Repair(mesh.triangles, {}, options);
  PrintReport(r, 0, options.max_rounds);

  std::string error;
  Expect(ValidateMesh(mesh.triangles, error), "组4: 网格结构有效 " + error);
  Expect(r.converged, "组4: 收敛");
  Expect(r.flips_applied == 1, "组4: 恰好翻转一次");
  Expect(r.min_angle_after_deg > r.min_angle_before_deg, "组4: 最小角提升");
  Expect(r.constraints_intact, "组4: 无约束时完整性恒为完好");
  Expect(p2t::IsDelaunay(mesh.triangles), "组4: 修复后满足空圆判据");
}

} // namespace

int main()
{
  ScenarioSingleConstraint();
  ScenarioMultiConstraint();
  ScenarioConflictingConstraints();
  ScenarioNoConstraint();

  if (g_failures == 0) {
    std::printf("全部四组用例通过, 报告可信。\n");
    return 0;
  }
  std::printf("有 %d 项预期未满足, 结果不可采信。\n", g_failures);
  return 1;
}

#!/usr/bin/env python3
"""Watertightness inspection and repair for triangle meshes (OBJ I/O).

Pipeline
--------
1. Optional vertex weld: merge vertices closer than epsilon. This closes
   hairline seams caused by duplicated vertices -- the classic "invisible
   gap" that only shows up when the slicer/printer complains.
2. Edge pairing: every undirected edge of a watertight oriented manifold
   is referenced by exactly two triangles, in opposite directions. Edges
   referenced once are boundary edges; more than twice is non-manifold;
   twice in the same direction is an orientation conflict.
3. Boundary edges are chained into closed loops. For every loop we compute
   the winding (Newell normal + projected signed area), perimeter, span and
   vector area, and decide whether the loop is closable at all
   (non-degenerate area) and small enough to auto-fill (span <= threshold).
4. Small gaps are filled with a centroid fan whose triangles are oriented
   opposite to the directed boundary edges, so the repaired mesh stays
   consistently oriented.
5. Loops that cannot be repaired (real holes, degenerate loops) and
   non-manifold / orientation-conflict edges are reported separately and
   never silently patched.
6. After repair the full watertightness check runs again; the result is
   part of the report.

CLI
---
    python3 watertight.py check  mesh.obj [--weld-eps E] [--report r.json]
    python3 watertight.py repair mesh.obj -o fixed.obj [--weld-eps E]
                                  [--gap-threshold T] [--report r.json]

Exit code is 0 when the (repaired) mesh is watertight, 1 otherwise.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field

Vec3 = tuple[float, float, float]
Face = tuple[int, int, int]


# ---------------------------------------------------------------------------
# small vec3 helpers
# ---------------------------------------------------------------------------

def v_sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def v_cross(a: Vec3, b: Vec3) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def v_len(a: Vec3) -> float:
    return math.sqrt(v_dot(a, a))


# ---------------------------------------------------------------------------
# mesh container + OBJ I/O
# ---------------------------------------------------------------------------

@dataclass
class Mesh:
    vertices: list[Vec3]
    faces: list[Face]

    def bbox_diagonal(self) -> float:
        if not self.vertices:
            return 0.0
        lo = [min(p[i] for p in self.vertices) for i in range(3)]
        hi = [max(p[i] for p in self.vertices) for i in range(3)]
        return v_len((hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]))


def load_obj(path: str) -> Mesh:
    """Load a Wavefront OBJ file. Polygons are fan-triangulated."""
    vertices: list[Vec3] = []
    faces: list[Face] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "v" and len(parts) >= 4:
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif parts[0] == "f" and len(parts) >= 4:
                idx = []
                for tok in parts[1:]:
                    i = int(tok.split("/")[0])
                    idx.append(i - 1 if i > 0 else len(vertices) + i)
                for k in range(1, len(idx) - 1):  # fan triangulation
                    faces.append((idx[0], idx[k], idx[k + 1]))
    return Mesh(vertices, faces)


def save_obj(mesh: Mesh, path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for p in mesh.vertices:
            fh.write(f"v {p[0]:.9g} {p[1]:.9g} {p[2]:.9g}\n")
        for f in mesh.faces:
            fh.write(f"f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n")


# ---------------------------------------------------------------------------
# step 1: vertex welding (hairline seam closure)
# ---------------------------------------------------------------------------

def weld_vertices(mesh: Mesh, eps: float) -> tuple[Mesh, dict]:
    """Merge vertices within `eps` of each other (first occurrence wins).

    Returns the welded mesh and a report dict. Degenerate or duplicate
    faces produced by the merge are dropped and counted.
    """
    report = {"epsilon": eps, "merged_vertices": 0, "removed_faces": 0}
    if eps <= 0.0 or not mesh.vertices:
        return mesh, report

    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    remap = list(range(len(mesh.vertices)))
    for i, p in enumerate(mesh.vertices):
        key = tuple(int(math.floor(c / eps)) for c in p)
        rep = None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        if v_len(v_sub(mesh.vertices[j], p)) <= eps:
                            rep = j
                            break
                    if rep is not None:
                        break
                if rep is not None:
                    break
            if rep is not None:
                break
        if rep is None:
            grid[key].append(i)
        else:
            remap[i] = rep
            report["merged_vertices"] += 1

    used = sorted({remap[v] for f in mesh.faces for v in f})
    new_index = {old: n for n, old in enumerate(used)}
    new_vertices = [mesh.vertices[old] for old in used]

    new_faces: list[Face] = []
    seen: set[Face] = set()
    for f in mesh.faces:
        g = (new_index[remap[f[0]]], new_index[remap[f[1]]], new_index[remap[f[2]]])
        if len(set(g)) != 3 or g in seen:
            report["removed_faces"] += 1
            continue
        seen.add(g)
        new_faces.append(g)
    return Mesh(new_vertices, new_faces), report


# ---------------------------------------------------------------------------
# step 2+3: edge pairing and boundary loop extraction
# ---------------------------------------------------------------------------

@dataclass
class BoundaryLoop:
    verts: list[int]          # vertex indices in boundary order (closed)
    perimeter: float
    span: float               # bbox diagonal of the loop vertices
    area: float               | None = None
    normal: Vec3 = (0.0, 0.0, 0.0)
    winding: str = "unknown"  # 'ccw'/'cw' seen from the side `normal` points to


@dataclass
class Analysis:
    boundary_loops: list[BoundaryLoop] = field(default_factory=list)
    nonmanifold_edges: list[tuple[int, int, int]] = field(default_factory=list)  # (a, b, count)
    orientation_conflicts: list[tuple[int, int]] = field(default_factory=list)
    degenerate_faces: int = 0
    warnings: list[str] = field(default_factory=list)

    @property
    def watertight(self) -> bool:
        return (not self.boundary_loops
                and not self.nonmanifold_edges
                and not self.orientation_conflicts
                and self.degenerate_faces == 0)


def _newell_normal(pts: list[Vec3]) -> Vec3:
    nx = ny = nz = 0.0
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        nx += (a[1] - b[1]) * (a[2] + b[2])
        ny += (a[2] - b[2]) * (a[0] + b[0])
        nz += (a[0] - b[0]) * (a[1] + b[1])
    return (nx, ny, nz)


def _projected_signed_area(pts: list[Vec3], axis: int) -> float:
    """Shoelace area after projecting onto the plane normal to `axis`.

    Axis-conforming coordinate order keeps the sign consistent with the
    right-hand rule about that axis.
    """
    if axis == 0:
        coords = [(p[1], p[2]) for p in pts]
    elif axis == 1:
        coords = [(p[2], p[0]) for p in pts]
    else:
        coords = [(p[0], p[1]) for p in pts]
    s = 0.0
    for i in range(len(coords)):
        x0, y0 = coords[i]
        x1, y1 = coords[(i + 1) % len(coords)]
        s += x0 * y1 - x1 * y0
    return 0.5 * s


def _loop_metrics(mesh: Mesh, verts: list[int]) -> BoundaryLoop:
    pts = [mesh.vertices[v] for v in verts]
    n = len(pts)
    perimeter = sum(v_len(v_sub(pts[(i + 1) % n], pts[i])) for i in range(n))
    lo = [min(p[i] for p in pts) for i in range(3)]
    hi = [max(p[i] for p in pts) for i in range(3)]
    span = v_len((hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]))

    normal = _newell_normal(pts)
    norm_len = v_len(normal)
    area = 0.5 * norm_len
    if norm_len > 0.0:
        axis = max(range(3), key=lambda i: abs(normal[i]))
        winding = "ccw" if _projected_signed_area(pts, axis) > 0 else "cw"
        unit = (normal[0] / norm_len, normal[1] / norm_len, normal[2] / norm_len)
    else:
        winding = "degenerate"
        unit = (0.0, 0.0, 0.0)
    return BoundaryLoop(verts=list(verts), perimeter=perimeter, span=span,
                        area=area, normal=unit, winding=winding)


def _extract_loops(boundary_edges: list[tuple[int, int]],
                   warnings: list[str]) -> list[list[int]]:
    """Chain directed boundary edges into closed loops of vertex indices."""
    unused = set(boundary_edges)
    out_edges: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for e in boundary_edges:
        out_edges[e[0]].append(e)

    loops: list[list[int]] = []
    while unused:
        edge = next(iter(unused))
        loop = [edge[0]]
        closed = False
        while True:
            unused.discard(edge)
            nxt = edge[1]
            if nxt == loop[0]:
                closed = True
                break
            loop.append(nxt)
            candidates = [e for e in out_edges[nxt] if e in unused]
            if not candidates:
                warnings.append(
                    f"boundary chain ends at open vertex {nxt}; chain discarded")
                break
            if len(candidates) > 1:
                warnings.append(
                    f"non-manifold boundary vertex {nxt}: "
                    f"{len(candidates)} outgoing boundary edges")
            edge = candidates[0]
        if closed and len(loop) >= 3:
            loops.append(loop)
    return loops


def analyze(mesh: Mesh) -> Analysis:
    """Pair edges and extract boundary loops. Pure check, no mutation."""
    result = Analysis()
    edge_dirs: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for f in mesh.faces:
        if len(set(f)) != 3:
            result.degenerate_faces += 1
            continue
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            edge_dirs[(min(a, b), max(a, b))].append((a, b))

    boundary: list[tuple[int, int]] = []
    for (lo, hi), dirs in sorted(edge_dirs.items()):
        if len(dirs) == 1:
            boundary.append(dirs[0])
        elif len(dirs) == 2:
            # consistent orientation <=> the two directed edges are opposite
            if dirs[0][0] == dirs[1][0]:
                result.orientation_conflicts.append((lo, hi))
        else:
            result.nonmanifold_edges.append((lo, hi, len(dirs)))

    raw_loops = _extract_loops(boundary, result.warnings)
    result.boundary_loops = [_loop_metrics(mesh, lp) for lp in raw_loops]
    return result


# ---------------------------------------------------------------------------
# step 4+5: classification and repair
# ---------------------------------------------------------------------------

def classify_loop(loop: BoundaryLoop, gap_threshold: float,
                  min_area: float) -> str:
    """'gap' (auto-fillable), 'hole' (too big) or 'degenerate' (not closable)."""
    if len(loop.verts) < 3 or loop.area is None or loop.area <= min_area:
        return "degenerate"
    if loop.span <= gap_threshold:
        return "gap"
    return "hole"


def _loop_to_dict(loop: BoundaryLoop, classification: str, action: str) -> dict:
    return {
        "vertices": loop.verts,
        "vertex_count": len(loop.verts),
        "perimeter": loop.perimeter,
        "span": loop.span,
        "area": loop.area,
        "normal": list(loop.normal),
        "winding": loop.winding,
        "classification": classification,
        "action": action,
    }


def repair(mesh: Mesh, gap_threshold: float,
           min_area: float = 1e-12) -> tuple[Mesh, dict]:
    """Fill boundary loops whose span is below `gap_threshold`.

    Fill triangles are wound opposite to the directed boundary edges so the
    repaired mesh stays consistently oriented. Unrepairable loops are only
    reported, never touched.
    """
    analysis = analyze(mesh)
    new_vertices = list(mesh.vertices)
    new_faces = list(mesh.faces)

    loops_report = []
    filled = reported_unrepairable = 0
    for loop in analysis.boundary_loops:
        cls = classify_loop(loop, gap_threshold, min_area)
        if cls == "gap":
            pts = [new_vertices[v] for v in loop.verts]
            centroid = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
            ci = len(new_vertices)
            new_vertices.append(centroid)
            k = len(loop.verts)
            for i in range(k):
                # boundary edge runs verts[i] -> verts[i+1]; the fill triangle
                # must traverse that edge in the opposite direction
                new_faces.append((ci, loop.verts[(i + 1) % k], loop.verts[i]))
            filled += 1
            loops_report.append(_loop_to_dict(loop, cls, "filled"))
        else:
            reported_unrepairable += 1
            loops_report.append(_loop_to_dict(loop, cls, "reported"))

    repaired = Mesh(new_vertices, new_faces)
    final = analyze(repaired)
    report = {
        "gap_threshold": gap_threshold,
        "loops": loops_report,
        "filled_loops": filled,
        "unrepairable_loops": reported_unrepairable,
        "nonmanifold_edges": [list(e) for e in analysis.nonmanifold_edges],
        "orientation_conflicts": [list(e) for e in analysis.orientation_conflicts],
        "degenerate_faces": analysis.degenerate_faces,
        "warnings": analysis.warnings,
        "final": {
            "watertight": final.watertight,
            "remaining_boundary_loops": len(final.boundary_loops),
            "remaining_nonmanifold_edges": len(final.nonmanifold_edges),
            "remaining_orientation_conflicts": len(final.orientation_conflicts),
            "vertices": len(repaired.vertices),
            "faces": len(repaired.faces),
        },
    }
    return repaired, report


def check(mesh: Mesh) -> dict:
    """Watertightness report without any modification."""
    analysis = analyze(mesh)
    return {
        "watertight": analysis.watertight,
        "boundary_loops": [
            _loop_to_dict(loop, "unclassified", "none")
            for loop in analysis.boundary_loops
        ],
        "nonmanifold_edges": [list(e) for e in analysis.nonmanifold_edges],
        "orientation_conflicts": [list(e) for e in analysis.orientation_conflicts],
        "degenerate_faces": analysis.degenerate_faces,
        "warnings": analysis.warnings,
        "vertices": len(mesh.vertices),
        "faces": len(mesh.faces),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_summary(report: dict, final: dict) -> None:
    print(f"vertices={final.get('vertices', '?')} faces={final.get('faces', '?')}")
    if "loops" in report:
        for i, lp in enumerate(report["loops"]):
            print(f"  loop {i}: verts={lp['vertex_count']} span={lp['span']:.6g} "
                  f"area={lp['area']:.6g} perimeter={lp['perimeter']:.6g} "
                  f"winding={lp['winding']} -> {lp['classification']} ({lp['action']})")
        print(f"filled={report['filled_loops']} "
              f"unrepairable={report['unrepairable_loops']}")
    elif "boundary_loops" in report:
        for i, lp in enumerate(report["boundary_loops"]):
            print(f"  loop {i}: verts={lp['vertex_count']} span={lp['span']:.6g} "
                  f"area={lp['area']:.6g} winding={lp['winding']}")
    watertight = final.get("watertight", report.get("watertight"))
    print(f"watertight={watertight}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "repair"):
        p = sub.add_parser(name)
        p.add_argument("mesh")
        p.add_argument("--weld-eps", type=float, default=0.0,
                       help="merge vertices closer than this before checking")
        p.add_argument("--report", help="write full JSON report here")
        if name == "repair":
            p.add_argument("-o", "--output", required=True)
            p.add_argument("--gap-threshold", type=float, default=None,
                           help="max loop span to auto-fill "
                                "(default: 5%% of the bbox diagonal)")
    args = parser.parse_args(argv)

    mesh = load_obj(args.mesh)
    if args.weld_eps > 0.0:
        mesh, weld_report = weld_vertices(mesh, args.weld_eps)
        print(f"weld: merged={weld_report['merged_vertices']} "
              f"removed_faces={weld_report['removed_faces']}")

    if args.command == "check":
        report = check(mesh)
        _print_summary(report, report)
        ok = report["watertight"]
    else:
        threshold = args.gap_threshold
        if threshold is None:
            threshold = 0.05 * mesh.bbox_diagonal()
        mesh, report = repair(mesh, threshold)
        save_obj(mesh, args.output)
        _print_summary(report, report["final"])
        ok = report["final"]["watertight"]

    if args.report:
        with open(args.report, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

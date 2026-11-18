#!/usr/bin/env python3
"""Acceptance tests for the watertightness checker/repairer.

Four required mesh groups plus one hairline-seam case:

  1. closed      - subdivided box, must be reported watertight with zero
                   findings (zero false positives)
  2. single gap  - box with one interior triangle removed -> one 3-vertex
                   loop, auto-filled, watertight afterwards
  3. multi gap   - box with three non-adjacent triangles removed -> three
                   3-vertex loops, all filled, watertight afterwards
  4. real hole   - box with a 2x2 quad region removed -> one 8-vertex loop
                   above the fill threshold, reported but NOT filled
  5. hairline    - box whose +x face vertices are duplicated with a 1e-7
                   offset (invisible crack); vertex welding must close it

All expected loop sizes, perimeters, areas and spans below were computed
by hand from the mesh construction and serve as the manual cross-check.
"""

import math
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from watertight import (Mesh, analyze, check, load_obj, repair, save_obj,
                        weld_vertices)

SIZE = 1.0
N = 4                      # grid cells per box face
H = SIZE / N               # grid cell edge length = 0.25
GAP_THRESHOLD = 0.5        # between gap span (~0.354) and hole span (~0.707)


def make_box(n=N, size=SIZE):
    """Closed box, each of the 6 faces an n x n grid of quads (2 tris each).

    Faces are emitted in this order: +z, -z, +x, -x, +y, -y.
    Triangle index of quad (i, j) on face f:
        f * 2*n*n + (i * n + j) * 2   (first triangle), +1 for the second.
    Quad (i, j) on the +z face covers [i*h,(i+1)*h] x [j*h,(j+1)*h] at z=size.
    """
    axes = [
        ((0.0, 0.0, size), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),  # +z
        ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),   # -z
        ((size, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),  # +x
        ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),   # -x
        ((0.0, size, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),  # +y
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),   # -y
    ]
    vid = {}
    vertices = []
    faces = []
    h = size / n

    def index(p):
        key = (round(p[0], 9), round(p[1], 9), round(p[2], 9))
        if key not in vid:
            vid[key] = len(vertices)
            vertices.append(key)
        return vid[key]

    for origin, eu, ev in axes:
        for i in range(n):
            for j in range(n):
                def pt(di, dj):
                    return tuple(origin[k] + (i + di) * h * eu[k]
                                 + (j + dj) * h * ev[k] for k in range(3))
                p00, p10 = index(pt(0, 0)), index(pt(1, 0))
                p11, p01 = index(pt(1, 1)), index(pt(0, 1))
                faces.append((p00, p10, p11))
                faces.append((p00, p11, p01))
    return Mesh(vertices, faces)


def tri_index(face, i, j, second=False, n=N):
    return face * 2 * n * n + (i * n + j) * 2 + (1 if second else 0)


def remove_faces(mesh, indices):
    keep = [f for k, f in enumerate(mesh.faces) if k not in set(indices)]
    return Mesh(list(mesh.vertices), keep)


def approx(a, b, tol=1e-9):
    return abs(a - b) <= tol


# hand-computed expectations for a single removed triangle (legs h, h,
# hypotenuse h*sqrt(2)):
TRI_LOOP_PERIMETER = H * (2.0 + math.sqrt(2.0))   # ~0.853553
TRI_LOOP_AREA = H * H / 2.0                       # 0.03125
TRI_LOOP_SPAN = H * math.sqrt(2.0)                # ~0.353553
# hand-computed expectations for the 2x2 quad hole (side 2h):
HOLE_LOOP_PERIMETER = 8.0 * H                     # 2.0
HOLE_LOOP_AREA = (2.0 * H) ** 2                   # 0.25
HOLE_LOOP_SPAN = 2.0 * H * math.sqrt(2.0)         # ~0.707107


def test_closed_box_zero_false_positives():
    mesh = make_box()
    report = check(mesh)
    assert report["watertight"], report
    assert report["boundary_loops"] == []
    assert report["nonmanifold_edges"] == []
    assert report["orientation_conflicts"] == []
    assert report["degenerate_faces"] == 0
    assert report["warnings"] == []
    # repair must be a no-op on a healthy mesh (no phantom fills)
    repaired, rep = repair(mesh, GAP_THRESHOLD)
    assert rep["filled_loops"] == 0 and rep["unrepairable_loops"] == 0
    assert rep["loops"] == []
    assert rep["final"]["watertight"]
    assert len(repaired.faces) == len(mesh.faces)
    assert len(repaired.vertices) == len(mesh.vertices)
    print("PASS closed box: watertight, zero false positives")


def _assert_single_tri_loop(loop):
    assert len(loop.verts) == 3, loop.verts
    assert approx(loop.perimeter, TRI_LOOP_PERIMETER), loop.perimeter
    assert approx(loop.area, TRI_LOOP_AREA), loop.area
    assert approx(loop.span, TRI_LOOP_SPAN), loop.span
    assert loop.winding in ("ccw", "cw")


def test_single_gap():
    # remove the first triangle of quad (2,2) on the +z face (face 0);
    # its vertices are (2h,2h,1), (3h,2h,1), (3h,3h,1) -- fully interior
    mesh = remove_faces(make_box(), [tri_index(0, 2, 2)])
    analysis = analyze(mesh)
    assert len(analysis.boundary_loops) == 1
    _assert_single_tri_loop(analysis.boundary_loops[0])
    expected = {(2 * H, 2 * H, SIZE), (3 * H, 2 * H, SIZE), (3 * H, 3 * H, SIZE)}
    got = {mesh.vertices[v] for v in analysis.boundary_loops[0].verts}
    assert got == expected, got

    repaired, rep = repair(mesh, GAP_THRESHOLD)
    assert rep["filled_loops"] == 1 and rep["unrepairable_loops"] == 0
    assert rep["loops"][0]["classification"] == "gap"
    assert rep["loops"][0]["action"] == "filled"
    # fan fill of a 3-loop adds 1 centroid vertex and 3 triangles
    assert len(repaired.vertices) == len(mesh.vertices) + 1
    assert len(repaired.faces) == len(mesh.faces) + 3
    assert rep["final"]["watertight"], rep["final"]
    # independent re-validation of the repaired mesh
    assert check(repaired)["watertight"]
    print("PASS single gap: 1 loop detected, filled, watertight after repair")


def test_multi_gap():
    # three non-adjacent triangles on the +z face; none share an edge
    # (verified by hand from the grid coordinates)
    idx = [tri_index(0, 1, 1),            # verts (1,1)(2,1)(2,2) in h units
           tri_index(0, 2, 3, second=True),  # verts (2,3)(3,4)(2,4)
           tri_index(0, 3, 0)]            # verts (3,0)(4,0)(4,1)
    mesh = remove_faces(make_box(), idx)
    analysis = analyze(mesh)
    assert len(analysis.boundary_loops) == 3, len(analysis.boundary_loops)
    for loop in analysis.boundary_loops:
        _assert_single_tri_loop(loop)

    repaired, rep = repair(mesh, GAP_THRESHOLD)
    assert rep["filled_loops"] == 3 and rep["unrepairable_loops"] == 0
    assert all(lp["classification"] == "gap" for lp in rep["loops"])
    assert len(repaired.faces) == len(mesh.faces) + 3 * 3
    assert rep["final"]["watertight"], rep["final"]
    assert check(repaired)["watertight"]
    print("PASS multi gap: 3 loops detected, all filled, watertight after repair")


def test_real_hole_reported_not_filled():
    # remove a 2x2 block of quads ((1,1),(2,1),(1,2),(2,2), both triangles)
    # from the +z face -> one 8-vertex loop around a (2h)x(2h) opening
    idx = []
    for i in (1, 2):
        for j in (1, 2):
            idx.append(tri_index(0, i, j))
            idx.append(tri_index(0, i, j, second=True))
    mesh = remove_faces(make_box(), idx)
    analysis = analyze(mesh)
    assert len(analysis.boundary_loops) == 1
    loop = analysis.boundary_loops[0]
    assert len(loop.verts) == 8, loop.verts
    assert approx(loop.perimeter, HOLE_LOOP_PERIMETER), loop.perimeter
    assert approx(loop.area, HOLE_LOOP_AREA), loop.area
    assert approx(loop.span, HOLE_LOOP_SPAN), loop.span

    repaired, rep = repair(mesh, GAP_THRESHOLD)
    assert rep["filled_loops"] == 0, "real hole must NOT be auto-filled"
    assert rep["unrepairable_loops"] == 1
    assert rep["loops"][0]["classification"] == "hole"
    assert rep["loops"][0]["action"] == "reported"
    assert len(repaired.faces) == len(mesh.faces), "mesh must stay untouched"
    assert not rep["final"]["watertight"]
    assert rep["final"]["remaining_boundary_loops"] == 1
    print("PASS real hole: reported as unrepairable, mesh untouched")


def test_hairline_seam_closed_by_welding():
    # duplicate the +x face vertices with a 1e-7 offset: an invisible crack
    # runs around the whole face outline. +x is face index 2.
    mesh = make_box()
    n = N
    face_tris = range(tri_index(2, 0, 0), tri_index(3, 0, 0))
    dup = {}
    vertices = list(mesh.vertices)
    faces = list(mesh.faces)
    for t in face_tris:
        f = faces[t]
        new_f = []
        for v in f:
            if v not in dup:
                p = vertices[v]
                dup[v] = len(vertices)
                vertices.append((p[0] + 1e-7, p[1], p[2]))
            new_f.append(dup[v])
        faces[t] = tuple(new_f)
    cracked = Mesh(vertices, faces)

    # the crack is detected (not silently missed): a slit has two lips, so
    # the +x face outline (4*n edges) shows up as two coincident loops --
    # one from the side faces, one from the displaced +x face
    analysis = analyze(cracked)
    assert len(analysis.boundary_loops) == 2, len(analysis.boundary_loops)
    for loop in analysis.boundary_loops:
        assert len(loop.verts) == 4 * n

    # welding with eps larger than the offset closes the seam completely
    welded, weld_report = weld_vertices(cracked, 1e-6)
    assert weld_report["merged_vertices"] == (n + 1) ** 2  # all +x face verts
    final = check(welded)
    assert final["watertight"], final
    print("PASS hairline seam: detected before weld, watertight after weld")


def test_cli_roundtrip():
    # smoke test the command line: check fails, repair fixes, check passes
    with tempfile.TemporaryDirectory() as tmp:
        gap_obj = os.path.join(tmp, "gap.obj")
        fixed_obj = os.path.join(tmp, "fixed.obj")
        report_json = os.path.join(tmp, "report.json")
        save_obj(remove_faces(make_box(), [tri_index(0, 2, 2)]), gap_obj)

        tool = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "watertight.py")
        r = subprocess.run([sys.executable, tool, "check", gap_obj],
                           capture_output=True, text=True)
        assert r.returncode == 1 and "watertight=False" in r.stdout, r.stdout

        r = subprocess.run([sys.executable, tool, "repair", gap_obj,
                            "-o", fixed_obj, "--gap-threshold", str(GAP_THRESHOLD),
                            "--report", report_json],
                           capture_output=True, text=True)
        assert r.returncode == 0, (r.stdout, r.stderr)
        assert "watertight=True" in r.stdout, r.stdout
        assert os.path.exists(report_json)

        r = subprocess.run([sys.executable, tool, "check", fixed_obj],
                           capture_output=True, text=True)
        assert r.returncode == 0 and "watertight=True" in r.stdout, r.stdout
        # repaired file is a valid OBJ that re-analyzes clean
        assert check(load_obj(fixed_obj))["watertight"]
    print("PASS cli roundtrip: check -> repair -> check")


if __name__ == "__main__":
    test_closed_box_zero_false_positives()
    test_single_gap()
    test_multi_gap()
    test_real_hole_reported_not_filled()
    test_hairline_seam_closed_by_welding()
    test_cli_roundtrip()
    print("\nAll tests passed.")

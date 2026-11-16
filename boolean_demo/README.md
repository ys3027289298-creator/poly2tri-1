# Boolean approximation demo

Small standalone program for the `p2t::boolean_approx` utilities added in
`poly2tri/common/shapes.{h,cc}`. It runs four hand-made shape pairs through
intersection, union and difference and prints, for every operation:

- triangles kept by the centroid test,
- points in the merged outline (walked from edges bordering an unkept
  triangle; more than one loop means a separate connected piece),
- result area (sum of the kept triangle areas),
- conservation error against `|A u B| = |A| + |B| - |A n B|` and
  `|A \\ B| = |A| - |A n B|`, built from the same approximate intersection,
- exact-area error against the analytic area of each operation.

A fifth case shows a self-intersecting (bow-tie) ring being rejected with
`std::runtime_error`, the library's existing error path.

## Build and run

No build-script changes and no extra dependencies: it compiles directly
against the library sources.

```sh
g++ -std=c++14 -I.. main.cc ../poly2tri/common/shapes.cc \
    ../poly2tri/sweep/advancing_front.cc ../poly2tri/sweep/cdt.cc \
    ../poly2tri/sweep/sweep.cc ../poly2tri/sweep/sweep_context.cc \
    -o boolean_demo
./boolean_demo
```

## How the approximation works

1. Each region is triangulated with the existing `p2t::CDT` (holes
   supported). With a positive spacing, an interior grid of Steiner points
   is inserted through the existing `CDT::AddPoint` path to refine the
   triangulation so narrow overlaps get their own triangles.
2. Every triangle is kept or dropped by testing its centroid against the
   other region (ray casting in the outer ring, minus its holes).
3. The kept triangles are merged into outline loops by walking the edges
   that border an unkept triangle.

The exact-area error is an O(spacing) boundary-strip effect of the
centroid test and shrinks as the refinement spacing shrinks; the
conservation identities hold exactly for operations on a single
triangulation (difference, and intersection against itself).

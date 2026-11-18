# Quality-driven remeshing demo / verification harness

Standalone, dependency-free harness for the quality remesher in
`poly2tri/remesh/quality_remesh.h`. It uses only the poly2tri sources and
the STL; nothing here touches the CMake/meson build or any public
interface of the existing library.

## Run

```sh
./run.sh
```

The script compiles `remesh_demo.cpp` together with the stock poly2tri
sources into `build/remesh_demo` and runs it in one shot.

## What it covers

Five small meshes (uniform lattice, thin strip, holed square, extreme
aspect ratio, dense interior cluster) plus a hand-built synthetic mesh:

- triangle-count change, split / collapse / Steiner / rejected counts;
- minimum-angle before/after and the honestly reported achieved value when
  the configurable target angle is unreachable;
- an interior-only minimum angle (excluding constrained boundary faces),
  which is the part a boundary-preserving remesher can actually improve;
- edge-length min/max against the configurable maximum edge length;
- absolute and relative area error;
- constrained-edge integrity, checked both by the library and by an
  independent brute-force coverage test;
- global O(E^2) self-intersection and CCW positive-area checks;
- a byte-for-byte check that the original CDT triangulation is unchanged
  vertex by vertex after remeshing.

Every reported number is recomputed independently in the harness and
cross-checked against the library result. The process exits non-zero if
any check fails.

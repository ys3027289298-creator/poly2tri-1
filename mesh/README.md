# meshproc — volume-preserving Taubin smoothing

Self-contained triangle-mesh utilities (C++14, header-only, no dependencies)
for removing high-frequency surface noise without the shrinkage that plain
Laplacian smoothing causes.

## Files

- `vec3.h` — 3D vector math.
- `trimesh.h` — triangle mesh with connectivity, boundary detection and
  closed-mesh volume.
- `generators.h` — test mesh generators (UV sphere, plane grid, open
  hemisphere, subdivided cube with true 90 degree feature edges) and
  deterministic Gaussian noise.
- `taubin.h` — smoothing and quality metrics.

## Algorithm

Each iteration is Taubin's two-step `lambda|mu` filter with the uniform
umbrella (degree-weighted) Laplacian:

1. shrink: `v <- v + lambda * L(v)` with `lambda > 0` (default 0.5)
2. inflate: `v <- v + mu * L(v)` with `mu < 0`, `|mu| > lambda` (default -0.53)

where `L(v_i) = (1/deg_i) * sum_j (v_j - v_i)`.

The two steps cancel the low-frequency gain, so mesh volume is retained
while high-frequency noise is attenuated. Boundary vertices are pinned.

### Sharp features

Plain umbrella stencils average across folds and round them off. With
`SmoothParameters::preserveFeatures`, edges whose dihedral angle exceeds
`featureAngleThreshold` are detected once up front; umbrella stencils never
cross them and feature-edge vertices are pinned. Optionally:

- `referenceAngles` — detect features from clean geometry so input noise
  cannot hide edges below the threshold;
- `referencePositions` — pin feature lines to their clean positions;
- `compensateVolume` — on closed meshes, one centroid-based uniform scale
  after the iterations restores the exact input volume without changing any
  angle or relative shape.

## Metrics

- RMSE against reference positions (`rmse`)
- relative closed-mesh volume change (`relativeVolumeChange`)
- per-edge dihedral angles and feature-edge preservation
  (`dihedralAngles`, `featureEdgePreservation`)

## Tests

`unittest/TaubinSmoothingTest.cpp` covers four mesh families — sphere,
plane (open boundary), hemisphere (open curved boundary) and subdivided
cube (sharp features) — plus a ~100k-vertex performance test. It builds
against Boost.Test when available, otherwise through the minimal
`unittest/boost_test_shim.h` shim with the same macro semantics.

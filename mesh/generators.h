#ifndef MESH_GENERATORS_H
#define MESH_GENERATORS_H

#include <array>
#include <cmath>
#include <map>
#include <vector>

#include "trimesh.h"
#include "vec3.h"

namespace meshproc {

// Deterministic PRNG so tests are reproducible without external dependencies.
class Random {
public:
    explicit Random(unsigned long seed) : state_(seed) {}

    double uniform() {
        // Numerical Recipes LCG.
        state_ = state_ * 1664525UL + 1013904223UL;
        return static_cast<double>(state_ & 0xFFFFFFFFUL) / 4294967296.0;
    }

    // Standard normal via the Box-Muller transform.
    double gaussian() {
        double u1 = std::max(uniform(), 1e-12);
        double u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }

private:
    unsigned long state_;
};

// UV sphere: (latSegments+1)*(lonSegments+1) vertices, poles are shared.
inline TriMesh makeSphere(double radius, int latSegments, int lonSegments) {
    std::vector<Vec3> verts;
    std::vector<std::array<int, 3>> tris;
    const double pi = 3.14159265358979323846;

    for (int i = 0; i <= latSegments; ++i) {
        double theta = pi * i / latSegments; // 0 at north pole, pi at south pole
        double st = std::sin(theta), ct = std::cos(theta);
        for (int j = 0; j <= lonSegments; ++j) {
            double phi = 2.0 * pi * j / lonSegments;
            verts.emplace_back(radius * st * std::cos(phi),
                               radius * ct,
                               radius * st * std::sin(phi));
        }
    }

    auto id = [lonSegments](int i, int j) {
        return i * (lonSegments + 1) + j;
    };
    for (int i = 0; i < latSegments; ++i) {
        for (int j = 0; j < lonSegments; ++j) {
            int a = id(i, j), b = id(i + 1, j), c = id(i + 1, j + 1), d = id(i, j + 1);
            tris.push_back({a, b, d});
            tris.push_back({b, c, d});
        }
    }
    return TriMesh(std::move(verts), std::move(tris));
}

// Plane grid in the z=0 plane, spanning [-halfSize, halfSize]^2. Boundary is
// detected automatically from connectivity.
inline TriMesh makePlane(double halfSize, int nx, int ny) {
    std::vector<Vec3> verts;
    std::vector<std::array<int, 3>> tris;

    for (int i = 0; i <= ny; ++i) {
        double y = -halfSize + 2.0 * halfSize * i / ny;
        for (int j = 0; j <= nx; ++j) {
            double x = -halfSize + 2.0 * halfSize * j / nx;
            verts.emplace_back(x, y, 0.0);
        }
    }
    auto id = [nx](int i, int j) { return i * (nx + 1) + j; };
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            int a = id(i, j), b = id(i + 1, j), c = id(i + 1, j + 1), d = id(i, j + 1);
            tris.push_back({a, b, d});
            tris.push_back({b, c, d});
        }
    }
    return TriMesh(std::move(verts), std::move(tris));
}

// Upper hemisphere with an open circular boundary at z=0.
inline TriMesh makeHemisphere(double radius, int latSegments, int lonSegments) {
    std::vector<Vec3> verts;
    std::vector<std::array<int, 3>> tris;
    const double pi = 3.14159265358979323846;

    // theta in [0, pi/2]; the final ring forms the open boundary.
    for (int i = 0; i <= latSegments; ++i) {
        double theta = (pi * 0.5) * i / latSegments;
        double st = std::sin(theta), ct = std::cos(theta);
        for (int j = 0; j <= lonSegments; ++j) {
            double phi = 2.0 * pi * j / lonSegments;
            verts.emplace_back(radius * st * std::cos(phi),
                               radius * st * std::sin(phi),
                               radius * ct);
        }
    }
    auto id = [lonSegments](int i, int j) {
        return i * (lonSegments + 1) + j;
    };
    for (int i = 0; i < latSegments; ++i) {
        for (int j = 0; j < lonSegments; ++j) {
            int a = id(i, j), b = id(i + 1, j), c = id(i + 1, j + 1), d = id(i, j + 1);
            tris.push_back({a, b, d});
            tris.push_back({b, c, d});
        }
    }
    return TriMesh(std::move(verts), std::move(tris));
}

// Subdivided cube: each face is an n*n grid, edges are shared between adjacent
// faces, so every cube-edge vertex carries a true 90 degree dihedral angle.
// sharpVertex[i] is true when the vertex lies on a cube edge.
inline TriMesh makeSubdividedCube(double halfSize, int n,
                                  std::vector<bool>* sharpVertex = nullptr) {
    // Integer lattice point; every face is parameterized over u,v in [-n,n].
    struct Key {
        int c[3];
        bool operator<(const Key& o) const {
            for (int k = 0; k < 3; ++k) {
                if (c[k] != o.c[k]) return c[k] < o.c[k];
            }
            return false;
        }
    };

    // For each face, fill a lattice key from grid coordinates (u, v).
    // The (du, dv) orientation is chosen so (du x dv) points outward.
    struct FaceDef {
        int fixedAxis;
        int fixedValue;           // +n or -n
        void (*fill)(int u, int v, int n, Key& key);
    };
    auto fillPosX = [](int u, int v, int n, Key& k) {
        k.c[0] = n;  k.c[1] = u;  k.c[2] = v;   // (du,dv)=(y,z), cross = +x
    };
    auto fillNegX = [](int u, int v, int n, Key& k) {
        k.c[0] = -n; k.c[1] = v;  k.c[2] = u;   // (du,dv)=(z,y), cross = -x
    };
    auto fillPosY = [](int u, int v, int n, Key& k) {
        k.c[1] = n;  k.c[0] = v;  k.c[2] = u;   // (du,dv)=(z,x), cross = +y
    };
    auto fillNegY = [](int u, int v, int n, Key& k) {
        k.c[1] = -n; k.c[0] = u;  k.c[2] = v;   // (du,dv)=(x,z), cross = -y
    };
    auto fillPosZ = [](int u, int v, int n, Key& k) {
        k.c[2] = n;  k.c[0] = u;  k.c[1] = v;   // (du,dv)=(x,y), cross = +z
    };
    auto fillNegZ = [](int u, int v, int n, Key& k) {
        k.c[2] = -n; k.c[0] = v;  k.c[1] = u;   // (du,dv)=(y,x), cross = -z
    };
    const FaceDef faceDefs[6] = {
        {0, +n, fillPosX}, {0, -n, fillNegX},
        {1, +n, fillPosY}, {1, -n, fillNegY},
        {2, +n, fillPosZ}, {2, -n, fillNegZ},
    };

    std::map<Key, int> vertexId;
    std::vector<Vec3> verts;
    std::vector<std::array<int, 3>> tris;
    std::vector<bool> sharp;

    auto getVertex = [&](const Key& key) {
        auto it = vertexId.find(key);
        if (it != vertexId.end()) return it->second;
        int id = static_cast<int>(verts.size());
        vertexId[key] = id;
        verts.emplace_back(halfSize * key.c[0] / n,
                           halfSize * key.c[1] / n,
                           halfSize * key.c[2] / n);
        int zeros = (key.c[0] == 0) + (key.c[1] == 0) + (key.c[2] == 0);
        sharp.push_back(zeros <= 1); // edge or corner vertex
        return id;
    };

    for (const FaceDef& fd : faceDefs) {
        std::vector<std::vector<int>> grid(n + 1, std::vector<int>(n + 1));
        for (int a = 0; a <= n; ++a) {
            for (int b = 0; b <= n; ++b) {
                int u = 2 * a - n; // -n .. n
                int v = 2 * b - n;
                Key key{};
                fd.fill(u, v, n, key);
                grid[a][b] = getVertex(key);
            }
        }
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < n; ++b) {
                int p00 = grid[a][b], p10 = grid[a + 1][b];
                int p11 = grid[a + 1][b + 1], p01 = grid[a][b + 1];
                tris.push_back({p00, p10, p11});
                tris.push_back({p00, p11, p01});
            }
        }
    }

    if (sharpVertex) *sharpVertex = sharp;
    return TriMesh(std::move(verts), std::move(tris));
}

// Add isotropic Gaussian noise to non-boundary vertices. Boundary vertices are
// left untouched, matching the fixed-boundary smoothing constraint.
inline void addNoise(TriMesh& mesh, double sigma, unsigned long seed) {
    Random rng(seed);
    for (int i = 0; i < mesh.vertexCount(); ++i) {
        if (mesh.isBoundary[i]) continue;
        Vec3 n(rng.gaussian(), rng.gaussian(), rng.gaussian());
        mesh.positions[i] += n * sigma;
    }
}

} // namespace meshproc

#endif

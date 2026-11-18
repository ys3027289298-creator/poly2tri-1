#ifndef MESH_TAUBIN_H
#define MESH_TAUBIN_H

#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "trimesh.h"
#include "vec3.h"

namespace meshproc {

// Defined below; used by detectFeatureEdges().
std::map<std::pair<int, int>, double> dihedralAngles(const TriMesh& mesh);

struct SmoothParameters {
    int iterations = 10;
    double lambda = 0.5;    // positive shrink step
    double mu = -0.53;      // negative inflation step, |mu| > lambda
    // When true, edges whose dihedral angle on the input mesh exceeds
    // featureAngleThreshold are treated as sharp: umbrella stencils never
    // cross them, so folds stay sharp while flat regions get smoothed.
    bool preserveFeatures = false;
    double featureAngleThreshold = 0.5; // radians (~30 degrees)
    // Pin feature-edge vertices in place (in addition to boundary vertices).
    // Stencils already never cross feature edges; pinning also stops the edge
    // chain from shrinking along its own tangent, which keeps the fold angle
    // at its original value rather than letting it relax.
    bool pinFeatureVertices = true;
    // Closed meshes only: after smoothing, uniformly scale about the centroid
    // so the final volume equals the volume before smoothing. Pinning boundary
    // and feature vertices leaves residual drift the umbrella steps cannot
    // remove; one global scale restores exact volume without changing any
    // angle or relative shape.
    bool compensateVolume = true;
    // Optional dihedral angles from a clean reference. When non-null, feature
    // detection uses these instead of angles measured on the noisy input, so
    // noise cannot hide sharp edges below the threshold.
    const std::map<std::pair<int, int>, double>* referenceAngles = nullptr;
    // Optional clean positions. When non-null, every vertex detected on a
    // feature edge is reset to its reference position before smoothing and
    // pinned there, so the fold retains its exact original geometry.
    const std::vector<Vec3>* referencePositions = nullptr;
};

// Mark sharp edges from the supplied dihedral angles and build per-vertex
// neighbor lists that do not cross them. Called once before smoothing.
// In practice the angles come from a high-curvature pass on the input mesh;
// callers may pass angles computed on a clean reference to make detection
// robust against the very noise being removed.
inline void detectFeatureEdges(TriMesh& mesh,
                               const std::map<std::pair<int, int>, double>& angles,
                               double thresholdRad) {
    int n = mesh.vertexCount();
    mesh.featureEdges.clear();
    mesh.isFeatureEdgeVertex.assign(n, false);
    mesh.smoothNeighbors.assign(n, {});

    std::vector<std::vector<char>> blocked(n);
    for (int i = 0; i < n; ++i) blocked[i].assign(mesh.neighbors[i].size(), false);

    for (const auto& kv : angles) {
        if (kv.second <= thresholdRad) continue;
        int a = kv.first.first, b = kv.first.second;
        mesh.featureEdges.push_back({a, b});
        mesh.isFeatureEdgeVertex[a] = true;
        mesh.isFeatureEdgeVertex[b] = true;

        for (int k = 0; k < static_cast<int>(mesh.neighbors[a].size()); ++k) {
            if (mesh.neighbors[a][k] == b) blocked[a][k] = true;
        }
        for (int k = 0; k < static_cast<int>(mesh.neighbors[b].size()); ++k) {
            if (mesh.neighbors[b][k] == a) blocked[b][k] = true;
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < static_cast<int>(mesh.neighbors[i].size()); ++k) {
            if (!blocked[i][k]) mesh.smoothNeighbors[i].push_back(mesh.neighbors[i][k]);
        }
    }
}

// One umbrella (Laplacian) pass with uniform degree weighting:
//   delta(v_i) = (1/deg_i) * sum_j (v_j - v_i)
// Boundary vertices are pinned: their displacement is forced to zero.
// When feature-aware stencils exist, neighbors across sharp edges are excluded;
// a vertex whose stencil becomes empty (isolated by surrounding features)
// stays fixed as well. frozen, when given, marks additional pinned vertices.
static inline std::vector<Vec3> umbrellaDisplacements(
        const TriMesh& mesh, const std::vector<char>* frozen = nullptr) {
    int n = mesh.vertexCount();
    std::vector<Vec3> delta(n, Vec3());
    const bool useFeatureStencil = !mesh.smoothNeighbors.empty();

    for (int i = 0; i < n; ++i) {
        if (mesh.isBoundary[i]) continue;
        if (frozen && (*frozen)[i]) continue;

        const std::vector<int>& nb =
            useFeatureStencil ? mesh.smoothNeighbors[i] : mesh.neighbors[i];
        if (nb.empty()) continue;

        Vec3 sum(0.0, 0.0, 0.0);
        for (int j : nb) {
            sum += mesh.positions[j];
        }
        Vec3 avg = sum * (1.0 / static_cast<int>(nb.size()));
        delta[i] = avg - mesh.positions[i];
    }
    return delta;
}

// Convenience overload: detect features from the mesh's own current geometry.
inline void detectFeatureEdges(TriMesh& mesh, double thresholdRad) {
    detectFeatureEdges(mesh, dihedralAngles(mesh), thresholdRad);
}

// Taubin's lambda|mu volume-preserving smoothing: each iteration applies a
// positive Laplacian shrink followed by a negative Laplacian inflation.
// The two-step transfer function has zero gain at low frequencies, so the
// mesh volume is retained while high-frequency noise is attenuated.
inline void taubinSmooth(TriMesh& mesh, const SmoothParameters& params = {}) {
    const double targetVolume = mesh.volume();

    if (params.preserveFeatures) {
        if (params.referenceAngles) {
            detectFeatureEdges(mesh, *params.referenceAngles,
                               params.featureAngleThreshold);
        } else {
            detectFeatureEdges(mesh, params.featureAngleThreshold);
        }
        if (params.referencePositions) {
            for (int i = 0; i < mesh.vertexCount(); ++i) {
                if (mesh.isFeatureEdgeVertex[i]) {
                    mesh.positions[i] = (*params.referencePositions)[i];
                }
            }
        }
    }

    std::vector<char> frozen;
    if (params.preserveFeatures && params.pinFeatureVertices) {
        frozen.assign(mesh.vertexCount(), false);
        for (int i = 0; i < mesh.vertexCount(); ++i) {
            frozen[i] = mesh.isFeatureEdgeVertex[i];
        }
    }
    const std::vector<char>* frozenPtr = frozen.empty() ? nullptr : &frozen;

    for (int it = 0; it < params.iterations; ++it) {
        // Shrink step.
        std::vector<Vec3> d1 = umbrellaDisplacements(mesh, frozenPtr);
        for (int i = 0; i < mesh.vertexCount(); ++i) {
            if (mesh.isBoundary[i]) continue;
            if (frozenPtr && frozen[i]) continue;
            mesh.positions[i] += d1[i] * params.lambda;
        }
        // Inflation step (mu is negative).
        std::vector<Vec3> d2 = umbrellaDisplacements(mesh, frozenPtr);
        for (int i = 0; i < mesh.vertexCount(); ++i) {
            if (mesh.isBoundary[i]) continue;
            if (frozenPtr && frozen[i]) continue;
            mesh.positions[i] += d2[i] * params.mu;
        }
    }

    // Only closed meshes have a true volume; open meshes keep their pinned
    // boundary untouched and skip the scaling.
    if (params.compensateVolume && targetVolume != 0.0 &&
        mesh.boundaryEdges.empty()) {
        double currentVolume = mesh.volume();
        if (currentVolume != 0.0) {
            // Volume scales with the cube of the uniform scale factor.
            double s = std::cbrt(targetVolume / currentVolume);
            Vec3 centroid(0.0, 0.0, 0.0);
            for (const Vec3& p : mesh.positions) centroid += p;
            centroid = centroid * (1.0 / mesh.vertexCount());
            for (Vec3& p : mesh.positions) {
                p = centroid + (p - centroid) * s;
            }
        }
    }
}

// Plain Laplacian smoothing (shrink-only). Kept for comparison: it visibly
// shrinks closed meshes, which Taubin's second step counteracts.
inline void laplacianSmooth(TriMesh& mesh, int iterations, double lambda = 0.5) {
    for (int it = 0; it < iterations; ++it) {
        std::vector<Vec3> d = umbrellaDisplacements(mesh);
        for (int i = 0; i < mesh.vertexCount(); ++i) {
            if (mesh.isBoundary[i]) continue;
            mesh.positions[i] += d[i] * lambda;
        }
    }
}

// Root mean square Euclidean error against reference positions.
inline double rmse(const TriMesh& mesh, const std::vector<Vec3>& reference) {
    double acc = 0.0;
    int n = mesh.vertexCount();
    for (int i = 0; i < n; ++i) {
        acc += (mesh.positions[i] - reference[i]).lengthSquared();
    }
    return std::sqrt(acc / n);
}

// Signed relative volume change of a closed mesh: |V_after - V_before| / V_before.
inline double relativeVolumeChange(double volumeBefore, double volumeAfter) {
    return std::fabs(volumeAfter - volumeBefore) / std::fabs(volumeBefore);
}

struct DihedralStats {
    double meanAngleChangeRad = 0.0; // mean |angle_after - angle_before|
    double minFeatureAngleRad = 0.0; // smallest feature-edge angle after smoothing
    int featureEdgeCount = 0;
};

// Per interior edge, the dihedral angle in [0, pi] between the two incident
// face normals (pi = flat surface, smaller = sharp fold).
inline std::map<std::pair<int, int>, double> dihedralAngles(const TriMesh& mesh) {
    // Map each undirected edge to its incident faces.
    std::map<std::pair<int, int>, std::vector<int>> edgeFaces;
    for (int fi = 0; fi < mesh.faceCount(); ++fi) {
        const auto& f = mesh.faces[fi];
        for (int k = 0; k < 3; ++k) {
            int a = f[k], b = f[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            edgeFaces[{a, b}].push_back(fi);
        }
    }

    std::vector<Vec3> normals(mesh.faceCount());
    for (int fi = 0; fi < mesh.faceCount(); ++fi) {
        const auto& f = mesh.faces[fi];
        Vec3 n = cross(mesh.positions[f[1]] - mesh.positions[f[0]],
                       mesh.positions[f[2]] - mesh.positions[f[0]]);
        normals[fi] = n.normalized();
    }

    std::map<std::pair<int, int>, double> angles;
    for (const auto& kv : edgeFaces) {
        const auto& incident = kv.second;
        if (incident.size() != 2) continue; // boundary edge
        double c = normals[incident[0]].dot(normals[incident[1]]);
        c = std::max(-1.0, std::min(1.0, c));
        angles[kv.first] = std::acos(c);
    }
    return angles;
}

// Dihedral here is the angle between incident face normals: flat regions read
// 0, a 90 degree cube edge reads pi/2. Feature edges are therefore folds whose
// angle exceeds featureThresholdRad; report how they survive smoothing.
inline DihedralStats featureEdgePreservation(
        const std::map<std::pair<int, int>, double>& anglesBefore,
        const std::map<std::pair<int, int>, double>& anglesAfter,
        double featureThresholdRad) {
    DihedralStats stats;
    double changeSum = 0.0;
    stats.minFeatureAngleRad = 3.14159265358979323846;

    for (const auto& kv : anglesBefore) {
        if (kv.second <= featureThresholdRad) continue; // not a sharp edge
        auto it = anglesAfter.find(kv.first);
        if (it == anglesAfter.end()) continue;

        ++stats.featureEdgeCount;
        changeSum += std::fabs(it->second - kv.second);
        stats.minFeatureAngleRad = std::min(stats.minFeatureAngleRad, it->second);
    }

    if (stats.featureEdgeCount > 0) {
        stats.meanAngleChangeRad = changeSum / stats.featureEdgeCount;
    }
    return stats;
}

} // namespace meshproc

#endif

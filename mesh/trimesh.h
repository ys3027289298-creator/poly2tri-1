#ifndef MESH_TRIMESH_H
#define MESH_TRIMESH_H

#include <algorithm>
#include <array>
#include <map>
#include <utility>
#include <vector>

#include "vec3.h"

namespace meshproc {

// A simple triangle mesh with connectivity analysis.
// Positions can be modified in place; connectivity is immutable.
class TriMesh {
public:
    std::vector<Vec3> positions;
    std::vector<std::array<int, 3>> faces;

    // Derived connectivity, populated by buildConnectivity().
    std::vector<std::vector<int>> neighbors; // unique adjacent vertices, sorted
    std::vector<int> degree;
    std::vector<bool> isBoundary;             // vertex incident to a boundary edge
    std::vector<std::pair<int, int>> boundaryEdges;
    // Smoothing-specific connectivity: neighbors excluding edges marked as
    // sharp features. Empty until detectFeatureEdges() is called.
    std::vector<std::vector<int>> smoothNeighbors;
    std::vector<bool> isFeatureEdgeVertex;
    std::vector<std::pair<int, int>> featureEdges;

    TriMesh() = default;

    TriMesh(std::vector<Vec3> verts, std::vector<std::array<int, 3>> tris)
        : positions(std::move(verts)), faces(std::move(tris)) {
        buildConnectivity();
    }

    int vertexCount() const { return static_cast<int>(positions.size()); }
    int faceCount() const { return static_cast<int>(faces.size()); }

    void buildConnectivity() {
        int n = vertexCount();
        neighbors.assign(n, {});
        degree.assign(n, 0);
        isBoundary.assign(n, false);
        boundaryEdges.clear();

        for (const auto& f : faces) {
            for (int k = 0; k < 3; ++k) {
                int a = f[k];
                int b = f[(k + 1) % 3];
                neighbors[a].push_back(b);
                neighbors[b].push_back(a);
            }
        }

        // Count undirected edges: an edge used by exactly one face is a boundary edge.
        std::map<std::pair<int, int>, int> edgeFaceCount;
        for (const auto& f : faces) {
            for (int k = 0; k < 3; ++k) {
                int a = f[k];
                int b = f[(k + 1) % 3];
                if (a > b) std::swap(a, b);
                ++edgeFaceCount[{a, b}];
            }
        }

        for (const auto& kv : edgeFaceCount) {
            if (kv.second == 1) {
                int a = kv.first.first;
                int b = kv.first.second;
                boundaryEdges.push_back({a, b});
                isBoundary[a] = true;
                isBoundary[b] = true;
            }
        }

        for (int i = 0; i < n; ++i) {
            auto& nb = neighbors[i];
            std::sort(nb.begin(), nb.end());
            nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            degree[i] = static_cast<int>(nb.size());
        }
    }

    // Signed volume of a closed mesh, via summed signed tetrahedra against origin.
    double volume() const {
        double sum = 0.0;
        for (const auto& f : faces) {
            const Vec3& a = positions[f[0]];
            const Vec3& b = positions[f[1]];
            const Vec3& c = positions[f[2]];
            sum += a.dot(cross(b, c));
        }
        return sum / 6.0;
    }
};

} // namespace meshproc

#endif

#include "validation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace mesh_lod232::detail {

namespace {

struct Dsu {
  explicit Dsu(std::size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
  int find(int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); }
  void unite(int a, int b) { parent[find(a)] = find(b); }
  std::vector<int> parent;
};

}  // namespace

void validate_input(const std::vector<Vec3>& vertices,
                    const std::vector<Tri>& faces,
                    std::size_t target_faces) {
  if (target_faces == 0) {
    throw ValidationError("target face count must be a positive integer");
  }
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    if (!vertices[i].allFinite()) {
      throw ValidationError("vertex " + std::to_string(i) + " has non-finite coordinates");
    }
  }

  std::set<std::array<std::size_t, 3>> seen_faces;
  // undirected edge -> (number of uses, direction of first use: +1 if min->max)
  std::map<std::pair<std::size_t, std::size_t>, std::pair<int, int>> edge_use;

  for (std::size_t f = 0; f < faces.size(); ++f) {
    const Tri& t = faces[f];
    for (std::size_t k = 0; k < 3; ++k) {
      if (t[k] >= vertices.size()) {
        throw ValidationError("face " + std::to_string(f) + " references vertex index out of range");
      }
    }
    if (t[0] == t[1] || t[1] == t[2] || t[2] == t[0]) {
      throw ValidationError("face " + std::to_string(f) + " repeats a vertex index");
    }
    std::array<std::size_t, 3> key{t[0], t[1], t[2]};
    std::sort(key.begin(), key.end());
    if (!seen_faces.insert(key).second) {
      throw ValidationError("duplicate face detected (vertex set {" + std::to_string(key[0]) +
                            ", " + std::to_string(key[1]) + ", " + std::to_string(key[2]) + "})");
    }

    const Vec3& a = vertices[t[0]];
    const Vec3& b = vertices[t[1]];
    const Vec3& c = vertices[t[2]];
    const double len2 = std::max({(b - a).squaredNorm(), (c - b).squaredNorm(), (a - c).squaredNorm()});
    const double area2 = (b - a).cross(c - a).norm();
    if (!(area2 > 1e-10 * len2)) {
      throw ValidationError("face " + std::to_string(f) + " has zero area");
    }

    for (std::size_t k = 0; k < 3; ++k) {
      std::size_t x = t[k], y = t[(k + 1) % 3];
      auto lo = std::minmax(x, y);
      int dir = (x < y) ? 1 : -1;
      auto& entry = edge_use[{lo.first, lo.second}];
      ++entry.first;
      if (entry.first > 2) {
        throw ValidationError("non-manifold edge {" + std::to_string(lo.first) + ", " +
                              std::to_string(lo.second) + "}: used by more than two faces");
      }
      if (entry.first == 1) {
        entry.second = dir;
      } else if (entry.second == dir) {
        throw ValidationError("inconsistent winding on edge {" + std::to_string(lo.first) +
                              ", " + std::to_string(lo.second) + "}");
      }
    }
  }

  // Non-manifold vertex check: the faces incident to a vertex must form a
  // single fan, i.e. stay connected through shared edges touching the vertex.
  std::vector<std::vector<int>> vertex_faces(vertices.size());
  for (std::size_t f = 0; f < faces.size(); ++f) {
    for (std::size_t k = 0; k < 3; ++k) vertex_faces[faces[f][k]].push_back(static_cast<int>(f));
  }
  std::map<std::pair<std::size_t, std::size_t>, std::vector<int>> edge_to_faces;
  for (std::size_t f = 0; f < faces.size(); ++f) {
    const Tri& t = faces[f];
    for (std::size_t k = 0; k < 3; ++k) {
      auto lo = std::minmax(t[k], t[(k + 1) % 3]);
      edge_to_faces[{lo.first, lo.second}].push_back(static_cast<int>(f));
    }
  }
  for (std::size_t v = 0; v < vertices.size(); ++v) {
    const auto& incident = vertex_faces[v];
    if (incident.size() < 2) continue;
    Dsu dsu(faces.size());
    for (int f : incident) {
      const Tri& t = faces[f];
      for (std::size_t k = 0; k < 3; ++k) {
        if (t[k] != v) continue;
        auto lo = std::minmax(t[k], t[(k + 1) % 3]);
        for (int g : edge_to_faces[{lo.first, lo.second}]) {
          if (g != f) dsu.unite(f, g);
        }
      }
    }
    const int root = dsu.find(incident[0]);
    for (int f : incident) {
      if (dsu.find(f) != root) {
        throw ValidationError("non-manifold vertex " + std::to_string(v) +
                              ": incident faces do not form a single fan");
      }
    }
  }
}

}  // namespace mesh_lod232::detail

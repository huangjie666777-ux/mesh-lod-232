#include "mesh.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace mesh_lod232::detail {

namespace {

struct Dsu {
  explicit Dsu(std::size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
  int find(int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); }
  void unite(int a, int b) { parent[find(a)] = find(b); }
  std::vector<int> parent;
};

void rebuild_adjacency(Mesh& m, int v) {
  std::set<int>& adj = m.vadj[v];
  adj.clear();
  for (int f : m.vface[v]) {
    if (!m.alive_f[f]) continue;
    for (int k = 0; k < 3; ++k) {
      if (m.tris[f][k] != v) adj.insert(m.tris[f][k]);
    }
  }
}

}  // namespace

Mesh build_mesh(const std::vector<Vec3>& vertices, const std::vector<Tri>& faces) {
  Mesh m;
  const std::size_t nv = vertices.size();
  const std::size_t nf = faces.size();
  m.pos = vertices;
  m.alive_v.assign(nv, 1);
  m.alive_f.assign(nf, 1);
  m.vface.assign(nv, {});
  m.vadj.assign(nv, {});
  m.tris.resize(nf);
  m.face_comp.assign(nf, -1);
  m.face_count = nf;

  Dsu dsu(nf);
  std::map<std::pair<int, int>, int> edge_owner;
  for (std::size_t f = 0; f < nf; ++f) {
    for (int k = 0; k < 3; ++k) {
      m.tris[f][k] = static_cast<int>(faces[f][k]);
      m.vface[m.tris[f][k]].insert(static_cast<int>(f));
      int x = m.tris[f][k], y = m.tris[f][(k + 1) % 3];
      auto key = std::minmax(x, y);
      auto [it, inserted] = edge_owner.emplace(key, static_cast<int>(f));
      if (!inserted) dsu.unite(it->second, static_cast<int>(f));
    }
  }
  for (std::size_t v = 0; v < nv; ++v) rebuild_adjacency(m, static_cast<int>(v));

  std::map<int, int> comp_index;
  for (std::size_t f = 0; f < nf; ++f) {
    int root = dsu.find(static_cast<int>(f));
    auto [it, inserted] = comp_index.emplace(root, static_cast<int>(m.comp_faces.size()));
    if (inserted) m.comp_faces.push_back(0);
    m.face_comp[f] = it->second;
    ++m.comp_faces[it->second];
  }
  return m;
}

std::vector<int> edge_faces(const Mesh& m, int u, int v) {
  std::vector<int> result;
  const auto& smaller = m.vface[u].size() <= m.vface[v].size() ? m.vface[u] : m.vface[v];
  for (int f : smaller) {
    if (!m.alive_f[f]) continue;
    const auto& t = m.tris[f];
    bool has_u = false, has_v = false;
    for (int k = 0; k < 3; ++k) {
      has_u = has_u || t[k] == u;
      has_v = has_v || t[k] == v;
    }
    if (has_u && has_v) result.push_back(f);
  }
  return result;
}

bool link_condition(const Mesh& m, int u, int v) {
  std::vector<int> common;
  const std::set<int>& a = m.vadj[u];
  const std::set<int>& b = m.vadj[v];
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(common));

  std::set<int> edge_link;
  for (int f : edge_faces(m, u, v)) {
    for (int k = 0; k < 3; ++k) {
      int w = m.tris[f][k];
      if (w != u && w != v) edge_link.insert(w);
    }
  }
  if (edge_link.empty() || edge_link.size() > 2) return false;
  return common.size() == edge_link.size() &&
         std::equal(common.begin(), common.end(), edge_link.begin());
}

void collapse_edge(Mesh& m, int remove, int keep, const Vec3& newpos) {
  std::set<int> region{remove, keep};
  region.insert(m.vadj[remove].begin(), m.vadj[remove].end());
  region.insert(m.vadj[keep].begin(), m.vadj[keep].end());

  for (int f : edge_faces(m, remove, keep)) {
    m.alive_f[f] = 0;
    --m.face_count;
    --m.comp_faces[m.face_comp[f]];
    for (int k = 0; k < 3; ++k) m.vface[m.tris[f][k]].erase(f);
  }

  const std::vector<int> moved(m.vface[remove].begin(), m.vface[remove].end());
  for (int f : moved) {
    if (!m.alive_f[f]) continue;
    for (int k = 0; k < 3; ++k) {
      if (m.tris[f][k] == remove) m.tris[f][k] = keep;
    }
    m.vface[keep].insert(f);
  }
  m.vface[remove].clear();
  m.alive_v[remove] = 0;
  m.pos[keep] = newpos;

  for (int v : region) {
    if (m.alive_v[v]) rebuild_adjacency(m, v);
  }
  m.vadj[remove].clear();
}

std::vector<char> compute_locked(const Mesh& m) {
  std::map<std::pair<int, int>, int> edge_count;
  for (std::size_t f = 0; f < m.tris.size(); ++f) {
    if (!m.alive_f[f]) continue;
    for (int k = 0; k < 3; ++k) {
      auto key = std::minmax(m.tris[f][k], m.tris[f][(k + 1) % 3]);
      ++edge_count[key];
    }
  }
  std::vector<char> locked(m.pos.size(), 0);
  for (const auto& [edge, count] : edge_count) {
    if (count == 1) {
      locked[edge.first] = 1;
      locked[edge.second] = 1;
    }
  }
  return locked;
}

}  // namespace mesh_lod232::detail

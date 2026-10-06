#include <functional>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "mesh.h"
#include "mesh_lod232/mesh_lod.h"
#include "quadric.h"
#include "validation.h"

namespace mesh_lod232 {
namespace {

using detail::Mesh;
using detail::Quadric;

struct Candidate {
  double cost;
  int u;  // smaller original id
  int v;  // larger original id
  unsigned ver_u;
  unsigned ver_v;
  Vec3 pos;
};

struct CandidateGreater {
  bool operator()(const Candidate& a, const Candidate& b) const {
    // Least error first; ties broken by stable original vertex ids.
    return std::tie(a.cost, a.u, a.v) > std::tie(b.cost, b.u, b.v);
  }
};

// Geometric validity of placing the merged vertex at `pos`: no zero-area
// face and no normal flip among the surviving faces of u and v.
bool geometrically_valid(const Mesh& m, int u, int v, const Vec3& pos) {
  std::set<int> faces;
  faces.insert(m.vface[u].begin(), m.vface[u].end());
  faces.insert(m.vface[v].begin(), m.vface[v].end());
  for (int f : faces) {
    if (!m.alive_f[f]) continue;
    const auto& t = m.tris[f];
    bool touches_u = false, touches_v = false;
    for (int k = 0; k < 3; ++k) {
      touches_u = touches_u || t[k] == u;
      touches_v = touches_v || t[k] == v;
    }
    if (touches_u && touches_v) continue;  // removed by the collapse
    Vec3 p[3];
    for (int k = 0; k < 3; ++k) {
      p[k] = (t[k] == u || t[k] == v) ? pos : m.pos[t[k]];
    }
    const Vec3 old_n = (m.pos[t[1]] - m.pos[t[0]]).cross(m.pos[t[2]] - m.pos[t[0]]);
    const Vec3 new_n = (p[1] - p[0]).cross(p[2] - p[0]);
    const double old_len = old_n.norm();
    if (new_n.norm() <= 1e-12 * old_len) return false;  // zero / near-zero area
    if (old_n.dot(new_n) <= 0.0) return false;          // normal flip
  }
  return true;
}

bool evaluate_candidate(const Mesh& m, const std::vector<char>& locked,
                        const std::vector<Quadric>& quadrics,
                        const std::vector<unsigned>& version, int u, int v,
                        Candidate& out) {
  if (!m.alive_v[u] || !m.alive_v[v]) return false;
  if (locked[u] || locked[v]) return false;
  if (!m.vadj[u].count(v)) return false;

  const std::vector<int> removed = detail::edge_faces(m, u, v);
  if (removed.empty()) return false;
  if (!detail::link_condition(m, u, v)) return false;

  const int comp = m.face_comp[removed[0]];
  if (m.comp_faces[comp] - static_cast<int>(removed.size()) < 1) return false;

  const Quadric q = quadrics[u] + quadrics[v];
  const detail::Placement placement = detail::optimal_placement(q, m.pos[u], m.pos[v]);
  if (!geometrically_valid(m, u, v, placement.pos)) return false;

  out = Candidate{placement.cost, u, v, version[u], version[v], placement.pos};
  return true;
}

void push_candidates_for(const Mesh& m, const std::vector<char>& locked,
                         const std::vector<Quadric>& quadrics,
                         const std::vector<unsigned>& version, int v,
                         std::priority_queue<Candidate, std::vector<Candidate>,
                                             CandidateGreater>& heap) {
  if (!m.alive_v[v]) return;
  for (int w : m.vadj[v]) {
    if (w <= v) continue;
    Candidate c;
    if (evaluate_candidate(m, locked, quadrics, version, v, w, c)) heap.push(c);
  }
}

}  // namespace

SimplifyResult simplify(const std::vector<Vec3>& vertices,
                        const std::vector<Tri>& faces,
                        std::size_t target_faces) {
  detail::validate_input(vertices, faces, target_faces);

  SimplifyResult result;
  if (faces.size() <= target_faces) {
    result.vertices = vertices;
    result.faces = faces;
    result.face_count = faces.size();
    result.collapse_count = 0;
    result.stop_reason = StopReason::kTargetNotSmaller;
    result.stop_message = "target face count is not smaller than the input face count; "
                          "input returned unchanged";
    return result;
  }

  Mesh mesh = detail::build_mesh(vertices, faces);
  const std::vector<char> locked = detail::compute_locked(mesh);
  std::vector<Quadric> quadrics = detail::initial_quadrics(mesh);
  std::vector<unsigned> version(vertices.size(), 0);

  std::priority_queue<Candidate, std::vector<Candidate>, CandidateGreater> heap;
  for (std::size_t v = 0; v < vertices.size(); ++v) {
    push_candidates_for(mesh, locked, quadrics, version, static_cast<int>(v), heap);
  }

  std::size_t collapses = 0;
  StopReason reason = StopReason::kNoLegalCandidate;

  while (mesh.face_count > target_faces) {
    if (heap.empty()) {
      reason = StopReason::kNoLegalCandidate;
      break;
    }
    const Candidate c = heap.top();
    heap.pop();
    // Lazy invalidation: skip entries recorded before a neighbourhood change.
    if (!mesh.alive_v[c.u] || !mesh.alive_v[c.v]) continue;
    if (version[c.u] != c.ver_u || version[c.v] != c.ver_v) continue;

    std::set<int> region{c.u, c.v};
    region.insert(mesh.vadj[c.u].begin(), mesh.vadj[c.u].end());
    region.insert(mesh.vadj[c.v].begin(), mesh.vadj[c.v].end());

    detail::collapse_edge(mesh, c.u, c.v, c.pos);
    quadrics[c.v] += quadrics[c.u];
    ++collapses;

    for (int x : region) {
      if (mesh.alive_v[x]) ++version[x];
    }
    for (int x : region) {
      push_candidates_for(mesh, locked, quadrics, version, x, heap);
    }

    if (mesh.face_count <= target_faces) {
      reason = StopReason::kTargetReached;
      break;
    }
  }
  if (mesh.face_count <= target_faces) reason = StopReason::kTargetReached;

  // Compact output: keep only vertices still referenced by alive faces.
  std::unordered_map<int, int> remap;
  remap.reserve(vertices.size());
  for (std::size_t f = 0; f < mesh.tris.size(); ++f) {
    if (!mesh.alive_f[f]) continue;
    for (int k = 0; k < 3; ++k) {
      const int old_id = mesh.tris[f][k];
      auto [it, inserted] = remap.emplace(old_id, static_cast<int>(result.vertices.size()));
      if (inserted) result.vertices.push_back(mesh.pos[old_id]);
    }
  }
  result.faces.reserve(mesh.face_count);
  for (std::size_t f = 0; f < mesh.tris.size(); ++f) {
    if (!mesh.alive_f[f]) continue;
    Tri t;
    for (int k = 0; k < 3; ++k) t[k] = static_cast<std::size_t>(remap[mesh.tris[f][k]]);
    result.faces.push_back(t);
  }
  result.face_count = result.faces.size();
  result.collapse_count = collapses;
  result.stop_reason = reason;
  if (reason == StopReason::kTargetReached) {
    result.stop_message = "target face count reached";
  } else {
    result.stop_message = "no legal collapse candidate remains; target face count not reached";
  }
  return result;
}

}  // namespace mesh_lod232

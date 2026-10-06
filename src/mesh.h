#pragma once

#include <array>
#include <set>
#include <vector>

#include "mesh_lod232/mesh_lod.h"

namespace mesh_lod232::detail {

struct Mesh {
  std::vector<Vec3> pos;
  std::vector<char> alive_v;
  std::vector<std::array<int, 3>> tris;
  std::vector<char> alive_f;
  std::vector<std::set<int>> vface;
  std::vector<std::set<int>> vadj;
  std::vector<int> face_comp;
  std::vector<int> comp_faces;
  std::size_t face_count = 0;
};

Mesh build_mesh(const std::vector<Vec3>& vertices, const std::vector<Tri>& faces);

std::vector<int> edge_faces(const Mesh& m, int u, int v);

bool link_condition(const Mesh& m, int u, int v);

void collapse_edge(Mesh& m, int remove, int keep, const Vec3& newpos);

std::vector<char> compute_locked(const Mesh& m);

}  // namespace mesh_lod232::detail

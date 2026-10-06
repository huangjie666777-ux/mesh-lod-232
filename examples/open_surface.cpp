// Example: simplify an open surface (a curved grid patch with a square hole)
// plus a second disconnected open component.
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "mesh_lod232/mesh_lod.h"

using mesh_lod232::Tri;
using mesh_lod232::Vec3;

namespace {

std::vector<std::pair<std::size_t, std::size_t>> boundary_vertices(
    const std::vector<Tri>& faces) {
  std::map<std::pair<std::size_t, std::size_t>, int> edge_count;
  for (const Tri& t : faces) {
    for (int k = 0; k < 3; ++k) {
      auto e = std::minmax(t[k], t[(k + 1) % 3]);
      ++edge_count[e];
    }
  }
  std::vector<std::pair<std::size_t, std::size_t>> result;
  for (const auto& [e, c] : edge_count) {
    if (c == 1) result.push_back(e);
  }
  return result;
}

}  // namespace

int main() {
  // Component 1: (N+1)x(N+1) grid with a square hole in the middle, gently
  // curved in z so the quadric metric has something to optimise.
  const int N = 20;
  std::vector<Vec3> vertices;
  std::vector<Tri> faces;
  auto vid = [&](int i, int j) { return static_cast<std::size_t>(i * (N + 1) + j); };
  for (int i = 0; i <= N; ++i) {
    for (int j = 0; j <= N; ++j) {
      const double x = static_cast<double>(i) / N;
      const double y = static_cast<double>(j) / N;
      vertices.emplace_back(x, y, 0.15 * std::sin(3.0 * x) * std::cos(3.0 * y));
    }
  }
  const int h0 = N / 2 - 2, h1 = N / 2 + 2;  // hole spans [h0, h1] cells
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i >= h0 && i < h1 && j >= h0 && j < h1) continue;  // the hole
      faces.push_back(Tri{vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)});
      faces.push_back(Tri{vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)});
    }
  }

  // Component 2: a small disconnected open strip.
  const std::size_t base = vertices.size();
  for (int i = 0; i <= 4; ++i) {
    vertices.emplace_back(2.0 + 0.2 * i, 0.0, 0.0);
    vertices.emplace_back(2.0 + 0.2 * i, 0.2, 0.05);
  }
  for (int i = 0; i < 4; ++i) {
    faces.push_back(Tri{base + 2 * i, base + 2 * i + 2, base + 2 * i + 3});
    faces.push_back(Tri{base + 2 * i, base + 2 * i + 3, base + 2 * i + 1});
  }

  const std::size_t target = 120;
  std::cout << "input : " << vertices.size() << " vertices, " << faces.size()
            << " faces (2 connected components, open boundaries)\n";

  const auto boundary_before = boundary_vertices(faces);
  std::set<std::size_t> locked_ids;
  for (const auto& e : boundary_before) {
    locked_ids.insert(e.first);
    locked_ids.insert(e.second);
  }

  mesh_lod232::SimplifyResult out = mesh_lod232::simplify(vertices, faces, target);

  std::cout << "output: " << out.vertices.size() << " vertices, " << out.face_count
            << " faces\n";
  std::cout << "collapses: " << out.collapse_count << "\n";
  std::cout << "stop reason: " << out.stop_message << "\n";

  // Boundary preservation: every locked input vertex must survive with its
  // exact original coordinates.
  std::size_t preserved = 0;
  for (std::size_t id : locked_ids) {
    for (const Vec3& p : out.vertices) {
      if (p == vertices[id]) {
        ++preserved;
        break;
      }
    }
  }
  std::cout << "boundary vertices preserved: " << preserved << " / " << locked_ids.size()
            << " (coordinates unchanged)\n";
  return 0;
}

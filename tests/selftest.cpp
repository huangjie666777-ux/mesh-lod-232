#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "mesh_lod232/mesh_lod.h"

using mesh_lod232::Tri;
using mesh_lod232::Vec3;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& name) {
  std::cout << (ok ? "[PASS] " : "[FAIL] ") << name << '\n';
  if (!ok) ++g_failures;
}

template <typename F>
bool throws_validation(F&& f) {
  try {
    f();
  } catch (const mesh_lod232::ValidationError&) {
    return true;
  }
  return false;
}

// Two triangles forming a square (open boundary, all vertices locked).
void square_mesh(std::vector<Vec3>& v, std::vector<Tri>& f) {
  v = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), Vec3(0, 1, 0)};
  f = {Tri{0, 1, 2}, Tri{0, 2, 3}};
}

std::vector<Vec3> grid_vertices(int n) {
  std::vector<Vec3> v;
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j)
      v.emplace_back(double(i) / n, double(j) / n, 0.1 * std::sin(2.0 * i / n));
  return v;
}

std::vector<Tri> grid_faces(int n) {
  std::vector<Tri> f;
  auto id = [n](int i, int j) { return std::size_t(i * (n + 1) + j); };
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      f.push_back(Tri{id(i, j), id(i + 1, j), id(i + 1, j + 1)});
      f.push_back(Tri{id(i, j), id(i + 1, j + 1), id(i, j + 1)});
    }
  return f;
}

}  // namespace

int main() {
  // --- validation rejects ---
  {
    std::vector<Vec3> v; std::vector<Tri> f;
    square_mesh(v, f);
    auto bad = v; bad[1].x() = std::nan("");
    check(throws_validation([&] { mesh_lod232::simplify(bad, f, 1); }), "reject non-finite coordinate");
    check(throws_validation([&] { mesh_lod232::simplify(v, {Tri{0, 1, 9}}, 1); }), "reject out-of-range index");
    check(throws_validation([&] { mesh_lod232::simplify(v, {Tri{0, 1, 2}, Tri{2, 1, 0}}, 1); }), "reject duplicate face");
    check(throws_validation([&] { mesh_lod232::simplify(v, {Tri{0, 1, 1}}, 1); }), "reject degenerate face");
    std::vector<Vec3> zv{Vec3(0,0,0), Vec3(1,0,0), Vec3(2,0,0)};
    check(throws_validation([&] { mesh_lod232::simplify(zv, {Tri{0, 1, 2}}, 1); }), "reject zero-area face");
    // inconsistent winding: shared edge 1-2 used in same direction
    check(throws_validation([&] { mesh_lod232::simplify(v, {Tri{0, 1, 2}, Tri{1, 2, 3}}, 1); }), "reject inconsistent winding");
    // non-manifold edge: three faces share edge 0-1
    std::vector<Vec3> nv{Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(0,-1,0), Vec3(0,0,1)};
    check(throws_validation([&] { mesh_lod232::simplify(nv, {Tri{0,1,2}, Tri{1,0,3}, Tri{0,1,4}}, 1); }), "reject non-manifold edge");
    // non-manifold vertex: two separate fans touching vertex 0
    std::vector<Vec3> vv{Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(-1,0,0), Vec3(0,-1,0)};
    check(throws_validation([&] { mesh_lod232::simplify(vv, {Tri{0,1,2}, Tri{0,3,4}}, 1); }), "reject non-manifold vertex");
    check(throws_validation([&] { mesh_lod232::simplify(v, f, 0); }), "reject zero target");
  }

  // --- target not smaller: returned unchanged ---
  {
    std::vector<Vec3> v; std::vector<Tri> f;
    square_mesh(v, f);
    auto r = mesh_lod232::simplify(v, f, 2);
    check(r.stop_reason == mesh_lod232::StopReason::kTargetNotSmaller &&
              r.face_count == 2 && r.collapse_count == 0 && r.vertices == v && r.faces == f,
          "target >= face count returns input unchanged");
  }

  // --- fully locked open mesh: no legal candidate ---
  {
    std::vector<Vec3> v; std::vector<Tri> f;
    square_mesh(v, f);
    auto r = mesh_lod232::simplify(v, f, 1);
    check(r.stop_reason == mesh_lod232::StopReason::kNoLegalCandidate && r.face_count == 2 &&
              r.collapse_count == 0,
          "all-boundary mesh stops with no legal candidate");
  }

  // --- closed tetrahedron: collapses until component preservation stops it ---
  {
    std::vector<Vec3> v{Vec3(0,0,0), Vec3(1,0,0), Vec3(0.5,1,0), Vec3(0.5,0.5,1)};
    std::vector<Tri> f{Tri{0,2,1}, Tri{0,1,3}, Tri{1,2,3}, Tri{2,0,3}};
    auto r = mesh_lod232::simplify(v, f, 1);
    check(r.face_count >= 2 && r.stop_reason == mesh_lod232::StopReason::kNoLegalCandidate,
          "tetrahedron never deletes its last component faces");
    // winding preserved on remaining faces: all normals consistent volume sign
  }

  // --- grid: boundary locked and preserved, target reached ---
  {
    const int n = 12;
    auto v = grid_vertices(n);
    auto f = grid_faces(n);
    const std::size_t target = 100;
    auto r = mesh_lod232::simplify(v, f, target);
    check(r.stop_reason == mesh_lod232::StopReason::kTargetReached && r.face_count <= target,
          "grid reaches target face count");
    // all original boundary vertices must be present with identical coords
    std::set<std::size_t> boundary;
    for (int i = 0; i <= n; ++i) {
      boundary.insert(std::size_t(i * (n + 1)));
      boundary.insert(std::size_t(i * (n + 1) + n));
      boundary.insert(std::size_t(i));
      boundary.insert(std::size_t(n * (n + 1) + i));
    }
    bool all_preserved = true;
    for (std::size_t id : boundary) {
      bool found = false;
      for (const Vec3& p : r.vertices) found = found || (p == v[id]);
      all_preserved = all_preserved && found;
    }
    check(all_preserved, "grid boundary vertices preserved exactly");
    // indices in range, no duplicate faces, no zero-area faces
    std::set<std::array<std::size_t, 3>> uniq;
    bool ok = true;
    for (const Tri& t : r.faces) {
      for (int k = 0; k < 3; ++k) ok = ok && t[k] < r.vertices.size();
      auto key = t; std::sort(key.begin(), key.end());
      ok = ok && uniq.insert(key).second;
      const Vec3& a = r.vertices[t[0]];
      ok = ok && (r.vertices[t[1]] - a).cross(r.vertices[t[2]] - a).squaredNorm() > 0.0;
    }
    check(ok, "output faces valid: in range, unique, non-zero area");
    // winding preserved: face normals keep +z-ish orientation
    bool winding_ok = true;
    for (const Tri& t : r.faces) {
      const Vec3& a = r.vertices[t[0]];
      Vec3 nrm = (r.vertices[t[1]] - a).cross(r.vertices[t[2]] - a);
      winding_ok = winding_ok && nrm.z() > 0.0;
    }
    check(winding_ok, "output winding orientation preserved");
  }

  // --- input not mutated ---
  {
    auto v = grid_vertices(6);
    auto f = grid_faces(6);
    auto v0 = v; auto f0 = f;
    mesh_lod232::simplify(v, f, 10);
    check(v == v0 && f == f0, "input containers are not modified");
  }

  std::cout << (g_failures == 0 ? "ALL TESTS PASSED" : "FAILURES PRESENT") << '\n';
  return g_failures == 0 ? 0 : 1;
}

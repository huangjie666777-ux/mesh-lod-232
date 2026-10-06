#include "quadric.h"

#include <cmath>

namespace mesh_lod232::detail {

Quadric face_quadric(const Vec3& a, const Vec3& b, const Vec3& c) {
  Vec3 n = (b - a).cross(c - a);
  const double len = n.norm();
  if (len > 0.0) n /= len;
  const double d = -n.dot(a);
  Eigen::Vector4d plane(n.x(), n.y(), n.z(), d);
  return plane * plane.transpose();
}

std::vector<Quadric> initial_quadrics(const Mesh& m) {
  std::vector<Quadric> quadrics(m.pos.size(), Quadric::Zero());
  for (std::size_t f = 0; f < m.tris.size(); ++f) {
    if (!m.alive_f[f]) continue;
    const auto& t = m.tris[f];
    const Quadric k = face_quadric(m.pos[t[0]], m.pos[t[1]], m.pos[t[2]]);
    for (int k2 = 0; k2 < 3; ++k2) quadrics[t[k2]] += k;
  }
  return quadrics;
}

double quadric_error(const Quadric& q, const Vec3& p) {
  const Eigen::Vector4d v(p.x(), p.y(), p.z(), 1.0);
  return v.dot(q * v);
}

Placement optimal_placement(const Quadric& q, const Vec3& a, const Vec3& b) {
  const Eigen::Matrix3d A = q.topLeftCorner<3, 3>();
  const Eigen::Vector3d rhs = -q.topRightCorner<3, 1>();

  const Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
  if (lu.isInvertible()) {
    const Vec3 p = lu.solve(rhs);
    if (p.allFinite()) {
      return {p, quadric_error(q, p)};
    }
  }

  // Singular system: fall back to the endpoints and the midpoint.
  Placement best{a, quadric_error(q, a)};
  const Placement at_b{b, quadric_error(q, b)};
  if (at_b.cost < best.cost) best = at_b;
  const Vec3 mid = 0.5 * (a + b);
  const Placement at_mid{mid, quadric_error(q, mid)};
  if (at_mid.cost < best.cost) best = at_mid;
  return best;
}

}  // namespace mesh_lod232::detail

#pragma once

#include <Eigen/Dense>

#include "mesh.h"

namespace mesh_lod232::detail {

using Quadric = Eigen::Matrix4d;

Quadric face_quadric(const Vec3& a, const Vec3& b, const Vec3& c);

std::vector<Quadric> initial_quadrics(const Mesh& m);

struct Placement {
  Vec3 pos;
  double cost = 0.0;
};

Placement optimal_placement(const Quadric& q, const Vec3& a, const Vec3& b);

double quadric_error(const Quadric& q, const Vec3& p);

}  // namespace mesh_lod232::detail

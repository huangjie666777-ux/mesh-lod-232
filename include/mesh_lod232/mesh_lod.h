#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace mesh_lod232 {

using Vec3 = Eigen::Vector3d;
using Tri = std::array<std::size_t, 3>;

enum class StopReason {
  kTargetNotSmaller,
  kTargetReached,
  kNoLegalCandidate,
};

struct SimplifyResult {
  std::vector<Vec3> vertices;
  std::vector<Tri> faces;
  std::size_t face_count = 0;
  std::size_t collapse_count = 0;
  StopReason stop_reason = StopReason::kTargetReached;
  std::string stop_message;
};

class ValidationError : public std::invalid_argument {
 public:
  explicit ValidationError(const std::string& what) : std::invalid_argument(what) {}
};

SimplifyResult simplify(const std::vector<Vec3>& vertices,
                        const std::vector<Tri>& faces,
                        std::size_t target_faces);

}  // namespace mesh_lod232

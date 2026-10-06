#pragma once

#include "mesh_lod232/mesh_lod.h"

namespace mesh_lod232::detail {

void validate_input(const std::vector<Vec3>& vertices,
                    const std::vector<Tri>& faces,
                    std::size_t target_faces);

}  // namespace mesh_lod232::detail

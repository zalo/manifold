// Copyright 2024 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "public.h"

namespace manifold {

/**
 * @ingroup Connections
 *
 * Result of a tetrahedralization operation, containing vertices and tetrahedra.
 */
struct TetMesh {
  /// The X-Y-Z positions of all vertices (may include Steiner points for
  /// constrained tetrahedralization)
  std::vector<glm::vec3> vertPos;
  /// The vertex indices of the four tetrahedron corners, for each tetrahedron.
  /// Ordering follows positive orientation (first three vertices form a face
  /// with outward normal pointing away from fourth vertex)
  std::vector<std::array<uint32_t, 4>> tetVerts;

  /// Number of tetrahedra
  size_t NumTet() const { return tetVerts.size(); }
  /// Number of vertices
  size_t NumVert() const { return vertPos.size(); }
};

/**
 * @ingroup Core
 *
 * Performs unconstrained Delaunay tetrahedralization on a set of 3D points.
 *
 * @param points Vector of 3D points to tetrahedralize
 * @param minQuality Minimum tetrahedron quality threshold (0.0 to 1.0).
 *                   Tetrahedra with quality below this are filtered out.
 *                   Default is 0.0 (no filtering).
 * @return TetMesh containing vertices and tetrahedra indices
 */
TetMesh DelaunayTetrahedralization(const std::vector<glm::vec3>& points,
                                   float minQuality = 0.0f);

class Manifold;

/**
 * @ingroup Core
 *
 * Performs constrained Delaunay tetrahedralization on a manifold mesh.
 * Uses an intersection-based approach to trim tetrahedra to the manifold
 * boundary:
 * 1. Performs unconstrained Delaunay tetrahedralization of manifold vertices
 * 2. Intersects each tetrahedron with the original manifold
 * 3. Discards tetrahedra with empty intersection (fully outside)
 * 4. Keeps tetrahedra fully inside the manifold
 * 5. Accumulates partial intersections, unions them, and re-tetrahedralizes
 * 6. Repeats until convergence
 *
 * @param manifold The input manifold mesh to tetrahedralize
 * @param minQuality Minimum tetrahedron quality threshold (0.0 to 1.0)
 * @param maxIterations Maximum number of refinement iterations. Default is 10.
 * @return TetMesh containing vertices and tetrahedra indices for the interior
 *         tetrahedralization of the manifold.
 */
TetMesh ConstrainedDelaunayTetrahedralization(const Manifold& manifold,
                                              float minQuality = 0.0f,
                                              int maxIterations = 10);

}  // namespace manifold

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

#include "tetrahedralization.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <set>
#include <unordered_set>
#include <vector>

#include "manifold.h"

namespace manifold {

namespace {

// Face indices for each tetrahedron face (vertices forming face opposite to
// vertex i)
const int kTetFaces[4][3] = {{2, 1, 0}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};

// Generate small random perturbation to avoid degenerate configurations
inline float RandEps() {
  constexpr float eps = 1e-6f;
  return -eps + 2.0f * (static_cast<float>(rand()) / RAND_MAX) * eps;
}

// Compute circumcenter of tetrahedron defined by four points
glm::vec3 GetCircumCenter(const glm::vec3& p0, const glm::vec3& p1,
                          const glm::vec3& p2, const glm::vec3& p3) {
  glm::vec3 b = p1 - p0;
  glm::vec3 c = p2 - p0;
  glm::vec3 d = p3 - p0;

  float det = 2.0f * (b.x * (c.y * d.z - c.z * d.y) -
                      b.y * (c.x * d.z - c.z * d.x) +
                      b.z * (c.x * d.y - c.y * d.x));
  if (det == 0.0f) return p0;

  glm::vec3 v = glm::cross(c, d) * glm::dot(b, b) +
                glm::cross(d, b) * glm::dot(c, c) +
                glm::cross(b, c) * glm::dot(d, d);
  return p0 + v / det;
}

// Compute tetrahedron quality metric (ratio of volume to mean squared edge
// length)
float TetQuality(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2,
                 const glm::vec3& p3) {
  glm::vec3 d0 = p1 - p0, d1 = p2 - p0, d2 = p3 - p0;
  glm::vec3 d3 = p2 - p1, d4 = p3 - p2, d5 = p1 - p3;

  float s0 = glm::length(d0), s1 = glm::length(d1), s2 = glm::length(d2);
  float s3 = glm::length(d3), s4 = glm::length(d4), s5 = glm::length(d5);

  float ms = (s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4 + s5 * s5) / 6.0f;
  float rms = std::sqrt(ms);

  if (rms == 0.0f) return 0.0f;

  float s = 12.0f / std::sqrt(2.0f);
  float vol = glm::dot(d0, glm::cross(d1, d2)) / 6.0f;
  return s * vol / (rms * rms * rms);
}

// Edge structure for tracking tetrahedron connectivity
struct Edge {
  int id0, id1, tetNr, faceNr;
  bool operator<(const Edge& o) const {
    return id0 < o.id0 || (id0 == o.id0 && id1 < o.id1);
  }
  bool operator==(const Edge& o) const {
    return id0 == o.id0 && id1 == o.id1;
  }
};

// Core incremental Delaunay tetrahedralization algorithm
std::vector<uint32_t> CreateTetIds(std::vector<glm::vec3>& verts,
                                   float minQuality) {
  std::vector<int> tetIds;
  std::vector<int> neighbors;
  std::vector<int> tetMarks;
  int tetMark = 0;
  int firstFreeTet = -1;

  std::vector<glm::vec3> planesN;
  std::vector<float> planesD;

  int firstBig = static_cast<int>(verts.size()) - 4;

  // Initialize with big bounding tetrahedron
  tetIds.push_back(firstBig);
  tetIds.push_back(firstBig + 1);
  tetIds.push_back(firstBig + 2);
  tetIds.push_back(firstBig + 3);
  tetMarks.push_back(0);

  for (int i = 0; i < 4; i++) {
    neighbors.push_back(-1);
    glm::vec3 fp0 = verts[firstBig + kTetFaces[i][0]];
    glm::vec3 fp1 = verts[firstBig + kTetFaces[i][1]];
    glm::vec3 fp2 = verts[firstBig + kTetFaces[i][2]];
    glm::vec3 n = glm::cross(fp1 - fp0, fp2 - fp0);
    n = glm::normalize(n);
    planesN.push_back(n);
    planesD.push_back(glm::dot(fp0, n));
  }

  // Insert each point incrementally
  for (int i = 0; i < firstBig; i++) {
    glm::vec3 p = verts[i];

    // Find non-deleted tet
    int tetNr = 0;
    while (tetIds[4 * tetNr] < 0) tetNr++;

    // Find containing tetrahedron using plane half-space tests
    tetMark++;
    bool found = false;

    while (!found) {
      if (tetNr < 0 || tetMarks[tetNr] == tetMark) break;
      tetMarks[tetNr] = tetMark;

      int id0 = tetIds[4 * tetNr];
      int id1 = tetIds[4 * tetNr + 1];
      int id2 = tetIds[4 * tetNr + 2];
      int id3 = tetIds[4 * tetNr + 3];

      glm::vec3 center =
          (verts[id0] + verts[id1] + verts[id2] + verts[id3]) * 0.25f;

      float minT = std::numeric_limits<float>::infinity();
      int minFaceNr = -1;

      for (int j = 0; j < 4; j++) {
        glm::vec3 n = planesN[4 * tetNr + j];
        float d = planesD[4 * tetNr + j];

        float hp = glm::dot(n, p) - d;
        float hc = glm::dot(n, center) - d;
        float t = hp - hc;
        if (t == 0) continue;

        t = -hc / t;
        if (t >= 0.0f && t < minT) {
          minT = t;
          minFaceNr = j;
        }
      }

      if (minT >= 1.0f)
        found = true;
      else
        tetNr = neighbors[4 * tetNr + minFaceNr];
    }

    if (!found) continue;

    // Find all tetrahedra violating Delaunay condition (circumsphere test)
    tetMark++;
    std::vector<int> violatingTets;
    std::vector<int> stack = {tetNr};

    while (!stack.empty()) {
      tetNr = stack.back();
      stack.pop_back();
      if (tetMarks[tetNr] == tetMark) continue;
      tetMarks[tetNr] = tetMark;
      violatingTets.push_back(tetNr);

      for (int j = 0; j < 4; j++) {
        int n = neighbors[4 * tetNr + j];
        if (n < 0 || tetMarks[n] == tetMark) continue;

        int nid0 = tetIds[4 * n];
        int nid1 = tetIds[4 * n + 1];
        int nid2 = tetIds[4 * n + 2];
        int nid3 = tetIds[4 * n + 3];

        glm::vec3 c =
            GetCircumCenter(verts[nid0], verts[nid1], verts[nid2], verts[nid3]);
        float r = glm::length(verts[nid0] - c);
        if (glm::length(p - c) < r) stack.push_back(n);
      }
    }

    // Remove violating tetrahedra and create new ones connecting to inserted
    // point
    std::vector<Edge> edges;

    for (int j = 0; j < static_cast<int>(violatingTets.size()); j++) {
      tetNr = violatingTets[j];

      int ids[4], ns[4];
      for (int k = 0; k < 4; k++) {
        ids[k] = tetIds[4 * tetNr + k];
        ns[k] = neighbors[4 * tetNr + k];
      }

      // Mark tet as deleted
      tetIds[4 * tetNr] = -1;
      tetIds[4 * tetNr + 1] = firstFreeTet;
      firstFreeTet = tetNr;

      for (int k = 0; k < 4; k++) {
        int neighborTet = ns[k];
        if (neighborTet >= 0 && tetMarks[neighborTet] == tetMark) continue;

        // Create new tetrahedron
        int newTetNr = firstFreeTet;
        if (newTetNr >= 0) {
          firstFreeTet = tetIds[4 * firstFreeTet + 1];
        } else {
          newTetNr = static_cast<int>(tetIds.size()) / 4;
          tetMarks.push_back(0);
          for (int l = 0; l < 4; l++) {
            tetIds.push_back(-1);
            neighbors.push_back(-1);
            planesN.push_back(glm::vec3());
            planesD.push_back(0.0f);
          }
        }

        int fid0 = ids[kTetFaces[k][2]];
        int fid1 = ids[kTetFaces[k][1]];
        int fid2 = ids[kTetFaces[k][0]];

        tetIds[4 * newTetNr] = fid0;
        tetIds[4 * newTetNr + 1] = fid1;
        tetIds[4 * newTetNr + 2] = fid2;
        tetIds[4 * newTetNr + 3] = i;

        neighbors[4 * newTetNr] = neighborTet;
        if (neighborTet >= 0) {
          for (int l = 0; l < 4; l++) {
            if (neighbors[4 * neighborTet + l] == tetNr)
              neighbors[4 * neighborTet + l] = newTetNr;
          }
        }

        neighbors[4 * newTetNr + 1] = -1;
        neighbors[4 * newTetNr + 2] = -1;
        neighbors[4 * newTetNr + 3] = -1;

        // Update plane equations for new tetrahedron
        for (int l = 0; l < 4; l++) {
          glm::vec3 fp0 = verts[tetIds[4 * newTetNr + kTetFaces[l][0]]];
          glm::vec3 fp1 = verts[tetIds[4 * newTetNr + kTetFaces[l][1]]];
          glm::vec3 fp2 = verts[tetIds[4 * newTetNr + kTetFaces[l][2]]];
          glm::vec3 newN = glm::cross(fp1 - fp0, fp2 - fp0);
          newN = glm::normalize(newN);
          planesN[4 * newTetNr + l] = newN;
          planesD[4 * newTetNr + l] = glm::dot(newN, fp0);
        }

        // Track edges for neighbor fixing
        Edge e;
        e.tetNr = newTetNr;
        if (fid0 < fid1) {
          e.id0 = fid0;
          e.id1 = fid1;
          e.faceNr = 1;
        } else {
          e.id0 = fid1;
          e.id1 = fid0;
          e.faceNr = 1;
        }
        edges.push_back(e);

        if (fid1 < fid2) {
          e.id0 = fid1;
          e.id1 = fid2;
          e.faceNr = 2;
        } else {
          e.id0 = fid2;
          e.id1 = fid1;
          e.faceNr = 2;
        }
        edges.push_back(e);

        if (fid2 < fid0) {
          e.id0 = fid2;
          e.id1 = fid0;
          e.faceNr = 3;
        } else {
          e.id0 = fid0;
          e.id1 = fid2;
          e.faceNr = 3;
        }
        edges.push_back(e);
      }
    }

    // Fix neighbor relationships by matching edges
    std::sort(edges.begin(), edges.end());

    size_t nr = 0;
    while (nr < edges.size()) {
      Edge e0 = edges[nr++];
      if (nr < edges.size() && edges[nr] == e0) {
        Edge e1 = edges[nr++];
        neighbors[4 * e0.tetNr + e0.faceNr] = e1.tetNr;
        neighbors[4 * e1.tetNr + e1.faceNr] = e0.tetNr;
      }
    }
  }

  // Filter results: remove deleted, boundary, and low-quality tetrahedra
  int numTets = static_cast<int>(tetIds.size()) / 4;
  std::vector<uint32_t> result;

  for (int i = 0; i < numTets; i++) {
    int id0 = tetIds[4 * i];
    int id1 = tetIds[4 * i + 1];
    int id2 = tetIds[4 * i + 2];
    int id3 = tetIds[4 * i + 3];

    // Skip deleted tets and those containing boundary vertices
    if (id0 < 0 || id0 >= firstBig || id1 >= firstBig || id2 >= firstBig ||
        id3 >= firstBig)
      continue;

    glm::vec3 p0 = verts[id0], p1 = verts[id1], p2 = verts[id2], p3 = verts[id3];

    float quality = TetQuality(p0, p1, p2, p3);
    if (quality < minQuality) continue;

    result.push_back(static_cast<uint32_t>(id0));
    result.push_back(static_cast<uint32_t>(id1));
    result.push_back(static_cast<uint32_t>(id2));
    result.push_back(static_cast<uint32_t>(id3));
  }

  return result;
}

// Hash function for triangle face (unordered vertex triple)
struct TriangleHash {
  size_t operator()(const std::array<uint32_t, 3>& tri) const {
    // Sort vertices to make hash order-independent
    uint32_t a = tri[0], b = tri[1], c = tri[2];
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return std::hash<uint64_t>{}(
        (static_cast<uint64_t>(a) << 42) | (static_cast<uint64_t>(b) << 21) | c);
  }
};

struct TriangleEqual {
  bool operator()(const std::array<uint32_t, 3>& a,
                  const std::array<uint32_t, 3>& b) const {
    std::array<uint32_t, 3> sa = a, sb = b;
    std::sort(sa.begin(), sa.end());
    std::sort(sb.begin(), sb.end());
    return sa == sb;
  }
};

}  // namespace

TetMesh DelaunayTetrahedralization(const std::vector<glm::vec3>& points,
                                   float minQuality) {
  TetMesh result;

  if (points.size() < 4) {
    return result;  // Not enough points for tetrahedralization
  }

  // Copy points and add small random perturbations to avoid degeneracies
  std::vector<glm::vec3> tetVerts;
  tetVerts.reserve(points.size() + 4);
  for (const auto& p : points) {
    tetVerts.push_back(
        glm::vec3(p.x + RandEps(), p.y + RandEps(), p.z + RandEps()));
  }

  // Compute bounding sphere
  glm::vec3 center(0, 0, 0);
  for (const auto& p : tetVerts) {
    center += p;
  }
  center /= static_cast<float>(tetVerts.size());

  float radius = 0.0f;
  for (const auto& p : tetVerts) {
    float d = glm::length(p - center);
    radius = std::max(radius, d);
  }

  // Add four corner vertices for bounding tetrahedron
  float s = 5.0f * radius;
  tetVerts.push_back(center + glm::vec3(-s, 0.0f, -s));
  tetVerts.push_back(center + glm::vec3(s, 0.0f, -s));
  tetVerts.push_back(center + glm::vec3(0.0f, s, s));
  tetVerts.push_back(center + glm::vec3(0.0f, -s, s));

  // Run core algorithm
  std::vector<uint32_t> tetIndices = CreateTetIds(tetVerts, minQuality);

  // Package results
  result.vertPos = points;  // Use original points (without perturbation)
  int numTets = static_cast<int>(tetIndices.size()) / 4;
  result.tetVerts.reserve(numTets);
  for (int i = 0; i < numTets; i++) {
    result.tetVerts.push_back({tetIndices[4 * i], tetIndices[4 * i + 1],
                               tetIndices[4 * i + 2], tetIndices[4 * i + 3]});
  }

  return result;
}

// Create a tetrahedron manifold from 4 vertices
Manifold CreateTetrahedronManifold(const glm::vec3& p0, const glm::vec3& p1,
                                   const glm::vec3& p2, const glm::vec3& p3) {
  // Build mesh for a single tetrahedron with 4 triangular faces
  Mesh tetMesh;
  tetMesh.vertPos = {p0, p1, p2, p3};

  // Compute orientation - ensure outward-facing normals
  glm::vec3 center = (p0 + p1 + p2 + p3) * 0.25f;

  // Helper to check if triangle faces outward from center
  auto facesOutward = [&center](const glm::vec3& a, const glm::vec3& b,
                                const glm::vec3& c) {
    glm::vec3 faceCenter = (a + b + c) / 3.0f;
    glm::vec3 normal = glm::cross(b - a, c - a);
    glm::vec3 toFace = faceCenter - center;
    return glm::dot(normal, toFace) > 0;
  };

  // Face 0: vertices 0, 2, 1 (opposite to vertex 3)
  if (facesOutward(p0, p2, p1)) {
    tetMesh.triVerts.push_back({0, 2, 1});
  } else {
    tetMesh.triVerts.push_back({0, 1, 2});
  }

  // Face 1: vertices 0, 1, 3 (opposite to vertex 2)
  if (facesOutward(p0, p1, p3)) {
    tetMesh.triVerts.push_back({0, 1, 3});
  } else {
    tetMesh.triVerts.push_back({0, 3, 1});
  }

  // Face 2: vertices 1, 2, 3 (opposite to vertex 0)
  if (facesOutward(p1, p2, p3)) {
    tetMesh.triVerts.push_back({1, 2, 3});
  } else {
    tetMesh.triVerts.push_back({1, 3, 2});
  }

  // Face 3: vertices 0, 3, 2 (opposite to vertex 1)
  if (facesOutward(p0, p3, p2)) {
    tetMesh.triVerts.push_back({0, 3, 2});
  } else {
    tetMesh.triVerts.push_back({0, 2, 3});
  }

  return Manifold(tetMesh);
}

// Compute signed volume of tetrahedron
float TetrahedronVolume(const glm::vec3& p0, const glm::vec3& p1,
                        const glm::vec3& p2, const glm::vec3& p3) {
  glm::vec3 d0 = p1 - p0, d1 = p2 - p0, d2 = p3 - p0;
  return std::abs(glm::dot(d0, glm::cross(d1, d2)) / 6.0f);
}

TetMesh ConstrainedDelaunayTetrahedralization(const Manifold& manifold,
                                              float minQuality,
                                              int /*maxIterations*/) {
  TetMesh result;

  if (manifold.IsEmpty()) {
    return result;
  }

  // Get mesh data from manifold
  Mesh mesh = manifold.GetMesh();
  if (mesh.vertPos.size() < 4 || mesh.triVerts.empty()) {
    return result;
  }

  const float volumeTolerance = 1e-6f;

  // Reset random seed for reproducibility
  srand(12345);

  // Perform unconstrained Delaunay tetrahedralization of surface vertices
  TetMesh tetMesh = DelaunayTetrahedralization(mesh.vertPos, minQuality);

  if (tetMesh.tetVerts.empty()) {
    return result;
  }

  // Filter tetrahedra to keep only those inside the manifold
  std::vector<std::array<uint32_t, 4>> insideTets;

  for (const auto& tet : tetMesh.tetVerts) {
    const glm::vec3& p0 = tetMesh.vertPos[tet[0]];
    const glm::vec3& p1 = tetMesh.vertPos[tet[1]];
    const glm::vec3& p2 = tetMesh.vertPos[tet[2]];
    const glm::vec3& p3 = tetMesh.vertPos[tet[3]];

    float tetVolume = TetrahedronVolume(p0, p1, p2, p3);
    if (tetVolume < volumeTolerance) {
      continue;
    }

    // Create manifold for this tetrahedron
    Manifold tetManifold = CreateTetrahedronManifold(p0, p1, p2, p3);
    if (tetManifold.IsEmpty() || tetManifold.Status() != Manifold::Error::NoError) {
      continue;
    }

    // Intersect with original manifold to check if inside
    Manifold intersection = tetManifold ^ manifold;
    if (intersection.IsEmpty()) {
      continue;  // Tetrahedron is outside
    }

    float intersectionVolume = intersection.GetProperties().volume;

    // Keep tetrahedra that are mostly inside (>50% overlap)
    if (intersectionVolume > 0.5f * tetVolume) {
      insideTets.push_back(tet);
    }
  }

  result.vertPos = tetMesh.vertPos;
  result.tetVerts = insideTets;
  return result;
}

}  // namespace manifold

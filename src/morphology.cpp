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

// Surface-only morphological opening and closing, after
//   Sellán, Kesten, Sheng & Jacobson, "Opening and Closing Surfaces",
//   ACM TOG (SIGGRAPH Asia) 2020.
//
// Closing by a ball of radius r is realized as a curvature-bounded
// minimum-curvature flow: vertices where a ball of radius r already fits
// (signed min principal curvature >= -1/r) are frozen and left bit-for-bit
// identical; the remaining concave vertices are moved by a semi-implicit
// minimum-curvature flow step until they too become ball-reachable. Opening
// is the exact dual (freeze where max principal curvature <= 1/r, flow the
// convex regions inward). No new third-party dependencies are introduced: the
// sparse semi-implicit solve uses a matrix-free Conjugate Gradient
// (morphology.h) and the per-face principal-curvature directions come from a
// closed-form 2x2 second-fundamental-form fit (Rusinkiewicz 2004).

#include <algorithm>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "impl.h"
#include "morphology.h"
#include "parallel.h"
#include "utils.h"

namespace {
using namespace manifold;

// Symmetric per-triangle 3x3 stiffness, stored as the 6 unique entries:
// {s00, s01, s02, s11, s12, s22}.
struct FaceStiff {
  double s[6];
};

// Per-vertex curvature accumulators, identical to
// properties.cpp::CurvatureAngles (and sharing its valence-correcting
// normalization). After normalization by `degree/(6*area)` these give
// Manifold's mean curvature (the SUM of the two principal curvatures) and
// Gaussian curvature (their PRODUCT); the principal curvatures then follow in
// closed form. Vertex area `area` (1/3 sum of incident triangle areas) doubles
// as the lumped FEM mass.
struct RawCurvature {
  VecView<double> H;       // integrated mean curvature
  VecView<double> K;       // angle defect (init to 2pi)
  VecView<double> area;    // vertex area
  VecView<double> degree;  // incident triangle-corner count (valence)
  const Halfedges& halfedge;
  VecView<const vec3> vertPos;
  VecView<const vec3> triNormal;

  void operator()(size_t tri) {
    vec3 edge[3];
    vec3 edgeLength(0.0);
    for (int i : {0, 1, 2}) {
      const int edgeIdx = 3 * tri + i;
      const int startVert = halfedge.Start(edgeIdx);
      const int endVert = halfedge.End(edgeIdx);
      edge[i] = vertPos[endVert] - vertPos[startVert];
      edgeLength[i] = la::length(edge[i]);
      edge[i] /= edgeLength[i];
      const int neighborTri = halfedge.Pair(edgeIdx) / 3;
      double s = la::dot(la::cross(triNormal[tri], triNormal[neighborTri]),
                         edge[i]);
      s = std::clamp(s, -1.0, 1.0);
      const double dihedral = 0.25 * edgeLength[i] * std::asin(s);
      AtomicAdd(H[startVert], dihedral);
      AtomicAdd(H[endVert], dihedral);
      AtomicAdd(degree[startVert], 1.0);
    }

    vec3 phi;
    phi[0] = std::acos(std::clamp(-la::dot(edge[2], edge[0]), -1.0, 1.0));
    phi[1] = std::acos(std::clamp(-la::dot(edge[0], edge[1]), -1.0, 1.0));
    phi[2] = kPi - phi[0] - phi[1];
    const double area3 =
        edgeLength[0] * edgeLength[1] * la::length(la::cross(edge[0], edge[1])) /
        6;

    for (int i : {0, 1, 2}) {
      const int vert = halfedge.Start(3 * tri + i);
      AtomicAdd(K[vert], -phi[i]);
      AtomicAdd(area[vert], area3);
    }
  }
};

// Triangle area from its three vertices.
inline double TriArea(vec3 p0, vec3 p1, vec3 p2) {
  return 0.5 * la::length(la::cross(p1 - p0, p2 - p0));
}

// Gradients of the three linear hat functions of a triangle (each a 3D vector
// lying in the triangle plane). grad_i = (n x e_i) / (2A) where e_i is the
// edge opposite vertex i. Their sum is zero, so the resulting stiffness has
// the constant field in its null space (as required for a Laplacian).
inline void HatGradients(vec3 p0, vec3 p1, vec3 p2, vec3 n, double area,
                         vec3& g0, vec3& g1, vec3& g2) {
  const double inv2A = 1.0 / (2.0 * area);
  g0 = la::cross(n, p2 - p1) * inv2A;
  g1 = la::cross(n, p0 - p2) * inv2A;
  g2 = la::cross(n, p1 - p0) * inv2A;
}

// Build the per-triangle stiffness S = A * G^T P G, where G maps the three
// vertex scalars to the in-plane gradient and P is the projection. For the
// isotropic (mean-curvature) flow P = I, giving the classic cotangent matrix.
// For the anisotropic (minimum-curvature) flow P = d d^T projects onto the
// per-face principal-curvature direction d, penalizing only stretching along
// that direction (paper Eq. 19).
inline FaceStiff BuildStiffness(vec3 g0, vec3 g1, vec3 g2, double area,
                                bool anisotropic, vec3 dir) {
  FaceStiff f;
  if (anisotropic) {
    const double a0 = la::dot(g0, dir);
    const double a1 = la::dot(g1, dir);
    const double a2 = la::dot(g2, dir);
    f.s[0] = area * a0 * a0;
    f.s[1] = area * a0 * a1;
    f.s[2] = area * a0 * a2;
    f.s[3] = area * a1 * a1;
    f.s[4] = area * a1 * a2;
    f.s[5] = area * a2 * a2;
  } else {
    f.s[0] = area * la::dot(g0, g0);
    f.s[1] = area * la::dot(g0, g1);
    f.s[2] = area * la::dot(g0, g2);
    f.s[3] = area * la::dot(g1, g1);
    f.s[4] = area * la::dot(g1, g2);
    f.s[5] = area * la::dot(g2, g2);
  }
  return f;
}

// Per-face principal-curvature direction, following the paper: fit a quadric
// height field to the triangle's three corners plus the three "flap" vertices
// (the opposite corner of each edge-neighbor triangle), in the face tangent
// frame (Cazals & Pouget 2005). The Hessian of that quadric is the second
// fundamental form; its eigenvector for the smaller (close) / larger (open)
// principal curvature is the returned direction. Fitting positions over this
// wider 6-point stencil is far smoother than a per-face normal difference,
// which is what keeps the min-curvature flow from rippling. The sign of the
// result is irrelevant: it is only used as the rank-1 projector d d^T.
//
// `flap[i] < 0` marks a missing neighbor (left out of the fit).
vec3 FacePrincipalDir(const vec3* pts, int nPts, vec3 centroid, vec3 faceNormal,
                      bool wantMin) {
  vec3 u = pts[1] - pts[0];
  u -= la::dot(u, faceNormal) * faceNormal;
  const double uLen = la::length(u);
  if (uLen < 1e-12) return vec3(1, 0, 0);
  u /= uLen;
  const vec3 v = la::cross(faceNormal, u);

  // Least-squares fit of z = a x^2 + b xy + c y^2 over the local (x,y,z) of
  // every point, where z is the height above the face's tangent plane.
  double N00 = 0, N01 = 0, N02 = 0, N11 = 0, N12 = 0, N22 = 0;
  double r0 = 0, r1 = 0, r2 = 0;
  for (int i = 0; i < nPts; ++i) {
    const vec3 d = pts[i] - centroid;
    const double x = la::dot(d, u), y = la::dot(d, v), z = la::dot(d, faceNormal);
    const double bx = x * x, by = x * y, bz = y * y;  // basis row
    N00 += bx * bx;
    N01 += bx * by;
    N02 += bx * bz;
    N11 += by * by;
    N12 += by * bz;
    N22 += bz * bz;
    r0 += bx * z;
    r1 += by * z;
    r2 += bz * z;
  }
  const mat3 M = mat3({N00, N01, N02}, {N01, N11, N12}, {N02, N12, N22});
  if (std::abs(la::determinant(M)) < 1e-24) return u;  // flat / degenerate
  const vec3 abc = la::mul(la::inverse(M), vec3(r0, r1, r2));
  // Second fundamental form is the Hessian [[2a, b], [b, 2c]].
  const double A = 2 * abc[0], B = abc[1], C = 2 * abc[2];

  const double tr = 0.5 * (A + C);
  const double disc = std::sqrt(std::max(0.25 * (A - C) * (A - C) + B * B, 0.0));
  // The quadric's height is measured along the OUTWARD normal, so its curvature
  // sign is opposite to our convention (convex positive). Hence closing's
  // minimum (most concave) principal direction is the quadric's LARGER
  // eigenvalue eigenvector, and opening's is the smaller.
  const double k = wantMin ? tr + disc : tr - disc;
  double cu = B, cv = k - A;
  if (std::abs(cu) + std::abs(cv) < 1e-12) {
    cu = k - C;
    cv = B;
  }
  const double len = std::sqrt(cu * cu + cv * cv);
  if (len < 1e-12) return u;
  return la::normalize((cu / len) * u + (cv / len) * v);
}

// Uniform spatial-hash grid over a fixed set of anchor points (the original
// concavity/convexity). Used to bound the region a morphological operation may
// ever change to a fixed distance from those anchors, so the moving front
// cannot creep into far-away flat regions. The cell size equals the query
// radius, so a 3x3x3 neighborhood search is exhaustive.
struct AnchorGrid {
  double cell = 1;
  std::unordered_map<int64_t, std::vector<vec3>> cells;
  static int64_t Key(int x, int y, int z) {
    return (static_cast<int64_t>(x + (1 << 20)) << 42) |
           (static_cast<int64_t>(y + (1 << 20)) << 21) |
           static_cast<int64_t>(z + (1 << 20));
  }
  void Insert(vec3 p) {
    cells[Key(static_cast<int>(std::floor(p[0] / cell)),
              static_cast<int>(std::floor(p[1] / cell)),
              static_cast<int>(std::floor(p[2] / cell)))]
        .push_back(p);
  }
  bool Near(vec3 p, double d) const {
    const double d2 = d * d;
    const int ix = static_cast<int>(std::floor(p[0] / cell));
    const int iy = static_cast<int>(std::floor(p[1] / cell));
    const int iz = static_cast<int>(std::floor(p[2] / cell));
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          const auto it = cells.find(Key(ix + dx, iy + dy, iz + dz));
          if (it == cells.end()) continue;
          for (const vec3& q : it->second)
            if (la::dot(p - q, p - q) <= d2) return true;
        }
    return false;
  }
};

// ---------------------------------------------------------------------------
// Local, frozen-aware isotropic remesher.
//
// Operates on a plain triangle soup extracted from the Impl, restricting every
// topological edit to edges whose BOTH endpoints are active (moving). Frozen
// vertices are never moved and frozen-incident faces are never split, collapsed
// or flipped, so unchanged regions of the surface keep their exact geometry.
// The result is validated for manifoldness before being committed back.
// ---------------------------------------------------------------------------
struct Soup {
  std::vector<vec3> pos;
  std::vector<char> frozen;
  std::vector<std::array<int, 3>> tri;  // {-1,-1,-1} marks a deleted triangle
};

inline int64_t EdgeKey(int u, int v) {
  if (u > v) std::swap(u, v);
  return (static_cast<int64_t>(u) << 32) | static_cast<uint32_t>(v);
}
inline int OppVert(const std::array<int, 3>& t, int a, int b) {
  for (int i = 0; i < 3; ++i)
    if (t[i] != a && t[i] != b) return t[i];
  return -1;
}
inline bool HasDirected(const std::array<int, 3>& t, int a, int b) {
  for (int i = 0; i < 3; ++i)
    if (t[i] == a && t[(i + 1) % 3] == b) return true;
  return false;
}
inline vec3 RawNormal(const Soup& m, const std::array<int, 3>& t) {
  return la::cross(m.pos[t[1]] - m.pos[t[0]], m.pos[t[2]] - m.pos[t[0]]);
}
inline int sq(int x) { return x * x; }

struct Adj {
  std::unordered_map<int64_t, std::vector<int>> edge;  // undirected edge -> tris
  std::vector<std::vector<int>> vtri;                  // vertex -> incident tris
};
Adj BuildAdj(const Soup& m) {
  Adj a;
  a.vtri.assign(m.pos.size(), {});
  for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
    if (m.tri[t][0] < 0) continue;
    for (int i = 0; i < 3; ++i) {
      a.edge[EdgeKey(m.tri[t][i], m.tri[t][(i + 1) % 3])].push_back(t);
      a.vtri[m.tri[t][i]].push_back(t);
    }
  }
  return a;
}

// Split triangle t's edge {a,b} by inserting midpoint vertex mIdx, preserving
// winding: triangle t is replaced and one new triangle is appended.
void SplitTri(Soup& m, int t, int a, int b, int mIdx) {
  const std::array<int, 3> T = m.tri[t];
  int i = 0;
  for (; i < 3; ++i) {
    const int x = T[i], y = T[(i + 1) % 3];
    if ((x == a && y == b) || (x == b && y == a)) break;
  }
  const int vi = T[i], vj = T[(i + 1) % 3], vk = T[(i + 2) % 3];
  m.tri[t] = {vi, mIdx, vk};
  m.tri.push_back({mIdx, vj, vk});
}

bool SplitPass(Soup& m, double high) {
  const Adj a = BuildAdj(m);
  std::vector<char> dirty(m.tri.size(), 0);
  const double h2 = high * high;
  bool any = false;
  for (const auto& kv : a.edge) {
    if (kv.second.size() != 2) continue;
    const int t0 = kv.second[0], t1 = kv.second[1];
    if (dirty[t0] || dirty[t1]) continue;
    const int u = static_cast<int>(kv.first >> 32);
    const int v = static_cast<int>(static_cast<uint32_t>(kv.first));
    // Edit only edges strictly interior to the active region: both endpoints
    // and both opposite verts active, so no frozen-incident face is altered.
    const int sc = OppVert(m.tri[t0], u, v), sd = OppVert(m.tri[t1], u, v);
    if (m.frozen[u] || m.frozen[v] || sc < 0 || sd < 0 || m.frozen[sc] ||
        m.frozen[sd])
      continue;
    if (la::dot(m.pos[u] - m.pos[v], m.pos[u] - m.pos[v]) <= h2) continue;
    const int mIdx = m.pos.size();
    m.pos.push_back(0.5 * (m.pos[u] + m.pos[v]));
    m.frozen.push_back(0);
    SplitTri(m, t0, u, v, mIdx);
    SplitTri(m, t1, u, v, mIdx);
    dirty[t0] = dirty[t1] = 1;
    any = true;
  }
  return any;
}

bool CollapsePass(Soup& m, double low) {
  const Adj a = BuildAdj(m);
  std::vector<char> dirtyV(m.pos.size(), 0);
  const double l2 = low * low;
  bool any = false;
  for (const auto& kv : a.edge) {
    if (kv.second.size() != 2) continue;
    const int u = static_cast<int>(kv.first >> 32);
    const int v = static_cast<int>(static_cast<uint32_t>(kv.first));
    if (m.frozen[u] || m.frozen[v] || dirtyV[u] || dirtyV[v]) continue;
    if (la::dot(m.pos[u] - m.pos[v], m.pos[u] - m.pos[v]) > l2) continue;
    const int t0 = kv.second[0], t1 = kv.second[1];
    const int c = OppVert(m.tri[t0], u, v), d = OppVert(m.tri[t1], u, v);
    if (c < 0 || d < 0) continue;

    // Link condition: the only shared neighbors of u and v are c and d,
    // otherwise the collapse would create non-manifold topology.
    std::unordered_set<int> nu;
    for (const int t : a.vtri[u])
      for (const int w : m.tri[t])
        if (w != u) nu.insert(w);
    std::vector<int> common;
    std::unordered_set<int> nv;
    for (const int t : a.vtri[v])
      for (const int w : m.tri[t])
        if (w != v && nv.insert(w).second && nu.count(w)) common.push_back(w);
    if (common.size() != 2) continue;
    if (!((common[0] == c && common[1] == d) ||
          (common[0] == d && common[1] == c)))
      continue;

    // Only collapse when the whole one-ring of both endpoints is active, so no
    // frozen-incident face is reshaped.
    bool frozenNbr = false;
    for (const int w : nu)
      if (m.frozen[w]) frozenNbr = true;
    for (const int w : nv)
      if (m.frozen[w]) frozenNbr = true;
    if (frozenNbr) continue;

    // Collapse v onto u, keeping u at its existing position. Placing the
    // merged vertex at the chord midpoint instead would push it off a curved
    // surface (outward on a concave fill), leaving a visible spike that
    // tangential relaxation cannot remove.
    const vec3 mid = m.pos[u];
    bool ok = true;
    auto check = [&](int t) {
      if (!ok || t == t0 || t == t1 || m.tri[t][0] < 0) return;
      const vec3 oldN = RawNormal(m, m.tri[t]);
      vec3 p[3];
      for (int i = 0; i < 3; ++i) {
        const int w = m.tri[t][i];
        p[i] = (w == u || w == v) ? mid : m.pos[w];
      }
      const vec3 nn = la::cross(p[1] - p[0], p[2] - p[0]);
      if (la::dot(nn, oldN) <= 0 || la::length(nn) < 1e-12) ok = false;
    };
    for (const int t : a.vtri[u]) check(t);
    for (const int t : a.vtri[v]) check(t);
    if (!ok) continue;

    m.pos[u] = mid;
    m.tri[t0] = {-1, -1, -1};
    m.tri[t1] = {-1, -1, -1};
    for (const int t : a.vtri[v])
      if (m.tri[t][0] >= 0)
        for (int& w : m.tri[t])
          if (w == v) w = u;
    dirtyV[u] = dirtyV[v] = 1;
    for (const int w : nu) dirtyV[w] = 1;
    for (const int w : nv) dirtyV[w] = 1;
    any = true;
  }
  return any;
}

bool FlipPass(Soup& m) {
  const Adj a = BuildAdj(m);
  std::vector<int> deg(m.pos.size());
  for (size_t i = 0; i < m.pos.size(); ++i) deg[i] = a.vtri[i].size();
  std::vector<char> dirty(m.tri.size(), 0);
  bool any = false;
  for (const auto& kv : a.edge) {
    if (kv.second.size() != 2) continue;
    const int t0 = kv.second[0], t1 = kv.second[1];
    if (dirty[t0] || dirty[t1]) continue;
    const int u = static_cast<int>(kv.first >> 32);
    const int v = static_cast<int>(static_cast<uint32_t>(kv.first));
    if (m.frozen[u] || m.frozen[v]) continue;
    const int ta = HasDirected(m.tri[t0], u, v) ? t0 : t1;
    const int tb = ta == t0 ? t1 : t0;
    const int c = OppVert(m.tri[ta], u, v), d = OppVert(m.tri[tb], u, v);
    if (c < 0 || d < 0 || c == d) continue;
    if (m.frozen[c] || m.frozen[d]) continue;  // keep flips interior to active
    if (a.edge.count(EdgeKey(c, d))) continue;  // flip would duplicate an edge

    // Flip only if it brings the four vertices' valences closer to 6.
    const int before = sq(deg[u] - 6) + sq(deg[v] - 6) + sq(deg[c] - 6) +
                       sq(deg[d] - 6);
    const int after = sq(deg[u] - 1 - 6) + sq(deg[v] - 1 - 6) +
                      sq(deg[c] + 1 - 6) + sq(deg[d] + 1 - 6);
    if (after >= before) continue;

    const std::array<int, 3> na = {c, u, d}, nb = {c, d, v};
    const vec3 oldN = RawNormal(m, m.tri[ta]) + RawNormal(m, m.tri[tb]);
    const vec3 nna = la::cross(m.pos[na[1]] - m.pos[na[0]],
                               m.pos[na[2]] - m.pos[na[0]]);
    const vec3 nnb = la::cross(m.pos[nb[1]] - m.pos[nb[0]],
                               m.pos[nb[2]] - m.pos[nb[0]]);
    if (la::dot(nna, oldN) <= 0 || la::dot(nnb, oldN) <= 0) continue;
    if (la::length(nna) < 1e-12 || la::length(nnb) < 1e-12) continue;

    m.tri[ta] = na;
    m.tri[tb] = nb;
    dirty[t0] = dirty[t1] = 1;
    --deg[u];
    --deg[v];
    ++deg[c];
    ++deg[d];
    any = true;
  }
  return any;
}

void RelaxPass(Soup& m, int iters) {
  for (int it = 0; it < iters; ++it) {
    const Adj a = BuildAdj(m);
    std::vector<vec3> vn(m.pos.size(), vec3(0.0));
    for (const auto& t : m.tri)
      if (t[0] >= 0) {
        const vec3 n = RawNormal(m, t);
        for (int i = 0; i < 3; ++i) vn[t[i]] += n;
      }
    std::vector<vec3> np = m.pos;
    for (size_t i = 0; i < m.pos.size(); ++i) {
      if (m.frozen[i] || a.vtri[i].empty()) continue;
      std::unordered_set<int> nb;
      for (const int t : a.vtri[i])
        for (const int w : m.tri[t])
          if (w != static_cast<int>(i)) nb.insert(w);
      if (nb.empty()) continue;
      vec3 ctr(0.0);
      for (const int w : nb) ctr += m.pos[w];
      ctr /= static_cast<double>(nb.size());
      vec3 delta = ctr - m.pos[i];
      const double L = la::length(vn[i]);
      if (L > 0) {
        const vec3 n = vn[i] / L;
        delta -= la::dot(delta, n) * n;  // tangential only: keep on surface
      }
      np[i] = m.pos[i] + 0.5 * delta;
    }
    m.pos = np;
  }
}

// Extract the soup, remesh the active region, validate, and commit back to the
// Impl (rebuilding the halfedge structure). Returns false (leaving the Impl
// untouched) if the edited soup is not a clean manifold.
bool RemeshActive(Manifold::Impl& impl, double h, const Vec<char>& frozenIn) {
  Soup m;
  const size_t nv = impl.NumVert(), nt = impl.NumTri();
  m.pos.resize(nv);
  m.frozen.resize(nv);
  for (size_t i = 0; i < nv; ++i) {
    m.pos[i] = impl.vertPos_[i];
    m.frozen[i] = frozenIn[i];
  }
  m.tri.resize(nt);
  for (size_t t = 0; t < nt; ++t)
    m.tri[t] = {impl.halfedge_.Start(3 * t), impl.halfedge_.Start(3 * t + 1),
                impl.halfedge_.Start(3 * t + 2)};

  // One isotropic remeshing iteration per flow step (Botsch et al.), matching
  // the reference implementation's remesh_iterations = 1: split long edges,
  // collapse short ones, flip toward valence 6, then tangentially relax.
  SplitPass(m, 4.0 / 3.0 * h);
  CollapsePass(m, 0.8 * h);
  FlipPass(m);
  RelaxPass(m, 1);

  // Compact and validate.
  std::vector<std::array<int, 3>> live;
  live.reserve(m.tri.size());
  for (const auto& t : m.tri)
    if (t[0] >= 0) live.push_back(t);

  std::unordered_set<int64_t> directed;
  std::unordered_map<int64_t, int> undirected;
  bool valid = true;
  for (const auto& t : live) {
    for (int i = 0; i < 3 && valid; ++i) {
      const int u = t[i], v = t[(i + 1) % 3];
      if (u == v) valid = false;
      const int64_t dk = (static_cast<int64_t>(u) << 32) |
                         static_cast<uint32_t>(v);
      if (!directed.insert(dk).second) valid = false;  // duplicate halfedge
      undirected[EdgeKey(u, v)]++;
    }
    if (!valid) break;
  }
  if (valid)
    for (const auto& kv : undirected)
      if (kv.second != 2) {
        valid = false;
        break;
      }
  if (!valid) return false;

  // Remove unreferenced verts and commit.
  std::vector<int> remap(m.pos.size(), -1);
  int cnt = 0;
  for (const auto& t : live)
    for (const int v : t)
      if (remap[v] < 0) remap[v] = cnt++;
  Vec<vec3> newPos(cnt);
  for (size_t i = 0; i < m.pos.size(); ++i)
    if (remap[i] >= 0) newPos[remap[i]] = m.pos[i];
  Vec<ivec3> tris(live.size());
  for (size_t i = 0; i < live.size(); ++i)
    tris[i] = ivec3(remap[live[i][0]], remap[live[i][1]], remap[live[i][2]]);

  impl.vertPos_ = newPos;
  impl.properties_ = Vec<double>();
  impl.numProp_ = 0;
  impl.halfedgeTangent_ = Vec<vec4>();
  impl.CreateHalfedges(tris);
  impl.InitializeOriginal();
  impl.CalculateBBox();
  impl.SetEpsilon(-1);
  impl.SetNormalsAndCoplanar();
  impl.CalculateVertNormals();
  return true;
}
}  // namespace

namespace manifold {

/**
 * Surface-only morphological opening/closing by a ball of radius `radius`,
 * implemented as a curvature-bounded curvature flow. `close == true` performs
 * a closing (fills concavities a ball of the given radius cannot enter, the
 * output contains the input); `close == false` performs an opening (shaves
 * convex features a ball cannot reach from outside, the output is contained in
 * the input).
 *
 * @param radius        Structuring-element ball radius. Larger values change
 *                      more of the surface.
 * @param edgeLength    Target edge length for an optional uniform pre-refine
 *                      (0 = flow the input mesh as-is).
 * @param maxIterations Cap on semi-implicit flow steps.
 * @param anisotropic   true: faithful minimum-/maximum-curvature flow (paper).
 *                      false: isotropic mean-curvature flow (faster, less exact
 *                      on saddle regions).
 */
void Manifold::Impl::MorphologicalFlow(double radius, double edgeLength,
                                       int maxIterations, bool anisotropic,
                                       bool close) {
  if (IsEmpty() || radius <= 0) return;

  // Optional global pre-refine. Frozen verts stay pinned afterward, so the
  // shape of unchanged regions is preserved (only their tessellation differs).
  if (edgeLength > 0) {
    Refine(
        [edgeLength](vec3 edge, vec4, vec4) {
          return static_cast<int>(la::length(edge) / edgeLength);
        },
        false);
  }

  // Characteristic length: drives the timestep and is scale-invariant.
  double h = edgeLength;
  if (h <= 0) {
    double total = 0;
    size_t count = 0;
    for (size_t e = 0; e < halfedge_.size(); ++e) {
      if (!halfedge_.IsForward(e)) continue;
      total += la::length(vertPos_[halfedge_.End(e)] -
                          vertPos_[halfedge_.Start(e)]);
      ++count;
    }
    h = count > 0 ? total / count : 1.0;
  }
  const double tau = h * h;
  const double kBound = 1.0 / radius;
  const double dispTol = 1e-4 * h;

  // Pre-pass: anchor the changeable region to the original concavity (for
  // closing) / convexity (for opening). A vertex may move only if it is both
  // curvature-active AND within `capDist` of an anchor, so the moving front
  // cannot creep ring-by-ring into distant flat regions.
  const double capDist = radius + 4.0 * h;
  AnchorGrid grid;
  grid.cell = capDist;
  {
    SetNormalsAndCoplanar();
    CalculateVertNormals();
    const size_t nv = NumVert(), nt = NumTri();
    Vec<double> H(nv, 0.0), K(nv, kTwoPi), area(nv, 0.0), degree(nv, 0.0);
    for_each(autoPolicy(nt, 1e4), countAt(0_uz), countAt(nt),
             RawCurvature{H, K, area, degree, halfedge_, vertPos_, faceNormal_});
    for (size_t i = 0; i < nv; ++i) {
      if (area[i] <= 0) continue;
      const double factor = degree[i] / (6 * area[i]);
      const double meanC = H[i] * factor, gaussC = K[i] * factor;
      const double root = std::sqrt(std::max(meanC * meanC - 4 * gaussC, 0.0));
      const double kMin = 0.5 * (meanC - root), kMax = 0.5 * (meanC + root);
      if (close ? (kMin < -kBound) : (kMax > kBound)) grid.Insert(vertPos_[i]);
    }
  }

  for (int iter = 0; iter < maxIterations; ++iter) {
    const size_t numVert = NumVert();
    const size_t numTri = NumTri();
    SetNormalsAndCoplanar();
    CalculateVertNormals();

    // --- Per-vertex signed principal curvature and the frozen mask. ---
    Vec<double> H(numVert, 0.0), K(numVert, kTwoPi), area(numVert, 0.0),
        degree(numVert, 0.0);
    auto policy = autoPolicy(numTri, 1e4);
    for_each(
        policy, countAt(0_uz), countAt(numTri),
        RawCurvature{H, K, area, degree, halfedge_, vertPos_, faceNormal_});

    Vec<char> frozen(numVert, 1);
    Vec<double> mass(numVert, 0.0);
    size_t activeCount = 0;
    for (size_t i = 0; i < numVert; ++i) {
      const double m = area[i];
      mass[i] = m > 0 ? m : 1.0;
      if (m <= 0) continue;
      // Same valence-correcting normalization as Impl::CalculateCurvature:
      // meanC = k1 + k2 (sum), gaussC = k1 * k2 (product).
      const double factor = degree[i] / (6 * m);
      const double meanC = H[i] * factor;
      const double gaussC = K[i] * factor;
      const double root = std::sqrt(std::max(meanC * meanC - 4 * gaussC, 0.0));
      const double kMin = 0.5 * (meanC - root);
      const double kMax = 0.5 * (meanC + root);
      // Active (moving) vertices are the ones a ball of radius r cannot reach,
      // bounded to the neighborhood of the original concavity/convexity.
      const bool active = (close ? (kMin < -kBound) : (kMax > kBound)) &&
                          grid.Near(vertPos_[i], capDist);
      if (active) {
        frozen[i] = 0;
        ++activeCount;
      }
    }
    if (activeCount == 0) break;  // ball fits everywhere: converged.

    // --- Per-face principal-curvature direction (anisotropic only). ---
    Vec<vec3> faceDir;
    if (anisotropic) {
      faceDir.resize(numTri);
      for_each_n(policy, countAt(0_uz), numTri, [&](size_t tri) {
        const int v0 = halfedge_.Start(3 * tri);
        const int v1 = halfedge_.Start(3 * tri + 1);
        const int v2 = halfedge_.Start(3 * tri + 2);
        // Triangle corners plus the three flap verts (opposite corner of each
        // edge neighbor) form the quadric-fit stencil.
        vec3 pts[6];
        int nPts = 3;
        pts[0] = vertPos_[v0];
        pts[1] = vertPos_[v1];
        pts[2] = vertPos_[v2];
        for (int i = 0; i < 3; ++i) {
          const int pair = halfedge_.Pair(3 * tri + i);
          if (pair < 0) continue;
          const int nt = pair / 3;
          const int a = halfedge_.Start(3 * tri + i);
          const int b = halfedge_.End(3 * tri + i);
          for (int k = 0; k < 3; ++k) {
            const int w = halfedge_.Start(3 * nt + k);
            if (w != a && w != b) {
              pts[nPts++] = vertPos_[w];
              break;
            }
          }
        }
        const vec3 centroid = (pts[0] + pts[1] + pts[2]) / 3.0;
        faceDir[tri] =
            FacePrincipalDir(pts, nPts, centroid, faceNormal_[tri], close);
      });
    }

    // --- Assemble per-triangle stiffness and the system diagonal. ---
    Vec<FaceStiff> stiff(numTri);
    Vec<double> diag(numVert);
    for (size_t i = 0; i < numVert; ++i) diag[i] = mass[i];
    for (size_t tri = 0; tri < numTri; ++tri) {
      const int v0 = halfedge_.Start(3 * tri);
      const int v1 = halfedge_.Start(3 * tri + 1);
      const int v2 = halfedge_.Start(3 * tri + 2);
      const vec3 p0 = vertPos_[v0], p1 = vertPos_[v1], p2 = vertPos_[v2];
      const double a = TriArea(p0, p1, p2);
      if (a < 1e-18) {
        stiff[tri] = FaceStiff{{0, 0, 0, 0, 0, 0}};
        continue;
      }
      vec3 g0, g1, g2;
      HatGradients(p0, p1, p2, faceNormal_[tri], a, g0, g1, g2);
      stiff[tri] = BuildStiffness(g0, g1, g2, a, anisotropic,
                                  anisotropic ? faceDir[tri] : vec3(0.0));
      diag[v0] += tau * stiff[tri].s[0];
      diag[v1] += tau * stiff[tri].s[3];
      diag[v2] += tau * stiff[tri].s[5];
    }

    // Matrix-free application of A = M + tau * L to a full-length vector.
    auto applyFull = [&](VecView<const double> p, VecView<double> out) {
      for (size_t i = 0; i < numVert; ++i) out[i] = mass[i] * p[i];
      for (size_t tri = 0; tri < numTri; ++tri) {
        const int v0 = halfedge_.Start(3 * tri);
        const int v1 = halfedge_.Start(3 * tri + 1);
        const int v2 = halfedge_.Start(3 * tri + 2);
        const double x0 = p[v0], x1 = p[v1], x2 = p[v2];
        const FaceStiff& f = stiff[tri];
        out[v0] += tau * (f.s[0] * x0 + f.s[1] * x1 + f.s[2] * x2);
        out[v1] += tau * (f.s[1] * x0 + f.s[3] * x1 + f.s[4] * x2);
        out[v2] += tau * (f.s[2] * x0 + f.s[4] * x1 + f.s[5] * x2);
      }
    };
    // Reduced operator A_ff: frozen rows are zeroed (Dirichlet elimination).
    auto applyAff = [&](VecView<const double> p, VecView<double> out) {
      applyFull(p, out);
      for (size_t i = 0; i < numVert; ++i)
        if (frozen[i]) out[i] = 0.0;
    };

    Vec<double> diagInv(numVert);
    for (size_t i = 0; i < numVert; ++i)
      diagInv[i] = frozen[i] ? 1.0 : 1.0 / diag[i];

    // --- Solve each coordinate independently. ---
    Vec<vec3> newPos = vertPos_;
    Vec<double> gd(numVert), rhs(numVert), Agd(numVert), x(numVert);
    for (int c = 0; c < 3; ++c) {
      for (size_t i = 0; i < numVert; ++i) {
        const double pos = vertPos_[i][c];
        gd[i] = frozen[i] ? pos : 0.0;       // fixed Dirichlet values
        x[i] = frozen[i] ? 0.0 : pos;        // initial guess for free DOFs
      }
      applyFull(gd, Agd);                     // couple frozen values into RHS
      for (size_t i = 0; i < numVert; ++i)
        rhs[i] = frozen[i] ? 0.0 : (mass[i] * vertPos_[i][c] - Agd[i]);
      ConjugateGradient(applyAff, diagInv, rhs, x, /*maxIter=*/500,
                        /*relTol=*/1e-7);
      for (size_t i = 0; i < numVert; ++i)
        if (!frozen[i]) newPos[i][c] = x[i];
    }

    // --- Convergence test (frozen verts are bit-for-bit unchanged). ---
    double maxDisp = 0;
    for (size_t i = 0; i < numVert; ++i) {
      if (frozen[i]) continue;
      maxDisp = std::max(maxDisp, la::length(newPos[i] - vertPos_[i]));
    }
    vertPos_ = newPos;

    // Local adaptive remeshing of the moving region: split stretched edges,
    // collapse over-short ones, flip toward valence 6, and tangentially relax
    // -- all confined to active-active edges, so frozen (unchanged) regions
    // keep their exact geometry and triangulation.
    RemeshActive(*this, h, frozen);

    if (maxDisp < dispTol) break;
  }

  // Re-validate geometry exactly as the eager refine / constructor paths do.
  CalculateBBox();
  SetEpsilon(-1);
  SetNormalsAndCoplanar();
  RemoveDegenerates();
  CalculateVertNormals();
  SortGeometry();
  meshRelation_.originalID = -1;
}

}  // namespace manifold

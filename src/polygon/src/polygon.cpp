// Copyright 2021 The Manifold Authors.
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

#include "polygon.h"

#include <algorithm>
#if MANIFOLD_PAR == 'T' && __has_include(<pstl/glue_execution_defs.h>)
#include <execution>
#endif
#include <list>
#include <map>
#if __has_include(<memory_resource>)
#include <memory_resource>
#endif
#include <queue>
#include <set>
#include <stack>

#include "optional_assert.h"

namespace {
using namespace manifold;

static ExecutionParams params;

#ifdef MANIFOLD_DEBUG
struct PolyEdge {
  int startVert, endVert;
};

bool OverlapAssert(bool condition, const char *file, int line,
                   const std::string &cond, const std::string &msg) {
  if (!params.processOverlaps) {
    ASSERT(condition, geometryErr, msg);
  }
  return condition;
}

/**
 * Only use directly inside of the SweepForward() and SweepBack() functions! If
 * the asserted condition is false, it implies the monotone subdivision has
 * failed. This is most likely due to the input polygons being overlapped by
 * more than the input precision, but if not, then it indicates a bug. Either
 * way subdivision processing stops: if params.processOverlaps is false, then an
 * exception is thrown. Otherwise this returns true from the sweep function,
 * causing polygons to be left in their original state.
 *
 * The input polygons are then triangulated by the monotone triangulator, which
 * is robust enough to create a manifold triangulation for all input, but it
 * will not be geometrically-valid in this case. It may create inverted
 * triangles which are significantly larger than precision, but it depends on
 * the nature of the overlap.
 */
#define OVERLAP_ASSERT(condition, msg)                                \
  if (!OverlapAssert(condition, __FILE__, __LINE__, #condition, msg)) \
    return true;

#define PRINT(msg) \
  if (params.verbose) std::cout << msg << std::endl;

std::vector<PolyEdge> Polygons2Edges(const PolygonsIdx &polys) {
  std::vector<PolyEdge> halfedges;
  for (const auto &poly : polys) {
    for (int i = 1; i < poly.size(); ++i) {
      halfedges.push_back({poly[i - 1].idx, poly[i].idx});
    }
    halfedges.push_back({poly.back().idx, poly[0].idx});
  }
  return halfedges;
}

std::vector<PolyEdge> Triangles2Edges(
    const std::vector<glm::ivec3> &triangles) {
  std::vector<PolyEdge> halfedges;
  halfedges.reserve(triangles.size() * 3);
  for (const glm::ivec3 &tri : triangles) {
    halfedges.push_back({tri[0], tri[1]});
    halfedges.push_back({tri[1], tri[2]});
    halfedges.push_back({tri[2], tri[0]});
  }
  return halfedges;
}

void CheckTopology(const std::vector<PolyEdge> &halfedges) {
  ASSERT(halfedges.size() % 2 == 0, topologyErr, "Odd number of halfedges.");
  size_t n_edges = halfedges.size() / 2;
  std::vector<PolyEdge> forward(halfedges.size()), backward(halfedges.size());

  auto end = std::copy_if(halfedges.begin(), halfedges.end(), forward.begin(),
                          [](PolyEdge e) { return e.endVert > e.startVert; });
  ASSERT(std::distance(forward.begin(), end) == n_edges, topologyErr,
         "Half of halfedges should be forward.");
  forward.resize(n_edges);

  end = std::copy_if(halfedges.begin(), halfedges.end(), backward.begin(),
                     [](PolyEdge e) { return e.endVert < e.startVert; });
  ASSERT(std::distance(backward.begin(), end) == n_edges, topologyErr,
         "Half of halfedges should be backward.");
  backward.resize(n_edges);

  std::for_each(backward.begin(), backward.end(),
                [](PolyEdge &e) { std::swap(e.startVert, e.endVert); });
  auto cmp = [](const PolyEdge &a, const PolyEdge &b) {
    return a.startVert < b.startVert ||
           (a.startVert == b.startVert && a.endVert < b.endVert);
  };
  std::sort(forward.begin(), forward.end(), cmp);
  std::sort(backward.begin(), backward.end(), cmp);
  for (int i = 0; i < n_edges; ++i) {
    ASSERT(forward[i].startVert == backward[i].startVert &&
               forward[i].endVert == backward[i].endVert,
           topologyErr, "Forward and backward edge do not match.");
    if (i > 0) {
      ASSERT(forward[i - 1].startVert != forward[i].startVert ||
                 forward[i - 1].endVert != forward[i].endVert,
             topologyErr, "Not a 2-manifold.");
      ASSERT(backward[i - 1].startVert != backward[i].startVert ||
                 backward[i - 1].endVert != backward[i].endVert,
             topologyErr, "Not a 2-manifold.");
    }
  }
}

void CheckTopology(const std::vector<glm::ivec3> &triangles,
                   const PolygonsIdx &polys) {
  std::vector<PolyEdge> halfedges = Triangles2Edges(triangles);
  std::vector<PolyEdge> openEdges = Polygons2Edges(polys);
  for (PolyEdge e : openEdges) {
    halfedges.push_back({e.endVert, e.startVert});
  }
  CheckTopology(halfedges);
}

void CheckGeometry(const std::vector<glm::ivec3> &triangles,
                   const PolygonsIdx &polys, float precision) {
  std::unordered_map<int, glm::vec2> vertPos;
  for (const auto &poly : polys) {
    for (int i = 0; i < poly.size(); ++i) {
      vertPos[poly[i].idx] = poly[i].pos;
    }
  }
  ASSERT(std::all_of(triangles.begin(), triangles.end(),
                     [&vertPos, precision](const glm::ivec3 &tri) {
                       return CCW(vertPos[tri[0]], vertPos[tri[1]],
                                  vertPos[tri[2]], precision) >= 0;
                     }),
         geometryErr, "triangulation is not entirely CCW!");
}

void Dump(const PolygonsIdx &polys) {
  for (auto poly : polys) {
    std::cout << "polys.push_back({" << std::setprecision(9) << std::endl;
    for (auto v : poly) {
      std::cout << "    {" << v.pos.x << ", " << v.pos.y << "},  //"
                << std::endl;
    }
    std::cout << "});" << std::endl;
  }
  for (auto poly : polys) {
    std::cout << "array([" << std::endl;
    for (auto v : poly) {
      std::cout << "  [" << v.pos.x << ", " << v.pos.y << "]," << std::endl;
    }
    std::cout << "])" << std::endl;
  }
}

void PrintFailure(const std::exception &e, const PolygonsIdx &polys,
                  std::vector<glm::ivec3> &triangles, float precision) {
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Triangulation failed! Precision = " << precision << std::endl;
  std::cout << e.what() << std::endl;
  Dump(polys);
  std::cout << "produced this triangulation:" << std::endl;
  for (int j = 0; j < triangles.size(); ++j) {
    std::cout << triangles[j][0] << ", " << triangles[j][1] << ", "
              << triangles[j][2] << std::endl;
  }
}
#else
#define OVERLAP_ASSERT(condition, msg) \
  if (!(condition)) return true;
#define PRINT(msg)
#endif

/**
 * The class first turns input polygons into monotone polygons, then
 * triangulates them using the above class.
 */
class Monotones {
 public:
  Monotones(const PolygonsIdx &polys, float precision) : precision_(precision) {
    VertItr start, last, current;
    float bound = 0;
    for (const SimplePolygonIdx &poly : polys) {
      for (int i = 0; i < poly.size(); ++i) {
        monotones_.push_back({poly[i].pos,  //
                              poly[i].idx,  //
                              0, monotones_.end(), monotones_.end(),
                              activeEdges_.end()});
        bound = glm::max(
            bound, glm::max(glm::abs(poly[i].pos.x), glm::abs(poly[i].pos.y)));

        current = std::prev(monotones_.end());
        if (i == 0)
          start = current;
        else
          Link(last, current);
        last = current;
      }
      Link(current, start);
    }

    if (precision_ < 0) precision_ = bound * kTolerance;

    if (SweepForward()) return;
    Check();

    if (SweepBack()) return;
    Check();
  }

  void Triangulate(std::vector<glm::ivec3> &triangles) {
    // Save the sweep-line order in the vert to check further down.
    int i = 1;
    for (auto &vert : monotones_) {
      vert.index = i++;
    }
    int triangles_left = monotones_.size();
    VertItr start = monotones_.begin();
    while (start != monotones_.end()) {
      PRINT(start->mesh_idx);
      Triangulator triangulator(start, precision_);
      start->SetProcessed(true);
      VertItr vR = start->right;
      VertItr vL = start->left;
      while (vR != vL) {
        // Process the neighbor vert that is next in the sweep-line.
        if (vR->index < vL->index) {
          PRINT(vR->mesh_idx);
          triangulator.ProcessVert(vR, true, false, triangles);
          vR->SetProcessed(true);
          vR = vR->right;
        } else {
          PRINT(vL->mesh_idx);
          triangulator.ProcessVert(vL, false, false, triangles);
          vL->SetProcessed(true);
          vL = vL->left;
        }
      }
      PRINT(vR->mesh_idx);
      triangulator.ProcessVert(vR, true, true, triangles);
      vR->SetProcessed(true);
      // validation
      ASSERT(triangulator.NumTriangles() > 0, topologyErr,
             "Monotone produced no triangles.");
      triangles_left -= 2 + triangulator.NumTriangles();
      // Find next monotone
      start = std::find_if(monotones_.begin(), monotones_.end(),
                           [](const VertAdj &v) { return !v.Processed(); });
    }
    ASSERT(triangles_left == 0, topologyErr,
           "Triangulation produced wrong number of triangles.");
  }

  // A variety of sanity checks on the data structure. Expensive checks are only
  // performed if params.intermediateChecks = true.
  void Check() {
#ifdef MANIFOLD_DEBUG
    if (!params.intermediateChecks) return;
    std::vector<PolyEdge> edges;
    for (VertItr vert = monotones_.begin(); vert != monotones_.end(); vert++) {
      vert->SetProcessed(false);
      edges.push_back({vert->mesh_idx, vert->right->mesh_idx});
      ASSERT(vert->right->right != vert, topologyErr, "two-edge monotone!");
      ASSERT(vert->left->right == vert, topologyErr,
             "monotone vert neighbors don't agree!");
    }
    if (params.verbose) {
      VertItr start = monotones_.begin();
      while (start != monotones_.end()) {
        start->SetProcessed(true);
        PRINT("monotone start: " << start->mesh_idx << ", " << start->pos.y);
        VertItr v = start->right;
        while (v != start) {
          PRINT(v->mesh_idx << ", " << v->pos.y);
          v->SetProcessed(true);
          v = v->right;
        }
        PRINT("");
        start = std::find_if(monotones_.begin(), monotones_.end(),
                             [](const VertAdj &v) { return !v.Processed(); });
      }
    }
#endif
  }

  float GetPrecision() const { return precision_; }

 private:
  struct VertAdj;
  struct Edge;
  enum VertType { Start, Backward, Forward, Merge, End, Skip };
#if __has_include(<memory_resource>)
  typedef std::pmr::list<VertAdj>::iterator VertItr;
  typedef std::pmr::list<Edge>::iterator EdgeItr;

  std::pmr::monotonic_buffer_resource mbr;
  std::pmr::polymorphic_allocator<int> pa{&mbr};
  std::pmr::list<VertAdj> monotones_{pa};   // sweep-line list of verts
  std::pmr::list<Edge> activeEdges_{pa};    // west to east monotone edges
  std::pmr::list<Edge> inactiveEdges_{pa};  // completed monotones
#else
  typedef std::list<VertAdj>::iterator VertItr;
  typedef std::list<Edge>::iterator EdgeItr;

  std::list<VertAdj> monotones_;   // sweep-line list of verts
  std::list<Edge> activeEdges_;    // west to east monotone edges
  std::list<Edge> inactiveEdges_;  // completed monotones
#endif
  float precision_;  // a triangle of this height or less is degenerate

  /**
   * This is the data structure of the polygons themselves. They are stored as a
   * list in sweep-line order. The left and right pointers form the polygons,
   * while the mesh_idx describes the input indices that will be transferred to
   * the output triangulation. The edgeRight value represents an extra contraint
   * from the mesh Boolean algorithm.
   */
  struct VertAdj {
    glm::vec2 pos;
    int mesh_idx;  // This is a global index into the manifold.
    int index;
    VertItr left, right;
    EdgeItr edgeL, edgeR;

    bool Processed() const { return index < 0; }
    void SetSkip() { index = -2; }
    void SetProcessed(bool processed) {
      if (index == -2) return;
      index = processed ? -1 : 0;
    }
    bool IsStart() const {
      return (left->pos.y >= pos.y && right->pos.y > pos.y) ||
             (left->pos.y == pos.y && right->pos.y == pos.y &&
              left->pos.x <= pos.x && right->pos.x < pos.x);
    }
    bool IsPast(const VertItr other, float precision) const {
      return pos.y > other->pos.y + precision;
    }
    bool operator<(const VertAdj &other) const { return pos.y < other.pos.y; }
  };

  /**
   * The EdgePairs form the two active edges of a monotone polygon as they are
   * being constructed. The sweep-line is horizontal and moves from -y to +y, or
   * South to North. The West edge is a backwards edge while the East edge is
   * forwards, a topological constraint. If the polygon is geometrically valid,
   * then the West edge will also be to the -x side of the East edge, hence the
   * name.
   *
   * The purpose of the certainty booleans is to represent if we're sure the
   * pairs (or monotones) are in the right order. This is uncertain if they are
   * degenerate, for instance if several active edges are colinear (within
   * tolerance). If the order is uncertain, then as each vert is processed, if
   * it yields new information, it can cause the order to be updated until
   * certain.
   */

  struct Edge {
    VertItr south;         //, chain;
    EdgeItr linked, next;  //, bottomEast;
    bool forward, linked2east, flipped, eastCertain;
    float minDegenerateY = 0;

    VertItr North() const { return forward ? south->right : south->left; }

    int EastOf(VertItr vert, float precision) const {
      const VertItr north = North();
      if (south->pos.x - precision > vert->pos.x &&
          north->pos.x - precision > vert->pos.x)
        return 1;
      if (south->pos.x + precision < vert->pos.x &&
          north->pos.x + precision < vert->pos.x)
        return -1;
      return CCW(south->pos, north->pos, vert->pos, precision);
    }
  };

  /**
   * This class takes sequential verts of a monotone polygon and outputs a
   * geometrically valid triangulation, step by step.
   */
  class Triangulator {
   public:
    Triangulator(VertItr vert, float precision) : precision_(precision) {
      reflex_chain_.push_back(vert);
      other_side_ = vert;
    }
    int NumTriangles() const { return triangles_output_; }

    /**
     * The vert, vi, must attach to the free end (specified by onRight) of the
     * polygon that has been input so far. The verts must also be processed in
     * sweep-line order to get a geometrically valid result. If not, then the
     * polygon is not monotone, as the result should be topologically valid, but
     * not geometrically. The parameter, last, must be set true only for the
     * final point, as this ensures the last triangle is output.
     */
    void ProcessVert(const VertItr vi, bool onRight, bool last,
                     std::vector<glm::ivec3> &triangles) {
      VertItr v_top = reflex_chain_.back();
      if (reflex_chain_.size() < 2) {
        reflex_chain_.push_back(vi);
        onRight_ = onRight;
        return;
      }
      reflex_chain_.pop_back();
      VertItr vj = reflex_chain_.back();
      if (onRight_ == onRight && !last) {
        // This only creates enough triangles to ensure the reflex chain is
        // still reflex.
        PRINT("same chain");
        int ccw = CCW(vi->pos, vj->pos, v_top->pos, precision_);
        while (ccw == (onRight_ ? 1 : -1) || ccw == 0) {
          AddTriangle(triangles, vi, vj, v_top);
          v_top = vj;
          reflex_chain_.pop_back();
          if (reflex_chain_.empty()) break;
          vj = reflex_chain_.back();
          ccw = CCW(vi->pos, vj->pos, v_top->pos, precision_);
        }
        reflex_chain_.push_back(v_top);
        reflex_chain_.push_back(vi);
      } else {
        // This branch empties the reflex chain and switches sides. It must be
        // used for the last point, as it will output all the triangles
        // regardless of geometry.
        PRINT("different chain");
        onRight_ = !onRight_;
        VertItr v_last = v_top;
        while (!reflex_chain_.empty()) {
          vj = reflex_chain_.back();
          AddTriangle(triangles, vi, v_last, vj);
          v_last = vj;
          reflex_chain_.pop_back();
        }
        reflex_chain_.push_back(v_top);
        reflex_chain_.push_back(vi);
        other_side_ = v_top;
      }
    }

   private:
    std::vector<VertItr> reflex_chain_;
    VertItr other_side_;  // The end vertex across from the reflex chain
    bool onRight_;        // The side the reflex chain is on
    int triangles_output_ = 0;
    const float precision_;

    void AddTriangle(std::vector<glm::ivec3> &triangles, VertItr v0, VertItr v1,
                     VertItr v2) {
      if (!onRight_) std::swap(v1, v2);
      triangles.emplace_back(v0->mesh_idx, v1->mesh_idx, v2->mesh_idx);
      ++triangles_output_;
      PRINT(triangles.back());
    }
  };

  void Link(VertItr left, VertItr right) {
    left->right = right;
    right->left = left;
  }

  void UpdateEdge(EdgeItr edge, VertItr vert) {
    edge->south = vert;
    vert->edgeL = edge;
    vert->edgeR = edge;
  }

  void LinkEdges(EdgeItr edge1, EdgeItr edge2) {
    edge1->linked = edge2;
    edge2->linked = edge1;
  }

  void CloseEnd(VertItr vert) {
    EdgeItr edgeR = vert->right->edgeL;
    EdgeItr edgeL = vert->left->edgeR;
    edgeR->south = vert;
    edgeL->south = vert;
    vert->edgeR = edgeR;
    vert->edgeL = edgeL;
    edgeL->eastCertain = true;
    edgeR->eastCertain = true;
  }

  /**
   * This function is shared between the forward and backward sweeps and
   * determines the topology of the vertex relative to the sweep line.
   */
  VertType ProcessVert(VertItr vert) {
    if (vert->right->Processed()) {
      if (vert->left->Processed()) {
        const EdgeItr edgeR = vert->right->edgeL;
        const EdgeItr edgeL = vert->left->edgeR;
        // if (CCW(vert->pos, vert->right->pos, vert->left->pos, precision_) ==
        //     0) {
        //   if (edgeL->linked2east && !edgeR->linked2east) {
        //     activeEdges_.splice(edgeL, activeEdges_, edgeR->linked,
        //                         std::next(edgeR));
        //   } else if (!edgeL->linked2east && edgeR->linked2east) {
        //     activeEdges_.splice(edgeR, activeEdges_, edgeL->linked,
        //                         std::next(edgeL));
        //   }
        // }

        if (std::next(edgeR) != edgeL && std::next(edgeL) != edgeR) {
          PRINT("Skip");
          return Skip;
        }

        edgeR->south = vert;
        edgeL->south = vert;
        vert->edgeR = edgeR;
        vert->edgeL = edgeL;
        LinkEdges(edgeL->linked, edgeR->linked);

        if (std::next(edgeR) == edgeL) {  // facing in
          PRINT("End");
          return End;
        } else {  // facing out
          PRINT("Merge");
          return Merge;
        }
      } else {
        const EdgeItr bwdEdge = vert->right->edgeL;
        const EdgeItr fwdEdge = std::next(bwdEdge);
        if (!vert->IsPast(vert->right, precision_) &&
            !fwdEdge->south->right->IsPast(vert, precision_) &&
            vert->IsPast(fwdEdge->south, precision_) &&
            vert->pos.x > fwdEdge->south->right->pos.x + precision_) {
          PRINT("Skip backward edge");
          return Skip;
        }
        UpdateEdge(bwdEdge, vert);
        PRINT("Backward");
        return Backward;
      }
    } else {
      if (vert->left->Processed()) {
        const EdgeItr fwdEdge = vert->left->edgeR;
        const EdgeItr bwdEdge = std::prev(fwdEdge);
        if (!vert->IsPast(vert->left, precision_) &&
            !bwdEdge->south->left->IsPast(vert, precision_) &&
            vert->IsPast(bwdEdge->south, precision_) &&
            vert->pos.x < bwdEdge->south->left->pos.x - precision_) {
          PRINT("Skip forward edge");
          return Skip;
        }
        UpdateEdge(fwdEdge, vert);
        PRINT("Forward");
        return Forward;
      } else {
        PRINT("Start");
        return Start;
      }
    }
  }

  /**
   * Remove this edge and its pair to the east, but save them and mark the edge
   * they were next to. When the reverse sweep happens, it will be placed next
   * to its last neighbor instead of using geometry.
   */
  void RemovePair(EdgeItr westEdge) {
    EdgeItr eastEdge = std::next(westEdge);
    EdgeItr nextEast = std::next(eastEdge);
    westEdge->next = eastEdge->next = nextEast;
    inactiveEdges_.splice(inactiveEdges_.end(), activeEdges_, westEdge,
                          nextEast);
  }

  VertType PlaceStart(VertItr vert) {
    EdgeItr eastEdge = activeEdges_.begin();
    while (eastEdge != activeEdges_.end() && eastEdge->EastOf(vert, 0) <= 0) {
      ++eastEdge;
    }

    bool isHole = CCW(vert->left->pos, vert->pos, vert->right->pos, 0) < 0;
    const bool holeCertain =
        CCW(vert->left->pos, vert->pos, vert->right->pos, precision_) != 0;
    const bool shouldBeStart =
        eastEdge == activeEdges_.end() || !eastEdge->forward;

    if (isHole == shouldBeStart) {  // invalid
      if (!holeCertain) {
        isHole = !isHole;
      } else {  // shift to a valid position
        if (eastEdge != activeEdges_.end() &&
            eastEdge->EastOf(vert, precision_) <= 0) {
          ++eastEdge;
        } else if (eastEdge != activeEdges_.begin() &&
                   std::prev(eastEdge)->EastOf(vert, precision_) >= 0) {
          --eastEdge;
        } else {
          return Skip;
        }
      }
    }

    const EdgeItr noEdge = activeEdges_.end();
    const EdgeItr newEastEdge = activeEdges_.insert(
        eastEdge, {vert, noEdge, noEdge, !isHole, false, false,
                   eastEdge == activeEdges_.end() ||
                       eastEdge->EastOf(vert, precision_) > 0});
    const EdgeItr newWestEdge = activeEdges_.insert(
        newEastEdge, {vert, noEdge, noEdge, isHole, true, false, holeCertain});
    vert->edgeR = isHole ? newWestEdge : newEastEdge;
    vert->edgeL = isHole ? newEastEdge : newWestEdge;
    LinkEdges(newEastEdge, newWestEdge);
    return Start;
  }

  /**
   * This is the key function for handling east-west degeneracies, and is the
   * purpose of running the sweep-line forwards and backwards. If the ordering
   * of inputEdge is uncertain, this function uses the edge ahead of vert to
   * check if this new bit of geometric information is enough to place the pair
   * with certainty. It can also invert the pair if it is determined to be a
   * hole, in which case the inputEdge becomes the eastEdge while the pair it is
   * inside of becomes the westEdge.
   *
   * This function normally returns false, but will instead return true if the
   * certainties conflict, indicating this vertex is not yet geometrically valid
   * and must be skipped.
   */
  bool ShiftEast(const EdgeItr inputEdge, float precision) {
    if (inputEdge->eastCertain) return false;

    const VertItr vert = inputEdge->North();
    EdgeItr eastEdge = std::next(inputEdge);
    bool swap = false;
    bool past = false;
    while (eastEdge != activeEdges_.end() &&
           eastEdge->EastOf(vert, precision) <= 0) {
      if (eastEdge == inputEdge->linked) {
        swap = !swap;
        past = true;
      }
      ++eastEdge;
      swap = !swap;
    }

    if (swap) {
    } else {
    }

    return false;
  }

  /**
   * Identical to the above function, but swapped to search westward instead.
   */
  // bool ShiftWest(const EdgeItr inputEdge, float precision) {
  //   if (inputEdge == activeEdges_.begin() ||
  //   std::prev(inputEdge)->eastCertain)
  //     return false;

  //   EdgeItr westEdge = inputEdge;
  //   while (westEdge != activeEdges_.begin()) {
  //     --westEdge;
  //   }

  //   return false;
  // }

  /**
   * This function sweeps forward (South to North) keeping track of the
   * monotones and reordering degenerates (monotone ordering in the x-direction
   * and sweep line ordering in the y-direction). The input polygons
   * (monotones_) is not changed during this process.
   */
  bool SweepForward() {
    // Reversed so that minimum element is at queue.top() / vector.back().
    auto cmp = [](VertItr a, VertItr b) { return *b < *a; };
    std::priority_queue<VertItr, std::vector<VertItr>, decltype(cmp)>
        nextAttached(cmp);

    std::vector<VertItr> starts;
    for (VertItr v = monotones_.begin(); v != monotones_.end(); v++) {
      if (v->IsStart()) {
        starts.push_back(v);
      }
    }
#if MANIFOLD_PAR == 'T' && __has_include(<pstl/glue_execution_defs.h>)
    std::sort(std::execution::par_unseq, starts.begin(), starts.end(), cmp);
#else
    std::sort(starts.begin(), starts.end(), cmp);
#endif

    std::vector<VertItr> skipped;
    VertItr insertAt = monotones_.begin();

    while (insertAt != monotones_.end()) {
      // fallback for completely degenerate polygons that have no starts.
      VertItr vert = insertAt;
      if (!nextAttached.empty() &&
          (starts.empty() ||
           !nextAttached.top()->IsPast(starts.back(), precision_))) {
        // Prefer neighbors, which may process starts without needing a new
        // pair.
        vert = nextAttached.top();
        nextAttached.pop();
      } else if (!starts.empty()) {
        // Create a new pair with the next vert from the sorted list of starts.
        vert = starts.back();
        starts.pop_back();
      } else {
        ++insertAt;
      }

      if (vert->Processed()) continue;

      PRINT("mesh_idx = " << vert->mesh_idx);

      OVERLAP_ASSERT(
          skipped.empty() || !vert->IsPast(skipped.back(), precision_),
          "Not Geometrically Valid! None of the skipped verts is valid.");

      VertType type = ProcessVert(vert);

      if (type == Start) {
        type = PlaceStart(vert);
      }

      // if (type != Skip && ShiftEast(vert->edge, 0) &&
      //     ShiftEast(vert->edge, precision_))
      //   type = Skip;
      // if (type != Skip && ShiftWest(vert->edge, 0) &&
      //     ShiftWest(vert->edge, precision_))
      //   type = Skip;
      // if (type == Start) {
      //   if (ShiftEast(vert->edge->linked, 0) &&
      //       ShiftEast(vert->edge->linked, precision_))
      //     type = Skip;
      //   if (type != Skip && ShiftWest(vert->edge->linked, 0) &&
      //       ShiftWest(vert->edge->linked, precision_))
      //     type = Skip;
      // }

      if (type == Skip) {
        OVERLAP_ASSERT(std::next(insertAt) != monotones_.end(),
                       "Not Geometrically Valid! Tried to skip final vert.");
        OVERLAP_ASSERT(
            !nextAttached.empty() || !starts.empty(),
            "Not Geometrically Valid! Tried to skip last queued vert.");
        skipped.push_back(vert);
        PRINT("Skipping vert");
        continue;
      }

      if (vert == insertAt)
        ++insertAt;
      else
        monotones_.splice(insertAt, monotones_, vert);

      switch (type) {
        case Backward:
          nextAttached.push(vert->left);
          break;
        case Forward:
          nextAttached.push(vert->right);
          break;
        case Start:
          nextAttached.push(vert->left);
          nextAttached.push(vert->right);
          break;
        case Merge:
          RemovePair(vert->edgeL);
          break;
        case End:
          RemovePair(vert->edgeR);
          break;
        case Skip:
          break;
      }

      vert->SetProcessed(true);
      // Push skipped verts back into unprocessed queue.
      while (!skipped.empty()) {
        starts.push_back(skipped.back());
        skipped.pop_back();
      }

#ifdef MANIFOLD_DEBUG
      if (params.verbose) ListActive();
#endif
    }
    return false;
  }  // namespace

  /**
   * This is the only function that actually changes monotones_; all the rest is
   * bookkeeping. This divides polygons by connecting two verts. It duplicates
   * these verts to break the polygons, then attaches them across to each other
   * with two new edges.
   */
  VertItr SplitVerts(VertItr north, VertItr south) {
    // at split events, add duplicate vertices to end of list and reconnect
    PRINT("split from " << north->mesh_idx << " to " << south->mesh_idx);

    VertItr northEast = monotones_.insert(north, *north);
    Link(north->left, northEast);
    northEast->SetProcessed(true);

    VertItr southEast = monotones_.insert(std::next(south), *south);
    Link(southEast, south->right);
    southEast->SetProcessed(true);

    Link(south, north);
    Link(northEast, southEast);

    return northEast;
  }

  VertItr CheckSplit(VertItr vert, EdgeItr westEdge) {
    if (westEdge->next != activeEdges_.end()) {
      vert = SplitVerts(vert, westEdge->next->south);
      westEdge->next = activeEdges_.end();  // unmark merge
    }
    return vert;
  }

  /**
   * This function sweeps back, splitting the input polygons
   * into monotone polygons without doing a single geometric calculation.
   * Instead everything is based on the topology saved from the forward sweep,
   * primarily the relative ordering of new monotones. Even though the sweep is
   * going back, the polygon is considered rotated, so we still refer to
   * sweeping from South to North and the pairs as ordered from West to East
   * (though this is now the opposite order from the forward sweep).
   */
  bool SweepBack() {
    for (auto &vert : monotones_) vert.SetProcessed(false);

    VertItr vert = monotones_.end();
    while (vert != monotones_.begin()) {
      --vert;

      if (vert->Processed()) continue;

      PRINT("mesh_idx = " << vert->mesh_idx);

      VertType type = ProcessVert(vert);
      OVERLAP_ASSERT(type != Skip, "Skip should not happen on reverse sweep!");

      if (type == Merge) {
        vert = CheckSplit(vert, vert->edgeR);
        const EdgeItr westOf = std::prev(vert->edgeL);
        CheckSplit(vert, westOf);
        westOf->next = vert->edgeL;
      } else if (type == End) {
        CheckSplit(vert, vert->edgeR);
      }

      if (type == Merge || type == End) {
        inactiveEdges_.splice(inactiveEdges_.end(), activeEdges_, vert->edgeR);
        inactiveEdges_.splice(inactiveEdges_.end(), activeEdges_, vert->edgeL);
      } else if (type == Forward) {
        CheckSplit(vert, std::prev(vert->edgeL));
      } else if (type == Backward) {
        CheckSplit(vert, vert->edgeR);
      } else if (type == Start) {
        // Due to sweeping in the opposite direction, east and west are
        // swapped and what was the next pair is now the previous pair and
        // begin and end are swapped.
        EdgeItr westEdge = vert->edgeL;
        EdgeItr eastEdge = vert->edgeR;
        EdgeItr eastOf = westEdge->next;

        if (std::next(eastEdge) == westEdge) std::swap(eastEdge, westEdge);

        if (!westEdge->flipped) {
          std::swap(westEdge, eastEdge);
          eastOf = eastOf == activeEdges_.end() ? activeEdges_.begin()
                                                : std::next(eastOf);
        }

        activeEdges_.splice(eastOf, inactiveEdges_, eastEdge);
        activeEdges_.splice(eastEdge, inactiveEdges_, westEdge);
        westEdge->forward ^= true;
        eastEdge->forward ^= true;
        const bool isHole = westEdge->forward;

        if (isHole) {
          EdgeItr westOf = std::prev(westEdge);
          VertItr split =
              westOf->next != activeEdges_.end() ? westOf->next->south
              : westOf->south->pos.y < eastOf->south->pos.y ? eastOf->south
                                                            : westOf->south;
          VertItr eastVert = SplitVerts(vert, split);
          westOf->next = activeEdges_.end();
          UpdateEdge(eastEdge, eastVert);
          UpdateEdge(westEdge, vert);
        } else {  // Start
          vert->edgeL = westEdge;
          vert->edgeR = eastEdge;
        }
        westEdge->next = activeEdges_.end();
        eastEdge->next = activeEdges_.end();
      }

      vert->SetProcessed(true);

#ifdef MANIFOLD_DEBUG
      if (params.verbose) ListActive();
#endif
    }
    return false;
  }
#ifdef MANIFOLD_DEBUG
  void ListEdge(const EdgeItr edge) const {
    std::cout << (edge->forward ? "Fwd" : "Bwd");
    std::cout << ": S = " << edge->south->mesh_idx
              << ", N = " << edge->North()->mesh_idx;
    std::cout << (edge->next == activeEdges_.end() ? " none" : " next");
    std::cout << (edge->eastCertain ? " certain" : " uncertain") << std::endl;
    EdgeItr same = edge->forward ? edge->south->edgeR : edge->south->edgeL;
    if (same != edge) std::cout << "edgeR does not point back!" << std::endl;
  }

  void ListActive() {
    std::cout << "active edges:" << std::endl;
    for (auto edge = activeEdges_.begin(); edge != activeEdges_.end(); ++edge)
      ListEdge(edge);
  }
#endif
};
}  // namespace

namespace manifold {

/**
 * @brief Triangulates a set of &epsilon;-valid polygons. If the input is not
 * &epsilon;-valid, the triangulation may overlap, but will always return a
 * manifold result that matches the input edge directions.
 *
 * @param polys The set of polygons, wound CCW and representing multiple
 * polygons and/or holes. These have 2D-projected positions as well as
 * references back to the original vertices.
 * @param precision The value of &epsilon;, bounding the uncertainty of the
 * input.
 * @return std::vector<glm::ivec3> The triangles, referencing the original
 * vertex indicies.
 */
std::vector<glm::ivec3> TriangulateIdx(const PolygonsIdx &polys,
                                       float precision) {
  std::vector<glm::ivec3> triangles;
  try {
    Monotones monotones(polys, precision);
    monotones.Triangulate(triangles);
#ifdef MANIFOLD_DEBUG
    if (params.intermediateChecks) {
      CheckTopology(triangles, polys);
      if (!params.processOverlaps) {
        CheckGeometry(triangles, polys, 2 * monotones.GetPrecision());
      }
    }
  } catch (const geometryErr &e) {
    if (!params.suppressErrors) {
      PrintFailure(e, polys, triangles, precision);
    }
    throw;
  } catch (const std::exception &e) {
    PrintFailure(e, polys, triangles, precision);
    throw;
#else
  } catch (const std::exception &e) {
#endif
  }
  return triangles;
}

/**
 * @brief Triangulates a set of &epsilon;-valid polygons. If the input is not
 * &epsilon;-valid, the triangulation may overlap, but will always return a
 * manifold result that matches the input edge directions.
 *
 * @param polygons The set of polygons, wound CCW and representing multiple
 * polygons and/or holes.
 * @param precision The value of &epsilon;, bounding the uncertainty of the
 * input.
 * @return std::vector<glm::ivec3> The triangles, referencing the original
 * polygon points in order.
 */
std::vector<glm::ivec3> Triangulate(const Polygons &polygons, float precision) {
  int idx = 0;
  PolygonsIdx polygonsIndexed;
  for (const auto &poly : polygons) {
    SimplePolygonIdx simpleIndexed;
    for (const glm::vec2 &polyVert : poly) {
      simpleIndexed.push_back({polyVert, idx++});
    }
    polygonsIndexed.push_back(simpleIndexed);
  }
  return TriangulateIdx(polygonsIndexed, precision);
}

ExecutionParams &PolygonParams() { return params; }

}  // namespace manifold

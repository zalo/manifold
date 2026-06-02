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
#include <cmath>
#include <functional>

#include "parallel.h"
#include "vec.h"

namespace manifold {

// Small, dependency-free linear-algebra helpers operating on Vec<double> of
// arbitrary length. Manifold's linalg.h only covers fixed-size (<=4) vectors,
// so these fill the gap needed by the matrix-free Conjugate Gradient solver
// used for the semi-implicit curvature-flow step (see morphology.cpp).

// dot = a . b
inline double VecDot(VecView<const double> a, VecView<const double> b) {
  return transform_reduce(
      autoPolicy(a.size()), countAt(0_uz), countAt(a.size()), 0.0,
      std::plus<double>(), [&a, &b](size_t i) { return a[i] * b[i]; });
}

// y += alpha * x
inline void VecAxpy(double alpha, VecView<const double> x, VecView<double> y) {
  for_each_n(autoPolicy(x.size()), countAt(0_uz), x.size(),
             [alpha, &x, &y](size_t i) { y[i] += alpha * x[i]; });
}

// y = x + beta * y
inline void VecXpby(VecView<const double> x, double beta, VecView<double> y) {
  for_each_n(autoPolicy(x.size()), countAt(0_uz), x.size(),
             [beta, &x, &y](size_t i) { y[i] = x[i] + beta * y[i]; });
}

/**
 * Matrix-free preconditioned Conjugate Gradient for a symmetric
 * positive-definite system A x = b.
 *
 * The operator A is supplied as a callable `applyA(in, out)` that computes
 * `out = A * in`; A is never materialized. `diagInv` is the inverse diagonal
 * of A (the Jacobi preconditioner). `x` is used as the initial guess and
 * receives the solution. Iterates until the relative residual drops below
 * `relTol` or `maxIter` is reached. Returns the number of iterations used.
 *
 * For the semi-implicit flow this is called once per coordinate. Frozen
 * (Dirichlet) vertices are handled by the caller's operator/RHS: their rows
 * are identity, so the residual there starts and stays at zero and CG never
 * perturbs them.
 */
template <typename ApplyA>
int ConjugateGradient(ApplyA&& applyA, VecView<const double> diagInv,
                      VecView<const double> b, VecView<double> x, int maxIter,
                      double relTol) {
  const size_t n = b.size();
  if (n == 0) return 0;
  Vec<double> r(n), z(n), p(n), Ap(n);

  applyA(x, r);                      // r = A x
  VecView<double> rView = r;         // r = b - A x
  for_each_n(autoPolicy(n), countAt(0_uz), n,
             [&rView, &b](size_t i) { rView[i] = b[i] - rView[i]; });

  const double bNorm2 = VecDot(b, b);
  const double thresh2 = relTol * relTol * (bNorm2 > 0 ? bNorm2 : 1.0);
  if (VecDot(r, r) <= thresh2) return 0;

  // z = M^-1 r ; p = z
  for_each_n(autoPolicy(n), countAt(0_uz), n, [&z, &r, &p, &diagInv](size_t i) {
    z[i] = diagInv[i] * r[i];
    p[i] = z[i];
  });
  double rzOld = VecDot(r, z);

  int iter = 0;
  for (; iter < maxIter; ++iter) {
    applyA(p, Ap);
    const double pAp = VecDot(p, Ap);
    if (!(pAp > 0)) break;  // breakdown / non-positive curvature guard
    const double alpha = rzOld / pAp;
    VecAxpy(alpha, p, x);                                       // x += a p
    VecAxpy(-alpha, Ap, r);                                     // r -= a Ap
    if (VecDot(r, r) <= thresh2) {
      ++iter;
      break;
    }
    for_each_n(autoPolicy(n), countAt(0_uz), n,
               [&z, &r, &diagInv](size_t i) { z[i] = diagInv[i] * r[i]; });
    const double rzNew = VecDot(r, z);
    const double beta = rzNew / rzOld;
    VecXpby(z, beta, p);  // p = z + beta p
    rzOld = rzNew;
  }
  return iter;
}

}  // namespace manifold

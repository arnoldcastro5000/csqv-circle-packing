// scalable_polish: a memory-light center+radius polish (scalable_polish.py port).
//
// It maximizes sum(r) minus a quadratic penalty for overlap and wall violations, over a
// rising penalty weight (continuation), using a projected L-BFGS. After the optimize,
// radii_lp projects to the EXACTLY optimal radii for the resulting centers and the
// verifier scores them. The optimizer is a means, not the objective: the verifier is the
// score gate, so an L-BFGS that reaches a comparable penalty minimum is sufficient.
//
// Variable layout matches the Python packing.ravel(): v = [x0,y0,r0, x1,y1,r1, ...].
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "geometry.hpp"
#include "lp.hpp"

namespace csqv {

// Penalty objective f(v) and gradient g(v). Faithful port of scalable_polish.py::_obj_grad.
inline double penalty_obj_grad(const std::vector<double>& v, int n, const std::vector<int>& pi,
                               const std::vector<int>& pj, double lam, std::vector<double>& g) {
  std::fill(g.begin(), g.end(), 0.0);
  double f = 0.0;
  for (int i = 0; i < n; ++i) {
    f -= v[3 * i + 2];       // -sum(r)
    g[3 * i + 2] -= 1.0;
  }
  const size_t P = pi.size();
  for (size_t p = 0; p < P; ++p) {
    const int i = pi[p], j = pj[p];
    const double xi = v[3 * i], yi = v[3 * i + 1], ri = v[3 * i + 2];
    const double xj = v[3 * j], yj = v[3 * j + 1], rj = v[3 * j + 2];
    const double dx = xi - xj, dy = yi - yj;
    const double d = std::sqrt(dx * dx + dy * dy) + 1e-12;
    const double o = (ri + rj) - d;
    if (o > 0.0) {
      f += lam * o * o;
      const double w = 2.0 * lam * o;
      const double ux = dx / d, uy = dy / d;
      g[3 * i + 2] += w;
      g[3 * j + 2] += w;
      g[3 * i] += w * (-ux);
      g[3 * j] += w * (ux);
      g[3 * i + 1] += w * (-uy);
      g[3 * j + 1] += w * (uy);
    }
  }
  for (int i = 0; i < n; ++i) {
    const double x = v[3 * i], y = v[3 * i + 1], r = v[3 * i + 2];
    double vl = r - x, vr = x + r - SIDE, vb = r - y, vt = y + r - SIDE;
    if (vl > 0.0) { f += lam * vl * vl; g[3 * i + 2] += 2.0 * lam * vl; g[3 * i] += 2.0 * lam * vl * (-1.0); }
    if (vr > 0.0) { f += lam * vr * vr; g[3 * i + 2] += 2.0 * lam * vr; g[3 * i] += 2.0 * lam * vr * (1.0); }
    if (vb > 0.0) { f += lam * vb * vb; g[3 * i + 2] += 2.0 * lam * vb; g[3 * i + 1] += 2.0 * lam * vb * (-1.0); }
    if (vt > 0.0) { f += lam * vt * vt; g[3 * i + 2] += 2.0 * lam * vt; g[3 * i + 1] += 2.0 * lam * vt * (1.0); }
  }
  return f;
}

// Clamp v into the box: x,y in [0, SIDE], r >= 0.
inline void project_box(std::vector<double>& v, int n) {
  for (int i = 0; i < n; ++i) {
    v[3 * i] = std::min(std::max(v[3 * i], 0.0), SIDE);
    v[3 * i + 1] = std::min(std::max(v[3 * i + 1], 0.0), SIDE);
    if (v[3 * i + 2] < 0.0) v[3 * i + 2] = 0.0;
  }
}

// Projected L-BFGS on the penalty objective for a fixed lambda.
inline void lbfgs_penalty(std::vector<double>& v, int n, const std::vector<int>& pi,
                          const std::vector<int>& pj, double lam, int maxiter) {
  const int dim = 3 * n;
  const int hist = 10;
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  std::vector<double> g(dim), gnew(dim), q(dim), d(dim), vnew(dim);
  double f = penalty_obj_grad(v, n, pi, pj, lam, g);

  for (int it = 0; it < maxiter; ++it) {
    // Two-loop recursion: d = -H * g.
    q = g;
    const int hs = static_cast<int>(S.size());
    std::vector<double> alpha(hs);
    for (int k = hs - 1; k >= 0; --k) {
      double sq = 0.0;
      for (int t = 0; t < dim; ++t) sq += S[k][t] * q[t];
      alpha[k] = rho[k] * sq;
      for (int t = 0; t < dim; ++t) q[t] -= alpha[k] * Y[k][t];
    }
    double gamma = 1.0;
    if (hs > 0) {
      double sy = 0.0, yy = 0.0;
      for (int t = 0; t < dim; ++t) { sy += S[hs - 1][t] * Y[hs - 1][t]; yy += Y[hs - 1][t] * Y[hs - 1][t]; }
      if (yy > 1e-30) gamma = sy / yy;
    }
    for (int t = 0; t < dim; ++t) q[t] *= gamma;
    for (int k = 0; k < hs; ++k) {
      double yq = 0.0;
      for (int t = 0; t < dim; ++t) yq += Y[k][t] * q[t];
      const double beta = rho[k] * yq;
      for (int t = 0; t < dim; ++t) q[t] += S[k][t] * (alpha[k] - beta);
    }
    for (int t = 0; t < dim; ++t) d[t] = -q[t];

    // Directional derivative; fall back to steepest descent if not a descent direction.
    double gd = 0.0;
    for (int t = 0; t < dim; ++t) gd += g[t] * d[t];
    if (gd >= 0.0) {
      for (int t = 0; t < dim; ++t) d[t] = -g[t];
      gd = 0.0;
      for (int t = 0; t < dim; ++t) gd += g[t] * d[t];
    }

    // Projected backtracking Armijo line search.
    double step = 1.0;
    const double c1 = 1e-4;
    double fnew = f;
    bool ok = false;
    for (int ls = 0; ls < 30; ++ls) {
      for (int t = 0; t < dim; ++t) vnew[t] = v[t] + step * d[t];
      project_box(vnew, n);
      fnew = penalty_obj_grad(vnew, n, pi, pj, lam, gnew);
      if (fnew <= f + c1 * step * gd) { ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;  // no progress at this lambda

    std::vector<double> s(dim), yv(dim);
    for (int t = 0; t < dim; ++t) { s[t] = vnew[t] - v[t]; yv[t] = gnew[t] - g[t]; }
    double sy = 0.0;
    for (int t = 0; t < dim; ++t) sy += s[t] * yv[t];
    if (sy > 1e-12) {
      if (static_cast<int>(S.size()) == hist) { S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin()); }
      S.push_back(s);
      Y.push_back(yv);
      rho.push_back(1.0 / sy);
    }
    v = vnew;
    g = gnew;
    f = fnew;

    double gn = 0.0;
    for (int t = 0; t < dim; ++t) gn += g[t] * g[t];
    if (gn < 1e-20) break;
  }
}

// The exact scored objective for FIXED centers: the exact-LP radii, then the verifier
// score. radii_lp is already strictly feasible, so the score equals sum(radii_lp) here.
// This is what the search maximizes, so the center polish maximizes exactly it.
inline double lp_score(const std::vector<double>& x, const std::vector<double>& y, int n) {
  std::vector<double> r = radii_lp(x, y, n);
  double s;
  if (!verify_and_score(x, y, r, n, TOL, s)) return -std::numeric_limits<double>::infinity();
  return s;
}

// Exact-LP center polish: locally maximize lp_score(centers) over the CENTERS ONLY, the
// radii always the exact LP optimum for the fixed centers. This is the terminal squeeze on
// a champion. The penalty polish (scalable_polish) stops at a soft-penalty minimum that
// leaves ~1e-5 of sum unclaimed (the gap Packomania's own re-optimization closed on N=121);
// this recovers it by ascending the true objective.
//
// Method: projected gradient ascent. The gradient of lp_score comes from a central finite
// difference over each center coordinate (analytic duals are not exposed by the lazy LP),
// and a backtracking line search accepts only STRICTLY improving steps. So the scored sum
// is MONOTONE non-decreasing: it is always safe to run before comparing a candidate to the
// incumbent best.
//
// Cost is O(n) LP solves per gradient. The worker runs it on each NEW best (and once on a
// resumed champion) so the saved champion and its emitted .pck are always fully squeezed,
// never a stale penalty-polish value. New bests are rare once a run warms up, so this is not
// a hot path; at large N a cold start sees more early bests, where the cost is highest (see
// ticket 34 for optional gating).
//
// Updates x, y (centers) and r (the exact radii for the polished centers). Returns the
// polished score. On a non-finite start it leaves the centers untouched and returns it.
inline double center_polish(std::vector<double>& x, std::vector<double>& y,
                            std::vector<double>& r, int n, int max_iter = 300,
                            double eps = 1e-6, double step0 = 1e-4) {
  double base = lp_score(x, y, n);
  if (!std::isfinite(base)) {
    r = radii_lp(x, y, n);
    return base;
  }
  // Minimum accepted improvement per step. The C++ exact LP is deterministic and strictly
  // feasible (far more precise than scipy's ~1e-9), and the emitted .pck is re-verified at
  // zero tolerance downstream, so this can sit well below scipy's floor to capture real
  // sub-1e-9 gains near the optimum without banking float noise.
  const double ftol = 1e-11;
  const double step_max = 1e-2;
  std::vector<double> gx(n), gy(n), tx(n), ty(n);
  double step = step0;

  for (int it = 0; it < max_iter; ++it) {
    // Central-difference gradient of lp_score. Near the optimum the gain per step is ~1e-6,
    // so a one-sided difference is too biased/noisy and stalls early; the symmetric
    // difference cancels the first-order error and keeps the ascent moving. Two solves per
    // coordinate. The perturbation stays inside the box; the divisor uses the actual span.
    double gnorm2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double ox = x[i];
      const double xp = std::min(ox + eps, SIDE), xm = std::max(ox - eps, 0.0);
      x[i] = xp;
      const double fxp = lp_score(x, y, n);
      x[i] = xm;
      const double fxm = lp_score(x, y, n);
      x[i] = ox;
      gx[i] = (fxp - fxm) / (xp - xm);

      const double oy = y[i];
      const double yp = std::min(oy + eps, SIDE), ym = std::max(oy - eps, 0.0);
      y[i] = yp;
      const double fyp = lp_score(x, y, n);
      y[i] = ym;
      const double fym = lp_score(x, y, n);
      y[i] = oy;
      gy[i] = (fyp - fym) / (yp - ym);

      gnorm2 += gx[i] * gx[i] + gy[i] * gy[i];
    }
    if (gnorm2 <= 0.0) break;  // exactly flat gradient: nothing to ascend (and guards 1/0)
    const double inv = 1.0 / std::sqrt(gnorm2);

    // Backtracking ascent along the normalized gradient; accept the first strict gain.
    bool improved = false;
    double s = step;
    for (int ls = 0; ls < 40; ++ls) {
      for (int i = 0; i < n; ++i) {
        tx[i] = std::min(std::max(x[i] + s * gx[i] * inv, 0.0), SIDE);
        ty[i] = std::min(std::max(y[i] + s * gy[i] * inv, 0.0), SIDE);
      }
      const double f = lp_score(tx, ty, n);
      if (f > base + ftol) {
        x = tx;
        y = ty;
        base = f;
        step = std::min(s * 2.0, step_max);  // grow the trial step for the next iteration
        improved = true;
        break;
      }
      s *= 0.5;
    }
    if (!improved) break;  // no improving step at any scale: local optimum reached
  }

  r = radii_lp(x, y, n);
  double sc;
  if (!verify_and_score(x, y, r, n, TOL, sc)) return base;
  return sc;
}

// Full scalable polish: continuation over rising lambda, then exact radii via radii_lp.
// Returns the packing (unit frame) and its verifier score.
inline double scalable_polish(std::vector<double>& x, std::vector<double>& y,
                              std::vector<double>& r, int n) {
  // All upper-triangle pairs for the penalty (matches the dense Python objective).
  std::vector<int> pi, pj;
  pi.reserve(static_cast<size_t>(n) * (n - 1) / 2);
  pj.reserve(static_cast<size_t>(n) * (n - 1) / 2);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) { pi.push_back(i); pj.push_back(j); }

  std::vector<double> v(3 * n);
  for (int i = 0; i < n; ++i) { v[3 * i] = x[i]; v[3 * i + 1] = y[i]; v[3 * i + 2] = r[i]; }

  const double rounds[] = {1e2, 1e3, 1e4, 1e5};
  for (double lam : rounds) lbfgs_penalty(v, n, pi, pj, lam, 200);

  for (int i = 0; i < n; ++i) { x[i] = v[3 * i]; y[i] = v[3 * i + 1]; }
  r = radii_lp(x, y, n);
  double score;
  if (!verify_and_score(x, y, r, n, TOL, score)) score = -1.0;
  return score;
}

}  // namespace csqv

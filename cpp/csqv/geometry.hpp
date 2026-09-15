// CSQV geometry: frame, containment, the exact feasibility shrink, and the scorer.
//
// This is a faithful C++ port of the trusted Python primitives in
// problems/csqv/data.py and problems/csqv/problem.py. The Python verifier stays the
// source of truth: a champion emitted here is always re-checked by the Python path
// before any submission. This header exists so the C++ search scores candidates the
// same way the harness does, so the search optimizes the real objective.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace csqv {

constexpr double SIDE = 1.0;
constexpr double TOL = 1e-9;  // matches data.py TOL

// A packing in the UNIT-square frame [0,1]x[0,1]. Columns x, y, r.
struct Packing {
  int n = 0;
  std::vector<double> x, y, r;
  explicit Packing(int n_ = 0) : n(n_), x(n_), y(n_), r(n_) {}
};

// Distance from a center to the nearest wall of the square [0, SIDE]^2. A circle at
// (xi, yi) fits exactly when its radius is at most this value (data.py::wall_slack).
inline double wall_slack(double xi, double yi) {
  return std::min(std::min(xi, SIDE - xi), std::min(yi, SIDE - yi));
}

// Largest uniform factor s in [0, 1] that makes the packing strictly valid for FIXED
// centers (problem.py::_feasible_scale). Centers must already sit inside the square.
// Returns false only for two coincident centers, which no shrink can separate.
inline bool feasible_scale(const std::vector<double>& x, const std::vector<double>& y,
                           const std::vector<double>& r, int n, double tol,
                           double& s_out) {
  double s = 1.0;
  for (int i = 0; i < n; ++i) {
    const double slack = std::max(wall_slack(x[i], y[i]), 0.0);
    if (r[i] > 0.0) s = std::min(s, slack / r[i]);
  }
  if (n >= 2) {
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double dx = x[i] - x[j], dy = y[i] - y[j];
        const double d = std::sqrt(dx * dx + dy * dy);
        if (d < tol) {  // coincident centers, uncorrectable
          s_out = 0.0;
          return false;
        }
        const double ps = r[i] + r[j];
        if (ps > 0.0) s = std::min(s, d / ps);
      }
    }
  }
  s_out = std::max(0.0, std::min(1.0, s));
  return true;
}

// Recompute feasibility and score the shrunk sum of radii (problem.py::verify_and_score).
// Returns false with score -inf for an INVALID packing (non-finite, negative radius, a
// center outside the square by more than tol, or coincident centers). Otherwise the
// score is s * sum(r) on the strictly-valid shrunk packing, so it cannot be gamed.
inline bool verify_and_score(const std::vector<double>& px, const std::vector<double>& py,
                             const std::vector<double>& pr, int n, double tol,
                             double& score) {
  const double neg_inf = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(px[i]) || !std::isfinite(py[i]) || !std::isfinite(pr[i])) {
      score = neg_inf;
      return false;
    }
  }
  for (int i = 0; i < n; ++i) {
    if (pr[i] < -tol) {
      score = neg_inf;
      return false;
    }
  }
  for (int i = 0; i < n; ++i) {
    if (px[i] < -tol || px[i] > SIDE + tol || py[i] < -tol || py[i] > SIDE + tol) {
      score = neg_inf;
      return false;
    }
  }
  std::vector<double> x(n), y(n), r(n);
  for (int i = 0; i < n; ++i) {
    x[i] = std::clamp(px[i], 0.0, SIDE);
    y[i] = std::clamp(py[i], 0.0, SIDE);
    r[i] = std::max(pr[i], 0.0);
  }
  double s;
  if (!feasible_scale(x, y, r, n, tol, s)) {
    score = neg_inf;
    return false;
  }
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += r[i];
  score = s * sum;
  return true;
}

// The strictly-valid packing the verifier scores (problem.py::canonical_feasible): centers
// clipped into the square, every radius scaled by the uniform factor s. Use this to EMIT a
// packing. Returns false if INVALID.
inline bool canonical_feasible(const std::vector<double>& px, const std::vector<double>& py,
                               const std::vector<double>& pr, int n, double tol,
                               Packing& out) {
  const double neg_inf = -std::numeric_limits<double>::infinity();
  (void)neg_inf;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(px[i]) || !std::isfinite(py[i]) || !std::isfinite(pr[i])) return false;
  }
  for (int i = 0; i < n; ++i)
    if (pr[i] < -tol) return false;
  for (int i = 0; i < n; ++i)
    if (px[i] < -tol || px[i] > SIDE + tol || py[i] < -tol || py[i] > SIDE + tol) return false;
  std::vector<double> x(n), y(n), r(n);
  for (int i = 0; i < n; ++i) {
    x[i] = std::clamp(px[i], 0.0, SIDE);
    y[i] = std::clamp(py[i], 0.0, SIDE);
    r[i] = std::max(pr[i], 0.0);
  }
  double s;
  if (!feasible_scale(x, y, r, n, tol, s)) return false;
  out = Packing(n);
  for (int i = 0; i < n; ++i) {
    out.x[i] = x[i];
    out.y[i] = y[i];
    out.r[i] = s * r[i];
  }
  return true;
}

}  // namespace csqv

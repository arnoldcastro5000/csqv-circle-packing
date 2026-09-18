// jam_slp: the trust-region SLP jamming optimizer, a C++ port of
// problems/csqv/terminal_squeeze.py::jam_slp (and its helpers _active_set, _active_jacobian,
// _active_residual, donev_jamming_margin, _solve_flex_lp).
//
// It drives a not-yet-jammed champion to its OWN same-basin jammed optimum: each iteration
// re-linearizes the influence-sphere active set, solves the Donev flex LP for a center step
// inside a trust box, and uses the exact feasibility-shrunk sum (radii_lp + verify_and_score)
// as the trust-region control. Termination is on the box-independent Donev jamming margin.
// See docs/research/phase-2-csqv-jamming-optimizer.md and ADR 0003.
//
// This is the lever that closes the terminal gap on real champions (ticket 47: 100% of gap on
// N=142/143). It builds on the general LP (general_lp.hpp) for the flex LP; the tuned radii-LP
// (lp.hpp) supplies the exact objective. MONOTONE in the exact score: it keeps the best-scoring
// centers, so it never regresses; on a bad input it returns the input unchanged.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "general_lp.hpp"
#include "geometry.hpp"
#include "lp.hpp"

namespace csqv {

namespace detail {
enum { WALL_LEFT = 0, WALL_RIGHT = 1, WALL_BOTTOM = 2, WALL_TOP = 3 };
}  // namespace detail

// The active contacts and wall tangencies at a packing. `pairs[k] = {i, j}` is a contacting
// circle pair (surface gap below pair_tol); `walls[k] = {circle, wall_id}` a tangent wall
// (residual below wall_tol). Wall ids match terminal_squeeze.py: 0 left, 1 right, 2 bottom, 3 top.
struct ActiveSet {
  std::vector<std::array<int, 2>> pairs;
  std::vector<std::array<int, 2>> walls;
};

inline ActiveSet active_set(const std::vector<double>& x, const std::vector<double>& y,
                            const std::vector<double>& r, int n, double pair_tol, double wall_tol) {
  ActiveSet as;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = x[i] - x[j], dy = y[i] - y[j];
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d - (r[i] + r[j]) < pair_tol) as.pairs.push_back({i, j});
    }
  }
  for (int i = 0; i < n; ++i) {
    const double res[4] = {x[i] - r[i], SIDE - x[i] - r[i], y[i] - r[i], SIDE - y[i] - r[i]};
    for (int w = 0; w < 4; ++w)
      if (res[w] < wall_tol) as.walls.push_back({i, w});
  }
  return as;
}

// The active-constraint values g_a (M,), rows = pairs then walls. Zero at a jammed contact.
inline std::vector<double> active_residual(const std::vector<double>& x, const std::vector<double>& y,
                                           const std::vector<double>& r, const ActiveSet& as) {
  const int P = static_cast<int>(as.pairs.size());
  std::vector<double> g(P + as.walls.size());
  for (int k = 0; k < P; ++k) {
    const int i = as.pairs[k][0], j = as.pairs[k][1];
    g[k] = std::hypot(x[i] - x[j], y[i] - y[j]) - r[i] - r[j];
  }
  for (size_t w = 0; w < as.walls.size(); ++w) {
    const int i = as.walls[w][0], wid = as.walls[w][1];
    double v = 0.0;
    if (wid == detail::WALL_LEFT) v = x[i] - r[i];
    else if (wid == detail::WALL_RIGHT) v = SIDE - x[i] - r[i];
    else if (wid == detail::WALL_BOTTOM) v = y[i] - r[i];
    else v = SIDE - y[i] - r[i];
    g[P + w] = v;
  }
  return g;
}

// The dense active-constraint Jacobian A (M x 3N), rows = pairs then walls, columns
// [dx_0, dy_0, dr_0, dx_1, ...]. A pair row is grad(dist_ij - r_i - r_j); a wall row is the
// gradient of its tangent residual. Row-major, ready for general_lp.
inline std::vector<double> active_jacobian_dense(const std::vector<double>& x,
                                                 const std::vector<double>& y, const ActiveSet& as,
                                                 int n) {
  const int P = static_cast<int>(as.pairs.size());
  const int M = P + static_cast<int>(as.walls.size());
  const int ncol = 3 * n;
  std::vector<double> A(static_cast<size_t>(M) * ncol, 0.0);
  auto At = [&](int row, int col) -> double& { return A[static_cast<size_t>(row) * ncol + col]; };
  for (int k = 0; k < P; ++k) {
    const int i = as.pairs[k][0], j = as.pairs[k][1];
    const double dx = x[i] - x[j], dy = y[i] - y[j];
    double d = std::hypot(dx, dy);
    if (d <= 0.0) d = 1.0;
    const double ux = dx / d, uy = dy / d;
    At(k, 3 * i) = ux;
    At(k, 3 * i + 1) = uy;
    At(k, 3 * j) = -ux;
    At(k, 3 * j + 1) = -uy;
    At(k, 3 * i + 2) = -1.0;
    At(k, 3 * j + 2) = -1.0;
  }
  for (int w = 0; w < static_cast<int>(as.walls.size()); ++w) {
    const int row = P + w;
    const int i = as.walls[w][0], wid = as.walls[w][1];
    if (wid == detail::WALL_LEFT) {
      At(row, 3 * i) = 1.0;
      At(row, 3 * i + 2) = -1.0;
    } else if (wid == detail::WALL_RIGHT) {
      At(row, 3 * i) = -1.0;
      At(row, 3 * i + 2) = -1.0;
    } else if (wid == detail::WALL_BOTTOM) {
      At(row, 3 * i + 1) = 1.0;
      At(row, 3 * i + 2) = -1.0;
    } else {
      At(row, 3 * i + 1) = -1.0;
      At(row, 3 * i + 2) = -1.0;
    }
  }
  return A;
}

// The Donev LP jamming margin: the best first-order sum-of-radii gain from an unjamming flex
// (maximize sum(dr) s.t. A.[dz;dr] >= 0, -1 <= all <= 1). ~0 means the contact set is rigid
// (collectively jammed). Returns +inf if the LP fails (treated as not certified).
inline double donev_jamming_margin(const std::vector<double>& x, const std::vector<double>& y,
                                   const std::vector<double>& r, const ActiveSet& as, int n) {
  (void)r;  // kept for signature parity: the Donev flex depends on centers, not the radii
  const int M = static_cast<int>(as.pairs.size() + as.walls.size());
  if (M == 0) return std::numeric_limits<double>::infinity();
  const int ncol = 3 * n;
  const std::vector<double> A = active_jacobian_dense(x, y, as, n);
  std::vector<double> c(ncol, 0.0);
  for (int i = 0; i < n; ++i) c[3 * i + 2] = 1.0;  // maximize sum(dr)
  const std::vector<int> sense(M, +1);             // A.delta >= 0
  const std::vector<double> rhs(M, 0.0);
  const std::vector<double> lo(ncol, -1.0), hi(ncol, 1.0);
  const LpResult res = general_lp(ncol, M, c, A, sense, rhs, lo, hi);
  if (res.status != LpStatus::Optimal) return std::numeric_limits<double>::infinity();
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += res.x[3 * i + 2];
  return s;
}

// The exact-LP feasibility-shrunk score for fixed centers (terminal_squeeze.py::_exact_score).
inline double exact_score(const std::vector<double>& x, const std::vector<double>& y, int n) {
  const std::vector<double> r = radii_lp(x, y, n);
  double s;
  verify_and_score(x, y, r, n, TOL, s);
  return s;
}

struct JamResult {
  std::vector<double> x, y, r;  // jammed centers and their radii_lp radii
  double score = 0.0;           // exact feasibility-shrunk sum at the output
};

// Drive a champion to its same-basin jammed optimum. Parameters and control flow match
// terminal_squeeze.py::jam_slp. `x0,y0` are unit-frame centers; returns the best-scoring
// jammed centers, their radii_lp radii, and the exact score. Monotone: never regresses.
inline JamResult jam_slp(const std::vector<double>& x0, const std::vector<double>& y0, int n,
                         int max_iter = 200, double box0 = -1.0, double gamma_tol = 1e-2,
                         double hand_tol = 1e-7, double jam_lo = 1e-8, double jam_check = 1e-6,
                         double rho_hi = 0.75, double shrink = 0.5, double grow = 2.0,
                         double box_cap = 8.0, double box_floor = 1e-12) {
  JamResult out;
  std::vector<double> bx = x0, by = y0;
  // Guard: need finite input and n >= 2 (matches the Python early return).
  bool ok = (n >= 2);
  for (int i = 0; ok && i < n; ++i)
    if (!std::isfinite(bx[i]) || !std::isfinite(by[i])) ok = false;
  std::vector<double> br = ok ? radii_lp(bx, by, n) : std::vector<double>(n, 0.0);
  double best_score = ok ? exact_score(bx, by, n) : 0.0;
  if (!ok) {
    out.x = bx;
    out.y = by;
    out.r = br;
    out.score = best_score;
    return out;
  }

  double box = (box0 > 0.0) ? box0 : 0.02 / std::sqrt(static_cast<double>(n));
  const double cap = box * box_cap;  // relative to the initial box, as in terminal_squeeze.py
  const int ncol = 3 * n;

  for (int iter = 0; iter < max_iter; ++iter) {
    const ActiveSet as = active_set(bx, by, br, n, gamma_tol, gamma_tol);
    const int M = static_cast<int>(as.pairs.size() + as.walls.size());
    if (M == 0) break;

    // Flex LP: maximize sum(dr) s.t. A.[dz;dr] >= -g, active vars in [-box, box], rattlers pinned.
    std::vector<char> active(n, 0);
    for (const auto& p : as.pairs) { active[p[0]] = 1; active[p[1]] = 1; }
    for (const auto& w : as.walls) active[w[0]] = 1;
    const std::vector<double> A = active_jacobian_dense(bx, by, as, n);
    const std::vector<double> g = active_residual(bx, by, br, as);
    std::vector<double> c(ncol, 0.0), lo(ncol, 0.0), hi(ncol, 0.0), rhs(M);
    for (int i = 0; i < n; ++i) {
      c[3 * i + 2] = 1.0;
      const double b = active[i] ? box : 0.0;
      for (int t = 0; t < 3; ++t) { lo[3 * i + t] = -b; hi[3 * i + t] = b; }
    }
    for (int k = 0; k < M; ++k) rhs[k] = -g[k];
    const std::vector<int> sense(M, +1);
    const LpResult res = general_lp(ncol, M, c, A, sense, rhs, lo, hi);

    // On any non-Optimal flex LP (which cannot happen for a feasible packing: dz=dr=0 always
    // satisfies A.0 >= -g since g >= 0) this degrades to a zero step, so the trust box shrinks
    // and the loop still terminates, gated by the box-independent Donev margin below. This
    // matches terminal_squeeze.py::_solve_flex_lp, which returns a zero no-op on LP failure.
    std::vector<double> dz_x(n, 0.0), dz_y(n, 0.0), dr(n, 0.0);
    double pred = 0.0;
    if (res.status == LpStatus::Optimal) {
      for (int i = 0; i < n; ++i) {
        dz_x[i] = res.x[3 * i];
        dz_y[i] = res.x[3 * i + 1];
        dr[i] = res.x[3 * i + 2];
        pred += dr[i];
      }
    }

    // Termination on the box-INDEPENDENT margin, gated by the cheap box-scaled pred.
    if (pred <= jam_check) {
      const ActiveSet tas = active_set(bx, by, br, n, hand_tol, hand_tol);
      if (tas.pairs.size() + tas.walls.size() > 0 &&
          donev_jamming_margin(bx, by, br, tas, n) <= jam_lo)
        break;
    }

    std::vector<double> cx(n), cy(n);
    for (int i = 0; i < n; ++i) {
      cx[i] = std::clamp(bx[i] + dz_x[i], 0.0, SIDE);
      cy[i] = std::clamp(by[i] + dz_y[i], 0.0, SIDE);
    }
    const std::vector<double> r_try = radii_lp(cx, cy, n);
    double score_try;
    verify_and_score(cx, cy, r_try, n, TOL, score_try);

    if (std::isfinite(score_try) && score_try > best_score) {
      double actual = 0.0, step_max = 0.0;
      for (int i = 0; i < n; ++i) {
        if (active[i]) actual += r_try[i] - br[i];
        step_max = std::max(step_max, std::max(std::fabs(dz_x[i]), std::fabs(dz_y[i])));
        step_max = std::max(step_max, std::fabs(dr[i]));
      }
      const bool hit_box = step_max >= 0.99 * box;
      bx = cx; by = cy; br = r_try; best_score = score_try;
      if (hit_box && actual > rho_hi * pred) box = std::min(box * grow, cap);
    } else {
      box *= shrink;
    }
    if (box < box_floor) break;
  }

  out.x = bx;
  out.y = by;
  out.r = br;
  out.score = best_score;
  return out;
}

}  // namespace csqv

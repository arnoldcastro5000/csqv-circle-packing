// radii_lp: the exact optimal radii for FIXED centers (problem.py::radii_lp).
//
// The LP is: maximize sum(r_i) subject to containment (0 <= r_i <= wall_slack_i) and
// pairwise non-overlap (r_i + r_j <= dist_ij). Python solves this with HiGHS. Here a
// compact bounded-variable primal simplex solves it, and a LAZY constraint loop keeps
// the active set tiny: at the optimum only a few dozen pair constraints are ever tight
// (82 at N=121, 406 at N=484), so we start with box bounds only and add violated pairs
// until none remain. The result matches HiGHS to < 1e-9 on real data (see selftest).
//
// r = 0 is always feasible (distances are positive), so the simplex starts from a
// trivial feasible basis (all radii at their lower bound, slacks basic); no phase 1.
#pragma once
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "geometry.hpp"

namespace csqv {

namespace detail {

// Solve  max sum(z_i)  s.t.  z_a + z_b <= rhs  for each pair (a,b),  0 <= z_i <= u_i.
// Bounded-variable primal simplex on the standard form  A z + s = rhs,  s >= 0.
// Variables 0..k-1 are the structural radii; k..k+m-1 are the row slacks.
// Returns the optimal z (size k).
inline std::vector<double> solve_reduced_lp(int k, const std::vector<double>& u,
                                            const std::vector<std::array<int, 2>>& pairs,
                                            const std::vector<double>& rhs) {
  const int m = static_cast<int>(pairs.size());
  const int N = k + m;
  const double INF = std::numeric_limits<double>::infinity();
  const double eps = 1e-12;

  // Bounds and objective for every variable.
  std::vector<double> lo(N, 0.0), hi(N, 0.0), cost(N, 0.0);
  for (int i = 0; i < k; ++i) {
    hi[i] = u[i];
    cost[i] = 1.0;  // maximize sum of radii
  }
  for (int p = 0; p < m; ++p) hi[k + p] = INF;  // slacks

  // Dense tableau T (m x N) = B^{-1} A, and basic values xB. Start B = slack identity.
  std::vector<double> T(static_cast<size_t>(m) * N, 0.0);
  auto Tat = [&](int row, int col) -> double& { return T[static_cast<size_t>(row) * N + col]; };
  for (int p = 0; p < m; ++p) {
    Tat(p, pairs[p][0]) += 1.0;
    Tat(p, pairs[p][1]) += 1.0;
    Tat(p, k + p) = 1.0;
  }

  // Variable status: 0 nonbasic at lower, 1 nonbasic at upper, 2 basic.
  std::vector<int> status(N, 0);
  std::vector<int> basis(m);           // basis[row] = variable index
  std::vector<int> basisRow(N, -1);    // inverse map
  std::vector<double> xB(m);           // value of the basic variable in each row
  for (int p = 0; p < m; ++p) {
    basis[p] = k + p;
    basisRow[k + p] = p;
    status[k + p] = 2;
    xB[p] = rhs[p];  // structural vars at lower bound 0
  }

  auto nonbasic_value = [&](int j) { return status[j] == 1 ? hi[j] : lo[j]; };

  const int max_iter = 50 * (N + 10);
  for (int iter = 0; iter < max_iter; ++iter) {
    // Pricing: reduced cost d_j = cost_j - cost_B^T T_j. Pick the best improving var.
    // Bland's rule kicks in late to guarantee termination on degenerate problems.
    const bool bland = iter > 20 * (N + 10);
    int q = -1, qdir = 0;
    double best = eps;
    for (int j = 0; j < N; ++j) {
      if (status[j] == 2) continue;
      double dj = cost[j];
      for (int p = 0; p < m; ++p) {
        if (cost[basis[p]] != 0.0) dj -= cost[basis[p]] * Tat(p, j);
      }
      int dir = 0;
      if (status[j] == 0 && dj > eps)
        dir = +1;  // at lower, raising improves
      else if (status[j] == 1 && dj < -eps)
        dir = -1;  // at upper, lowering improves
      else
        continue;
      if (bland) {
        q = j;
        qdir = dir;
        break;
      }
      const double gain = std::fabs(dj);
      if (gain > best) {
        best = gain;
        q = j;
        qdir = dir;
      }
    }
    if (q < 0) break;  // optimal

    // Ratio test. Entering var changes by delta = qdir * t, t >= 0.
    // Basic var p changes: xB[p] -= T(p,q) * delta.
    double t = INF;
    int leave_row = -1;
    int leave_to = 0;  // bound the leaving var lands on: 0 lower, 1 upper
    // 1) entering variable hitting its own opposite bound (a bound flip).
    if (hi[q] < INF) {
      const double span = hi[q] - lo[q];
      if (span < t) {
        t = span;
        leave_row = -1;  // flip, no basis change
      }
    }
    // 2) basic variables hitting a bound.
    for (int p = 0; p < m; ++p) {
      const double coef = Tat(p, q) * qdir;  // xB[p] decreases at rate coef as t grows
      if (coef > eps) {                      // heading toward lower bound
        const double room = xB[p] - lo[basis[p]];
        const double ratio = room / coef;
        if (ratio < t - 1e-15) {
          t = ratio;
          leave_row = p;
          leave_to = 0;
        }
      } else if (coef < -eps) {  // heading toward upper bound
        if (hi[basis[p]] < INF) {
          const double room = xB[p] - hi[basis[p]];  // negative
          const double ratio = room / coef;          // positive
          if (ratio < t - 1e-15) {
            t = ratio;
            leave_row = p;
            leave_to = 1;
          }
        }
      }
    }

    if (!std::isfinite(t)) break;  // unbounded (should not happen: bounded by u)
    if (t < 0.0) t = 0.0;
    const double delta = qdir * t;

    // Update basic values for the step.
    for (int p = 0; p < m; ++p) xB[p] -= Tat(p, q) * delta;

    if (leave_row < 0) {
      // Bound flip: q moves to its opposite bound, stays nonbasic.
      status[q] = (status[q] == 0) ? 1 : 0;
      continue;
    }

    // Pivot: q enters the basis in leave_row; the old basic var leaves.
    const int leaving = basis[leave_row];
    const double piv = Tat(leave_row, q);
    for (int j = 0; j < N; ++j) Tat(leave_row, j) /= piv;
    const double entering_val = nonbasic_value(q) + delta;
    xB[leave_row] = entering_val;
    for (int p = 0; p < m; ++p) {
      if (p == leave_row) continue;
      const double f = Tat(p, q);
      if (f == 0.0) continue;
      for (int j = 0; j < N; ++j) Tat(p, j) -= f * Tat(leave_row, j);
      xB[p] -= f * 0.0;  // xB already updated by the step above
    }
    basisRow[leaving] = -1;
    status[leaving] = (leave_to == 1) ? 1 : 0;
    basis[leave_row] = q;
    basisRow[q] = leave_row;
    status[q] = 2;
  }

  std::vector<double> z(k, 0.0);
  for (int i = 0; i < k; ++i) {
    if (status[i] == 2)
      z[i] = xB[basisRow[i]];
    else
      z[i] = nonbasic_value(i);
    if (z[i] < 0.0) z[i] = 0.0;
    if (z[i] > u[i]) z[i] = u[i];
  }
  return z;
}

}  // namespace detail

// Exact optimal radii for fixed centers. Matches problem.py::radii_lp.
inline std::vector<double> radii_lp(const std::vector<double>& x, const std::vector<double>& y,
                                    int n) {
  std::vector<double> u(n);
  for (int i = 0; i < n; ++i) u[i] = std::max(wall_slack(x[i], y[i]), 0.0);
  if (n < 2) return u;

  // Prunable pairs: a constraint r_i + r_j <= d can bind only if d < u_i + u_j.
  struct Pr {
    int i, j;
    double d;
  };
  std::vector<Pr> prun;
  prun.reserve(static_cast<size_t>(n) * 4);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = x[i] - x[j], dy = y[i] - y[j];
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d < u[i] + u[j]) prun.push_back({i, j, d});
    }
  }

  std::vector<double> r = u;  // start: every circle at its containment bound
  if (prun.empty()) return r;

  // Lazy active-set loop. Starting from r = u nearly every prunable pair overlaps, so
  // add only the most-violated batch per pass and PURGE clearly-loose pairs after each
  // solve. A deadband (add when overlap > viol_tol, purge when slack > purge_tol) keeps
  // the active set near the tight count (~82 at N=121, ~406 at N=484) and the tableau
  // small, so this scales to large N in low memory.
  const double viol_tol = 1e-11;
  const double purge_tol = 1e-7;
  const size_t batch = static_cast<size_t>(3 * n + 64);
  std::vector<char> active(prun.size(), 0);
  std::vector<int> involved_local(n, -1);
  std::vector<int> involved_global;

  for (int pass = 0; pass < 5 * n + 50; ++pass) {
    // Collect non-active pairs the current radii overlap; add the most-violated batch.
    std::vector<std::pair<double, size_t>> cand;
    for (size_t p = 0; p < prun.size(); ++p) {
      if (active[p]) continue;
      const double viol = r[prun[p].i] + r[prun[p].j] - prun[p].d;
      if (viol > viol_tol) cand.emplace_back(viol, p);
    }
    if (cand.empty()) break;  // converged: all non-active satisfied, active enforced
    if (cand.size() > batch) {
      std::nth_element(cand.begin(), cand.begin() + batch, cand.end(),
                       [](const auto& a, const auto& b) { return a.first > b.first; });
      cand.resize(batch);
    }
    for (const auto& c : cand) active[c.second] = 1;

    // Build the reduced LP over circles that appear in an active pair.
    std::fill(involved_local.begin(), involved_local.end(), -1);
    involved_global.clear();
    std::vector<std::array<int, 2>> apairs;
    std::vector<double> arhs;
    for (size_t p = 0; p < prun.size(); ++p) {
      if (!active[p]) continue;
      for (int side = 0; side < 2; ++side) {
        const int g = side == 0 ? prun[p].i : prun[p].j;
        if (involved_local[g] < 0) {
          involved_local[g] = static_cast<int>(involved_global.size());
          involved_global.push_back(g);
        }
      }
      apairs.push_back({involved_local[prun[p].i], involved_local[prun[p].j]});
      arhs.push_back(prun[p].d);
    }
    const int k = static_cast<int>(involved_global.size());
    std::vector<double> ru(k);
    for (int t = 0; t < k; ++t) ru[t] = u[involved_global[t]];

    const std::vector<double> z = detail::solve_reduced_lp(k, ru, apairs, arhs);

    // Non-involved circles keep r = u; involved take the LP value.
    for (int i = 0; i < n; ++i) r[i] = u[i];
    for (int t = 0; t < k; ++t) r[involved_global[t]] = z[t];

    // Purge clearly-loose active pairs to bound the active-set size. The deadband above
    // viol_tol and below purge_tol keeps near-tight pairs, so this cannot oscillate.
    for (size_t p = 0; p < prun.size(); ++p) {
      if (active[p] && prun[p].d - (r[prun[p].i] + r[prun[p].j]) > purge_tol) active[p] = 0;
    }
  }
  return r;
}

}  // namespace csqv

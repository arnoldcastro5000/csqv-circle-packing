// general_lp: a self-contained general linear program solver for the terminal squeeze.
//
// Maximize c.x subject to, for each row i, (A_i . x) sense_i rhs_i, with lo_j <= x_j <= hi_j.
// sense: -1 means '<=', 0 means '==', +1 means '>='. Bounds may be infinite.
//
// The tuned radii-LP simplex in lp.hpp is hard-specialized (fixed cost 1, rows of the form
// z_a + z_b <= rhs, structural lower bound 0, no phase 1), so it cannot express the flex LP the
// jamming optimizer needs (arbitrary rows, negative bounds, a general objective). This is the
// separate general solver that jam_slp and the Donev jamming certificate build on; the radii-LP
// stays untouched so the hot path never regresses (ADR 0003).
//
// Method: a two-phase bounded-variable primal simplex. Each structural variable is shifted to a
// nonnegative variable at its finite bound (a fully free variable is split into two nonnegative
// parts), inequalities get a slack or surplus, and every row gets an artificial. Phase 1 drives
// the artificials out (or proves the LP infeasible); phase 2 optimizes the real objective (or
// proves it unbounded). The pivot, ratio test, and Bland anti-cycling mirror lp.hpp.
#pragma once
#include <cmath>
#include <limits>
#include <vector>

namespace csqv {

enum class LpStatus { Optimal, Infeasible, Unbounded };

struct LpResult {
  LpStatus status = LpStatus::Infeasible;
  std::vector<double> x;    // size n; meaningful only when status == Optimal
  double objective = 0.0;   // c.x at the optimum
};

namespace detail {

// Bounded-variable primal simplex core on the standard form  T x = xB (as a basis view),
// 0 <= x_j <= hi_j, maximize cost.x. `T` is the current tableau B^{-1}A (m x N, row-major),
// `xB` the basic values, `basis`/`basisRow`/`status` the basis bookkeeping (status: 0 nonbasic
// at lower 0, 1 nonbasic at upper hi, 2 basic). Returns 0 if it reached an optimum, 1 if the LP
// is unbounded in the current objective, 2 if it hit the iteration cap without converging (the
// caller must NOT trust the basis in that case). All lower bounds are 0 (the caller shifts).
inline int glp_simplex_core(int N, int m, std::vector<double>& T, std::vector<double>& xB,
                            std::vector<int>& basis, std::vector<int>& basisRow,
                            std::vector<int>& status, const std::vector<double>& cost,
                            const std::vector<double>& hi) {
  const double INF = std::numeric_limits<double>::infinity();
  const double eps = 1e-12;
  auto Tat = [&](int row, int col) -> double& { return T[static_cast<size_t>(row) * N + col]; };
  auto nonbasic_value = [&](int j) { return status[j] == 1 ? hi[j] : 0.0; };

  // Bland's rule guarantees termination, so the cap is only a pathology backstop; keep it
  // generous and switch to Bland early enough that the anti-cycling exit fires first.
  const int max_iter = 200 * (N + 10);
  for (int iter = 0; iter < max_iter; ++iter) {
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
        dir = +1;
      else if (status[j] == 1 && dj < -eps)
        dir = -1;
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
    if (q < 0) return 0;  // optimal

    double t = INF;
    int leave_row = -1;
    int leave_to = 0;
    if (hi[q] < INF) {  // entering var flips to its own opposite bound
      const double span = hi[q];  // lo is 0
      if (span < t) {
        t = span;
        leave_row = -1;
      }
    }
    for (int p = 0; p < m; ++p) {
      const double coef = Tat(p, q) * qdir;  // xB[p] decreases at rate coef as t grows
      if (coef > eps) {                       // toward lower bound 0
        const double room = xB[p] - 0.0;
        const double ratio = room / coef;
        if (ratio < t - 1e-15) {
          t = ratio;
          leave_row = p;
          leave_to = 0;
        }
      } else if (coef < -eps) {  // toward upper bound
        if (hi[basis[p]] < INF) {
          const double room = xB[p] - hi[basis[p]];  // <= 0
          const double ratio = room / coef;          // >= 0
          if (ratio < t - 1e-15) {
            t = ratio;
            leave_row = p;
            leave_to = 1;
          }
        }
      }
    }

    if (!std::isfinite(t)) return 1;  // unbounded
    if (t < 0.0) t = 0.0;
    const double delta = qdir * t;

    for (int p = 0; p < m; ++p) xB[p] -= Tat(p, q) * delta;

    if (leave_row < 0) {  // bound flip, q stays nonbasic
      status[q] = (status[q] == 0) ? 1 : 0;
      continue;
    }

    const int leaving = basis[leave_row];
    const double piv = Tat(leave_row, q);
    for (int j = 0; j < N; ++j) Tat(leave_row, j) /= piv;
    xB[leave_row] = nonbasic_value(q) + delta;
    for (int p = 0; p < m; ++p) {
      if (p == leave_row) continue;
      const double f = Tat(p, q);
      if (f == 0.0) continue;
      for (int j = 0; j < N; ++j) Tat(p, j) -= f * Tat(leave_row, j);
    }
    basisRow[leaving] = -1;
    status[leaving] = (leave_to == 1) ? 1 : 0;
    basis[leave_row] = q;
    basisRow[q] = leave_row;
    status[q] = 2;
  }
  return 2;  // iteration cap hit without an optimality/unbounded exit: do not trust the basis
}

}  // namespace detail

inline LpResult general_lp(int n, int m, const std::vector<double>& c, const std::vector<double>& A,
                           const std::vector<int>& sense, const std::vector<double>& rhs,
                           const std::vector<double>& lo, const std::vector<double>& hi) {
  const double INF = std::numeric_limits<double>::infinity();
  LpResult out;

  // 1) Shift each structural variable to a nonnegative variable (or two, if fully free).
  //    x_j = base_j + col1_sign * y[col1]  [ - y[col2] when free ]
  //    kind: 0 => x = lo + y (upper y = hi-lo);  1 => x = hi - y (y in [0, inf));
  //          2 => free, x = y+ - y- (both in [0, inf)).
  std::vector<int> kind(n), col1(n), col2(n, -1);
  std::vector<double> ycost;   // objective per y column (maximize)
  std::vector<double> yhi;     // upper bound per y column (0 lower)
  auto add_y = [&](double cost_j, double upper) {
    ycost.push_back(cost_j);
    yhi.push_back(upper);
    return static_cast<int>(ycost.size()) - 1;
  };
  std::vector<double> base(n, 0.0);
  for (int j = 0; j < n; ++j) {
    const bool loFin = std::isfinite(lo[j]);
    const bool hiFin = std::isfinite(hi[j]);
    if (loFin) {
      kind[j] = 0;
      base[j] = lo[j];
      col1[j] = add_y(c[j], hiFin ? (hi[j] - lo[j]) : INF);
    } else if (hiFin) {
      kind[j] = 1;
      base[j] = hi[j];
      col1[j] = add_y(-c[j], INF);
    } else {
      kind[j] = 2;
      base[j] = 0.0;
      col1[j] = add_y(c[j], INF);
      col2[j] = add_y(-c[j], INF);
    }
  }

  // 2) Build the equality rows over the y columns, folding the shift constants into the rhs.
  //    Then add one slack/surplus per inequality row.
  std::vector<double> beq(m);
  // Accumulate each row's structural + slack coefficients in a compact per-row list.
  std::vector<std::vector<std::pair<int, double>>> sparse(m);
  for (int i = 0; i < m; ++i) {
    double b = rhs[i];
    for (int j = 0; j < n; ++j) {
      const double aij = A[static_cast<size_t>(i) * n + j];
      if (aij == 0.0) continue;
      b -= aij * base[j];
      if (kind[j] == 0) {
        sparse[i].emplace_back(col1[j], aij);
      } else if (kind[j] == 1) {
        sparse[i].emplace_back(col1[j], -aij);
      } else {
        sparse[i].emplace_back(col1[j], aij);
        sparse[i].emplace_back(col2[j], -aij);
      }
    }
    beq[i] = b;
  }
  // Slacks: '<=' gets +slack, '>=' gets -slack, '==' none. All slacks are y columns in [0, inf).
  for (int i = 0; i < m; ++i) {
    if (sense[i] == -1) {
      const int s = add_y(0.0, INF);
      sparse[i].emplace_back(s, 1.0);
    } else if (sense[i] == +1) {
      const int s = add_y(0.0, INF);
      sparse[i].emplace_back(s, -1.0);
    }
  }

  const int nY = static_cast<int>(ycost.size());  // structural + slack columns
  const int N = nY + m;                            // + one artificial per row
  // 3) Assemble the dense tableau. Artificials occupy columns nY..nY+m-1 as an identity block.
  std::vector<double> T(static_cast<size_t>(m) * N, 0.0);
  auto Tat = [&](int r, int col) -> double& { return T[static_cast<size_t>(r) * N + col]; };
  std::vector<double> xB(m);
  for (int i = 0; i < m; ++i) {
    double sign = 1.0;
    if (beq[i] < 0.0) sign = -1.0;  // keep the artificial's basic value nonnegative
    for (const auto& [col, val] : sparse[i]) Tat(i, col) += sign * val;
    Tat(i, nY + i) = 1.0;  // artificial
    xB[i] = sign * beq[i];
  }

  std::vector<int> status(N, 0), basis(m), basisRow(N, -1);
  for (int i = 0; i < m; ++i) {
    basis[i] = nY + i;
    basisRow[nY + i] = i;
    status[nY + i] = 2;
  }

  std::vector<double> hiA(N, INF);
  for (int j = 0; j < nY; ++j) hiA[j] = yhi[j];
  // artificials keep upper INF during phase 1.

  // 4) Phase 1: minimize the sum of artificials (maximize its negation).
  std::vector<double> cost1(N, 0.0);
  for (int i = 0; i < m; ++i) cost1[nY + i] = -1.0;
  const int rc1 = detail::glp_simplex_core(N, m, T, xB, basis, basisRow, status, cost1, hiA);

  // Declare Infeasible only from a CONVERGED phase 1 (rc1 == 0): if the cap was hit (rc1 == 2)
  // the artificial sum is not certified minimal, so a false Infeasible must not be reported.
  double art_sum = 0.0;
  for (int i = 0; i < m; ++i)
    if (basis[i] >= nY) art_sum += std::fabs(xB[i]);
  if (rc1 == 0 && art_sum > 1e-7) {
    out.status = LpStatus::Infeasible;
    return out;
  }

  // 5) Freeze artificials at 0 and optimize the real objective. A degenerate artificial may
  // remain basic at a value in (0, 1e-7]; snap it to exactly 0 first so that after hi -> 0 no
  // basic value sits above its upper bound (which would give a negative phase-2 ratio).
  for (int i = 0; i < m; ++i)
    if (basis[i] >= nY) xB[i] = 0.0;
  for (int i = 0; i < m; ++i) hiA[nY + i] = 0.0;
  std::vector<double> cost2(N, 0.0);
  for (int j = 0; j < nY; ++j) cost2[j] = ycost[j];
  const int rc = detail::glp_simplex_core(N, m, T, xB, basis, basisRow, status, cost2, hiA);
  if (rc == 1) {
    out.status = LpStatus::Unbounded;
    return out;
  }

  // 6) Recover the y values, then map back to the structural x and score the objective.
  std::vector<double> y(nY, 0.0);
  for (int j = 0; j < nY; ++j) {
    if (status[j] == 2)
      y[j] = xB[basisRow[j]];
    else
      y[j] = (status[j] == 1) ? hiA[j] : 0.0;
    if (y[j] < 0.0) y[j] = 0.0;
  }
  out.x.assign(n, 0.0);
  for (int j = 0; j < n; ++j) {
    if (kind[j] == 0)
      out.x[j] = base[j] + y[col1[j]];
    else if (kind[j] == 1)
      out.x[j] = base[j] - y[col1[j]];
    else
      out.x[j] = y[col1[j]] - y[col2[j]];
  }
  double obj = 0.0;
  for (int j = 0; j < n; ++j) obj += c[j] * out.x[j];
  out.objective = obj;
  out.status = LpStatus::Optimal;
  return out;
}

}  // namespace csqv

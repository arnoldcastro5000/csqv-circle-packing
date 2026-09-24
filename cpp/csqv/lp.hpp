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
//
// Pricing is INCREMENTAL. A basis of this LP holds only edge-vertex incidence columns and
// slack columns. So the tableau B^-1 [A I] is half-integral: each entry is a multiple of
// 1/2 (Nemhauser-Trotter 1975; Balinski 1965). A measurement of 1e10 entries at N=90
// and N=121 found |entry| <= 2. So each pivot is +-1/2,
// +-1 or +-2, and each tableau and reduced-cost operation is exact in IEEE-754.
// Because of this, each pivot updates the reduced-cost row (d_j -= d_q * T(leave, j)). A
// full recompute from the tableau is not necessary. The update gives the same bits as the
// recompute (both start from cost and subtract; see "Signed zeros" below), so the radii
// and the duals stay bit-identical.
// A guard checks each pivot and each updated d_j. It does not check every tableau entry.
// If the guard fails, the solve uses the full recompute (the original pricing) from that
// pivot to its end.
//
// Signed zeros. The pivot row is sparse, so a basis change divides and eliminates only over
// its nonzero columns. A dense pass also computes x - f * (+-0) at each zero column. That
// changes x only if x is an exact zero, and then only its sign. So the sparse pass leaves
// other zero signs in the tableau T than the dense pass. These signs cannot reach the
// output, because d and xB never hold -0:
// - In round-to-nearest, x - y is -0 only if x is -0 (and y is +0), and x + y is -0 only
//   if x and y are both -0. A subtraction from a value that is not -0 never gives -0.
// - d starts as cost (+0 or 1) and changes only by subtraction. xB starts as rhs + 0.0 and
//   changes only by subtraction, or takes a bound plus the step: lo = +0 or hi = u + 0.0,
//   so neither is -0. The dual reduced costs start from cost and subtract too, and each
//   pair dual is clamped to at least +0.
// Two values that differ only in the sign of zero compare equal, so the control flow is the
// same. Each division has a nonzero divisor (the pivot, the ratio-test coefficient), so a
// zero sign never becomes the sign of an infinity. Thus the radii, the duals and d are the
// same bits as with the dense pass, and none of them is -0.
//
// Column q. The entering column q of the tableau is sparse too (about 5-8% nonzero at
// N=90..143). The tableau is row-major, so a pass down column q reads one double per row
// with a stride of N. Each iteration reads column q once and keeps its nonzero rows, in
// increasing row order. The ratio test, the xB step and the elimination use only this list:
// - The ratio test skips a zero coefficient (it is neither > eps nor < -eps), and the list
//   keeps the row order, so the leaving row and its tie-break are the same.
// - At a zero row, the dense xB step computes xB[p] - (+-0). This is xB[p], because xB
//   never holds -0 (see above) and the step is finite.
// - The division of the pivot row changes only the leaving row, and the elimination skips
//   that row. So the listed values of the other rows are still column q at the elimination.
//
// Entering pick. The pricing scan reads the status and d_j of all N columns in each
// iteration. But a basis change updates d only at the nonzero columns of the pivot row, and
// only q and the leaving variable change status (only q at a bound flip). So the incremental
// pricing keeps one key per column: the gain |d_j| of an eligible column, else 0. Each
// iteration updates the keys of the changed columns only. The pick is the largest key, the
// smallest j on a tie. This is the column of the scan, because the scan also takes the
// largest gain above eps and changes its pick only for a strictly larger gain.
// A key holds the bits of the gain as an int64. For doubles >= +0, the int64 order is the
// order of the values, and two equal values have equal bits: std::fabs never gives -0, and d
// holds no NaN (the guard rejects NaN). The compiler vectorizes an int64 max. It does not
// vectorize a double max without fast-math, because that max must keep the order of a NaN.
// The Bland phase, the full recompute and the CSQV_LP_VERIFY_PRICING build use the scan.
//
// Pivot row. A basis change divides the nonzeros of the leaving row and keeps their columns in
// nz, in increasing order. A dense pass over the row reads all N columns, but the row has only
// about 7-8 nonzeros. So each tableau row keeps a block mask: one bit per block of
// kRowBlock columns. Each nonzero of the row is in a marked block. A marked block can hold only
// zeros. The mask of the leaving row gives nz from its marked blocks only:
// - Each row starts with 3 nonzeros (its two radii and its slack), and its mask marks their
//   blocks.
// - The division changes only the nonzeros of the leaving row. The elimination changes a
//   row only at the columns of nz, and each of them is in a marked block of the leaving
//   row. So the elimination adds the mask of the leaving row to the mask of each changed
//   row. A bound flip does not change the tableau.
// - The pass reads the marked blocks in increasing order and each block in increasing
//   order, and it keeps a column if its entry is != 0 (as the dense pass does). So nz is the
//   same list as from the dense pass. After the pass, the mask of the leaving row marks
//   only the blocks that hold a nonzero.
//
// Tableau buffer. The dense tableau is large (up to about 2.5 MB at N=121), and the worker
// solves more than a thousand LPs per second. A new heap block per solve costs a page fault
// per page on a heap that returns large freed blocks to the OS: the Windows heap does, glibc
// does not by default. So each thread keeps one tableau buffer, and each solve resets it to
// zeros with assign. The values are the same as in a new block. The buffer keeps the largest
// size that its thread used (about 44 MB at N=484). The buffer is thread_local. MinGW
// emulates thread-local storage with a function call per access, so a solve reads the buffer
// once. A solve must not start another solve (none does), because the two would share the
// buffer.
#pragma once
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "geometry.hpp"
#include "stage_timer.hpp"

namespace csqv {

// The number of solves where the half-integral guard failed and the pricing changed to the
// full recompute. For this LP it must stay 0. A count above 0 means the LP form changed.
inline std::atomic<long>& lp_pricing_fallbacks() {
  static std::atomic<long> count{0};
  return count;
}

namespace detail {

// Recompute is the original pricing: a full dot product per column per iteration. It is
// the guard's fallback, and the unit test compares the incremental pricing against it.
enum class Pricing { Incremental, Recompute };

// Test hook for the fallback path. If K > 0, the guard fails at the K-th basis change of
// each solve, so the solve continues with the full recompute from that pivot. The result
// must stay bit-identical. lp_test.cpp passes K directly; the bench fallback build sets the
// default with -DCSQV_LP_FORCE_FALLBACK_AT=K. 0 (the default) disables the hook.
#ifndef CSQV_LP_FORCE_FALLBACK_AT
#define CSQV_LP_FORCE_FALLBACK_AT 0
#endif

#ifdef CSQV_LP_VERIFY_PRICING
// Bench-only check (cpp/bench/pricing_check.cpp). At each pricing, the solve also computes
// the full recompute of each incremental d_j and counts the compares and the mismatches.
// It is slow; never define CSQV_LP_VERIFY_PRICING in a worker build.
inline std::atomic<long>& lp_pricing_compares() {
  static std::atomic<long> count{0};
  return count;
}
inline std::atomic<long>& lp_pricing_mismatches() {
  static std::atomic<long> count{0};
  return count;
}
// The same check for the entering pick: at each pricing that can use the keys, the solve
// also computes the pick from the keys and counts the picks and the mismatches with the scan.
inline std::atomic<long>& lp_pricing_picks() {
  static std::atomic<long> count{0};
  return count;
}
inline std::atomic<long>& lp_pricing_pick_mismatches() {
  static std::atomic<long> count{0};
  return count;
}
#endif

// True if v is an exact multiple of 1/2 and |2v| < 2^51. NaN and infinity give false.
// The test rounds 2v to an integer: adding and then subtracting 1.5 * 2^52 does this in
// exact IEEE-754 arithmetic. If the rounding changes nothing, 2v is an integer.
// The test has no branch, so the update loop that calls it vectorizes. (std::floor is a
// library call at the x86-64 baseline target, CSQV_ARCH=; a guard with it took about 3x
// more time. Branch-free, the test also stays fast in a build for that target.)
inline bool is_half_integral(double v) {
  const double h = 2.0 * v;
  const double kRound = 6755399441055744.0;              // 1.5 * 2^52
  const bool small = std::fabs(h) < 2251799813685248.0;  // 2^51
  const bool whole = (h + kRound) - kRound == h;
  return small & whole;  // bitwise, not &&: a branch would stop the vectorization
}

// True if a pivot keeps the tableau exact: division by +-1/2, +-1 or +-2 is exact.
inline bool is_half_integral_pivot(double piv) {
  const double a = std::fabs(piv);
  return a == 0.5 || a == 1.0 || a == 2.0;
}

// The incremental pricing step for a basis change: d_j -= d_q * row_j for each j in nz, the
// nonzero columns of row. row is the pivot row after the division by piv. A zero column
// changes nothing, because d never holds -0 (see "Signed zeros" above). The guard checks
// only the updated d_j; each other d_j passed it before. Returns false if the half-integral
// guard fails (piv or an updated d_j). Then the caller must stop the use of d.
inline bool update_reduced_costs(std::vector<double>& d, const double* row,
                                 const std::vector<int>& nz, int q, double piv) {
  const double dq = d[q];
  int inexact = 0;
  for (const int j : nz) {
    d[j] -= dq * row[j];
    inexact |= !is_half_integral(d[j]);
  }
  return !inexact && is_half_integral_pivot(piv);
}

// The pricing key of a column (see "Entering pick" above): the gain |d_j| as int64 bits if
// the column is eligible, else 0. A column is eligible if a move off its bound improves the
// objective by more than eps: d_j > eps at the lower bound (status 0), d_j < -eps at the
// upper bound (status 1). A basic column (status 2) is not eligible.
inline int64_t pricing_key(int status, double dj, double eps) {
  const bool eligible = (status == 0 && dj > eps) || (status == 1 && dj < -eps);
  const double gain = eligible ? std::fabs(dj) : 0.0;
  int64_t key;
  std::memcpy(&key, &gain, sizeof key);
  return key;
}

// The entering column from the keys: the largest key, the smallest j on a tie. Returns -1 if
// all keys are 0 (no column is eligible). The first pass finds the largest key, the second
// pass finds its first column. Both loops vectorize.
inline int pick_entering(const int64_t* key, int n) {
  int64_t top = 0;
  for (int j = 0; j < n; ++j) top = key[j] > top ? key[j] : top;
  if (top == 0) return -1;
  for (int j = 0; j < n; ++j)
    if (key[j] == top) return j;
  return -1;  // not reached: some key equals top
}

// The block masks of the tableau rows (see "Pivot row" above). Bit b of word w marks the
// columns (64 * w + b) * kRowBlock to that + kRowBlock - 1.
inline constexpr int kRowBlock = 8;

// The number of mask words per row for n columns.
inline int row_mask_words(int n) { return (n + 64 * kRowBlock - 1) / (64 * kRowBlock); }

// Marks the block of column j in a row mask.
inline void mark_column(uint64_t* mask, int j) {
  mask[j / (64 * kRowBlock)] |= uint64_t{1} << (j / kRowBlock % 64);
}

// The index of the lowest set bit of v (v != 0).
inline int lowest_bit(uint64_t v) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(v);
#else
  int b = 0;
  for (; (v & 1) == 0; v >>= 1) ++b;
  return b;
#endif
}

// The nonzero columns of row (n columns), in increasing order, into nz. The pass reads only
// the blocks that mask marks, so each nonzero must be in a marked block. After the pass, mask
// marks only the blocks that hold a nonzero. buf is scratch space for n columns: each column of
// a block is written, and the count moves on only at a nonzero, so the loop has no branch.
inline void row_nonzeros(const double* row, int n, uint64_t* mask, int* buf,
                         std::vector<int>& nz) {
  int c = 0;
  for (int w = 0; w < row_mask_words(n); ++w) {
    uint64_t bits = mask[w], exact = 0;
    while (bits != 0) {
      const int b = lowest_bit(bits);
      bits &= bits - 1;
      const int j0 = (64 * w + b) * kRowBlock;
      const int j1 = j0 + kRowBlock < n ? j0 + kRowBlock : n;
      const int before = c;
      for (int j = j0; j < j1; ++j) {
        buf[c] = j;
        c += row[j] != 0.0;
      }
      exact |= static_cast<uint64_t>(c != before) << b;
    }
    mask[w] = exact;
  }
  nz.assign(buf, buf + c);
}

// The tableau buffer of the calling thread (see "Tableau buffer" above).
inline std::vector<double>& tableau_buffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}

// Solve  max sum(z_i)  s.t.  z_a + z_b <= rhs  for each pair (a,b),  0 <= z_i <= u_i.
// Bounded-variable primal simplex on the standard form  A z + s = rhs,  s >= 0.
// Variables 0..k-1 are the structural radii; k..k+m-1 are the row slacks.
// Returns the optimal z (size k).
//
// Optional DUALS for the envelope-theorem center gradient (ticket 40). When `pair_dual` is
// non-null it is filled (size m) with the shadow price y_p >= 0 of each pair constraint:
// y_p = d(optimal sum)/d(rhs_p), which for a max LP is minus the reduced cost of that row's
// slack. When `struct_rc` is non-null it is filled (size k) with the reduced cost of each
// structural radius; for a radius pinned at its upper bound u_i this is the wall shadow price
// w_i = d(optimal sum)/d(u_i) >= 0, and it is ~0 for an interior (basic) radius. Both come
// from the final optimal basis, so they are exact for the LP the caller passed.
inline std::vector<double> solve_reduced_lp(int k, const std::vector<double>& u,
                                            const std::vector<std::array<int, 2>>& pairs,
                                            const std::vector<double>& rhs,
                                            std::vector<double>* pair_dual = nullptr,
                                            std::vector<double>* struct_rc = nullptr,
                                            Pricing pricing = Pricing::Incremental,
                                            int force_fallback_at = CSQV_LP_FORCE_FALLBACK_AT) {
  const int m = static_cast<int>(pairs.size());
  const int N = k + m;
  const double INF = std::numeric_limits<double>::infinity();
  const double eps = 1e-12;

  // Bounds and objective for every variable.
  std::vector<double> lo(N, 0.0), hi(N, 0.0), cost(N, 0.0);
  for (int i = 0; i < k; ++i) {
    hi[i] = u[i] + 0.0;  // + 0.0 changes a -0 bound to +0 (see "Signed zeros" above)
    cost[i] = 1.0;       // maximize sum of radii
  }
  for (int p = 0; p < m; ++p) hi[k + p] = INF;  // slacks

  // Dense tableau T (m x N) = B^{-1} A, and basic values xB. Start B = slack identity.
  std::vector<double>& T = tableau_buffer();
  T.assign(static_cast<size_t>(m) * N, 0.0);
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
    xB[p] = rhs[p] + 0.0;  // structural vars at lower bound 0; + 0.0 as for hi
  }

  auto nonbasic_value = [&](int j) { return status[j] == 1 ? hi[j] : lo[j]; };

  // Reduced costs d_j = cost_j - cost_B^T T_j. The slack basis has cost_B = 0, so d = cost.
  bool incremental = pricing == Pricing::Incremental;
  std::vector<double> d = cost;
  // The pricing keys (see "Entering pick" above). Only the incremental pricing uses them.
  std::vector<int64_t> key(N, 0);
  if (incremental)
    for (int j = 0; j < N; ++j) key[j] = pricing_key(status[j], d[j], eps);
  std::vector<int> nz;  // the nonzero columns of the current pivot row
  nz.reserve(N);
  std::vector<int> nz_buf(N);  // scratch space for row_nonzeros
  // The block mask of each tableau row (see "Pivot row" above): W words per row.
  const int W = row_mask_words(N);
  std::vector<uint64_t> row_mask(static_cast<size_t>(m) * W, 0);
  for (int p = 0; p < m; ++p) {
    uint64_t* const mask = &row_mask[static_cast<size_t>(p) * W];
    mark_column(mask, pairs[p][0]);
    mark_column(mask, pairs[p][1]);
    mark_column(mask, k + p);
  }
  struct ColEntry {
    int row;
    double val;
  };
  std::vector<ColEntry> col;  // the nonzero rows of column q (see "Column q" above)
  col.reserve(m);
  int basis_changes = 0;  // counts pivots for the force_fallback_at test hook
#ifdef CSQV_LP_VERIFY_PRICING
  long compares = 0, mismatches = 0, picks = 0, pick_mismatches = 0;
#endif

  const int max_iter = 50 * (N + 10);
  for (int iter = 0; iter < max_iter; ++iter) {
    // Pricing: reduced cost d_j = cost_j - cost_B^T T_j. Pick the best improving var.
    // Bland's rule kicks in late to guarantee termination on degenerate problems.
    const bool bland = iter > 20 * (N + 10);
    int q = -1, qdir = 0;
    double best = eps;
#ifdef CSQV_LP_VERIFY_PRICING
    const bool use_keys = false;  // the scan picks; the check below compares the key pick
#else
    const bool use_keys = incremental && !bland;
#endif
    if (use_keys) {
      q = pick_entering(key.data(), N);
      if (q >= 0) qdir = status[q] == 0 ? +1 : -1;
    } else {
      for (int j = 0; j < N; ++j) {
        if (status[j] == 2) continue;
        double dj;
        if (incremental) {
          dj = d[j];
#ifdef CSQV_LP_VERIFY_PRICING
          double rj = cost[j];
          for (int p = 0; p < m; ++p) {
            if (cost[basis[p]] != 0.0) rj -= cost[basis[p]] * Tat(p, j);
          }
          ++compares;
          if (!(rj == dj)) ++mismatches;  // == ignores the sign of an exact zero
#endif
        } else {
          dj = cost[j];
          for (int p = 0; p < m; ++p) {
            if (cost[basis[p]] != 0.0) dj -= cost[basis[p]] * Tat(p, j);
          }
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
    }
#ifdef CSQV_LP_VERIFY_PRICING
    if (incremental && !bland) {
      ++picks;
      if (pick_entering(key.data(), N) != q) ++pick_mismatches;
    }
#endif
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
    // 2) basic variables hitting a bound. Only a nonzero row of column q can reach one.
    col.clear();
    for (int p = 0; p < m; ++p) {
      const double v = Tat(p, q);
      if (v != 0.0) col.push_back({p, v});
    }
    for (const ColEntry& e : col) {
      const int p = e.row;
      const double coef = e.val * qdir;  // xB[p] decreases at rate coef as t grows
      if (coef > eps) {                  // heading toward lower bound
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
    for (const ColEntry& e : col) xB[e.row] -= e.val * delta;

    if (leave_row < 0) {
      // Bound flip: q moves to its opposite bound, stays nonbasic.
      status[q] = (status[q] == 0) ? 1 : 0;
      if (incremental) key[q] = pricing_key(status[q], d[q], eps);
      continue;
    }

    // Pivot: q enters the basis in leave_row; the old basic var leaves. The pivot row is
    // sparse (about 3% nonzero at N=90), so divide only its nonzeros and keep their columns
    // in nz. The row mask gives nz without a pass over all N columns (see "Pivot row"
    // above). The d update and the elimination below skip the zero columns.
    const int leaving = basis[leave_row];
    const double piv = Tat(leave_row, q);
    double* const lrow = &Tat(leave_row, 0);
    uint64_t* const lmask = &row_mask[static_cast<size_t>(leave_row) * W];
    row_nonzeros(lrow, N, lmask, nz_buf.data(), nz);
    for (const int j : nz) lrow[j] /= piv;
    // A bound flip (above) does not change d. A basis change updates d with the divided
    // pivot row. If the guard fails, the rest of this solve uses the full recompute.
    ++basis_changes;
    if (incremental && (!update_reduced_costs(d, lrow, nz, q, piv) ||
                        basis_changes == force_fallback_at)) {
      incremental = false;
      lp_pricing_fallbacks().fetch_add(1, std::memory_order_relaxed);
    }
    const double entering_val = nonbasic_value(q) + delta;
    xB[leave_row] = entering_val;
    for (const ColEntry& e : col) {
      if (e.row == leave_row) continue;
      const double f = e.val;
      double* const prow = &Tat(e.row, 0);
      for (const int j : nz) prow[j] -= f * lrow[j];
      uint64_t* const emask = &row_mask[static_cast<size_t>(e.row) * W];
      for (int w = 0; w < W; ++w) emask[w] |= lmask[w];
    }
    basisRow[leaving] = -1;
    status[leaving] = (leave_to == 1) ? 1 : 0;
    basis[leave_row] = q;
    basisRow[q] = leave_row;
    status[q] = 2;
    // The changed d_j and the two changed statuses are all in nz: the pivot row is piv != 0
    // at q, and 1 at the leaving variable (a basic column is a unit column).
    if (incremental)
      for (const int j : nz) key[j] = pricing_key(status[j], d[j], eps);
  }

  std::vector<double> z(k, 0.0);
  for (int i = 0; i < k; ++i) {
    if (status[i] == 2)
      z[i] = xB[basisRow[i]];
    else
      z[i] = nonbasic_value(i);
    if (z[i] < 0.0) z[i] = 0.0;
    if (z[i] > hi[i]) z[i] = hi[i];
  }
#ifdef CSQV_LP_VERIFY_PRICING
  lp_pricing_compares().fetch_add(compares, std::memory_order_relaxed);
  lp_pricing_mismatches().fetch_add(mismatches, std::memory_order_relaxed);
  lp_pricing_picks().fetch_add(picks, std::memory_order_relaxed);
  lp_pricing_pick_mismatches().fetch_add(pick_mismatches, std::memory_order_relaxed);
#endif

  // Optimal-basis duals (ticket 40). Reduced cost of column j is cost_j - cost_B^T (B^-1 A)_j;
  // the tableau row for a basic variable already holds B^-1 A, so this is a single dot product
  // over the rows whose basic variable carries a non-zero objective (only the radii do).
  if (pair_dual != nullptr || struct_rc != nullptr) {
    auto rcost = [&](int j) {
      double d = cost[j];
      for (int p = 0; p < m; ++p) {
        if (cost[basis[p]] != 0.0) d -= cost[basis[p]] * Tat(p, j);
      }
      return d;
    };
    if (pair_dual != nullptr) {
      pair_dual->assign(m, 0.0);
      for (int p = 0; p < m; ++p) {
        double y = -rcost(k + p);  // shadow price of the pair row = -(slack reduced cost)
        (*pair_dual)[p] = y > 0.0 ? y : 0.0;  // clamp tiny negatives from round-off
      }
    }
    if (struct_rc != nullptr) {
      struct_rc->assign(k, 0.0);
      for (int i = 0; i < k; ++i) (*struct_rc)[i] = rcost(i);
    }
  }
  return z;
}

}  // namespace detail

// Exact optimal radii for fixed centers. Matches problem.py::radii_lp.
inline std::vector<double> radii_lp(const std::vector<double>& x, const std::vector<double>& y,
                                    int n) {
  CSQV_STAGE_LP_SCOPE();  // the stage worker's LP clock (stage_timer.hpp); off in the worker
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

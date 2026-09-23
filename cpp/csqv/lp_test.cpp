// Unit test for lp.hpp::detail::solve_reduced_lp incremental reduced-cost pricing.
//
// The radii-LP tableau is half-integral (every entry is k/2, Nemhauser-Trotter), so the
// reduced costs are exact half-integers in IEEE-754 and the incremental row update
// d_j -= d_q * T(leave, j) gives the SAME doubles as the full recompute. The locked
// invariants:
//   1) BIT-IDENTICAL: incremental pricing and the full recompute return the same radii and
//      the same duals, bit for bit, on real radii LPs (bench and growpush configurations).
//   2) NO FALLBACK: the half-integral guard never trips on those LPs.
//   3) GUARD PREDICATES: is_half_integral accepts k/2 and rejects anything else;
//      is_half_integral_pivot accepts only +-1/2, +-1, +-2.
//   4) GUARD STEP: update_reduced_costs applies d_j -= d_q * row_j, and it reports a
//      failure for a non-half-integral row or pivot (the signal for the fallback).
// Build: see cpp/CMakeLists.txt (target csqv_lp_test) or cpp/Makefile.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include "benchmark.hpp"
#include "lp.hpp"
#include "search.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                             \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                  \
    }                                                                \
  } while (0)

static bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() &&
         (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

// The full LP over every prunable pair of fixed centers (no lazy loop): a large, degenerate
// instance with many pivots, which is what stresses the pricing.
struct FullLp {
  int k = 0;
  std::vector<double> u, rhs;
  std::vector<std::array<int, 2>> pairs;
};

static FullLp full_lp(const std::vector<double>& x, const std::vector<double>& y, int n) {
  FullLp lp;
  lp.k = n;
  lp.u.resize(n);
  for (int i = 0; i < n; ++i) lp.u[i] = std::max(csqv::wall_slack(x[i], y[i]), 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = x[i] - x[j], dy = y[i] - y[j];
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d < lp.u[i] + lp.u[j]) {
        lp.pairs.push_back({i, j});
        lp.rhs.push_back(d);
      }
    }
  }
  return lp;
}

static void check_modes_agree(const FullLp& lp, const char* what) {
  using csqv::detail::Pricing;
  std::vector<double> dual_inc, rc_inc, dual_rec, rc_rec;
  const long before = csqv::lp_pricing_fallbacks().load();
  const std::vector<double> z_inc = csqv::detail::solve_reduced_lp(
      lp.k, lp.u, lp.pairs, lp.rhs, &dual_inc, &rc_inc, Pricing::Incremental);
  const std::vector<double> z_rec = csqv::detail::solve_reduced_lp(
      lp.k, lp.u, lp.pairs, lp.rhs, &dual_rec, &rc_rec, Pricing::Recompute);
  if (!same_bits(z_inc, z_rec)) std::fprintf(stderr, "  radii differ: %s\n", what);
  if (!same_bits(dual_inc, dual_rec)) std::fprintf(stderr, "  pair duals differ: %s\n", what);
  if (!same_bits(rc_inc, rc_rec)) std::fprintf(stderr, "  wall duals differ: %s\n", what);
  CHECK(same_bits(z_inc, z_rec), "incremental radii are bit-identical to the recompute");
  CHECK(same_bits(dual_inc, dual_rec), "incremental pair duals are bit-identical");
  CHECK(same_bits(rc_inc, rc_rec), "incremental structural reduced costs are bit-identical");
  CHECK(csqv::lp_pricing_fallbacks().load() == before, "the half-integral guard never trips");
}

int main() {
  // 1 + 2 on the startup benchmark configurations (random centers, N up to 60).
  {
    const std::pair<int, uint64_t> cases[] = {{5, 1}, {13, 2}, {30, 3}, {60, 4}};
    int checked = 0;
    for (auto [n, seed] : cases) {
      std::vector<double> x, y;
      csqv::bench_centers(n, seed, x, y);
      const FullLp lp = full_lp(x, y, n);
      CHECK(!lp.pairs.empty(), "bench configuration has prunable pairs");
      check_modes_agree(lp, "bench");
      ++checked;
    }
    CHECK(checked == 4, "all bench cases ran");
  }

  // 1 + 2 on grown, near-jammed configurations: the worker's real LP regime, where many
  // pair constraints are tight at once and the simplex is highly degenerate.
  {
    const char* kinds[] = {"grid", "jitter", "hex", "random"};
    for (int n : {20, 45, 70}) {
      for (int s = 0; s < 4; ++s) {
        csqv::Rng rng(static_cast<uint64_t>(100 + s) * 0x9e3779b97f4a7c15ULL + 1);
        std::vector<double> x, y;
        csqv::construct(n, kinds[s], rng, x, y);
        csqv::growpush(x, y, n, rng);
        check_modes_agree(full_lp(x, y, n), kinds[s]);
      }
    }
  }

  // 3: the guard predicates.
  {
    using csqv::detail::is_half_integral;
    using csqv::detail::is_half_integral_pivot;
    for (double v : {0.0, -0.0, 0.5, -0.5, 1.0, -1.5, 2.0, 7.5, -12.0})
      CHECK(is_half_integral(v), "k/2 is half-integral");
    for (double v : {0.25, -0.75, 1.0 / 3.0, 1e-12, 0.5 + 1e-15,
                     std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity()})
      CHECK(!is_half_integral(v), "a non-multiple of 1/2 is rejected");
    for (double p : {0.5, -0.5, 1.0, -1.0, 2.0, -2.0})
      CHECK(is_half_integral_pivot(p), "+-1/2, +-1, +-2 are exact pivots");
    for (double p : {0.0, 0.25, 1.5, 4.0, -3.0, 1.0 + 1e-15})
      CHECK(!is_half_integral_pivot(p), "any other pivot is rejected");
  }

  // 4: the guarded update step, on hand-made rows.
  {
    using csqv::detail::update_reduced_costs;
    // Exact case: d_q = 1, a half-integral row with a unit at q. d_q becomes 0.
    std::vector<double> d = {1.0, 0.5, 0.0, -1.0};
    const std::vector<double> row = {1.0, 0.5, -0.5, 2.0};
    CHECK(update_reduced_costs(d, row.data(), 0, 2.0), "a half-integral step passes the guard");
    CHECK(d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.5 && d[3] == -3.0,
          "the step applies d_j -= d_q * row_j");
    // A row entry of 1/3 makes d_j inexact: the guard must fail.
    std::vector<double> d2 = {1.0, 0.5};
    const std::vector<double> row2 = {1.0, 1.0 / 3.0};
    CHECK(!update_reduced_costs(d2, row2.data(), 0, 1.0), "an inexact d_j fails the guard");
    // A half-integral row with a pivot of 3 (inexact division): the guard must fail.
    std::vector<double> d3 = {1.0, 0.5};
    const std::vector<double> row3 = {1.0, 0.5};
    CHECK(!update_reduced_costs(d3, row3.data(), 0, 3.0), "an inexact pivot fails the guard");
  }

  if (g_failures == 0) std::printf("lp_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}

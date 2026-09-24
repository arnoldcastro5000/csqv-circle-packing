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
//   4) GUARD STEP: update_reduced_costs applies d_j -= d_q * row_j over the nonzero columns
//      of the pivot row, and it reports a failure for a non-half-integral row or pivot (the
//      signal for the fallback). The skipped zero columns give the same bits as the dense
//      update, because d never holds -0.
//   5) FALLBACK MID-SOLVE: a guard failure forced at the K-th basis change (the test hook)
//      switches the solve to the full recompute, counts one fallback, and still returns
//      the same radii and duals, bit for bit. Normal runs never reach this path.
//   6) NO NEGATIVE ZERO: the radii and the duals never hold -0, also when a bound u_i is -0
//      (a center at x = -0). The sparse elimination leaves other zero signs in the tableau
//      than the dense one; this invariant is why those signs cannot reach the output.
//   7) GOLDEN BITS: the radii and the duals of a fixed set of LPs hash to a recorded value.
//      Invariants 1 and 5 compare two paths of the same solver, so a change to a step that
//      both paths share (the ratio test, its tie-break, the elimination) passes them. This
//      check catches it. A bit-identical speed change must keep the value. A deliberate
//      change of the LP output records the new value that the failure message prints.
//   8) ENTERING PICK: the incremental pricing keeps one key per column (the gain of an
//      eligible column, else 0) and picks the column with the largest key, the smallest j on
//      a tie. On random half-integral reduced costs with many ties, this pick is the column
//      of the full scan: the largest |d_j| > eps of an eligible column, the smallest j on a
//      tie.
// Build: see cpp/CMakeLists.txt (target csqv_lp_test) or cpp/Makefile.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
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

// Forces the guard to fail at the K-th basis change for each K in kForceAt. Returns the
// number of solves where the forced failure happened (a solve with fewer than K basis
// changes does not reach it).
static int check_forced_fallback(const FullLp& lp, const char* what) {
  using csqv::detail::Pricing;
  std::vector<double> dual_rec, rc_rec;
  const std::vector<double> z_rec = csqv::detail::solve_reduced_lp(
      lp.k, lp.u, lp.pairs, lp.rhs, &dual_rec, &rc_rec, Pricing::Recompute);
  int tripped = 0;
  for (int force_at : {1, 2, 3, 7, 25}) {
    std::vector<double> dual_fb, rc_fb;
    const long before = csqv::lp_pricing_fallbacks().load();
    const std::vector<double> z_fb = csqv::detail::solve_reduced_lp(
        lp.k, lp.u, lp.pairs, lp.rhs, &dual_fb, &rc_fb, Pricing::Incremental, force_at);
    const long trips = csqv::lp_pricing_fallbacks().load() - before;
    CHECK(trips == 0 || trips == 1, "a solve counts at most one fallback");
    tripped += static_cast<int>(trips);
    if (!same_bits(z_fb, z_rec) || !same_bits(dual_fb, dual_rec) || !same_bits(rc_fb, rc_rec))
      std::fprintf(stderr, "  fallback at basis change %d differs: %s\n", force_at, what);
    CHECK(same_bits(z_fb, z_rec), "fallback radii are bit-identical to the recompute");
    CHECK(same_bits(dual_fb, dual_rec), "fallback pair duals are bit-identical");
    CHECK(same_bits(rc_fb, rc_rec), "fallback structural reduced costs are bit-identical");
  }
  return tripped;
}

static bool has_negative_zero(const std::vector<double>& v) {
  return std::any_of(v.begin(), v.end(), [](double a) { return a == 0.0 && std::signbit(a); });
}

// Solves in both pricing modes and checks the radii and the duals for -0. Returns 1.
static int check_no_negative_zero(const FullLp& lp) {
  using csqv::detail::Pricing;
  for (Pricing pricing : {Pricing::Incremental, Pricing::Recompute}) {
    std::vector<double> dual, rc;
    const std::vector<double> z =
        csqv::detail::solve_reduced_lp(lp.k, lp.u, lp.pairs, lp.rhs, &dual, &rc, pricing);
    CHECK(!has_negative_zero(z), "no radius is -0");
    CHECK(!has_negative_zero(dual), "no pair dual is -0");
    CHECK(!has_negative_zero(rc), "no structural reduced cost is -0");
  }
  return 1;
}

// Grid centers for invariant 7: a near-square grid of n points with exact coordinates
// (c + 0.5) / cols, plus a jitter of up to +-jitter/2 per coordinate. The jitter comes from
// bench_centers, the raw mt19937_64 stream with a fixed map. So every conforming compiler
// and C library gives the same centers. A grid with no jitter has many equal distances, so
// its LP has many ties in the pricing and in the ratio test.
static void grid_centers(int n, uint64_t seed, double jitter, std::vector<double>& x,
                         std::vector<double>& y) {
  std::vector<double> jx, jy;
  csqv::bench_centers(n, seed, jx, jy);
  const int cols = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n))));
  x.resize(n);
  y.resize(n);
  for (int i = 0; i < n; ++i) {
    x[i] = ((i % cols) + 0.5) / cols + jitter * (jx[i] - 0.5);
    y[i] = ((i / cols) + 0.5) / cols + jitter * (jy[i] - 0.5);
  }
}

// The golden LP set of invariant 7. Only portable inputs: the startup benchmark centers and
// the grid centers. (construct and growpush use std::uniform_real_distribution and
// std::normal_distribution, which differ between standard libraries.) For each, the radii
// and the duals of the full LP, then the radii of radii_lp (the lazy active-set loop).
static uint64_t golden_hash() {
  std::vector<double> out;
  auto add = [&](const std::vector<double>& x, const std::vector<double>& y, int n) {
    const FullLp lp = full_lp(x, y, n);
    std::vector<double> dual, rc;
    const std::vector<double> z =
        csqv::detail::solve_reduced_lp(lp.k, lp.u, lp.pairs, lp.rhs, &dual, &rc);
    const std::vector<double> r = csqv::radii_lp(x, y, n);
    for (const std::vector<double>& v : {z, dual, rc, r}) out.insert(out.end(), v.begin(), v.end());
  };
  const std::pair<int, uint64_t> cases[] = {{5, 1}, {13, 2}, {30, 3}, {60, 4}, {90, 5}, {121, 6}};
  for (auto [n, seed] : cases) {
    std::vector<double> x, y;
    csqv::bench_centers(n, seed, x, y);
    add(x, y, n);
  }
  for (int n : {20, 45, 70, 90, 121}) {
    for (double jitter : {0.0, 1e-6, 1e-3, 3e-2}) {
      std::vector<double> x, y;
      grid_centers(n, static_cast<uint64_t>(1000 + n), jitter, x, y);
      add(x, y, n);
    }
  }
  return csqv::hash_doubles(out);
}

// The entering column of the full pricing scan (the rule of solve_reduced_lp): the eligible
// column with the largest |d_j| > eps; a later column wins only with a strictly larger gain.
static int scan_pick(const std::vector<int>& status, const std::vector<double>& d, double eps) {
  int q = -1;
  double best = eps;
  for (size_t j = 0; j < d.size(); ++j) {
    const bool eligible = (status[j] == 0 && d[j] > eps) || (status[j] == 1 && d[j] < -eps);
    if (eligible && std::fabs(d[j]) > best) {
      best = std::fabs(d[j]);
      q = static_cast<int>(j);
    }
  }
  return q;
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

  // 5 on the same two families. The forced failure must happen in most solves, or the hook
  // has no effect and the check proves nothing. (A small LP can end before the K-th basis
  // change.)
  {
    int tripped = 0, solves = 0;
    const std::pair<int, uint64_t> cases[] = {{5, 1}, {13, 2}, {30, 3}, {60, 4}};
    for (auto [n, seed] : cases) {
      std::vector<double> x, y;
      csqv::bench_centers(n, seed, x, y);
      tripped += check_forced_fallback(full_lp(x, y, n), "bench");
      solves += 5;
    }
    const char* kinds[] = {"grid", "jitter", "hex", "random"};
    for (int n : {20, 45, 70}) {
      for (int s = 0; s < 4; ++s) {
        csqv::Rng rng(static_cast<uint64_t>(100 + s) * 0x9e3779b97f4a7c15ULL + 1);
        std::vector<double> x, y;
        csqv::construct(n, kinds[s], rng, x, y);
        csqv::growpush(x, y, n, rng);
        tripped += check_forced_fallback(full_lp(x, y, n), kinds[s]);
        solves += 5;
      }
    }
    std::printf("lp_test: forced fallback tripped in %d of %d solves\n", tripped, solves);
    CHECK(tripped >= solves * 3 / 4, "the forced fallback trips in most solves");
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

  // 4: the guarded update step, on hand-made rows. nz lists the nonzero columns of the row.
  {
    using csqv::detail::update_reduced_costs;
    // Exact case: d_q = 1, a half-integral row with a unit at q. d_q becomes 0.
    std::vector<double> d = {1.0, 0.5, 0.0, -1.0};
    const std::vector<double> row = {1.0, 0.5, -0.5, 2.0};
    CHECK(update_reduced_costs(d, row.data(), {0, 1, 2, 3}, 0, 2.0),
          "a half-integral step passes the guard");
    CHECK(d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.5 && d[3] == -3.0,
          "the step applies d_j -= d_q * row_j");
    // A row with zero columns of both signs. The sparse step skips them and must give the
    // same bits as the dense step d_j -= d_q * row_j over every column.
    const std::vector<double> d0 = {0.0, -1.5, 1.0, 0.0, 2.0, -0.5};
    const std::vector<double> zrow = {-0.0, 1.0, 0.0, 0.5, -0.0, -2.0};
    for (double dq : {1.0, -1.0, 2.0}) {
      std::vector<double> dense = d0, sparse = d0;
      dense[1] = sparse[1] = dq;
      const double q_val = dense[1];
      for (size_t j = 0; j < dense.size(); ++j) dense[j] -= q_val * zrow[j];
      CHECK(update_reduced_costs(sparse, zrow.data(), {1, 3, 5}, 1, -1.0),
            "a sparse half-integral step passes the guard");
      CHECK(same_bits(sparse, dense), "the sparse step gives the bits of the dense step");
    }
    // A row entry of 1/3 makes d_j inexact: the guard must fail.
    std::vector<double> d2 = {1.0, 0.5};
    const std::vector<double> row2 = {1.0, 1.0 / 3.0};
    CHECK(!update_reduced_costs(d2, row2.data(), {0, 1}, 0, 1.0),
          "an inexact d_j fails the guard");
    // A half-integral row with a pivot of 3 (inexact division): the guard must fail.
    std::vector<double> d3 = {1.0, 0.5};
    const std::vector<double> row3 = {1.0, 0.5};
    CHECK(!update_reduced_costs(d3, row3.data(), {0, 1}, 0, 3.0),
          "an inexact pivot fails the guard");
  }

  // 6: no -0 in the output. A hand-made LP first: circle 0 has the bound u_0 = -0 (its
  // center is at x = -0). The simplex flips it to its upper bound, so the radius is the
  // bound itself. Then the real LP families.
  {
    const std::vector<double> u = {-0.0, 0.5};
    const std::vector<std::array<int, 2>> pairs = {{0, 1}};
    const std::vector<double> rhs = {1.0};
    std::vector<double> dual, rc;
    const std::vector<double> z = csqv::detail::solve_reduced_lp(2, u, pairs, rhs, &dual, &rc);
    CHECK(z.size() == 2 && z[0] == 0.0 && z[1] == 0.5, "the hand-made LP solves");
    CHECK(!has_negative_zero(z), "a -0 bound gives a +0 radius");
    CHECK(!has_negative_zero(dual) && !has_negative_zero(rc), "a -0 bound gives no -0 dual");

    int checked = 0;
    const std::pair<int, uint64_t> cases[] = {{5, 1}, {13, 2}, {30, 3}, {60, 4}};
    for (auto [n, seed] : cases) {
      std::vector<double> x, y;
      csqv::bench_centers(n, seed, x, y);
      checked += check_no_negative_zero(full_lp(x, y, n));
    }
    const char* kinds[] = {"grid", "jitter", "hex", "random"};
    for (int n : {20, 45, 70}) {
      for (int s = 0; s < 4; ++s) {
        csqv::Rng rng(static_cast<uint64_t>(100 + s) * 0x9e3779b97f4a7c15ULL + 1);
        std::vector<double> x, y;
        csqv::construct(n, kinds[s], rng, x, y);
        csqv::growpush(x, y, n, rng);
        checked += check_no_negative_zero(full_lp(x, y, n));
      }
    }
    CHECK(checked == 16, "all no-negative-zero cases ran");
  }

  // 7: the golden bits.
  {
    const uint64_t kGolden = 0x0a4c1ce804e3e72cULL;
    const uint64_t h = golden_hash();
    if (h != kGolden)
      std::fprintf(stderr, "  golden hash 0x%016llxULL, expected 0x%016llxULL\n",
                   static_cast<unsigned long long>(h), static_cast<unsigned long long>(kGolden));
    CHECK(h == kGolden, "the LP outputs keep the recorded bits");
  }

  // 8: the entering pick from the keys. Hand-made keys first, then random columns.
  {
    using csqv::detail::pick_entering;
    using csqv::detail::pricing_key;
    const double eps = 1e-12;
    CHECK(pricing_key(0, 1.5, eps) == pricing_key(1, -1.5, eps), "the key is the gain");
    CHECK(pricing_key(0, -1.0, eps) == 0 && pricing_key(1, 1.0, eps) == 0,
          "a column that moves the wrong way has key 0");
    CHECK(pricing_key(2, 1.0, eps) == 0 && pricing_key(2, -1.0, eps) == 0,
          "a basic column has key 0");
    CHECK(pricing_key(0, 0.0, eps) == 0 && pricing_key(1, -0.0, eps) == 0 &&
              pricing_key(0, 1e-13, eps) == 0 && pricing_key(1, -1e-13, eps) == 0,
          "a gain of eps or less has key 0");
    CHECK(pricing_key(0, 0.5, eps) < pricing_key(0, 1.0, eps) &&
              pricing_key(0, 1.0, eps) < pricing_key(1, -2.0, eps) &&
              pricing_key(1, -2.0, eps) < pricing_key(0, 1e300, eps),
          "the keys have the order of the gains");
    const std::vector<int64_t> none = {0, 0, 0};
    CHECK(pick_entering(none.data(), 3) == -1, "no eligible column gives -1");
    CHECK(pick_entering(none.data(), 0) == -1, "no column gives -1");
    const int64_t one = pricing_key(0, 1.0, eps), two = pricing_key(0, 2.0, eps);
    const std::vector<int64_t> tie = {0, one, 0, one, one};
    CHECK(pick_entering(tie.data(), 5) == 1, "a tie picks the smallest column");
    const std::vector<int64_t> later = {one, 0, one, two, two, one};
    CHECK(pick_entering(later.data(), 6) == 3, "a larger gain wins; its tie picks the first");

    // Random columns: status 0, 1 or 2 and reduced costs k/2 with |k| <= 6, so most gains
    // tie. The raw mt19937_64 stream with a fixed map gives the same cases everywhere.
    std::mt19937_64 rng(20260924);
    int agree = 0, cases = 0;
    for (int c = 0; c < 20000; ++c) {
      const int n = static_cast<int>(rng() % 300);
      std::vector<int> status(n);
      std::vector<double> d(n);
      std::vector<int64_t> keys(n);
      for (int j = 0; j < n; ++j) {
        status[j] = static_cast<int>(rng() % 3);
        d[j] = (static_cast<int>(rng() % 13) - 6) * 0.5;
        keys[j] = pricing_key(status[j], d[j], eps);
      }
      agree += pick_entering(keys.data(), n) == scan_pick(status, d, eps);
      ++cases;
    }
    if (agree != cases) std::fprintf(stderr, "  pick differs in %d of %d cases\n", cases - agree, cases);
    CHECK(agree == cases, "the key pick is the column of the full scan");
  }

  if (g_failures == 0) std::printf("lp_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}

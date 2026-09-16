// Unit test for polish.hpp::center_polish (exact-LP center ascent).
//
// center_polish locally maximizes sum(radii_lp(centers)) over the centers, with the radii
// always the exact LP optimum for the fixed centers. It is the terminal "squeeze" the
// worker runs on a champion: it recovers the last ~1e-5 of sum that the penalty polish
// leaves on the table (the gap that Packomania's own re-optimization closed on N=121).
//
// The locked invariants:
//   1) MONOTONE: it never lowers the scored sum (safe to run before comparing to best).
//   2) IMPROVES a configuration with obvious slack (a clear ascent, not a no-op).
//   3) FEASIBLE and IN-BOX: the result passes verify_and_score with centers in [0, SIDE].
//   4) CONSISTENT: the returned score equals verify_and_score on the returned packing.
// Build: see cpp/CMakeLists.txt (target csqv_center_polish_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>
#include <vector>

#include "polish.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                              \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                  \
    }                                                                \
  } while (0)

static double score_of(const std::vector<double>& x, const std::vector<double>& y, int n) {
  std::vector<double> r = csqv::radii_lp(x, y, n);
  double s;
  if (!csqv::verify_and_score(x, y, r, n, csqv::TOL, s)) return -1.0;
  return s;
}

static bool in_box(const std::vector<double>& x, const std::vector<double>& y, int n) {
  for (int i = 0; i < n; ++i)
    if (x[i] < 0.0 || x[i] > csqv::SIDE || y[i] < 0.0 || y[i] > csqv::SIDE) return false;
  return true;
}

int main() {
  // 1 + 2 + 3 + 4: a suboptimal N=2 configuration with obvious slack. Two circles crowded
  // near the middle can only share the tiny gap between them; pushing them apart lets both
  // grow. The reachable local optimum (a vertical stack, both radius 1/4) sums to 0.5, well
  // above the crowded start (~0.2). center_polish must climb clearly and stay feasible.
  {
    const int n = 2;
    std::vector<double> x = {0.5, 0.5}, y = {0.4, 0.6}, r;
    const double before = score_of(x, y, n);
    const double after = csqv::center_polish(x, y, r, n);
    CHECK(after > before + 0.1, "N=2 crowded start must improve clearly");
    CHECK(after <= 0.585787 + 1e-6, "N=2 sum cannot exceed the two-circle optimum");
    CHECK(in_box(x, y, n), "N=2 centers stay in the box");
    double sc;
    CHECK(csqv::verify_and_score(x, y, r, n, csqv::TOL, sc), "N=2 result is feasible");
    CHECK(std::fabs(sc - after) < 1e-12, "N=2 returned score matches verify_and_score");
    CHECK(std::fabs(after - score_of(x, y, n)) < 1e-12, "N=2 radii are the exact LP optimum");
  }

  // 1 (monotone) on a near-optimal input: the symmetric two-circle diagonal optimum
  // (r = 1/(2+sqrt2)) is already a local maximum. center_polish must NOT regress it.
  {
    const int n = 2;
    const double r0 = 1.0 / (2.0 + std::sqrt(2.0));
    std::vector<double> x = {r0, 1.0 - r0}, y = {r0, 1.0 - r0}, r;
    const double before = score_of(x, y, n);
    const double after = csqv::center_polish(x, y, r, n);
    CHECK(after >= before - 1e-10, "near-optimal input must not regress (monotone)");
    CHECK(in_box(x, y, n), "monotone case: centers stay in the box");
  }

  // 1 + 3 on a larger case: a 3x3 grid has slack (grid packings are not optimal for the
  // sum objective). center_polish must not regress and must stay feasible.
  {
    const int n = 9;
    std::vector<double> x(n), y(n), r;
    const double g[3] = {1.0 / 6.0, 3.0 / 6.0, 5.0 / 6.0};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        x[a * 3 + b] = g[b];
        y[a * 3 + b] = g[a];
      }
    const double before = score_of(x, y, n);
    const double after = csqv::center_polish(x, y, r, n);
    CHECK(after >= before - 1e-10, "N=9 grid must not regress");
    CHECK(in_box(x, y, n), "N=9 centers stay in the box");
    double sc;
    CHECK(csqv::verify_and_score(x, y, r, n, csqv::TOL, sc), "N=9 result is feasible");
  }

  // Analytic dual gradient (ticket 40). On a GENERIC config (deterministic pseudo-random
  // centers, away from the walls and the x=y / x=1-y tie lines) the value function is smooth,
  // so the analytic gradient from the LP duals must agree with the central finite-difference
  // gradient of lp_score and must be a genuine ascent direction.
  {
    const int n = 8;
    std::vector<double> x(n), y(n);
    unsigned s = 2463534242u;
    auto rnd = [&]() {  // xorshift, deterministic; returns [0,1)
      s ^= s << 13; s ^= s >> 17; s ^= s << 5;
      return (s >> 8) / 16777216.0;
    };
    for (int i = 0; i < n; ++i) { x[i] = 0.15 + 0.70 * rnd(); y[i] = 0.15 + 0.70 * rnd(); }

    std::vector<double> gx, gy;
    csqv::center_gradient(x, y, n, gx, gy);

    // Central finite-difference gradient of lp_score, same eps as center_polish.
    const double eps = 1e-6;
    std::vector<double> fx(n), fy(n);
    for (int i = 0; i < n; ++i) {
      const double ox = x[i];
      const double xp = std::min(ox + eps, csqv::SIDE), xm = std::max(ox - eps, 0.0);
      x[i] = xp; const double fxp = csqv::lp_score(x, y, n);
      x[i] = xm; const double fxm = csqv::lp_score(x, y, n);
      x[i] = ox; fx[i] = (fxp - fxm) / (xp - xm);
      const double oy = y[i];
      const double yp = std::min(oy + eps, csqv::SIDE), ym = std::max(oy - eps, 0.0);
      y[i] = yp; const double fyp = csqv::lp_score(x, y, n);
      y[i] = ym; const double fym = csqv::lp_score(x, y, n);
      y[i] = oy; fy[i] = (fyp - fym) / (yp - ym);
    }

    double dot = 0.0, na = 0.0, nf = 0.0;
    for (int i = 0; i < n; ++i) {
      dot += gx[i] * fx[i] + gy[i] * fy[i];
      na += gx[i] * gx[i] + gy[i] * gy[i];
      nf += fx[i] * fx[i] + fy[i] * fy[i];
    }
    CHECK(na > 1e-12 && nf > 1e-12, "both gradients are non-trivial on a generic config");
    const double cosine = dot / (std::sqrt(na) * std::sqrt(nf) + 1e-30);
    CHECK(cosine > 0.999, "analytic dual gradient matches the finite-difference gradient");

    // The analytic gradient is a real ascent direction: a small step along it raises lp_score.
    const double base = csqv::lp_score(x, y, n);
    const double inv = 1.0 / std::sqrt(na);
    std::vector<double> ax(n), ay(n);
    for (int i = 0; i < n; ++i) {
      ax[i] = std::min(std::max(x[i] + 1e-5 * gx[i] * inv, 0.0), csqv::SIDE);
      ay[i] = std::min(std::max(y[i] + 1e-5 * gy[i] * inv, 0.0), csqv::SIDE);
    }
    CHECK(csqv::lp_score(ax, ay, n) > base, "a step along the analytic gradient increases the score");
  }

  // center_polish_analytic parity: same crowded N=2 start as case 1. It must climb clearly,
  // stay feasible, and reach the same two-circle optimum as the finite-difference polish.
  {
    const int n = 2;
    std::vector<double> x = {0.5, 0.5}, y = {0.4, 0.6}, r;
    const double before = score_of(x, y, n);
    const double after = csqv::center_polish_analytic(x, y, r, n);
    CHECK(after > before + 0.1, "analytic polish: crowded start must improve clearly");
    CHECK(after <= 0.585787 + 1e-6, "analytic polish: cannot exceed the two-circle optimum");
    CHECK(in_box(x, y, n), "analytic polish: centers stay in the box");
    double sc;
    CHECK(csqv::verify_and_score(x, y, r, n, csqv::TOL, sc), "analytic polish: result is feasible");
    CHECK(std::fabs(sc - after) < 1e-12, "analytic polish: returned score matches verify_and_score");
  }

  // Monotone on a near-optimal input: the symmetric two-circle diagonal optimum must not regress.
  {
    const int n = 2;
    const double r0 = 1.0 / (2.0 + std::sqrt(2.0));
    std::vector<double> x = {r0, 1.0 - r0}, y = {r0, 1.0 - r0}, r;
    const double before = score_of(x, y, n);
    const double after = csqv::center_polish_analytic(x, y, r, n);
    CHECK(after >= before - 1e-10, "analytic polish: near-optimal input must not regress");
  }

  if (g_failures == 0) {
    std::printf("center_polish_test: all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "center_polish_test: %d check(s) FAILED\n", g_failures);
  return 1;
}

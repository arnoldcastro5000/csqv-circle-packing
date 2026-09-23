// Unit test for jam_slp.hpp (the trust-region SLP jamming optimizer and its helpers).
//
// Locked invariants (ticket 02):
//   1) MONOTONE + FEASIBLE: jam_slp never lowers the exact score and emits a feasible packing;
//      its output radii equal radii_lp(output centers).
//   2) NO-OP ON A JAMMED INPUT: the 2x2 grid is already collectively jammed (Donev margin ~0),
//      so jam_slp holds it at the optimum sum 1.0.
//   3) ASCENDS AN UNJAMMED INPUT: a grid nudged off the walls is not jammed; jam_slp drives it
//      back to the jammed optimum, raising the score and dropping the Donev margin.
//   4) n < 2 returns the input unchanged.
// Build: see cpp/CMakeLists.txt (target csqv_jam_slp_test) or cpp/Makefile.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "jam_slp.hpp"
#include "lp.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                             \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                  \
    }                                                                \
  } while (0)

static bool approx(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

int main() {
  using namespace csqv;

  // A 2x2 grid: centers at 0.25/0.75. radii_lp gives r = 0.25 (touching walls and neighbors).
  const std::vector<double> gx = {0.25, 0.75, 0.25, 0.75};
  const std::vector<double> gy = {0.25, 0.25, 0.75, 0.75};
  const int n = 4;

  // 2) jammed grid: jam_slp holds the optimum; output radii == radii_lp(output).
  {
    const double before = exact_score(gx, gy, n);
    CHECK(approx(before, 1.0, 1e-6), "jammed grid score ~ 1.0 before");
    JamResult r = jam_slp(gx, gy, n);
    CHECK(r.score >= before - 1e-12, "jam_slp monotone on jammed grid");
    CHECK(approx(r.score, 1.0, 1e-6), "jam_slp holds optimum 1.0");
    const std::vector<double> rlp = radii_lp(r.x, r.y, n);
    double maxdiff = 0.0;
    for (int i = 0; i < n; ++i) maxdiff = std::max(maxdiff, std::fabs(rlp[i] - r.r[i]));
    CHECK(maxdiff < 1e-12, "output radii == radii_lp(output centers)");
    // feasible at the verifier tolerance
    double sc;
    const bool feas = verify_and_score(r.x, r.y, r.r, n, TOL, sc);
    CHECK(feas && std::isfinite(sc), "jammed grid output feasible");
    // Donev margin at the jammed optimum is ~0.
    const ActiveSet as = active_set(r.x, r.y, r.r, n, 1e-7, 1e-7);
    const double margin = donev_jamming_margin(r.x, r.y, r.r, as, n);
    CHECK(margin <= 1e-6, "Donev margin ~ 0 at the jammed optimum");
  }

  // 3) unjammed grid nudged inward off the walls: jam_slp must ascend back toward 1.0.
  {
    const std::vector<double> px = {0.27, 0.73, 0.27, 0.73};
    const std::vector<double> py = {0.27, 0.27, 0.73, 0.73};
    const double before = exact_score(px, py, n);
    CHECK(before < 0.95, "unjammed grid starts sub-optimal");
    // Donev margin before: only contacts are active (walls are >gamma off), so an unjamming
    // flex exists -> margin clearly positive.
    const ActiveSet as0 = active_set(px, py, radii_lp(px, py, n), n, 1e-2, 1e-2);
    const double margin0 = donev_jamming_margin(px, py, radii_lp(px, py, n), as0, n);
    CHECK(margin0 > 1e-3, "unjammed grid has a positive Donev margin before");

    JamResult r = jam_slp(px, py, n);
    CHECK(r.score > before + 1e-6, "jam_slp strictly ascends the unjammed grid");
    CHECK(r.score > 0.99, "jam_slp reaches near the jammed optimum");
    const ActiveSet as1 = active_set(r.x, r.y, r.r, n, 1e-7, 1e-7);
    const double margin1 = donev_jamming_margin(r.x, r.y, r.r, as1, n);
    CHECK(margin1 < margin0, "Donev margin drops after jamming");
    double sc;
    const bool feas = verify_and_score(r.x, r.y, r.r, n, TOL, sc);
    CHECK(feas && std::isfinite(sc), "unjammed grid output feasible");
  }

  // 4) n < 2 returns the input unchanged.
  {
    const std::vector<double> sx = {0.5}, sy = {0.5};
    JamResult r = jam_slp(sx, sy, 1);
    CHECK(r.x.size() == 1 && approx(r.x[0], 0.5, 0.0), "n<2 returns input x");
  }

  if (g_failures == 0)
    std::printf("jam_slp_test: OK\n");
  else
    std::printf("jam_slp_test: %d FAILURE(S)\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}

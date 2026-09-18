// Unit test for general_lp.hpp::general_lp (the general two-phase bounded-variable simplex).
//
// Locked invariants (ticket 01):
//   1) '<=' LP with a trivial feasible start solves without phase 1.
//   2) '>=' and '==' rows (no feasible slack start) solve via phase 1.
//   3) negative lower bounds work, including the no-constraint (box-only) case.
//   4) a fully free variable (both bounds infinite) solves via the split.
//   5) an infeasible LP is reported Infeasible; an unbounded objective is reported Unbounded.
//   6) a flex-LP-shaped problem (the shape jam_slp builds) reaches the hand-checked optimum.
// Build: see cpp/CMakeLists.txt (target csqv_general_lp_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "general_lp.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                              \
  do {                                                                \
    if (!(cond)) {                                                    \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__);  \
      ++g_failures;                                                   \
    }                                                                 \
  } while (0)

static const double INF = std::numeric_limits<double>::infinity();

static bool approx(double a, double b, double tol = 1e-7) { return std::fabs(a - b) <= tol; }

int main() {
  using csqv::general_lp;
  using csqv::LpStatus;

  // 1) max x+y s.t. x+y <= 4, with 0 <= x,y <= 3. Optimum sum = 4.
  {
    auto r = general_lp(2, 1, {1, 1}, {1, 1}, {-1}, {4}, {0, 0}, {3, 3});
    CHECK(r.status == LpStatus::Optimal, "case1 optimal");
    CHECK(approx(r.objective, 4.0), "case1 obj == 4");
    CHECK(r.x[0] >= -1e-9 && r.x[0] <= 3 + 1e-9, "case1 x in box");
    CHECK(approx(r.x[0] + r.x[1], 4.0), "case1 x+y == 4");
  }

  // 2) max x1 s.t. x1+x2 == 3, x1-x2 >= 1, x >= 0. Phase 1 needed. Optimum x1 = 3, x2 = 0.
  {
    auto r = general_lp(2, 2, {1, 0}, {1, 1, 1, -1}, {0, +1}, {3, 1}, {0, 0}, {INF, INF});
    CHECK(r.status == LpStatus::Optimal, "case2 optimal");
    CHECK(approx(r.objective, 3.0), "case2 obj == 3");
    CHECK(approx(r.x[0], 3.0) && approx(r.x[1], 0.0), "case2 x == (3,0)");
    CHECK(approx(r.x[0] + r.x[1], 3.0), "case2 equality holds");
    CHECK(r.x[0] - r.x[1] >= 1 - 1e-7, "case2 inequality holds");
  }

  // 3) negative lower bound, no constraints (box only): max -x with -3 <= x <= 2. Optimum x = -3.
  {
    auto r = general_lp(1, 0, {-1}, {}, {}, {}, {-3}, {2});
    CHECK(r.status == LpStatus::Optimal, "case3 optimal");
    CHECK(approx(r.objective, 3.0), "case3 obj == 3");
    CHECK(approx(r.x[0], -3.0), "case3 x == -3");
  }

  // 4) free variable: max x s.t. x + y == 5, x free, 0 <= y <= 10. Optimum x = 5 (y = 0).
  {
    auto r = general_lp(2, 1, {1, 0}, {1, 1}, {0}, {5}, {-INF, 0}, {INF, 10});
    CHECK(r.status == LpStatus::Optimal, "case4a optimal");
    CHECK(approx(r.objective, 5.0), "case4a obj == 5");
    CHECK(approx(r.x[0], 5.0) && approx(r.x[1], 0.0), "case4a x == (5,0)");
  }
  // 4b) same but maximize -x: the free var must go negative. Optimum x = -5 (y = 10).
  {
    auto r = general_lp(2, 1, {-1, 0}, {1, 1}, {0}, {5}, {-INF, 0}, {INF, 10});
    CHECK(r.status == LpStatus::Optimal, "case4b optimal");
    CHECK(approx(r.objective, 5.0), "case4b obj == 5");
    CHECK(approx(r.x[0], -5.0) && approx(r.x[1], 10.0), "case4b x == (-5,10)");
  }

  // 5a) infeasible: x <= 2 and x >= 5, 0 <= x. Reported Infeasible.
  {
    auto r = general_lp(1, 2, {1}, {1, 1}, {-1, +1}, {2, 5}, {0}, {INF});
    CHECK(r.status == LpStatus::Infeasible, "case5a infeasible");
  }
  // 5b) unbounded: max x, 0 <= x, no upper constraint. Reported Unbounded.
  {
    auto r = general_lp(1, 0, {1}, {}, {}, {}, {0}, {INF});
    CHECK(r.status == LpStatus::Unbounded, "case5b unbounded");
  }

  // 6) flex-LP shape (the jam_slp form): max dr1+dr2 s.t. dr1 - dz <= 0.01, dr2 + dz <= 0.01,
  //    dz in [-0.05, 0.05], dr in [0, 0.03]. Hand optimum: dr1 = dr2 = 0.01 (dz = 0), sum = 0.02.
  //    Variables ordered [dz, dr1, dr2].
  {
    std::vector<double> A = {-1, 1, 0,   // dr1 - dz <= 0.01
                             +1, 0, 1};  // dr2 + dz <= 0.01
    auto r = general_lp(3, 2, {0, 1, 1}, A, {-1, -1}, {0.01, 0.01},
                        {-0.05, 0.0, 0.0}, {0.05, 0.03, 0.03});
    CHECK(r.status == LpStatus::Optimal, "case6 optimal");
    CHECK(approx(r.objective, 0.02), "case6 obj == 0.02");
    // feasibility of the returned point
    CHECK(r.x[1] - r.x[0] <= 0.01 + 1e-9, "case6 row1 feasible");
    CHECK(r.x[2] + r.x[0] <= 0.01 + 1e-9, "case6 row2 feasible");
    CHECK(r.x[0] >= -0.05 - 1e-9 && r.x[0] <= 0.05 + 1e-9, "case6 dz in box");
  }

  // 7) redundant equality rows: x+y == 2 and 2x+2y == 4 (the same constraint). Phase 1 leaves an
  //    artificial basic at 0 in the redundant row; the snap-before-freeze must keep phase 2 sound.
  //    max x, x,y >= 0. Optimum x = 2, y = 0.
  {
    auto r = general_lp(2, 2, {1, 0}, {1, 1, 2, 2}, {0, 0}, {2, 4}, {0, 0}, {INF, INF});
    CHECK(r.status == LpStatus::Optimal, "case7 optimal");
    CHECK(approx(r.objective, 2.0), "case7 obj == 2");
    CHECK(approx(r.x[0], 2.0) && approx(r.x[1], 0.0), "case7 x == (2,0)");
    CHECK(approx(r.x[0] + r.x[1], 2.0), "case7 equality holds");
  }

  if (g_failures == 0)
    std::printf("general_lp_test: OK\n");
  else
    std::printf("general_lp_test: %d FAILURE(S)\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}

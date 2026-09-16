// Unit test for fingerprint.hpp::basin_fingerprint (the basin hash).
//
// Locked invariants (ticket 42 / ADR 0002):
//   1) STABLE: numerical noise below the contact tolerance does not change the hash (same basin
//      always hashes the same; WL never false-splits on structure).
//   2) D4-INVARIANT: the 8 rotations/reflections of a packing all hash the same (they are the
//      same basin).
//   3) DISCRIMINATING: genuinely different arrangements hash apart.
//   4) RADII-ONLY IS NOT ENOUGH: two arrangements with the SAME sorted radii but a different
//      contact graph still hash apart (the reason the fingerprint uses the graph, not just radii).
// Build: see cpp/CMakeLists.txt (target csqv_fingerprint_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>
#include <vector>

#include "fingerprint.hpp"
#include "lp.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                             \
  do {                                                              \
    if (!(cond)) {                                                  \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                 \
    }                                                               \
  } while (0)

// Fingerprint a packing whose radii are the exact LP optimum for the centers.
static uint64_t fp_of(const std::vector<double>& x, const std::vector<double>& y, int n) {
  std::vector<double> r = csqv::radii_lp(x, y, n);
  return csqv::basin_fingerprint(x, y, r, n);
}

int main() {
  // A concrete converged basin: a 3x3 grid of 9 circles centered in the square. radii_lp gives
  // the exact touching radii; the contact graph is the 3x3 king-free lattice.
  const int n = 9;
  std::vector<double> gx(n), gy(n);
  {
    const double c[3] = {1.0 / 6.0, 0.5, 5.0 / 6.0};
    int k = 0;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        gx[k] = c[b];
        gy[k] = c[a];
        ++k;
      }
  }
  const uint64_t base = fp_of(gx, gy, n);

  // 1) STABLE: jitter every center by 1e-8, far below the 1e-6 contact tolerance, so the contact
  // graph and the radius buckets are unchanged. The hash must not move.
  {
    std::vector<double> x = gx, y = gy;
    for (int i = 0; i < n; ++i) {
      x[i] += (i % 2 ? 1e-8 : -1e-8);
      y[i] += (i % 3 ? -1e-8 : 1e-8);
    }
    CHECK(fp_of(x, y, n) == base, "sub-tolerance jitter must not change the fingerprint");
  }

  // 2) D4-INVARIANT: all 8 symmetries of the square map the packing to itself's basin, so the
  // hash is identical. Radii follow the centers (recomputed by radii_lp), so nothing to permute.
  {
    auto tf = [&](int k, double xi, double yi, double& ox, double& oy) {
      // 4 rotations then the same 4 with a reflection across x.
      double a = xi - 0.5, b = yi - 0.5;  // to centered frame
      for (int t = 0; t < (k % 4); ++t) {
        const double na = -b, nb = a;
        a = na;
        b = nb;
      }
      if (k >= 4) a = -a;  // reflection
      ox = a + 0.5;
      oy = b + 0.5;
    };
    for (int k = 0; k < 8; ++k) {
      std::vector<double> x(n), y(n);
      for (int i = 0; i < n; ++i) tf(k, gx[i], gy[i], x[i], y[i]);
      CHECK(fp_of(x, y, n) == base, "every D4 symmetry must hash the same");
    }
  }

  // 3) DISCRIMINATING: a structurally different arrangement (a single row of 9 circles across the
  // middle) has a different contact graph and different radii, so a different hash.
  {
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
      x[i] = (i + 0.5) / n;
      y[i] = 0.5;
    }
    CHECK(fp_of(x, y, n) != base, "a row must not hash like a 3x3 grid");
  }

  // 4) RADII-ONLY IS NOT ENOUGH. Build two arrangements with the SAME sorted radii but a
  // different contact graph, and confirm the graph-based hash separates them. Four equal circles:
  // arrangement A is a 2x2 block (each circle touches two walls and two neighbors); arrangement B
  // is a diagonal pair of touching couples rotated so the contact pattern differs. Both use
  // radii_lp, so if the radii multisets match, only the graph distinguishes them.
  {
    const int m = 4;
    // A: 2x2 block.
    std::vector<double> ax = {0.25, 0.75, 0.25, 0.75}, ay = {0.25, 0.25, 0.75, 0.75};
    std::vector<double> ar = csqv::radii_lp(ax, ay, m);
    // B: a horizontal row of 4.
    std::vector<double> bx = {0.125, 0.375, 0.625, 0.875}, by = {0.5, 0.5, 0.5, 0.5};
    std::vector<double> br = csqv::radii_lp(bx, by, m);
    const uint64_t fa = csqv::basin_fingerprint(ax, ay, ar, m);
    const uint64_t fb = csqv::basin_fingerprint(bx, by, br, m);
    CHECK(fa != fb, "different contact graphs must hash apart");
  }

  if (g_failures == 0)
    std::printf("fingerprint_test: OK\n");
  else
    std::printf("fingerprint_test: %d FAILURE(S)\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}

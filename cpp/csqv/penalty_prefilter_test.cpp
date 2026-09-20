// Unit test for polish.hpp::penalty_obj_grad after the squared-distance pre-filter.
//
// The pre-filter must be EXACT: it only skips the sqrt for pairs that cannot overlap, so
// the objective f and the gradient g must be BIT-FOR-BIT identical to the dense sweep that
// tests every pair with the eager sqrt. This test recomputes f and g with the original
// eager arithmetic (a local reference below) and asserts exact equality against the
// production penalty_obj_grad, over configurations that mix overlapping and far pairs.
// Build: see cpp/CMakeLists.txt (target csqv_penalty_prefilter_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>
#include <vector>

#include "polish.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                             \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                  \
    }                                                                \
  } while (0)

// Reference objective+gradient with the ORIGINAL eager arithmetic (no pre-filter): every
// pair gets the sqrt, then the o>0 branch. This mirrors polish.hpp before the change, so a
// bit-exact match proves the pre-filter changed nothing.
static double ref_obj_grad(const std::vector<double>& v, int n, const std::vector<int>& pi,
                           const std::vector<int>& pj, double lam, std::vector<double>& g) {
  std::fill(g.begin(), g.end(), 0.0);
  double f = 0.0;
  for (int i = 0; i < n; ++i) {
    f -= v[3 * i + 2];
    g[3 * i + 2] -= 1.0;
  }
  const size_t P = pi.size();
  for (size_t p = 0; p < P; ++p) {
    const int i = pi[p], j = pj[p];
    const double xi = v[3 * i], yi = v[3 * i + 1], ri = v[3 * i + 2];
    const double xj = v[3 * j], yj = v[3 * j + 1], rj = v[3 * j + 2];
    const double dx = xi - xj, dy = yi - yj;
    const double d = std::sqrt(dx * dx + dy * dy) + 1e-12;
    const double o = (ri + rj) - d;
    if (o > 0.0) {
      f += lam * o * o;
      const double w = 2.0 * lam * o;
      const double ux = dx / d, uy = dy / d;
      g[3 * i + 2] += w;
      g[3 * j + 2] += w;
      g[3 * i] += w * (-ux);
      g[3 * j] += w * (ux);
      g[3 * i + 1] += w * (-uy);
      g[3 * j + 1] += w * (uy);
    }
  }
  for (int i = 0; i < n; ++i) {
    const double x = v[3 * i], y = v[3 * i + 1], r = v[3 * i + 2];
    double vl = r - x, vr = x + r - csqv::SIDE, vb = r - y, vt = y + r - csqv::SIDE;
    if (vl > 0.0) { f += lam * vl * vl; g[3 * i + 2] += 2.0 * lam * vl; g[3 * i] += 2.0 * lam * vl * (-1.0); }
    if (vr > 0.0) { f += lam * vr * vr; g[3 * i + 2] += 2.0 * lam * vr; g[3 * i] += 2.0 * lam * vr * (1.0); }
    if (vb > 0.0) { f += lam * vb * vb; g[3 * i + 2] += 2.0 * lam * vb; g[3 * i + 1] += 2.0 * lam * vb * (-1.0); }
    if (vt > 0.0) { f += lam * vt * vt; g[3 * i + 2] += 2.0 * lam * vt; g[3 * i + 1] += 2.0 * lam * vt * (1.0); }
  }
  return f;
}

// Assert exact bit equality of f and every g entry against the reference, for one config.
static void assert_exact(const std::vector<double>& v, int n, const char* label) {
  std::vector<int> pi, pj;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) { pi.push_back(i); pj.push_back(j); }
  const double lam = 1e4;

  std::vector<double> g(3 * n), gref(3 * n);
  const double f = csqv::penalty_obj_grad(v, n, pi, pj, lam, g);
  const double fref = ref_obj_grad(v, n, pi, pj, lam, gref);

  CHECK(f == fref, label);  // bit-exact objective
  bool gmatch = true;
  for (int t = 0; t < 3 * n; ++t)
    if (g[t] != gref[t]) gmatch = false;
  CHECK(gmatch, label);  // bit-exact gradient
}

int main() {
  // Deterministic pseudo-random configs, radii >= 0 (as the polish always feeds after
  // project_box). Mix crowded (many overlaps) and spread (few overlaps) cases.
  unsigned s = 88172645u;
  auto rnd = [&]() {
    s ^= s << 13; s ^= s >> 17; s ^= s << 5;
    return (s >> 8) / 16777216.0;
  };

  // Crowded: small box, generous radii -> many overlapping pairs exercise the o>0 path.
  for (int trial = 0; trial < 5; ++trial) {
    const int n = 12 + trial * 6;
    std::vector<double> v(3 * n);
    for (int i = 0; i < n; ++i) {
      v[3 * i] = 0.2 + 0.6 * rnd();
      v[3 * i + 1] = 0.2 + 0.6 * rnd();
      v[3 * i + 2] = 0.02 + 0.10 * rnd();
    }
    assert_exact(v, n, "crowded config: f and g bit-exact vs dense reference");
  }

  // Spread: full box, tiny radii -> most pairs skip; the pre-filter's common path.
  for (int trial = 0; trial < 5; ++trial) {
    const int n = 20 + trial * 10;
    std::vector<double> v(3 * n);
    for (int i = 0; i < n; ++i) {
      v[3 * i] = rnd();
      v[3 * i + 1] = rnd();
      v[3 * i + 2] = 0.005 + 0.02 * rnd();
    }
    assert_exact(v, n, "spread config: f and g bit-exact vs dense reference");
  }

  // Edge: some radii exactly 0 (zero-radius pairs must never trigger the penalty).
  {
    const int n = 8;
    std::vector<double> v(3 * n);
    for (int i = 0; i < n; ++i) {
      v[3 * i] = 0.3 + 0.4 * rnd();
      v[3 * i + 1] = 0.3 + 0.4 * rnd();
      v[3 * i + 2] = (i % 2 == 0) ? 0.0 : 0.06;
    }
    assert_exact(v, n, "zero-radius edge: f and g bit-exact vs dense reference");
  }

  // Edge: near-touching pair (distance just under the radius sum) sits in the o>0 sliver.
  {
    const int n = 2;
    std::vector<double> v = {0.40, 0.5, 0.10, 0.59, 0.5, 0.10};  // d=0.19 < 0.20 = r+r
    assert_exact(v, n, "near-touching pair: f and g bit-exact vs dense reference");
  }

  if (g_failures == 0) {
    std::printf("penalty_prefilter_test: all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "penalty_prefilter_test: %d check(s) FAILED\n", g_failures);
  return 1;
}

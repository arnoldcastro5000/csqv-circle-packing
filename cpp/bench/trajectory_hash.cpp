// Trajectory hash (gate G2). One FNV-1a hash over the bytes of every output on the path
// construct -> growpush -> center_gradient (radii + dual gradient) -> center_polish_analytic.
// `gate` hashes only the startup benchmark radii (the G1 cases). A bit-identical change
// keeps both hashes.
//
// Usage: trajectory_hash <n> <seeds> <polish_iters>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include "benchmark.hpp"
#include "polish.hpp"
#include "search.hpp"

static uint64_t g_hash = 1469598103934665603ULL;

static void mix(const std::vector<double>& v) {
  for (double d : v) {
    uint64_t b;
    std::memcpy(&b, &d, sizeof b);
    for (int k = 0; k < 8; ++k) {
      g_hash ^= (b >> (8 * k)) & 0xff;
      g_hash *= 1099511628211ULL;
    }
  }
}

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 121;
  const int seeds = argc > 2 ? std::atoi(argv[2]) : 3;
  const int polish_iters = argc > 3 ? std::atoi(argv[3]) : 5;
  const std::pair<int, uint64_t> cases[] = {{5, 1}, {13, 2}, {30, 3}, {60, 4}, {100, 5}, {121, 6}};
  for (auto [cn, cs] : cases) {
    std::vector<double> x, y;
    csqv::bench_centers(cn, cs, x, y);
    mix(csqv::radii_lp(x, y, cn));
  }
  const uint64_t gate = g_hash;
  const auto t0 = std::chrono::steady_clock::now();
  const char* kinds[] = {"hex", "jitter", "random"};
  for (int s = 0; s < seeds; ++s) {
    csqv::Rng g(1000 + s);
    std::vector<double> x, y, gx, gy;
    csqv::construct(n, kinds[s % 3], g, x, y);
    csqv::growpush(x, y, n, g);
    mix(x);
    mix(y);
    std::vector<double> r = csqv::center_gradient(x, y, n, gx, gy);
    mix(r);
    mix(gx);
    mix(gy);
    const double score = csqv::center_polish_analytic(x, y, r, n, polish_iters);
    mix(x);
    mix(y);
    mix(r);
    mix(std::vector<double>{score});
  }
  const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("n=%d seeds=%d gate=%016llx all=%016llx time=%.3f\n", n, seeds,
              static_cast<unsigned long long>(gate), static_cast<unsigned long long>(g_hash), el);
  return 0;
}

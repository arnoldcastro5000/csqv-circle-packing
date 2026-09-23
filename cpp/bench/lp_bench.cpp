// LP-only timing: radii_lp on fixed construct and growpush configurations (8 seeds, 16
// configurations, 3 repeats). The acc value is the sum of every radius; a bit-identical
// change keeps it.
//
// Usage: lp_bench <n>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "lp.hpp"
#include "search.hpp"

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 90;
  const char* kinds[] = {"grid", "jitter", "hex", "random"};
  std::vector<std::vector<double>> xs, ys;
  for (int s = 0; s < 8; ++s) {
    csqv::Rng g(500 + s);
    std::vector<double> x, y;
    csqv::construct(n, kinds[s % 4], g, x, y);
    xs.push_back(x);
    ys.push_back(y);
    csqv::growpush(x, y, n, g);
    xs.push_back(x);
    ys.push_back(y);
  }
  double acc = 0.0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int rep = 0; rep < 3; ++rep)
    for (size_t c = 0; c < xs.size(); ++c)
      for (double r : csqv::radii_lp(xs[c], ys[c], n)) acc += r;
  const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("N=%d lp_calls=%zu time=%.3f acc=%a\n", n, 3 * xs.size(), el, acc);
  return 0;
}

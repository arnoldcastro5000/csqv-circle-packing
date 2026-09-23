// Perf driver: the worker's restart hot path (construct -> growpush -> radii_lp ->
// scalable_polish), deterministic for a seed, one thread. It prints restarts/s and the
// sealed checksum (the sum of the polished scores).
//
// Usage: perf_driver <n> <iters> <seed>
// Sealed checksums (seed 777): `perf_driver 90 20 777` = 98.374190007874,
// `perf_driver 143 8 777` = 49.910861299027. A bit-identical change keeps the hex value.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "geometry.hpp"
#include "lp.hpp"
#include "polish.hpp"
#include "search.hpp"

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 90;
  const int iters = argc > 2 ? std::atoi(argv[2]) : 30;
  const int seed = argc > 3 ? std::atoi(argv[3]) : 12345;
  csqv::Rng rng(static_cast<uint64_t>(seed) * 0x9e3779b97f4a7c15ULL + 1);
  const char* kinds[] = {"grid", "jitter", "hex", "random"};
  std::vector<double> x, y, r;
  double checksum = 0.0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int restart = 0; restart < iters; ++restart) {
    csqv::construct(n, kinds[restart % 4], rng, x, y);
    csqv::growpush(x, y, n, rng);
    r = csqv::radii_lp(x, y, n);
    checksum += csqv::scalable_polish(x, y, r, n);
  }
  const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("N=%d iters=%d time=%.3f rate=%.4f/s checksum=%.12f hex=%a\n", n, iters, el,
              iters / el, checksum, checksum);
  return 0;
}

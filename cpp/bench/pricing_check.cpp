// In-loop pricing check. Build with -DCSQV_LP_VERIFY_PRICING: then each pricing step in
// solve_reduced_lp compares the incremental reduced cost d_j with the full recompute. This
// driver runs the worker's LP calls (construct, growpush, center_gradient,
// center_polish_analytic) and prints the counts. The gate: mismatches=0 and fallbacks=0.
//
// Usage: pricing_check <n> <seeds>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "lp.hpp"
#include "polish.hpp"
#include "search.hpp"

#ifndef CSQV_LP_VERIFY_PRICING
#error "build pricing_check with -DCSQV_LP_VERIFY_PRICING"
#endif

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 90;
  const int seeds = argc > 2 ? std::atoi(argv[2]) : 4;
  const char* kinds[] = {"grid", "jitter", "hex", "random"};
  for (int s = 0; s < seeds; ++s) {
    csqv::Rng g(7000 + s);
    std::vector<double> x, y, gx, gy;
    csqv::construct(n, kinds[s % 4], g, x, y);
    csqv::radii_lp(x, y, n);
    csqv::growpush(x, y, n, g);
    csqv::radii_lp(x, y, n);
    csqv::center_gradient(x, y, n, gx, gy);
    std::vector<double> r = csqv::radii_lp(x, y, n);
    csqv::center_polish_analytic(x, y, r, n, 3);
  }
  const long mismatches = csqv::detail::lp_pricing_mismatches().load();
  const long fallbacks = csqv::lp_pricing_fallbacks().load();
  std::printf("n=%d seeds=%d compares=%ld mismatches=%ld fallbacks=%ld\n", n, seeds,
              csqv::detail::lp_pricing_compares().load(), mismatches, fallbacks);
  return mismatches == 0 && fallbacks == 0 ? 0 : 1;
}

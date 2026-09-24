// Unit test: the build does not contract a * b + c into a fused multiply-add (FMA).
//
// An FMA rounds once, a * b then + c rounds twice, so the two can give different doubles.
// The sealed checksums and the startup radii hashes assume the two roundings. Each -march
// that has FMA (x86-64-v3, znver3, native on a modern CPU, and every arm64 target) lets GCC
// contract, unless the build also gives -ffp-contract=off. This test fails in such a build.
//
// a = 1 + 2^-27 and b = 1 - 2^-27, so a * b = 1 - 2^-54 exactly. Rounded, a * b is 1.0 and
// a * b - 1 is 0. Fused, a * b - 1 is -2^-54. The inputs are volatile, so the compiler
// cannot fold the expression at compile time.
// Build: see cpp/CMakeLists.txt (target csqv_fp_contract_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>

int main() {
  volatile double va = 1.0 + std::ldexp(1.0, -27);
  volatile double vb = 1.0 - std::ldexp(1.0, -27);
  volatile double vc = -1.0;
  const double a = va, b = vb, c = vc;
  const double r = a * b + c;
  if (r != 0.0) {
    std::fprintf(stderr,
                 "FAIL: a * b + c = %a, not 0: the build contracts to FMA. Add -ffp-contract=off.\n",
                 r);
    return 1;
  }
  std::printf("fp_contract_test: all checks passed\n");
  return 0;
}

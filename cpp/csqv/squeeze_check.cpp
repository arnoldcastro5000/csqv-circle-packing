// squeeze_check: measure the first-order center_polish_analytic squeeze on a champion.
//
// This is the permanent home of the squeeze-audit harness (ticket 45). It loads a
// centered-frame champion (a Packomania .pck or an n<N>-best.txt), re-derives the exact
// radii with radii_lp, then applies center_polish_analytic (the first-order LP-dual
// squeeze, cpp/csqv/polish.hpp) and reports the residual gain. It is BOTH the squeeze
// audit ("is this .pck fully squeezed?") and the ANALYTIC BASELINE the terminal
// second-order post-processor (problems/csqv/terminal_squeeze.py) is measured against.
//
// A .pck stores "x y r" rows in the centered side-1 square with a header (line 1 = the
// largest radius, line 2 = the author). This reads every whitespace line with three
// numeric tokens as an (x, y, r) row and ignores the header, so it parses a .pck and an
// n-best.txt the same way. The r column is ignored: the radii are always re-derived.
//
// Build (or via `make -C cpp` / CMake, target csqv_squeeze_check):
//   g++ -O3 -std=c++17 -Icsqv csqv/squeeze_check.cpp -o build/csqv_squeeze_check
// Usage: squeeze_check <champion.pck> [conditioned_out.txt]
//   Prints:  the raw radii_lp sum, the analytic-squeezed sum, and the gain.
//   With conditioned_out.txt it also WRITES the analytic-conditioned centers (centered
//   frame, "x y r" at 15 dp), the warm start for the second-order squeeze.
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "lp.hpp"
#include "polish.hpp"

namespace {

// Read every line with exactly three numeric tokens as a centered-frame (x, y, r) row.
// Returns false if the file cannot be opened or holds no such row.
bool load_centered(const std::string& path, std::vector<double>& x, std::vector<double>& y) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream is(line);
    double a, b, c;
    if (is >> a >> b >> c) {
      std::string extra;
      if (is >> extra) continue;  // more than three tokens: not a coordinate row
      x.push_back(a + 0.5);        // centered -> unit
      y.push_back(b + 0.5);
    }
  }
  return !x.empty();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <champion.pck> [conditioned_out.txt]\n", argv[0]);
    return 2;
  }
  std::vector<double> x, y;
  if (!load_centered(argv[1], x, y)) {
    std::fprintf(stderr, "cannot read a centered-frame champion from %s\n", argv[1]);
    return 2;
  }
  const int n = static_cast<int>(x.size());

  std::vector<double> r = csqv::radii_lp(x, y, n);
  double raw;
  if (!csqv::verify_and_score(x, y, r, n, csqv::TOL, raw)) {
    std::fprintf(stderr, "the champion is INVALID (radii_lp did not verify)\n");
    return 1;
  }

  std::vector<double> px = x, py = y, pr = r;
  const double squeezed = csqv::center_polish_analytic(px, py, pr, n);

  std::printf("N=%d\n", n);
  std::printf("radii_lp sum        = %.12f\n", raw);
  std::printf("center_polish_analytic = %.12f\n", squeezed);
  std::printf("analytic gain       = %.3e\n", squeezed - raw);

  if (argc >= 3) {
    std::ofstream out(argv[2]);
    if (!out) {
      std::fprintf(stderr, "cannot write %s\n", argv[2]);
      return 1;
    }
    for (int i = 0; i < n; ++i)
      out << std::scientific << std::setprecision(15) << (px[i] - 0.5) << " " << (py[i] - 0.5)
          << " " << pr[i] << "\n";
    std::printf("wrote conditioned centers to %s\n", argv[2]);
  }
  return 0;
}

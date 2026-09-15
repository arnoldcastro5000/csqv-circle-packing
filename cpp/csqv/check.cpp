// Parity harness: read UNIT-frame centers, print radii_lp + verify_and_score.
// A Python oracle (tools/make_golden.py) produces the expected values from the trusted
// primitives; tools/run_parity.py diffs this program's output against them.
//
// Input:  line 1 = n, then n lines "x y" (unit square [0,1]).
// Output: line 1 = score (%.12f)
//         line 2 = max_pair_overlap max_wall_violation  (feasibility, want <= TOL)
//         then n lines = radius (%.12f).
#include <cmath>
#include <cstdio>
#include <vector>

#include "geometry.hpp"
#include "lp.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <centers_file>\n", argv[0]);
    return 2;
  }
  std::FILE* f = std::fopen(argv[1], "r");
  if (!f) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  int n = 0;
  if (std::fscanf(f, "%d", &n) != 1 || n <= 0) {
    std::fclose(f);
    return 2;
  }
  std::vector<double> x(n), y(n);
  for (int i = 0; i < n; ++i) {
    if (std::fscanf(f, "%lf %lf", &x[i], &y[i]) != 2) {
      std::fclose(f);
      return 2;
    }
  }
  std::fclose(f);

  std::vector<double> r = csqv::radii_lp(x, y, n);
  double score;
  csqv::verify_and_score(x, y, r, n, csqv::TOL, score);

  // Feasibility of the raw radii_lp output (before any shrink).
  double max_overlap = 0.0, max_wall = 0.0;
  for (int i = 0; i < n; ++i) {
    max_wall = std::max(max_wall, r[i] - std::max(csqv::wall_slack(x[i], y[i]), 0.0));
    for (int j = i + 1; j < n; ++j) {
      const double dx = x[i] - x[j], dy = y[i] - y[j];
      const double d = std::sqrt(dx * dx + dy * dy);
      max_overlap = std::max(max_overlap, (r[i] + r[j]) - d);
    }
  }
  std::printf("%.12f\n", score);
  std::printf("%.3e %.3e\n", max_overlap, max_wall);
  for (int i = 0; i < n; ++i) std::printf("%.12f\n", r[i]);
  return 0;
}

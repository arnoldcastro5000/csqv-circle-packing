// squeeze_check: measure the C++ terminal squeeze on a champion, in the C++ frame.
//
// Two modes on a centered-frame champion (a Packomania .pck or an n<N>-best.txt):
//
//   default        the first-order center_polish_analytic squeeze (cpp/csqv/polish.hpp). This
//                  is the squeeze audit ("is this .pck fully squeezed?") and the ANALYTIC
//                  BASELINE the second-order post-processor is measured against.
//   --jam          the trust-region SLP jammer (cpp/csqv/jam_slp.hpp), the shipped worker seal.
//                  Prints the raw radii_lp sum, the jammed sum, the gain, and the Donev jamming
//                  margin before and after (~0 = jammed). This is the C++-FRAME half of the
//                  acceptance gate; the acceptance harness re-scores the result with the trusted
//                  Python verifier (source of truth) and diffs it against Python jam_slp in one
//                  frame.
//
// A .pck stores "x y r" rows in the centered side-1 square with a header (line 1 = the largest
// radius, line 2 = the author). This reads every whitespace line with exactly three numeric
// tokens as an (x, y, r) row and ignores the header, so it parses a .pck and an n-best.txt the
// same way. The r column is ignored: the radii are always re-derived with radii_lp.
//
// Build (or via `make -C cpp` / CMake, target csqv_squeeze_check):
//   g++ -O3 -std=c++17 -Icsqv csqv/squeeze_check.cpp -o build/csqv_squeeze_check
// Usage:
//   squeeze_check <champion> [--jam] [--out-centers PATH] [--out-pck PATH]
//                 [--author NAME] [--dp N]
//   <champion>       a .pck or n<N>-best.txt (centered frame).
//   --out-centers    write the squeezed centers (centered "x y r" at 15 dp). A bare second
//                    positional path is accepted as this, for backward compatibility.
//   --out-pck        write a strictly-feasible native .pck of the squeezed packing (--jam only).
//   --author         author line for --out-pck (default "csqv").
//   --dp             .pck decimals for --out-pck (default 15).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "jam_slp.hpp"
#include "lp.hpp"
#include "pck.hpp"
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

// The Donev jamming margin at the hand-off active set (matches terminal_squeeze.py::_margin:
// pair/wall tol 1e-7). ~0 means the contact set is rigid (collectively jammed).
double jamming_margin(const std::vector<double>& x, const std::vector<double>& y,
                      const std::vector<double>& r, int n) {
  const csqv::ActiveSet as = csqv::active_set(x, y, r, n, 1e-7, 1e-7);
  return csqv::donev_jamming_margin(x, y, r, as, n);
}

// Write centered-frame "x y r" rows at 15 dp (scientific), the second-order warm start / the
// jammed centers for the Python parity re-score.
bool write_centers(const std::string& path, const std::vector<double>& x,
                   const std::vector<double>& y, const std::vector<double>& r, int n) {
  std::ofstream out(path);
  if (!out) return false;
  for (int i = 0; i < n; ++i)
    out << std::scientific << std::setprecision(15) << (x[i] - 0.5) << " " << (y[i] - 0.5) << " "
        << r[i] << "\n";
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  std::string champion, out_centers, out_pck, author = "csqv";
  int dp = 15;
  bool jam = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "%s needs a value\n", name);
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "--jam") jam = true;
    else if (a == "--out-centers") out_centers = next("--out-centers");
    else if (a == "--out-pck") out_pck = next("--out-pck");
    else if (a == "--author") author = next("--author");
    else if (a == "--dp") dp = std::atoi(next("--dp").c_str());
    else if (!a.empty() && a[0] == '-') {
      std::fprintf(stderr, "unknown option %s\n", a.c_str());
      return 2;
    } else if (champion.empty()) champion = a;
    else if (out_centers.empty()) out_centers = a;  // legacy bare second positional
    else {
      std::fprintf(stderr, "unexpected argument %s\n", a.c_str());
      return 2;
    }
  }
  if (champion.empty()) {
    std::fprintf(stderr,
                 "usage: %s <champion> [--jam] [--out-centers PATH] [--out-pck PATH] "
                 "[--author NAME] [--dp N]\n",
                 argv[0]);
    return 2;
  }
  if (!out_pck.empty() && !jam) {
    std::fprintf(stderr, "--out-pck requires --jam\n");
    return 2;
  }

  std::vector<double> x, y;
  if (!load_centered(champion, x, y)) {
    std::fprintf(stderr, "cannot read a centered-frame champion from %s\n", champion.c_str());
    return 2;
  }
  const int n = static_cast<int>(x.size());

  std::vector<double> r = csqv::radii_lp(x, y, n);
  double raw;
  if (!csqv::verify_and_score(x, y, r, n, csqv::TOL, raw)) {
    std::fprintf(stderr, "the champion is INVALID (radii_lp did not verify)\n");
    return 1;
  }
  std::printf("N=%d\n", n);
  std::printf("radii_lp sum        = %.12f\n", raw);

  std::vector<double> ox, oy, orr;  // squeezed centers/radii to (optionally) write out
  if (jam) {
    const double margin_before = jamming_margin(x, y, r, n);
    const csqv::JamResult jr = csqv::jam_slp(x, y, n);
    const double margin_after = jamming_margin(jr.x, jr.y, jr.r, n);
    std::printf("jam_slp sum         = %.12f\n", jr.score);
    std::printf("jam gain            = %.3e\n", jr.score - raw);
    std::printf("Donev margin        before=%.3e  after=%.3e\n", margin_before, margin_after);
    ox = jr.x;
    oy = jr.y;
    orr = jr.r;
  } else {
    std::vector<double> px = x, py = y, pr = r;
    const double squeezed = csqv::center_polish_analytic(px, py, pr, n);
    std::printf("center_polish_analytic = %.12f\n", squeezed);
    std::printf("analytic gain       = %.3e\n", squeezed - raw);
    ox = px;
    oy = py;
    orr = pr;
  }

  if (!out_centers.empty()) {
    if (!write_centers(out_centers, ox, oy, orr, n)) {
      std::fprintf(stderr, "cannot write %s\n", out_centers.c_str());
      return 1;
    }
    std::printf("wrote centers to %s\n", out_centers.c_str());
  }
  if (!out_pck.empty()) {
    csqv::Packing canon;
    if (!csqv::canonical_feasible(ox, oy, orr, n, csqv::TOL, canon)) {
      std::fprintf(stderr, "canonical_feasible failed on the jammed packing\n");
      return 1;
    }
    const double sum = csqv::emit_pck(out_pck, canon, author, dp);
    if (sum < 0.0) {
      std::fprintf(stderr, "emit_pck could not write a strictly feasible %s\n", out_pck.c_str());
      return 1;
    }
    std::printf("wrote .pck (sum %.12f) to %s\n", sum, out_pck.c_str());
  }
  return 0;
}

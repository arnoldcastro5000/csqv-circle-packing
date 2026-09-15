// Startup accuracy/precision self-benchmark.
//
// On every worker start this runs a FIXED, deterministic suite of the exact primitives
// (radii_lp + verify_and_score) and writes the results to the output folder. The point is
// cross-validation:
//   - ACCURACY: cpp/tools/cross_validate.py recomputes each case with the trusted Python
//     primitives from the embedded centers and diffs.
//   - PRECISION / REPRODUCIBILITY: the cases are generated from the RAW mt19937_64 sequence
//     (standard-specified, identical on every conforming compiler) mapped to [0,1) with a
//     fixed exact formula, NOT std::uniform_real_distribution (which is implementation
//     defined). So a Linux build and a Windows build see the SAME centers, and their
//     benchmark files diff directly: identical sums and radii hashes mean the numerics are
//     reproducible across the two OSes.
//
// Values are written at 17 significant digits (round-trip) plus a bitwise FNV-1a hash of
// the radii, so a single-ULP divergence is detectable.
#pragma once
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "lp.hpp"

namespace csqv {

// Deterministic centers from the raw mt19937_64 stream, identical across platforms.
inline void bench_centers(int n, uint64_t seed, std::vector<double>& x, std::vector<double>& y) {
  std::mt19937_64 g(seed);
  x.resize(n);
  y.resize(n);
  auto next01 = [&]() { return (g() >> 11) * (1.0 / 9007199254740992.0); };  // 53-bit -> [0,1)
  for (int i = 0; i < n; ++i) {
    x[i] = next01();
    y[i] = next01();
  }
}

// FNV-1a over the raw bytes of each double. x86 Linux and Windows are both little-endian,
// so identical doubles hash identically across the two.
inline uint64_t hash_doubles(const std::vector<double>& v) {
  uint64_t h = 1469598103934665603ULL;
  for (double d : v) {
    uint64_t bits;
    std::memcpy(&bits, &d, sizeof(bits));
    for (int b = 0; b < 8; ++b) {
      h ^= (bits >> (8 * b)) & 0xffULL;
      h *= 1099511628211ULL;
    }
  }
  return h;
}

inline std::string build_info() {
  std::string os =
#if defined(_WIN32)
      "windows";
#elif defined(__linux__)
      "linux";
#elif defined(__APPLE__)
      "macos";
#else
      "unknown";
#endif
  std::string compiler =
#if defined(__clang__)
      std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
      std::string("gcc ") + __VERSION__;
#else
      "unknown";
#endif
#if defined(__MINGW64__)
  compiler += " (mingw-w64)";
#elif defined(__MINGW32__)
  compiler += " (mingw32)";
#endif
  std::string ieee =
#if defined(__STDC_IEC_559__)
      "yes";
#else
      "unverified";
#endif
  return "platform: " + os + "  compiler: " + compiler + "  ieee754: " + ieee +
         "  sizeof(long double): " + std::to_string(sizeof(long double));
}

// Run the suite, write the report to path. Returns true if every case is strictly feasible
// (max overlap and wall violation both <= TOL), i.e. the primitives behave correctly here.
inline bool run_benchmark(const std::string& path, bool& all_feasible) {
  // Fixed (N, seed) cases. Sizes stay small so startup cost is well under a second.
  const std::pair<int, uint64_t> cases[] = {
      {5, 1}, {13, 2}, {30, 3}, {60, 4}, {100, 5}, {121, 6}};

  std::ofstream f(path);
  if (!f) return false;
  f.setf(std::ios::fmtflags(0));

  char when[32];
  std::time_t t = std::time(nullptr);
  std::strftime(when, sizeof(when), "%Y-%m-%d %H:%M:%S", std::localtime(&t));

  f << "# CSQV worker accuracy/precision benchmark v1\n";
  f << "# generated: " << when << "\n";
  f << "# " << build_info() << "\n";
  f << "# cross-validate:  PYTHONPATH=. python3 cpp/tools/cross_validate.py " << path << "\n";
  f << "# reproducibility: diff this file against another OS's benchmark (sums + hashes)\n";
  f << "# values at 17 significant digits; radii_hash = FNV-1a over the radii bit patterns\n";

  all_feasible = true;
  auto put = [&](const char* k, double v) {
    char b[64];
    std::snprintf(b, sizeof(b), "%s=%.17g", k, v);
    f << b;
  };

  for (auto [n, seed] : cases) {
    std::vector<double> x, y;
    bench_centers(n, seed, x, y);
    std::vector<double> r = radii_lp(x, y, n);
    double score;
    verify_and_score(x, y, r, n, TOL, score);

    double max_overlap = 0.0, max_wall = 0.0, sum = 0.0;
    for (int i = 0; i < n; ++i) {
      sum += r[i];
      max_wall = std::max(max_wall, r[i] - std::max(wall_slack(x[i], y[i]), 0.0));
      for (int j = i + 1; j < n; ++j) {
        const double dx = x[i] - x[j], dy = y[i] - y[j];
        max_overlap = std::max(max_overlap, (r[i] + r[j]) - std::sqrt(dx * dx + dy * dy));
      }
    }
    const bool feasible = max_overlap <= TOL && max_wall <= TOL;
    all_feasible = all_feasible && feasible;

    char hx[32];
    std::snprintf(hx, sizeof(hx), "0x%016llx", static_cast<unsigned long long>(hash_doubles(r)));

    f << "CASE N=" << n << " seed=" << seed << "\n";
    f << "SUMMARY ";
    put("sum", sum);
    f << " ";
    put("overlap", max_overlap);
    f << " ";
    put("wall", max_wall);
    f << " ";
    put("score", score);
    f << " feasible=" << (feasible ? 1 : 0) << " radii_hash=" << hx << "\n";
    f << "DATA\n";
    for (int i = 0; i < n; ++i) {
      char b[80];
      std::snprintf(b, sizeof(b), "%.17g %.17g %.17g\n", x[i], y[i], r[i]);
      f << b;
    }
    f << "END\n";
  }
  return true;
}

}  // namespace csqv

// Seed spreading: per-thread low-discrepancy cold constructions (ticket 42, spec section 1).
//
// The old worker builds a fresh construction from kinds[restarts % 4], a list SHARED across
// threads, so two threads fall into the same basin and coverage is < 3x for 3 threads. Seed
// spreading gives each thread its OWN disjoint, evenly-spread slice of a low-discrepancy sequence
// over the construction family and its parameters, so threads start in different regions. It
// steers WHERE a search starts, never how deep it digs (ADR 0002).
//
// Sequence choice. The spec named "Sobol". This uses the R2 (Roberts generalized-golden-ratio)
// low-discrepancy sequence instead: it has the same even-spread (low star-discrepancy) property,
// but it is a five-line dependency-free recurrence with no direction-number tables, so it fits the
// worker's self-contained no-dependency rule and is simpler to maintain. Swap in a full Sobol
// generator here if the exact named sequence is ever required; the call sites do not change.
#pragma once
#include <cmath>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "search.hpp"

namespace csqv {

// M_PI is not defined under MinGW's strict C++17, so use a local Pi constant (portable across the
// native and the Windows cross builds).
constexpr double kPi = 3.14159265358979323846;

// The R2 low-discrepancy sequence over the unit hypercube [0,1)^d. Point k, coordinate i, is
// frac(0.5 + alpha_i * k) with alpha_i = frac(1 / phi^(i+1)) and phi the generalized plastic
// constant solving phi^(d+1) = phi + 1. Successive points fill the cube evenly with no clumps.
struct LowDisc {
  double alpha[8] = {0};
  int dims = 0;
  explicit LowDisc(int d) : dims(d) {
    // Solve phi^(d+1) = phi + 1 by the contraction phi <- (1 + phi)^(1/(d+1)).
    double phi = 2.0;
    for (int it = 0; it < 64; ++it) phi = std::pow(1.0 + phi, 1.0 / (d + 1));
    double inv = 1.0;
    for (int i = 0; i < d && i < 8; ++i) {
      inv /= phi;                      // inv = phi^-(i+1)
      double a = inv - std::floor(inv);  // frac
      alpha[i] = a;
    }
  }
  double coord(uint64_t k, int i) const {
    const double v = 0.5 + alpha[i] * static_cast<double>(k);
    return v - std::floor(v);
  }
};

// Base centered-frame layout for a family, in the centered square [-0.5, 0.5]^2 before transform.
// grid/jitter/hex/random reuse the trusted construct(); the three new families build directly.
inline void spread_base(int n, int family, Rng& g, std::vector<double>& bx,
                        std::vector<double>& by) {
  bx.assign(n, 0.0);
  by.assign(n, 0.0);
  const double golden = kPi * (3.0 - std::sqrt(5.0));  // ~2.39996 rad, the golden angle
  if (family <= 3) {
    const char* kinds[] = {"grid", "jitter", "hex", "random"};
    std::vector<double> ux, uy;
    construct(n, kinds[family], g, ux, uy);  // unit-frame [0,1]
    for (int i = 0; i < n; ++i) {
      bx[i] = ux[i] - 0.5;
      by[i] = uy[i] - 0.5;
    }
    return;
  }
  if (family == 4) {  // phyllotaxis (Vogel sunflower)
    const double R = 0.47;
    for (int i = 0; i < n; ++i) {
      const double rr = R * std::sqrt((i + 0.5) / n);
      const double th = i * golden;
      bx[i] = rr * std::cos(th);
      by[i] = rr * std::sin(th);
    }
  } else if (family == 5) {  // Archimedean spiral
    const double turns = 3.0 + std::sqrt(static_cast<double>(n)) / 4.0;
    const double th_max = 2.0 * kPi * turns;
    const double R = 0.47;
    for (int i = 0; i < n; ++i) {
      const double f = (n == 1) ? 0.0 : static_cast<double>(i) / (n - 1);
      const double th = f * th_max;
      const double rr = R * f;
      bx[i] = rr * std::cos(th);
      by[i] = rr * std::sin(th);
    }
  } else {  // family == 6: concentric rings
    const int rings = std::max(1, static_cast<int>(std::round(std::sqrt(n / 3.0))));
    int placed = 0;
    for (int ring = 0; ring < rings && placed < n; ++ring) {
      const double rr = 0.47 * (ring + 1) / rings;
      const int remaining_rings = rings - ring;
      int cnt = (n - placed) / remaining_rings;
      if (cnt < 1) cnt = 1;
      if (ring == rings - 1) cnt = n - placed;  // last ring takes the rest
      for (int t = 0; t < cnt && placed < n; ++t) {
        const double th = 2.0 * kPi * t / cnt;
        bx[placed] = rr * std::cos(th);
        by[placed] = rr * std::sin(th);
        ++placed;
      }
    }
  }
}

// set to 5, 6 or 7 to enable experimental families. 4 is the baseline.
constexpr int kSpreadFamilies = 4;

// Build a spread cold construction. `k` is this thread's global low-discrepancy index (use
// disjoint per-thread slices, e.g. k = base + thread + draw * num_threads). The low-discrepancy
// point picks the family and the rotation/aspect/jitter, so successive constructions on one thread
// fan out across the space instead of repeating one kind. Output centers are in the unit frame,
// clamped inside the box, ready for growpush.
inline void construct_spread(int n, const LowDisc& ld, uint64_t k, Rng& g, std::vector<double>& x,
                            std::vector<double>& y) {
  const int family = std::min(kSpreadFamilies - 1,
                              static_cast<int>(ld.coord(k, 0) * kSpreadFamilies));
  const double rot = 2.0 * kPi * ld.coord(k, 1);
  const double aspect = 0.8 + 0.4 * ld.coord(k, 2);        // 0.8 .. 1.2
  const double jit = 0.004 + 0.03 * ld.coord(k, 3);        // jitter sd
  std::vector<double> bx, by;
  spread_base(n, family, g, bx, by);
  const double cs = std::cos(rot), sn = std::sin(rot);
  x.assign(n, 0.0);
  y.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double a = bx[i] * aspect, b = by[i] / aspect;  // anisotropic stretch about center
    double rx = a * cs - b * sn, ry = a * sn + b * cs;  // rotate
    rx += gauss(g, jit);
    ry += gauss(g, jit);
    x[i] = std::min(std::max(rx + 0.5, 0.001), 0.999);
    y[i] = std::min(std::max(ry + 0.5, 0.001), 0.999);
  }
}

// Escape-on-hit kick (ticket 42, spec section 5). A SMALL center perturbation collapses back to
// the same basin (the worker's N=60/N=90 evidence), so to hop to an ADJACENT UNSEEN basin the kick
// must change the CONTACT GRAPH. This removes the k weakest circles (smallest radius) and reinserts
// each at the emptiest spot (the point of largest clearance to all other circles and the walls),
// which is the C++ analogue of the Python defect_move primitive. The caller re-minimizes
// (growpush + radii_lp + polish) and re-fingerprints afterward.
inline void defect_kick(std::vector<double>& x, std::vector<double>& y, int n, Rng& g, int kremove) {
  if (n <= 1 || kremove <= 0) return;
  std::vector<double> r = radii_lp(x, y, n);
  std::vector<int> idx(n);
  for (int i = 0; i < n; ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return r[a] < r[b]; });
  const int k = std::min(kremove, n - 1);
  const int G = 24;  // clearance-scan grid resolution
  for (int t = 0; t < k; ++t) {
    const int rem = idx[t];
    double best_clear = -1e9, bxr = x[rem], byr = y[rem];
    for (int gi = 0; gi < G; ++gi) {
      for (int gj = 0; gj < G; ++gj) {
        const double cx = (gi + 0.5) / G, cy = (gj + 0.5) / G;
        double clear = std::min(std::min(cx, 1.0 - cx), std::min(cy, 1.0 - cy));  // to walls
        for (int i = 0; i < n; ++i) {
          if (i == rem) continue;
          const double dx = cx - x[i], dy = cy - y[i];
          clear = std::min(clear, std::sqrt(dx * dx + dy * dy) - r[i]);
        }
        if (clear > best_clear) { best_clear = clear; bxr = cx; byr = cy; }
      }
    }
    // Nudge off the exact grid node so re-descents from repeated kicks are not identical.
    x[rem] = std::min(std::max(bxr + gauss(g, 0.01), 0.001), 0.999);
    y[rem] = std::min(std::max(byr + gauss(g, 0.01), 0.001), 0.999);
  }
}

}  // namespace csqv

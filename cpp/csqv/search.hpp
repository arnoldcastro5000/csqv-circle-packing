// Construction seeds and the grow-push relaxation (run_scalable.py port).
//
// These are the stochastic front end of one restart: build a starting arrangement, then
// relax overlaps outward with an LP-scaled push, before the L-BFGS polish. The RNG is a
// std::mt19937_64, so a C++ run does not reproduce a numpy run bit-for-bit; that is fine,
// the search is stochastic and the verifier is the score gate.
#pragma once
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "lp.hpp"

namespace csqv {

using Rng = std::mt19937_64;

inline double uni(Rng& g) { return std::uniform_real_distribution<double>(0.0, 1.0)(g); }
inline double gauss(Rng& g, double sd) { return std::normal_distribution<double>(0.0, sd)(g); }

// Build n starting centers of a given kind. Ports run_scalable.py::construct.
inline void construct(int n, const std::string& kind, Rng& g, std::vector<double>& x,
                      std::vector<double>& y) {
  const int s = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n))));
  x.assign(n, 0.0);
  y.assign(n, 0.0);
  auto grid_line = [&](double m) {
    std::vector<double> xs(s);
    for (int t = 0; t < s; ++t) xs[t] = (s == 1) ? m : m + (1.0 - 2.0 * m) * t / (s - 1);
    int c = 0;
    for (int a = 0; a < s && c < n; ++a)
      for (int b = 0; b < s && c < n; ++b) {
        x[c] = xs[b];
        y[c] = xs[a];
        ++c;
      }
  };
  if (kind == "grid") {
    grid_line(0.02 + 0.03 * uni(g));
  } else if (kind == "jitter") {
    grid_line(0.03);  // linspace(0.03,0.97)
    for (int i = 0; i < n; ++i) {
      x[i] += gauss(g, 0.02);
      y[i] += gauss(g, 0.02);
    }
  } else if (kind == "hex") {
    const double dx = 1.0 / s;
    std::vector<double> px, py;
    double yy = dx * 0.6;
    int row = 0;
    while (yy < 1.0 && static_cast<int>(px.size()) < n * 2) {
      const double off = (row % 2) ? dx / 2 : 0.0;
      double xx = dx * 0.6 + off;
      while (xx < 1.0) {
        px.push_back(xx);
        py.push_back(yy);
        xx += dx;
      }
      yy += dx * 0.866;
      ++row;
    }
    if (static_cast<int>(px.size()) >= n) {
      for (int i = 0; i < n; ++i) {
        x[i] = px[i];
        y[i] = py[i];
      }
    } else {
      construct(n, "grid", g, x, y);
      return;
    }
  } else {  // random
    for (int i = 0; i < n; ++i) {
      x[i] = uni(g);
      y[i] = uni(g);
    }
  }
  for (int i = 0; i < n; ++i) {
    x[i] = std::min(std::max(x[i], 0.001), 0.999);
    y[i] = std::min(std::max(y[i], 0.001), 0.999);
  }
}

// LP-scaled overlap relaxation. Ports run_scalable.py::growpush.
inline void growpush(std::vector<double>& x, std::vector<double>& y, int n, Rng& g,
                     int outer = 12, int inner = 15) {
  std::vector<double> fx(n), fy(n);
  for (int o = 0; o < outer; ++o) {
    std::vector<double> r = radii_lp(x, y, n);
    const double scale = 1.0 + 0.05 + 0.05 * uni(g);
    for (int i = 0; i < n; ++i) r[i] *= scale;
    for (int it = 0; it < inner; ++it) {
      std::fill(fx.begin(), fx.end(), 0.0);
      std::fill(fy.begin(), fy.end(), 0.0);
      for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
          const double dx = x[i] - x[j], dy = y[i] - y[j];
          const double d = std::sqrt(dx * dx + dy * dy) + 1e-12;
          const double ov = std::max(0.0, (r[i] + r[j]) - d);
          if (ov <= 0.0) continue;
          const double ux = dx / d, uy = dy / d;
          fx[i] += 0.25 * ov * ux;
          fy[i] += 0.25 * ov * uy;
          fx[j] -= 0.25 * ov * ux;
          fy[j] -= 0.25 * ov * uy;
        }
      }
      for (int i = 0; i < n; ++i) {
        x[i] = std::min(std::max(x[i] + fx[i], r[i]), 1.0 - r[i]);
        y[i] = std::min(std::max(y[i] + fy[i], r[i]), 1.0 - r[i]);
      }
    }
  }
}

}  // namespace csqv

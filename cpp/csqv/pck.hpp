// pck: emit a Packomania .pck from a strictly-valid (shrunk) packing.
//
// This is the shared native emitter used by the worker end-of-run seal (worker.cpp) and by
// the squeeze acceptance harness (squeeze_check.cpp). It mirrors the reference Python emitter so
// a sealed record needs NO Python runtime: the file is strictly feasible AS WRITTEN, which a
// zero-tolerance re-check then confirms on the exact text.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "geometry.hpp"

namespace csqv {

// Emit a Packomania .pck from the strictly-valid (shrunk) packing, in the format from
// packomania.com/hints.html: line 1 = the largest radius (bare number), line 2 = the
// author(s, comma-separated), then one "x y r" line per circle sorted by INCREASING
// radius, at dp decimals, in the centered side-1 square [-0.5, 0.5]^2 (container center at 0,0).
//
// Packomania re-checks at zero tolerance, so the file must be strictly feasible AS
// WRITTEN. Rounding the canonical (already-shrunk) packing to dp decimals can reintroduce a
// ~1e-12 overlap, so this applies the SMALLEST safety shrink that keeps the ROUNDED
// coordinates strictly non-overlapping and contained. The cost to the sum is ~1e-9, far
// below any record margin. Returns the emitted sum of radii, or -1 if it cannot be made
// feasible (never expected).
inline double emit_pck(const std::string& path, const Packing& p, const std::string& author,
                       int dp) {
  const int n = p.n;
  // Format a value to dp decimals and PARSE IT BACK: the parsed value is exactly what a
  // re-check reads from the file, so feasibility must hold on these, not on the in-memory
  // doubles. This makes any precision (12, 15, ...) strictly feasible as written.
  auto fmt = [dp](double v, std::string& s) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", dp, v);
    s = buf;
    return std::strtod(buf, nullptr);
  };
  std::vector<std::string> sx(n), sy(n), sr(n);
  std::vector<double> cx(n), cy(n), cr(n);

  // Target a ROBUST feasibility gap, not merely <= 0. At 15 dp the natural gap sits at the
  // float noise floor (~1e-16), where a re-check with a different summation order could see
  // a positive overlap. A 1e-12 cushion survives that and still costs ~1e-12 on the sum,
  // far below any record margin.
  const double cushion = 1e-12;
  double margin = 0.0;
  const double base = std::pow(10.0, -dp);  // one ULP of the last written decimal
  bool feasible = false;
  for (int attempt = 0; attempt < 16 && !feasible; ++attempt) {
    for (int i = 0; i < n; ++i) {
      cx[i] = fmt(p.x[i] - 0.5, sx[i]);
      cy[i] = fmt(p.y[i] - 0.5, sy[i]);
      cr[i] = fmt(std::max(p.r[i] - margin, 0.0), sr[i]);
    }
    feasible = true;
    for (int i = 0; i < n && feasible; ++i) {
      if (std::fabs(cx[i]) + cr[i] > 0.5 - cushion || std::fabs(cy[i]) + cr[i] > 0.5 - cushion)
        feasible = false;
      for (int j = i + 1; j < n && feasible; ++j) {
        const double dx = cx[i] - cx[j], dy = cy[i] - cy[j];
        if (cr[i] + cr[j] - std::sqrt(dx * dx + dy * dy) > -cushion) feasible = false;
      }
    }
    if (!feasible) margin = (margin < base) ? base : margin * 10.0;
  }
  if (!feasible) return -1.0;

  std::vector<int> idx(n);
  for (int i = 0; i < n; ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return cr[a] < cr[b]; });
  double rmax = 0.0, sum = 0.0;
  for (int i = 0; i < n; ++i) { rmax = std::max(rmax, cr[i]); sum += cr[i]; }

  std::ofstream f(path);
  char rbuf[64];
  std::snprintf(rbuf, sizeof(rbuf), "%.*f", dp, rmax);
  f << rbuf << "\n" << author << "\n";
  for (int i : idx) f << sx[i] << " " << sy[i] << " " << sr[i] << "\n";
  return sum;
}

}  // namespace csqv

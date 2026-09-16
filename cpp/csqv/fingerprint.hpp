// Basin fingerprint: a hash that names the BASIN a converged packing sits in.
//
// Two searches reach the SAME basin when their converged CONTACT GRAPHS match, not when their
// float coordinates match (CONTEXT.md, term Basin). So this hashes the tolerance-thresholded,
// node-colored contact graph through Weisfeiler-Leman (WL) color refinement, never the raw
// floats. See ticket 42 / ADR 0002 and docs/research/phase-2-csqv-basin-dedup-execution.md.
//
// Properties (locked, tested in fingerprint_test.cpp):
//   - STABLE under numerical noise: two descents to one basin (coordinates equal to ~1e-9) hash
//     identically, because only the contact structure and coarse radius buckets enter the hash.
//   - D4-INVARIANT: a rotated or reflected packing has an isomorphic colored graph, so the same
//     hash. Achieved by (a) using the wall-touch COUNT (not which wall) as a node color, which is
//     reflection invariant, and (b) hashing the SORTED multiset of final node colors, which is
//     permutation invariant.
//   - DISCRIMINATING: genuinely different arrangements hash apart. WL never false-splits (same
//     basin always hashes the same); it can rarely false-merge two graphs it cannot separate,
//     which only skips a new basin (the safe direction, the incumbent is already saved).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "geometry.hpp"

namespace csqv {

// splitmix64 finalizer: a fast 64-bit avalanche mix used to combine sub-hashes.
inline uint64_t fp_mix(uint64_t z) {
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
inline uint64_t fp_combine(uint64_t h, uint64_t v) { return fp_mix(h ^ fp_mix(v)); }

// Hash of the basin the packing (x, y, r) sits in. `r` should be the exact-LP radii for the
// centers (radii_lp), so the contacts are the true ones. `tol` separates a real contact (gap
// ~1e-9..1e-12) from a real gap (orders larger); ~1e-6 sits cleanly in that window. `rbucket`
// quantizes radii into coarse buckets for the node color; two descents to one basin share the
// same buckets because their radii agree to ~1e-9.
inline uint64_t basin_fingerprint(const std::vector<double>& x, const std::vector<double>& y,
                                  const std::vector<double>& r, int n, double tol = 1e-6,
                                  double rbucket = 1e-4) {
  if (n <= 0) return 0;

  // Initial node color = (wall-touch count in {0..4}, quantized radius bucket). A circle touches
  // a wall when its gap to that wall (center distance minus radius) is within tol. The COUNT, not
  // which wall, keeps this reflection/rotation invariant.
  std::vector<uint64_t> color(n);
  for (int i = 0; i < n; ++i) {
    int walls = 0;
    if (x[i] - r[i] <= tol) ++walls;
    if (SIDE - x[i] - r[i] <= tol) ++walls;
    if (y[i] - r[i] <= tol) ++walls;
    if (SIDE - y[i] - r[i] <= tol) ++walls;
    const long bucket = std::lround(std::max(r[i], 0.0) / rbucket);
    color[i] = fp_combine(0x1234567ULL + static_cast<uint64_t>(walls),
                          static_cast<uint64_t>(bucket));
  }

  // Contact adjacency: edge (i, j) when dist(c_i, c_j) - (r_i + r_j) <= tol. O(n^2), in line with
  // the rest of the worker.
  std::vector<std::vector<int>> adj(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = x[i] - x[j], dy = y[i] - y[j];
      const double gap = std::sqrt(dx * dx + dy * dy) - (r[i] + r[j]);
      if (gap <= tol) {
        adj[i].push_back(j);
        adj[j].push_back(i);
      }
    }
  }

  // WL refinement: each round, a node's new color mixes its own color with the SORTED multiset of
  // its neighbors' colors (sorting makes it order independent). Iterate until the partition stops
  // refining (distinct-color count stable) or n rounds, whichever comes first.
  auto distinct_count = [](const std::vector<uint64_t>& c) {
    std::vector<uint64_t> t = c;
    std::sort(t.begin(), t.end());
    return static_cast<int>(std::unique(t.begin(), t.end()) - t.begin());
  };
  int prev_distinct = distinct_count(color);
  std::vector<uint64_t> next(n);
  std::vector<uint64_t> nbr;
  for (int round = 0; round < n; ++round) {
    for (int i = 0; i < n; ++i) {
      nbr.clear();
      for (int j : adj[i]) nbr.push_back(color[j]);
      std::sort(nbr.begin(), nbr.end());
      uint64_t h = fp_combine(0xABCDEFULL, color[i]);
      for (uint64_t v : nbr) h = fp_combine(h, v);
      next[i] = h;
    }
    color.swap(next);
    const int d = distinct_count(color);
    if (d == prev_distinct) break;  // stable partition; further rounds cannot refine it
    prev_distinct = d;
  }

  // Fingerprint = hash of the SORTED multiset of final node colors (permutation invariant, so it
  // does not depend on the node ordering and is D4 invariant).
  std::vector<uint64_t> final_colors = color;
  std::sort(final_colors.begin(), final_colors.end());
  uint64_t fp = fp_mix(static_cast<uint64_t>(n));
  for (uint64_t c : final_colors) fp = fp_combine(fp, c);
  return fp;
}

}  // namespace csqv

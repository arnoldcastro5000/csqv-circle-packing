// Unit test for the per-thread tableau buffer of lp.hpp::detail::solve_reduced_lp.
//
// Each solve needs a dense m x (k + m) tableau (up to about 2.5 MB at N=121). A new heap
// block per solve costs a page fault per page on a heap that returns large freed blocks to
// the OS (the Windows heap does; glibc does not by default). So each thread keeps one
// tableau buffer and each solve reuses it. The locked invariants:
//   1) NO LARGE ALLOCATION ON REPEAT: after one radii_lp on grown centers, a second
//      radii_lp on the same centers allocates no heap block of 128 KiB or more. The same
//      holds for a repeated solve with duals (the center polish path).
//   2) NO STALE STATE: a small LP that a thread solves after a large LP gives the same radii
//      and duals, bit for bit, as the same small LP on a fresh thread (a new, empty buffer).
//   3) ONE BUFFER PER THREAD: threads that solve different LPs at the same time each get
//      the bits of the same solve on a fresh thread.
// Build: see cpp/CMakeLists.txt (target csqv_lp_tableau_test) or cpp/Makefile.
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <thread>
#include <vector>

#include "lp.hpp"
#include "search.hpp"

// Counts the heap blocks of at least kLarge bytes. All vectors allocate through here.
static constexpr std::size_t kLarge = 128 * 1024;
static std::atomic<long> g_large_allocs{0};

void* operator new(std::size_t size) {
  if (size >= kLarge) g_large_allocs.fetch_add(1, std::memory_order_relaxed);
  if (void* p = std::malloc(size ? size : 1)) return p;
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

static int g_failures = 0;
#define CHECK(cond, msg)                                             \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                  \
    }                                                                \
  } while (0)

static bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() &&
         (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

// Grown, near-jammed centers: the worker's real LP regime.
static void grown_centers(int n, uint64_t seed, std::vector<double>& x,
                          std::vector<double>& y) {
  csqv::Rng rng(seed * 0x9e3779b97f4a7c15ULL + 1);
  csqv::construct(n, "jitter", rng, x, y);
  csqv::growpush(x, y, n, rng);
}

// The full LP over every prunable pair (no lazy loop): the largest tableau for the centers.
struct FullLp {
  int k = 0;
  std::vector<double> u, rhs;
  std::vector<std::array<int, 2>> pairs;
};

static FullLp full_lp(int n, uint64_t seed) {
  std::vector<double> x, y;
  grown_centers(n, seed, x, y);
  FullLp lp;
  lp.k = n;
  lp.u.resize(n);
  for (int i = 0; i < n; ++i) lp.u[i] = std::max(csqv::wall_slack(x[i], y[i]), 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = x[i] - x[j], dy = y[i] - y[j];
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d < lp.u[i] + lp.u[j]) {
        lp.pairs.push_back({i, j});
        lp.rhs.push_back(d);
      }
    }
  }
  return lp;
}

struct Solution {
  std::vector<double> z, dual, rc;
};

static Solution solve(const FullLp& lp) {
  Solution s;
  s.z = csqv::detail::solve_reduced_lp(lp.k, lp.u, lp.pairs, lp.rhs, &s.dual, &s.rc);
  return s;
}

// The solve on a new thread, so on a new, empty tableau buffer.
static Solution solve_fresh(const FullLp& lp) {
  Solution s;
  std::thread([&] { s = solve(lp); }).join();
  return s;
}

static bool same_solution(const Solution& a, const Solution& b) {
  return same_bits(a.z, b.z) && same_bits(a.dual, b.dual) && same_bits(a.rc, b.rc);
}

int main() {
  const FullLp big = full_lp(121, 7);
  const FullLp small = full_lp(20, 8);
  const FullLp mid = full_lp(60, 9);
  const size_t big_tableau = big.pairs.size() * (big.k + big.pairs.size()) * sizeof(double);
  CHECK(big_tableau >= 4 * kLarge, "the large LP needs a tableau far above the threshold");

  // 1: radii_lp (the lazy loop, several solves per call), then a solve with duals.
  {
    std::vector<double> x, y;
    grown_centers(121, 11, x, y);
    const std::vector<double> r1 = csqv::radii_lp(x, y, 121);
    const long before = g_large_allocs.load();
    const std::vector<double> r2 = csqv::radii_lp(x, y, 121);
    const long large = g_large_allocs.load() - before;
    if (large != 0) std::fprintf(stderr, "  repeated radii_lp: %ld large blocks\n", large);
    CHECK(large == 0, "a repeated radii_lp allocates no large block");
    CHECK(same_bits(r1, r2), "a repeated radii_lp gives the same radii");

    solve(big);
    const long before_dual = g_large_allocs.load();
    solve(big);
    CHECK(g_large_allocs.load() == before_dual,
          "a repeated solve with duals allocates no large block");
  }

  // 2: the buffer holds the large tableau from the solves above; the small solve must not
  // see its old entries.
  {
    const Solution fresh = solve_fresh(small);
    solve(big);
    CHECK(same_solution(solve(small), fresh), "a small LP after a large LP = a fresh thread");
    solve(mid);
    CHECK(same_solution(solve(small), fresh), "a small LP after a mid LP = a fresh thread");
  }

  // 3: threads solve different LPs at the same time, many times each.
  {
    const FullLp* lps[] = {&big, &small, &mid, &small};
    Solution want[4];
    for (int t = 0; t < 4; ++t) want[t] = solve_fresh(*lps[t]);
    std::atomic<int> mismatches{0};
    std::vector<std::thread> threads;
    for (int t = 0; t < 4; ++t) {
      threads.emplace_back([&, t] {
        for (int rep = 0; rep < 6; ++rep) {
          const FullLp& lp = *lps[(t + rep) % 4];
          if (!same_solution(solve(lp), want[(t + rep) % 4])) mismatches.fetch_add(1);
        }
      });
    }
    for (std::thread& th : threads) th.join();
    CHECK(mismatches.load() == 0, "concurrent solves = fresh-thread solves");
  }

  if (g_failures == 0) std::printf("lp_tableau_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}

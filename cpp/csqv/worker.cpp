// CSQV record worker (native, multi-threaded). Port of results/run_scalable.py.
//
// Per restart: varied construction -> grow-push -> scalable_polish (L-BFGS penalty) ->
// center_polish (exact-LP squeeze on a new best) -> keep the best, with perturb-the-best
// restarts for diversity. Each thread is an
// independent worker with its own seed; the best packing per seed persists atomically, so
// a long unattended run never loses a record, and a relaunch resumes from it. A cross of
// the live record prints *** RECORD BEATEN ***.
//
// The trusted Python verifier stays the source of truth: re-check any champion with
//   PYTHONPATH=. python3 cpp/tools/verify_champion.py <this-worker-output.txt>
// before treating it as a record or submitting it.
//
// Usage: csqv_worker <n> <budget_s> <base_seed> [threads] [out_dir]
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <random>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "benchmark.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "lp.hpp"
#include "polish.hpp"
#include "report.hpp"
#include "search.hpp"

namespace {

using csqv::beats_record;
using csqv::record_supplied;
using csqv::save_champion;

using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

std::mutex g_io;  // guards stdout and the global-best file

// Emit a Packomania .pck from the strictly-valid (shrunk) packing, in the format from
// packomania.com/hints.html: line 1 = the largest radius (bare number), line 2 = the
// author(s, comma-separated), then one "x y r" line per circle sorted by INCREASING
// radius, at 12 dp, in the centered side-1 square [-0.5, 0.5]^2 (container center at 0,0).
//
// Packomania re-checks at zero tolerance, so the file must be strictly feasible AS
// WRITTEN. Rounding the canonical (already-shrunk) packing to 12 dp can reintroduce a
// ~1e-12 overlap, so this applies the SMALLEST safety shrink that keeps the ROUNDED
// coordinates strictly non-overlapping and contained. The cost to the sum is ~1e-9, far
// below any record margin. Returns the emitted sum of radii, or -1 if it cannot be made
// feasible (never expected).
double emit_pck(const std::string& path, const csqv::Packing& p, const std::string& author,
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

// Try to resume from a saved centered-frame champion. Returns true and fills best on success.
bool resume(const std::string& path, int n, std::vector<double>& bx, std::vector<double>& by,
            double& best) {
  std::ifstream f(path);
  if (!f) return false;
  std::vector<double> x, y;
  std::string line;
  while (std::getline(f, line)) {
    size_t s = line.find_first_not_of(" \t");
    if (s == std::string::npos || line[s] == '#') continue;
    std::istringstream is(line);
    double cx, cy, cr;
    if (is >> cx >> cy >> cr) {
      x.push_back(cx + 0.5);  // centered -> unit
      y.push_back(cy + 0.5);
    }
  }
  if (static_cast<int>(x.size()) != n) return false;
  std::vector<double> r = csqv::radii_lp(x, y, n);
  double score;
  if (csqv::verify_and_score(x, y, r, n, csqv::TOL, score) && score > best) {
    best = score;
    bx = x;
    by = y;
    return true;
  }
  return false;
}

}  // namespace

struct Shared {
  int n;
  double budget;
  double record;
  std::string out_dir;
  std::string author;
  int pck_dp = 15;
  Clock::time_point t0;
  std::atomic<double> global_best{-1.0};
  std::mutex best_mtx;  // serializes the global-best file write within this process
};

void worker_loop(Shared& sh, int seed) {
  csqv::Rng rng(static_cast<uint64_t>(seed) * 0x9e3779b97f4a7c15ULL + 1);
  const int n = sh.n;
  const char* kinds[] = {"grid", "jitter", "hex", "random"};
  const std::string seed_path = sh.out_dir + "/n" + std::to_string(n) + "-s" +
                                std::to_string(seed) + "-best.txt";

  double best = -1.0;
  std::vector<double> bx, by;
  if (resume(seed_path, n, bx, by, best)) {
    {
      std::lock_guard<std::mutex> lk(g_io);
      if (record_supplied(sh.record))
        std::printf("  [s%d] RESUMED best %.9f gap %+.5f%%\n", seed, best,
                    csqv::gap_percent(sh.record, best));
      else
        std::printf("  [s%d] RESUMED best %.9f\n", seed, best);
      std::fflush(stdout);
    }
    // Refine the resumed champion with the exact-LP center squeeze (monotone; it can only
    // lift the score). This immediately banks the ~1e-5 the penalty polish left behind, so a
    // resumed run starts from a fully-squeezed incumbent.
    std::vector<double> rr;  // center_polish fills this with the exact radii for the result.
    double sc = csqv::center_polish(bx, by, rr, n);
    if (std::isfinite(sc) && sc > best) {
      best = sc;
      save_champion(seed_path, n, seed, bx, by, rr, sc);
    }
    // Seed the GLOBAL best from the resumed champion so n<N>-best.txt and the .pck reflect it
    // even if this run never improves (monotonic, so it only lifts the global, never lowers it).
    if (std::isfinite(sc)) {
      std::lock_guard<std::mutex> lk(sh.best_mtx);
      const std::string gpath = sh.out_dir + "/n" + std::to_string(n) + "-best.txt";
      csqv::Packing canon;
      if (csqv::canonical_feasible(bx, by, rr, n, csqv::TOL, canon) &&
          save_champion(gpath, n, seed, bx, by, rr, sc)) {
        emit_pck(sh.out_dir + "/csqv" + std::to_string(n) + ".pck", canon, sh.author, sh.pck_dp);
        if (sc > sh.global_best.load()) sh.global_best.store(sc);
      }
    }
  }

  std::vector<double> x, y, r;
  long restarts = 0;
  while (elapsed(sh.t0) < sh.budget) {
    if (!bx.empty() && csqv::uni(rng) < 0.6) {
      const double scale = (csqv::uni(rng) < 0.34) ? 0.006 : (csqv::uni(rng) < 0.5 ? 0.012 : 0.03);
      x = bx;
      y = by;
      for (int i = 0; i < n; ++i) {
        x[i] = std::min(std::max(x[i] + csqv::gauss(rng, scale), 0.001), 0.999);
        y[i] = std::min(std::max(y[i] + csqv::gauss(rng, scale), 0.001), 0.999);
      }
    } else {
      csqv::construct(n, kinds[restarts % 4], rng, x, y);
    }
    csqv::growpush(x, y, n, rng);
    r = csqv::radii_lp(x, y, n);
    double s = csqv::scalable_polish(x, y, r, n);
    // A new best earns the exact-LP center squeeze: it recovers the last ~1e-5 of sum the
    // penalty polish leaves unclaimed. It is monotone, so `s` only rises and stays > best.
    if (s > best) s = csqv::center_polish(x, y, r, n);
    ++restarts;
    if (s > best) {
      best = s;
      bx = x;
      by = y;
      save_champion(seed_path, n, seed, x, y, r, s);
      // Update the global best. The mutex serializes threads in THIS process; save_champion's
      // read-before-write makes the global files monotonic across processes too, and the .pck
      // is emitted ONLY when the .txt write actually won, so the two stay consistent.
      {
        std::lock_guard<std::mutex> lk(sh.best_mtx);
        const std::string gpath = sh.out_dir + "/n" + std::to_string(n) + "-best.txt";
        csqv::Packing canon;
        if (csqv::canonical_feasible(x, y, r, n, csqv::TOL, canon) &&
            save_champion(gpath, n, seed, x, y, r, s)) {
          emit_pck(sh.out_dir + "/csqv" + std::to_string(n) + ".pck", canon, sh.author, sh.pck_dp);
          if (s > sh.global_best.load()) sh.global_best.store(s);
          std::lock_guard<std::mutex> lkio(g_io);
          if (record_supplied(sh.record)) {
            const char* cross = beats_record(sh.record, s) ? "  *** RECORD BEATEN ***" : "";
            std::printf("  [%6.0fs s%d r%ld] best %.9f gap %+.5f%%%s\n", elapsed(sh.t0), seed,
                        restarts, s, csqv::gap_percent(sh.record, s), cross);
          } else {
            std::printf("  [%6.0fs s%d r%ld] best %.9f\n", elapsed(sh.t0), seed, restarts, s);
          }
          std::fflush(stdout);
        }
      }
    }
  }
  std::lock_guard<std::mutex> lk(g_io);
  std::printf("  [s%d] done best=%.9f restarts=%ld\n", seed, best, restarts);
  std::fflush(stdout);
}

int main(int argc, char** argv) {
  // The optional `--record <live_sum>` flag is advisory only (console gap + banner). Extract
  // it first, then parse the positionals from what remains.
  std::vector<std::string> args(argv + 1, argv + argc);
  const csqv::RecordArg rec = csqv::parse_record_flag(args);
  if (args.size() < 3) {
    std::fprintf(stderr,
                 "usage: %s <n> <budget_s> <base_seed> [threads] [out_dir] [author] [pck_dp] "
                 "[--record <live_sum>]\n",
                 argv[0]);
    return 2;
  }
  if (rec.malformed)
    std::fprintf(stderr,
                 "warning: --record needs a positive number; ignoring it, no gap will be shown\n");
  Shared sh;
  sh.n = std::atoi(args[0].c_str());
  sh.budget = std::atof(args[1].c_str());
  const int base_seed = std::atoi(args[2].c_str());
  const int threads = (args.size() > 3) ? std::atoi(args[3].c_str()) : 1;
  sh.out_dir = (args.size() > 4) ? args[4] : "results-cpp";
  sh.author = (args.size() > 5) ? args[5] : "Arnold Castro";
  sh.pck_dp = (args.size() > 6) ? std::atoi(args[6].c_str()) : 15;  // .pck decimals; feasible as written
  sh.record = rec.value;  // no_record() when --record is absent/malformed; the search never uses it
  sh.t0 = Clock::now();

  std::error_code ec;
  std::filesystem::create_directories(sh.out_dir, ec);

  {
    std::lock_guard<std::mutex> lk(g_io);
    char rec[48];
    if (record_supplied(sh.record))
      std::snprintf(rec, sizeof(rec), "%.9f", sh.record);
    else
      std::snprintf(rec, sizeof(rec), "not supplied (pass --record for a gap)");
    std::printf("CSQV worker: N=%d budget=%.0fs threads=%d base_seed=%d record=%s out=%s\n",
                sh.n, sh.budget, threads, base_seed, rec, sh.out_dir.c_str());
    std::fflush(stdout);
  }

  // Startup accuracy/precision self-benchmark, saved to the output folder. It lets you
  // cross-validate this build against the trusted Python primitives, and against a run on
  // the other OS, at every start.
  {
    char stamp[24];
    std::time_t nowt = std::time(nullptr);
    std::strftime(stamp, sizeof(stamp), "%Y%m%d-%H%M%S", std::localtime(&nowt));
    const std::string bpath = sh.out_dir + "/benchmark-" + stamp + ".txt";
    bool all_feasible = false;
    const bool wrote = csqv::run_benchmark(bpath, all_feasible);
    std::lock_guard<std::mutex> lk(g_io);
    if (!wrote)
      std::printf("benchmark: FAILED to write %s\n", bpath.c_str());
    else
      std::printf("benchmark: %s  primitives_feasible=%s\n", bpath.c_str(),
                  all_feasible ? "YES" : "NO (WARNING: numeric divergence on this build)");
    std::fflush(stdout);
  }

  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t) pool.emplace_back(worker_loop, std::ref(sh), base_seed + t);
  for (auto& th : pool) th.join();

  std::lock_guard<std::mutex> lk(g_io);
  const double gb = sh.global_best.load();
  if (record_supplied(sh.record))
    std::printf("\nN=%d: global_best=%.9f record=%.9f gap=%+.5f%%\n", sh.n, gb, sh.record,
                csqv::gap_percent(sh.record, gb));
  else
    std::printf("\nN=%d: global_best=%.9f (no --record supplied; verify with "
                "cpp/tools/verify_champion.py)\n",
                sh.n, gb);
  return 0;
}

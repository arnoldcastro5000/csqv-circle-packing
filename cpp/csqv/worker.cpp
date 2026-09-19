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
// for the trusted-frame score, and validate the emitted file as written with
//   python3 cpp/tools/validate_pck.py <this-worker-output.pck>
// before treating it as a record or submitting it.
//
// At the end of a run it prints the pre-squeeze best, runs a one-shot terminal jamming seal
// (jam_slp) on the single global best, writes the sealed packing to the .pck and n<N>-sealed.txt,
// and prints the final sealed value. Pass --no-squeeze to skip the seal.
//
// Usage: csqv_worker <n> <budget_s> <base_seed> [threads] [out_dir] [author] [pck_dp]
//        [--record <live_sum>] [--no-squeeze]
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
#include <unordered_set>
#include <vector>

#include "benchmark.hpp"
#include "fingerprint.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "jam_slp.hpp"
#include "lp.hpp"
#include "pck.hpp"
#include "polish.hpp"
#include "report.hpp"
#include "search.hpp"
#include "spread.hpp"

namespace {

using csqv::beats_record;
using csqv::emit_pck;
using csqv::record_supplied;
using csqv::save_champion;

using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

std::mutex g_io;  // guards stdout and the global-best file

// The native .pck emitter lives in csqv/pck.hpp (shared with the squeeze acceptance harness).

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

  // The unit-frame centers of the current global best, captured on each winning publish (under
  // best_mtx). The end-of-run terminal jamming squeeze (ticket 03) runs ONCE on these, after all
  // threads join, so it seals the single best packing rather than every per-thread champion.
  bool has_global = false;
  int gb_seed = -1;
  std::vector<double> gb_x, gb_y;  // radii are recomputed by jam_slp, so they are not stored
  bool no_squeeze = false;  // --no-squeeze: skip the terminal jamming seal (A/B and speed tests)

  // Basin-dedup + seed-spread mode (ticket 42 / ADR 0002). spread_mode = the new behavior (spread
  // seeds + seeds-only dedup gate + escape + basin counting). count_mode = basin COUNTING ONLY (no
  // spread, no gate, no escape); it exists so the A/B can instrument the current search as a fair
  // control without changing the shipped baseline. With NEITHER flag the worker runs the current
  // independent-parallel search byte-for-byte (no fingerprint, no lock on the hot path).
  bool spread_mode = false;
  bool count_mode = false;
  int base_seed = 0;
  int threads = 1;
  std::unordered_set<uint64_t> seen;  // basin fingerprints found across ALL threads
  std::mutex seen_mtx;                 // guards `seen`
  // Memory cap on `seen` (fingerprints are ~16 B each; the target is N<=100, tens of MB). When the
  // cap is hit the set stops growing and the gate degrades gracefully to plain multistart (no more
  // skips), so a multi-day run cannot climb without bound. See code-review F2.
  long seen_cap = 4000000;

  // Record a basin fingerprint. Returns true if it is NEW (should be kept, not escaped), false if
  // another search already found this basin. When the cap is hit it returns true without inserting,
  // so the gate stops interfering (graceful degradation). `seen.size()` is the distinct-basin count.
  bool record_basin(uint64_t fp) {
    std::lock_guard<std::mutex> lk(seen_mtx);
    if (static_cast<long>(seen.size()) >= seen_cap) return true;
    return seen.insert(fp).second;
  }
  long distinct_basins() {
    std::lock_guard<std::mutex> lk(seen_mtx);
    return static_cast<long>(seen.size());
  }
  bool instrumented() const { return spread_mode || count_mode; }
};

// Publish a champion to the global best file + its .pck, monotonically and mutex-guarded, and
// bump the global_best gauge. Returns true if the global write won. The per-seed file is saved
// by the caller: it is the resume anchor and must persist even when a global publish is skipped
// (canonical_feasible is the stricter deliverable gate; the per-seed champion is kept regardless).
bool publish_global(Shared& sh, int seed, const std::vector<double>& x,
                    const std::vector<double>& y, const std::vector<double>& r, double s) {
  std::lock_guard<std::mutex> lk(sh.best_mtx);
  const std::string gpath = sh.out_dir + "/n" + std::to_string(sh.n) + "-best.txt";
  csqv::Packing canon;
  if (csqv::canonical_feasible(x, y, r, sh.n, csqv::TOL, canon) &&
      save_champion(gpath, sh.n, seed, x, y, r, s)) {
    emit_pck(sh.out_dir + "/csqv" + std::to_string(sh.n) + ".pck", canon, sh.author, sh.pck_dp);
    if (s > sh.global_best.load()) sh.global_best.store(s);
    // Snapshot the winning centers so the end-of-run terminal squeeze can seal them (ticket 03).
    sh.gb_x = x;
    sh.gb_y = y;
    sh.gb_seed = seed;
    sh.has_global = true;
    return true;
  }
  return false;
}

// Console line for a newly published global best, with the record gap and the RECORD BEATEN
// banner when a record was supplied. `pos` is the source token ("r<restarts>" mid-run, "final"
// for the exit squeeze).
void announce(Shared& sh, int seed, const std::string& pos, double s) {
  std::lock_guard<std::mutex> lk(g_io);
  if (record_supplied(sh.record)) {
    const char* cross = beats_record(sh.record, s) ? "  *** RECORD BEATEN ***" : "";
    std::printf("  [%6.0fs s%d %s] best %.9f gap %+.5f%%%s\n", elapsed(sh.t0), seed, pos.c_str(),
                s, csqv::gap_percent(sh.record, s), cross);
  } else {
    std::printf("  [%6.0fs s%d %s] best %.9f\n", elapsed(sh.t0), seed, pos.c_str(), s);
  }
  std::fflush(stdout);
}

void worker_loop(Shared& sh, int seed) {
  csqv::Rng rng(static_cast<uint64_t>(seed) * 0x9e3779b97f4a7c15ULL + 1);
  const int n = sh.n;
  const char* kinds[] = {"grid", "jitter", "hex", "random"};
  const std::string seed_path = sh.out_dir + "/n" + std::to_string(n) + "-s" +
                                std::to_string(seed) + "-best.txt";

  double best = -1.0;
  bool squeezed = false;  // whether the current incumbent (bx, by) has had its center squeeze
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
    double sc = csqv::center_polish_analytic(bx, by, rr, n);
    if (std::isfinite(sc)) {
      squeezed = true;  // bx, by now sit at the local optimum: the incumbent is squeezed.
      if (sc > best) best = sc;
      save_champion(seed_path, n, seed, bx, by, rr, sc);
      // Seed the GLOBAL best from the resumed champion so n<N>-best.txt and the .pck reflect it
      // even if this run never improves (monotonic, so it only lifts the global, never lowers it).
      publish_global(sh, seed, bx, by, rr, sc);
    }
  }

  // Seed-spread state (spread_mode only). Each thread owns a disjoint, evenly-spread slice of the
  // low-discrepancy sequence (stride = thread count), so threads start cold searches in different
  // regions instead of sharing kinds[restarts % 4]. See ticket 42, spec section 1.
  const int t = seed - sh.base_seed;
  const int stride = std::max(1, sh.threads);
  csqv::LowDisc ld(4);
  uint64_t draw = 0;

  std::vector<double> x, y, r;
  long restarts = 0;

  // Bank a descent the instant it improves this thread's best, so no improvement is ever lost, even
  // one that the escape hop later overwrites (code-review F4). Monotone: saves + publishes only on a
  // strict gain. Sets squeezed=false so the terminal squeeze still lifts the champion.
  auto consider = [&](double s, const std::vector<double>& cx, const std::vector<double>& cy,
                      const std::vector<double>& cr) {
    if (s > best) {
      best = s;
      bx = cx;
      by = cy;
      squeezed = false;
      save_champion(seed_path, n, seed, cx, cy, cr, s);
      if (publish_global(sh, seed, cx, cy, cr, s))
        announce(sh, seed, "r" + std::to_string(restarts), s);
    }
  };
  while (elapsed(sh.t0) < sh.budget) {
    // A warm restart perturbs the thread's own best: this is DEPTH (the MBH walk), never gated by
    // the basin dedup. A cold restart builds a fresh construction: this is the SEED, which the
    // dedup gates in spread_mode.
    const bool warm = (!bx.empty() && csqv::uni(rng) < 0.6);
    if (warm) {
      const double scale = (csqv::uni(rng) < 0.34) ? 0.006 : (csqv::uni(rng) < 0.5 ? 0.012 : 0.03);
      x = bx;
      y = by;
      for (int i = 0; i < n; ++i) {
        x[i] = std::min(std::max(x[i] + csqv::gauss(rng, scale), 0.001), 0.999);
        y[i] = std::min(std::max(y[i] + csqv::gauss(rng, scale), 0.001), 0.999);
      }
    } else if (sh.spread_mode) {
      csqv::construct_spread(n, ld, static_cast<uint64_t>(t) + draw * stride, rng, x, y);
      ++draw;
    } else {
      csqv::construct(n, kinds[restarts % 4], rng, x, y);
    }
    csqv::growpush(x, y, n, rng);
    r = csqv::radii_lp(x, y, n);
    double s = csqv::scalable_polish(x, y, r, n);
    // NO mid-run center squeeze. center_polish costs up to max_iter x 4n exact-LP solves, and a
    // warm mid-N push produces frequent SMALL (refinement) new bests, so squeezing each one
    // stalls the restart rate (the user-observed slowdown at N=123/124). The squeeze is monotone,
    // and the final squeeze at loop exit plus the resume squeeze on relaunch fully squeeze the
    // delivered champion, so a mid-run squeeze would change only intermediate on-disk state, never
    // the delivered per-seed file or .pck. So defer the whole squeeze to the end. See ticket 35.
    ++restarts;

    // Bank any improvement from THIS descent before the dedup can overwrite it (F4).
    consider(s, x, y, r);

    // Basin bookkeeping, ONLY when instrumented (spread or count). The plain baseline (no flag)
    // stays byte-identical, with no fingerprint and no lock on the hot path (F1). Fingerprint every
    // converged optimum and record it, so `distinct_basins` is the coverage metric in both A/B arms
    // (--count vs --spread). In spread_mode ONLY, a COLD seed that re-descends an already-seen basin
    // escapes: hop with a combinatorial defect kick (a small perturbation collapses back) up to 3
    // times, then the next iteration draws a fresh spread seed. Warm (DEPTH) restarts are never
    // gated, and every hop is still banked by consider().
    if (sh.instrumented()) {
      uint64_t fp = csqv::basin_fingerprint(x, y, r, n);
      bool was_new = sh.record_basin(fp);
      if (sh.spread_mode && !warm && !was_new) {
        const int kremove = std::max(1, n / 25);
        for (int attempt = 0; attempt < 3 && !was_new; ++attempt) {
          csqv::defect_kick(x, y, n, rng, kremove);
          csqv::growpush(x, y, n, rng);
          r = csqv::radii_lp(x, y, n);
          s = csqv::scalable_polish(x, y, r, n);
          consider(s, x, y, r);  // a hop can find a new best too; never lose it
          fp = csqv::basin_fingerprint(x, y, r, n);
          was_new = sh.record_basin(fp);
        }
      }
    }
  }

  // Final squeeze on clean exit (budget reached). The worker never squeezes mid-run, so this
  // thread's champion is saved un-squeezed; squeeze it once here so the delivered per-seed file,
  // the global best, and the .pck are fully squeezed. It is monotone, so it only lifts. Skip it
  // when the incumbent is already squeezed (a resumed champion that saw no new best), which
  // spares a redundant O(4n) ascent. A hard-killed run (SIGTERM/taskkill) skips this, but the
  // resume squeeze recovers it on the next launch. See tickets 34 and 35.
  if (!bx.empty() && !squeezed) {
    std::vector<double> rr;
    const double sc = csqv::center_polish_analytic(bx, by, rr, n);
    if (std::isfinite(sc) && sc > best) {
      best = sc;
      save_champion(seed_path, n, seed, bx, by, rr, sc);
      if (publish_global(sh, seed, bx, by, rr, sc)) announce(sh, seed, "final", sc);
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
  // --spread turns on the basin-dedup + seed-spread path (ticket 42). --count turns on basin
  // COUNTING ONLY (the A/B control that instruments the current search without the gate/escape).
  // Absent both = the byte-identical production baseline. Parse and strip before the positionals.
  bool spread_mode = false, count_mode = false, no_squeeze = false;
  for (auto it = args.begin(); it != args.end();) {
    if (*it == "--spread") {
      spread_mode = true;
      it = args.erase(it);
    } else if (*it == "--count") {
      count_mode = true;
      it = args.erase(it);
    } else if (*it == "--no-squeeze") {
      no_squeeze = true;
      it = args.erase(it);
    } else {
      ++it;
    }
  }
  if (args.size() < 3) {
    std::fprintf(stderr,
                 "usage: %s <n> <budget_s> <base_seed> [threads] [out_dir] [author] [pck_dp] "
                 "[--record <live_sum>] [--no-squeeze]\n",
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
  sh.spread_mode = spread_mode;
  sh.count_mode = count_mode;
  sh.no_squeeze = no_squeeze;
  sh.base_seed = base_seed;
  sh.threads = threads;
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
    const char* mode = sh.spread_mode ? "spread" : (sh.count_mode ? "count" : "baseline");
    std::printf("CSQV worker: N=%d budget=%.0fs threads=%d base_seed=%d record=%s mode=%s out=%s\n",
                sh.n, sh.budget, threads, base_seed, rec, mode, sh.out_dir.c_str());
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

  // Pre-squeeze best: the value the search reached before the terminal seal.
  const double pre = sh.global_best.load();
  {
    std::lock_guard<std::mutex> lk(g_io);
    if (record_supplied(sh.record))
      std::printf("\nN=%d: pre-squeeze best=%.12f record=%.9f gap=%+.5f%%\n", sh.n, pre, sh.record,
                  csqv::gap_percent(sh.record, pre));
    else
      std::printf("\nN=%d: pre-squeeze best=%.12f\n", sh.n, pre);
    std::fflush(stdout);
  }

  // Terminal jamming seal (ticket 03): run the SLP jammer ONCE on the global best, after the
  // joins. The per-thread final step is only the first-order analytic polish; this drives the
  // single best packing to its jammed optimum. It is monotone (keep-better), so it never
  // regresses. The SEALED packing goes to the .pck (the submittable artifact) and to a new
  // n<N>-sealed.txt; the raw n<N>-best.txt is left as the search's own champion. The
  // Python verifier stays the source of truth: re-check the .pck before any submission.
  // final_val is the TOUCHING score (verify_and_score), the same frame the whole run reports
  // (announce / publish_global / the pre-squeeze best). The emitted .pck is shrunk below this by
  // emit's feasibility cushion (~1e-13 at 15 dp), so validate_pck.py reads a hair less; that is
  // the submittable value and the user re-checks it before any submit. Keeping final_val in the
  // touching frame keeps the final line consistent with every mid-run line.
  double final_val = pre;
  bool sealed = false;
  bool seal_write_failed = false;  // gain found, but the .pck could not be rewritten feasibly
  if (!sh.no_squeeze && sh.has_global && sh.n >= 2) {
    std::vector<double> gx, gy;
    int gseed;
    {
      std::lock_guard<std::mutex> lk(sh.best_mtx);
      gx = sh.gb_x;
      gy = sh.gb_y;
      gseed = sh.gb_seed;
    }
    const csqv::JamResult jr = csqv::jam_slp(gx, gy, sh.n);
    if (jr.score > pre) {
      csqv::Packing canon;
      const std::string pck = sh.out_dir + "/csqv" + std::to_string(sh.n) + ".pck";
      // Claim a seal ONLY if the .pck was actually rewritten strictly feasible. emit_pck returns
      // -1 (WITHOUT writing) when it cannot round the packing feasibly; ignoring that would leave
      // the lower pre-squeeze .pck on disk while the run reports (and may claim a record on) the
      // higher unsaved value. canonical_feasible failing is the same false-success case.
      if (csqv::canonical_feasible(jr.x, jr.y, jr.r, sh.n, csqv::TOL, canon) &&
          emit_pck(pck, canon, sh.author, sh.pck_dp) >= 0.0) {
        save_champion(sh.out_dir + "/n" + std::to_string(sh.n) + "-sealed.txt", sh.n, gseed, jr.x,
                      jr.y, jr.r, jr.score);
        sh.global_best.store(jr.score);
        final_val = jr.score;
        sealed = true;
      } else {
        seal_write_failed = true;
      }
    }
  }

  std::lock_guard<std::mutex> lk(g_io);
  if (sh.no_squeeze)
    std::printf("N=%d: terminal squeeze SKIPPED (--no-squeeze)\n", sh.n);
  else if (sealed)
    std::printf("N=%d: terminal squeeze sealed %.12f -> %.12f (delta %+.3e); .pck + "
                "n%d-sealed.txt updated\n",
                sh.n, pre, final_val, final_val - pre, sh.n);
  else if (seal_write_failed)
    std::printf("N=%d: terminal squeeze found gain but could NOT write a feasible .pck; the .pck "
                "still holds the pre-squeeze best %.12f (not sealed)\n",
                sh.n, pre);
  else
    std::printf("N=%d: terminal squeeze: no gain over %.12f (already jammed)\n", sh.n, pre);

  if (record_supplied(sh.record)) {
    const char* cross = beats_record(sh.record, final_val) ? "  *** RECORD BEATEN ***" : "";
    std::printf("N=%d: final=%.9f record=%.9f gap=%+.5f%%%s\n", sh.n, final_val, sh.record,
                csqv::gap_percent(sh.record, final_val), cross);
  } else {
    std::printf("N=%d: final=%.9f (no --record supplied; score with verify_champion.py, "
                "validate the .pck with validate_pck.py)\n",
                sh.n, final_val);
  }
  if (sh.instrumented())
    std::printf("N=%d: distinct_basins=%ld spread=%s (the A/B coverage metric)\n", sh.n,
                sh.distinct_basins(), sh.spread_mode ? "on" : "off");
  return 0;
}

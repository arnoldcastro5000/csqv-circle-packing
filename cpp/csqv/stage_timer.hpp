// Stage timers for the worker re-profile: where the time of a restart goes.
//
// A worker built with -DCSQV_STAGE_TIMERS (make -C cpp stages, which writes
// build/bench/csqv_worker_stages) adds the wall time of each stage to per-thread totals and
// prints a table at exit. Without the macro, each CSQV_STAGE_* macro expands to nothing, so
// the production worker compiles to the same code.
//
// A mark gives the time since the previous mark of the same thread to one stage. radii_lp
// adds its own time and one call to a per-thread LP clock, so a mark also gives the LP time
// and the LP calls of the interval to its stage. The table shows each restart stage as a part
// of the restart loop, and the radii_lp part of each stage.
//
// The clock calls cost about 50 ns per mark and per radii_lp call, far below a restart
// (milliseconds). The stage worker follows the same search path as the production worker for
// a seed (the timers do not touch the arithmetic), but it runs a little slower.
#pragma once
#include <array>
#include <cstdio>

#ifdef CSQV_STAGE_TIMERS
#include <chrono>
#include <mutex>
#endif

namespace csqv::stages {

// Startup, final polish and seal run once per thread or per run; the others run per restart.
enum Stage {
  kStartup,
  kSeed,
  kGrowpush,
  kRadiiLp,
  kPolish,
  kBank,
  kBasin,
  kEscape,
  kFinalPolish,
  kSeal,
  kCount
};

inline const char* stage_name(int s) {
  static const char* kNames[kCount] = {"startup", "seed",   "growpush", "radii_lp",     "polish",
                                       "bank",    "basin",  "escape",   "final_polish", "seal"};
  return kNames[s];
}

inline bool in_restart_loop(int s) { return s >= kSeed && s <= kEscape; }

// The totals of one thread or of a whole run, per stage.
struct Totals {
  std::array<double, kCount> sec{};     // wall time
  std::array<double, kCount> lp_sec{};  // the part of sec spent in radii_lp
  std::array<long, kCount> marks{};
  std::array<long, kCount> lp_calls{};

  void add(const Totals& o) {
    for (int s = 0; s < kCount; ++s) {
      sec[s] += o.sec[s];
      lp_sec[s] += o.lp_sec[s];
      marks[s] += o.marks[s];
      lp_calls[s] += o.lp_calls[s];
    }
  }
};

// Gives the time between marks to stages. The caller supplies the times (seconds from any
// fixed origin) and the LP clock readings, so this class does no clock calls itself.
class Recorder {
 public:
  void start(double now, double lp_now, long lp_calls_now) {
    started_ = true;
    last_ = now;
    last_lp_ = lp_now;
    last_calls_ = lp_calls_now;
  }
  void mark(Stage s, double now, double lp_now, long lp_calls_now) {
    if (!started_) return;
    t_.sec[s] += now - last_;
    t_.lp_sec[s] += lp_now - last_lp_;
    t_.lp_calls[s] += lp_calls_now - last_calls_;
    ++t_.marks[s];
    start(now, lp_now, lp_calls_now);
  }
  const Totals& totals() const { return t_; }

 private:
  bool started_ = false;
  double last_ = 0.0, last_lp_ = 0.0;
  long last_calls_ = 0;
  Totals t_;
};

// One row per stage that has time or marks. loop% is the part of the restart loop (the sum
// of the restart stages); a one-shot stage shows "-" for loop% and per restart.
inline void print_table(std::FILE* f, const Totals& t, long restarts, double budget) {
  double loop = 0.0;
  for (int s = 0; s < kCount; ++s)
    if (in_restart_loop(s)) loop += t.sec[s];
  std::fprintf(f, "stage times (sum over threads): restarts=%ld budget=%.3fs restarts/s=%.4f\n",
               restarts, budget, budget > 0.0 ? restarts / budget : 0.0);
  std::fprintf(f, "  %-12s %9s %7s %11s %10s %9s\n", "stage", "seconds", "loop%", "per restart",
               "radii_lp%", "lp/rst");
  for (int s = 0; s < kCount; ++s) {
    if (t.sec[s] <= 0.0 && t.marks[s] == 0) continue;
    char part[16] = "-", per[24] = "-", calls[16] = "-";
    if (in_restart_loop(s) && loop > 0.0)
      std::snprintf(part, sizeof(part), "%6.1f%%", 100.0 * t.sec[s] / loop);
    if (in_restart_loop(s) && restarts > 0) {
      std::snprintf(per, sizeof(per), "%8.3f ms", 1e3 * t.sec[s] / restarts);
      std::snprintf(calls, sizeof(calls), "%8.1f", static_cast<double>(t.lp_calls[s]) / restarts);
    }
    const double lp_part = t.sec[s] > 0.0 ? 100.0 * t.lp_sec[s] / t.sec[s] : 0.0;
    std::fprintf(f, "  %-12s %9.3f %7s %11s %9.1f%% %9s\n", stage_name(s), t.sec[s], part, per,
                 lp_part, calls);
  }
  std::fprintf(f, "  %-12s %9.3f %7s\n", "loop", loop, "100.0%");
}

#ifdef CSQV_STAGE_TIMERS

inline double now_sec() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// The radii_lp time and calls of this thread.
struct LpClock {
  double sec = 0.0;
  long calls = 0;
};
inline LpClock& lp_clock() {
  thread_local LpClock c;
  return c;
}

// Adds the time of one radii_lp call to the LP clock (CSQV_STAGE_LP_SCOPE in radii_lp).
class LpScope {
 public:
  LpScope() : t0_(now_sec()) {}
  ~LpScope() {
    LpClock& c = lp_clock();
    c.sec += now_sec() - t0_;
    ++c.calls;
  }
  LpScope(const LpScope&) = delete;
  LpScope& operator=(const LpScope&) = delete;

 private:
  double t0_;
};

inline Recorder& recorder() {
  thread_local Recorder r;
  return r;
}
inline Totals& run_totals() {
  static Totals t;
  return t;
}
inline std::mutex& run_mutex() {
  static std::mutex m;
  return m;
}

inline void start_now() {
  const LpClock& c = lp_clock();
  recorder().start(now_sec(), c.sec, c.calls);
}
inline void mark_now(Stage s) {
  const LpClock& c = lp_clock();
  recorder().mark(s, now_sec(), c.sec, c.calls);
}
// Adds this thread's totals to the run totals and clears them.
inline void flush() {
  std::lock_guard<std::mutex> lk(run_mutex());
  run_totals().add(recorder().totals());
  recorder() = Recorder{};
}
inline void report(double budget) {
  std::lock_guard<std::mutex> lk(run_mutex());
  print_table(stdout, run_totals(), run_totals().marks[kPolish], budget);
  std::fflush(stdout);
}

#define CSQV_STAGE_START() ::csqv::stages::start_now()
#define CSQV_STAGE_MARK(s) ::csqv::stages::mark_now(::csqv::stages::s)
#define CSQV_STAGE_FLUSH() ::csqv::stages::flush()
#define CSQV_STAGE_REPORT(budget) ::csqv::stages::report(budget)
#define CSQV_STAGE_LP_SCOPE() const ::csqv::stages::LpScope csqv_stage_lp_scope_

#else

#define CSQV_STAGE_START() ((void)0)
#define CSQV_STAGE_MARK(s) ((void)0)
#define CSQV_STAGE_FLUSH() ((void)0)
#define CSQV_STAGE_REPORT(budget) ((void)0)
#define CSQV_STAGE_LP_SCOPE() ((void)0)

#endif

}  // namespace csqv::stages

// Unit test for the stage timers of csqv/stage_timer.hpp (the worker re-profile).
//
// The locked invariants:
//   1) ATTRIBUTION: a mark gives the time since the previous mark of the recorder to its
//      stage, and the radii_lp time and calls in the same interval to the same stage.
//   2) START: a mark before start records nothing; start again skips the time between.
//   3) MERGE: Totals::add sums each stage; the per-thread recorders of a run merge into one
//      run total, with no lost marks.
//   4) TABLE: the report gives each restart stage its part of the restart loop, the time
//      per restart, and the radii_lp part of the stage.
//   5) LP CLOCK: each radii_lp call adds one call and a time >= 0 to the thread's LP clock.
// The test builds with -DCSQV_STAGE_TIMERS (cpp/Makefile), as the stage worker does.
// Build: see cpp/CMakeLists.txt (target csqv_stage_timer_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

#include "lp.hpp"
#include "search.hpp"
#include "stage_timer.hpp"

#ifndef CSQV_STAGE_TIMERS
#error "stage_timer_test needs -DCSQV_STAGE_TIMERS"
#endif

namespace st = csqv::stages;

static int g_failures = 0;
#define CHECK(cond, msg)                                             \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                  \
    }                                                                \
  } while (0)

static bool near(double a, double b) { return std::fabs(a - b) < 1e-12; }

static void test_attribution() {
  st::Recorder rec;
  rec.start(10.0, 1.0, 5);
  rec.mark(st::kSeed, 10.5, 1.0, 5);
  rec.mark(st::kGrowpush, 13.5, 3.0, 17);
  rec.mark(st::kRadiiLp, 13.75, 3.2, 18);
  rec.mark(st::kGrowpush, 14.75, 4.0, 30);
  const st::Totals& t = rec.totals();
  CHECK(near(t.sec[st::kSeed], 0.5), "seed time");
  CHECK(near(t.sec[st::kGrowpush], 4.0), "growpush time sums over two marks");
  CHECK(near(t.sec[st::kRadiiLp], 0.25), "radii_lp stage time");
  CHECK(t.marks[st::kGrowpush] == 2 && t.marks[st::kSeed] == 1, "mark counts");
  CHECK(near(t.lp_sec[st::kGrowpush], 2.8), "growpush lp time");
  CHECK(near(t.lp_sec[st::kRadiiLp], 0.2), "radii_lp stage lp time");
  CHECK(t.lp_calls[st::kGrowpush] == 24 && t.lp_calls[st::kRadiiLp] == 1, "lp calls");
  CHECK(t.lp_calls[st::kSeed] == 0 && near(t.lp_sec[st::kSeed], 0.0), "no lp in seed");
}

static void test_start() {
  st::Recorder rec;
  rec.mark(st::kSeed, 5.0, 0.0, 0);
  CHECK(rec.totals().marks[st::kSeed] == 0, "a mark before start records nothing");
  rec.start(6.0, 0.0, 0);
  rec.mark(st::kPolish, 7.0, 0.0, 0);
  rec.start(100.0, 50.0, 9);  // skip the wait (for example a join)
  rec.mark(st::kSeal, 101.0, 50.5, 10);
  const st::Totals& t = rec.totals();
  CHECK(near(t.sec[st::kPolish], 1.0), "polish time");
  CHECK(near(t.sec[st::kSeal], 1.0), "start again skips the gap");
  CHECK(near(t.lp_sec[st::kSeal], 0.5) && t.lp_calls[st::kSeal] == 1, "seal lp from restart");
}

static void test_merge() {
  st::Totals a, b;
  a.sec[st::kBank] = 1.5;
  a.marks[st::kBank] = 3;
  b.sec[st::kBank] = 0.5;
  b.marks[st::kBank] = 1;
  b.lp_calls[st::kBank] = 7;
  a.add(b);
  CHECK(near(a.sec[st::kBank], 2.0) && a.marks[st::kBank] == 4 && a.lp_calls[st::kBank] == 7,
        "Totals::add sums");

  // Real threads through the macro layer: each marks 1000 restarts, then flushes.
  st::run_totals() = st::Totals{};
  const int kThreads = 4, kMarks = 1000;
  std::vector<std::thread> pool;
  for (int t = 0; t < kThreads; ++t)
    pool.emplace_back([] {
      CSQV_STAGE_START();
      for (int i = 0; i < kMarks; ++i) {
        CSQV_STAGE_MARK(kGrowpush);
        CSQV_STAGE_MARK(kPolish);
      }
      CSQV_STAGE_FLUSH();
    });
  for (auto& th : pool) th.join();
  const st::Totals& r = st::run_totals();
  CHECK(r.marks[st::kGrowpush] == kThreads * kMarks, "no lost growpush marks");
  CHECK(r.marks[st::kPolish] == kThreads * kMarks, "no lost polish marks");
  CHECK(r.sec[st::kGrowpush] >= 0.0 && r.sec[st::kPolish] >= 0.0, "times are not negative");
}

static std::string table_text(const st::Totals& t, long restarts, double wall) {
  std::FILE* f = std::tmpfile();
  st::print_table(f, t, restarts, wall);
  std::rewind(f);
  std::string s;
  char buf[256];
  while (std::fgets(buf, sizeof(buf), f)) s += buf;
  std::fclose(f);
  return s;
}

static bool has_line(const std::string& text, const std::string& start, const std::string& part) {
  size_t p = 0;
  while ((p = text.find(start, p)) != std::string::npos) {
    const size_t e = text.find('\n', p);
    if (text.substr(p, e - p).find(part) != std::string::npos) return true;
    p = e;
  }
  return false;
}

static void test_table() {
  st::Totals t;
  t.sec[st::kSeed] = 0.5;
  t.sec[st::kGrowpush] = 6.0;
  t.lp_sec[st::kGrowpush] = 4.5;
  t.lp_calls[st::kGrowpush] = 1200;
  t.sec[st::kRadiiLp] = 0.5;
  t.sec[st::kPolish] = 3.0;
  t.marks[st::kPolish] = 100;
  t.sec[st::kSeal] = 2.0;  // a one-shot stage: not in the restart loop
  const std::string s = table_text(t, 100, 12.0);
  // Loop = 0.5 + 6 + 0.5 + 3 = 10 s over 100 restarts.
  CHECK(has_line(s, "  growpush", " 60.0%"), "growpush part of the loop");
  CHECK(has_line(s, "  growpush", " 60.000 ms"), "growpush ms per restart");
  CHECK(has_line(s, "  growpush", " 75.0%"), "radii_lp part of growpush");
  CHECK(has_line(s, "  growpush", " 12.0"), "radii_lp calls per restart");
  CHECK(has_line(s, "  seal", "   -"), "a one-shot stage has no loop part");
  CHECK(s.find("restarts=100") != std::string::npos, "restart count");
  CHECK(s.find("restarts/s=8.3333") != std::string::npos, "restarts per budget second");
  CHECK(s.find("  escape") == std::string::npos, "a stage with no time is not a row");
  if (g_failures) std::fputs(s.c_str(), stderr);
}

static void test_lp_clock() {
  std::vector<double> x, y;
  csqv::Rng rng(7);
  csqv::construct(30, "jitter", rng, x, y);
  const st::LpClock before = st::lp_clock();
  const std::vector<double> r = csqv::radii_lp(x, y, 30);
  const st::LpClock after = st::lp_clock();
  CHECK(r.size() == 30, "radii_lp result");
  CHECK(after.calls == before.calls + 1, "one radii_lp call adds one call");
  CHECK(after.sec >= before.sec, "radii_lp adds a time >= 0");
}

int main() {
  test_attribution();
  test_start();
  test_merge();
  test_table();
  test_lp_clock();
  if (g_failures) {
    std::fprintf(stderr, "stage_timer_test: %d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("stage_timer_test: all checks passed\n");
  return 0;
}

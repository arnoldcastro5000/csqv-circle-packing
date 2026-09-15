// Unit test for report.hpp (operator-supplied --record flag, advisory gap logic).
// Build: see cpp/CMakeLists.txt (target csqv_report_test) or cpp/Makefile (make test).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "report.hpp"

static int g_failures = 0;
#define CHECK(cond, msg)                                              \
  do {                                                               \
    if (!(cond)) {                                                   \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__); \
      ++g_failures;                                                 \
    }                                                               \
  } while (0)

static bool near(double a, double b) { return std::fabs(a - b) < 1e-9; }

int main() {
  // record_supplied: only a finite positive number counts.
  CHECK(!csqv::record_supplied(csqv::no_record()), "NaN is not a record");
  CHECK(!csqv::record_supplied(0.0), "zero is not a record");
  CHECK(!csqv::record_supplied(-1.0), "negative is not a record");
  CHECK(csqv::record_supplied(5.797468812113), "a positive value is a record");

  // parse_record_flag: `--record <value>` form is extracted, positionals preserved.
  {
    std::vector<std::string> a{"121", "100", "0", "--record", "5.797468812113"};
    const csqv::RecordArg r = csqv::parse_record_flag(a);
    CHECK(near(r.value, 5.797468812113) && !r.malformed, "record parsed from --record <value>");
    CHECK((a == std::vector<std::string>{"121", "100", "0"}), "positionals preserved");
  }
  // `--record=<value>` form, flag in the middle.
  {
    std::vector<std::string> a{"121", "--record=5.5", "100", "0"};
    const csqv::RecordArg r = csqv::parse_record_flag(a);
    CHECK(near(r.value, 5.5) && !r.malformed, "record parsed from --record=<value>");
    CHECK((a == std::vector<std::string>{"121", "100", "0"}), "middle flag removed cleanly");
  }
  // Absent flag: no record, not malformed, args untouched.
  {
    std::vector<std::string> a{"121", "100", "0"};
    const csqv::RecordArg r = csqv::parse_record_flag(a);
    CHECK(!csqv::record_supplied(r.value) && !r.malformed, "no flag means no record, not malformed");
    CHECK((a == std::vector<std::string>{"121", "100", "0"}), "args untouched without a flag");
  }
  // Trailing `--record` with no value: dropped, flagged malformed (a typo, not "absent").
  {
    std::vector<std::string> a{"121", "100", "0", "--record"};
    const csqv::RecordArg r = csqv::parse_record_flag(a);
    CHECK(!csqv::record_supplied(r.value) && r.malformed, "trailing --record is malformed");
    CHECK((a == std::vector<std::string>{"121", "100", "0"}), "trailing flag removed");
  }
  // Malformed values are flagged, never silently accepted or silently dropped.
  for (const char* bad : {"foo", "5.8x", "5,8", "-1", "0", "", "nan", "inf"}) {
    std::vector<std::string> a{"121", "100", "0", "--record", bad};
    const csqv::RecordArg r = csqv::parse_record_flag(a);
    CHECK(!csqv::record_supplied(r.value) && r.malformed,
          "malformed --record value is flagged malformed");
    CHECK((a == std::vector<std::string>{"121", "100", "0"}), "malformed flag+value removed");
  }
  // A valid value with surrounding whitespace still parses (strtod skips leading spaces).
  {
    std::vector<std::string> a{"121", "100", "0", "--record", "5.8"};
    const csqv::RecordArg r = csqv::parse_record_flag(a);
    CHECK(near(r.value, 5.8) && !r.malformed, "clean value parses");
  }

  // gap_percent: sign and magnitude.
  CHECK(near(csqv::gap_percent(10.0, 10.0), 0.0), "gap is 0 at the record");
  CHECK(near(csqv::gap_percent(10.0, 5.0), 50.0), "gap is +50% below the record");
  CHECK(near(csqv::gap_percent(10.0, 20.0), -100.0), "gap is negative past the record");

  // beats_record: needs a supplied record AND a clear margin.
  CHECK(csqv::beats_record(5.797468812113, 5.799074421369), "clears the record by a margin");
  CHECK(!csqv::beats_record(5.797468812113, 5.797468812113), "a tie does not beat");
  CHECK(!csqv::beats_record(5.797468812113, 5.797468812500), "within 1e-6 does not beat");
  CHECK(!csqv::beats_record(csqv::no_record(), 999.0), "no record means nothing beats it");

  if (g_failures == 0) std::printf("csqv_report_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}

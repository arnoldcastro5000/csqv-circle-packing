// Advisory record reporting for the CSQV worker.
//
// The search NEVER uses the record: it maximizes the sum of radii regardless. The record
// is only a reference line for the operator watching the console (the gap and the
// *** RECORD BEATEN *** banner). Packomania records move daily, so the worker does NOT
// hardcode one: the operator supplies the live value at launch with `--record <sum>`, and
// when it is absent the worker simply omits the gap rather than compare against a stale
// constant. The Python verifier (cpp/tools/verify_champion.py) stays the authority on
// whether a packing actually beats the live record.
#pragma once
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace csqv {

// Sentinel for "no record supplied". A supplied record must be finite and positive.
inline double no_record() { return std::numeric_limits<double>::quiet_NaN(); }
inline bool record_supplied(double record) { return std::isfinite(record) && record > 0.0; }

// Result of scanning args for the optional --record flag.
struct RecordArg {
  double value = no_record();  // the record to use, or no_record() when none is usable
  bool malformed = false;      // --record was present but its value was not a positive number
};

// Extract an optional `--record <value>` (or `--record=<value>`) flag from args, IN PLACE:
// the flag and its value are erased, and the positional arguments left in order. Returns the
// parsed record. When the flag is present but its value is empty, non-numeric, or not a
// positive finite number, `malformed` is set so the caller can WARN rather than silently
// ignore a typo (`--record 5,8` or `--record foo` must not masquerade as "no record").
inline RecordArg parse_record_flag(std::vector<std::string>& args) {
  RecordArg out;
  bool present = false;
  bool have_raw = false;
  std::string raw;
  for (std::size_t i = 0; i < args.size();) {
    const std::string& a = args[i];
    if (a == "--record") {
      present = true;
      if (i + 1 < args.size()) {
        raw = args[i + 1];
        have_raw = true;
        args.erase(args.begin() + i, args.begin() + i + 2);
      } else {
        args.erase(args.begin() + i);  // trailing flag with no value
      }
    } else if (a.rfind("--record=", 0) == 0) {
      present = true;
      raw = a.substr(9);
      have_raw = true;
      args.erase(args.begin() + i);
    } else {
      ++i;
    }
  }
  if (!present) return out;  // absent: value=no_record(), not malformed
  if (have_raw) {
    const char* s = raw.c_str();
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    // The WHOLE token must parse (no trailing junk) and be a positive finite number.
    if (end != s && *end == '\0' && record_supplied(v)) {
      out.value = v;
      return out;
    }
  }
  out.malformed = true;
  return out;
}

// Gap of `sum` below `record`, in percent (negative once the sum passes the record).
// Only meaningful when record_supplied(record).
inline double gap_percent(double record, double sum) {
  return 100.0 * (record - sum) / record;
}

// True only when a record was supplied AND `sum` clears it by more than the noise margin.
inline bool beats_record(double record, double sum) {
  return record_supplied(record) && sum > record + 1e-6;
}

}  // namespace csqv

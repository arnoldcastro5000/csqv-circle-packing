// Champion file I/O for the CSQV worker: atomic + monotonic saves.
//
// A champion file must never regress, even when several worker PROCESSES share one
// out_dir (the in-process best is invisible across processes). Each save reads the sum
// already on disk and writes only if it strictly beats it, into a per-writer temp that is
// then atomically REPLACED over the target. The temp name is unique per writer, so
// concurrent writers never share a `.tmp`.
//
// Cross-platform atomic replace: POSIX `rename` replaces the destination atomically, but
// the Windows CRT `rename`/`std::rename` FAILS when the destination exists. That freezes
// the champion file after its first write and leaks the temp. So Windows uses
// `MoveFileExA(..., MOVEFILE_REPLACE_EXISTING)`, which replaces atomically. Either way a
// failed replace removes the temp, so no `.tmp` is ever left behind.
#pragma once
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <system_error>
#include <vector>

#ifdef _WIN32
// Keep windows.h from defining min/max macros (they break std::min/std::max) and from
// pulling in rarely-needed subsystems. MinGW's libstdc++ predefines NOMINMAX, so guard to
// avoid a redefinition warning; MSVC does not, so the guard still protects that build.
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace csqv {

// Read the `sum=` recorded in a champion file's header, or -1 if absent/unreadable.
inline double read_saved_sum(const std::string& path) {
  std::ifstream f(path);
  std::string line;
  if (!f || !std::getline(f, line)) return -1.0;
  const auto p = line.find("sum=");
  if (p == std::string::npos) return -1.0;
  try {
    return std::stod(line.substr(p + 4));
  } catch (...) {
    return -1.0;
  }
}

// A process-and-thread-unique suffix for temp files, so two writers never share a .tmp.
inline uint64_t unique_tag() {
  static const uint64_t base =
      (static_cast<uint64_t>(std::random_device{}()) << 32) ^
      static_cast<uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count());
  static std::atomic<uint64_t> ctr{0};
  return base ^ ctr.fetch_add(1);
}

// Atomically REPLACE `path` with `tmp` (already written). Returns true on success. On
// failure the temp is removed, so a lost/failed write never leaks a `.tmp`.
//
// Windows note: `MoveFileExA` needs `MOVEFILE_REPLACE_EXISTING` to overwrite; without it
// (as plain `std::rename` does) the call fails whenever `path` exists. `MOVEFILE_WRITE_THROUGH`
// asks Windows to flush the replace to disk before returning. The POSIX branch does NOT
// fsync, so crash-durability is not symmetric across platforms; a follow-up ticket adds
// fsync parity if an unattended run must survive power loss. Both branches are RACE-free
// against a partial reader (the target is swapped in one filesystem operation).
inline bool atomic_replace(const std::string& tmp, const std::string& path) {
#ifdef _WIN32
  const bool ok = MoveFileExA(tmp.c_str(), path.c_str(),
                              MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) != 0;
#else
  std::error_code ec;
  std::filesystem::rename(tmp, path, ec);  // POSIX rename replaces atomically
  const bool ok = !ec;
#endif
  if (!ok) {
    std::error_code rmec;
    std::filesystem::remove(tmp, rmec);  // never leave a stray temp behind
  }
  return ok;
}

// Write a champion in the CENTERED frame [-0.5,0.5]^2, atomically and MONOTONICALLY: it
// writes only if `sum` strictly beats the value already on disk. WITHIN one process the
// caller's mutex serializes writers, so the file never regresses. ACROSS processes the
// read-then-replace has a narrow TOCTOU window (two writers can both pass the guard on the
// same stale on-disk sum), so a true cross-process guarantee needs OS file locking (a
// follow-up ticket); prefer ONE multi-threaded process per out_dir. Matches
// results/run_scalable.py::save_best so the Python tools re-read it unchanged. Returns true
// if it wrote.
//
// A partial/failed write must never replace a good champion, so the temp stream is checked
// after writing; on any stream error the temp is removed and no replace happens.
inline bool save_champion(const std::string& path, int n, int seed, const std::vector<double>& x,
                          const std::vector<double>& y, const std::vector<double>& r, double sum) {
  if (read_saved_sum(path) >= sum) return false;  // monotonic guard
  const std::string tmp = path + ".tmp." + std::to_string(unique_tag());
  bool ok;
  {
    std::ofstream f(tmp);
    f.setf(std::ios::fixed);
    f.precision(12);
    f << "# N=" << n << " seed=" << seed << " sum=" << sum << "\n";
    for (int i = 0; i < n; ++i)
      f << (x[i] - 0.5) << " " << (y[i] - 0.5) << " " << r[i] << "\n";
    f.flush();
    ok = static_cast<bool>(f);  // false if any write (disk full, I/O error) failed
  }
  if (!ok) {
    std::error_code ec;
    std::filesystem::remove(tmp, ec);  // never replace a good champion with a truncated file
    return false;
  }
  return atomic_replace(tmp, path);
}

}  // namespace csqv

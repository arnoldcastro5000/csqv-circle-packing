// Unit test for io.hpp (champion save: atomic, monotonic, no temp leak).
//
// Runs on any OS. The Windows-only freeze/leak bug (std::rename cannot replace an existing
// file) is fixed in atomic_replace; this test locks the OS-visible invariants: an improving
// save REPLACES the file, a worse save is a no-op, and NO `.tmp` is ever left behind, win,
// lose, or fail. Build: see cpp/CMakeLists.txt (target csqv_save_test) or cpp/Makefile.
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "io.hpp"

namespace fs = std::filesystem;

static int g_failures = 0;
#define CHECK(cond, msg)                                                    \
  do {                                                                      \
    if (!(cond)) {                                                          \
      std::fprintf(stderr, "FAIL: %s (line %d)\n", (msg), __LINE__);        \
      ++g_failures;                                                         \
    }                                                                       \
  } while (0)

// Count leftover `*.tmp.*` files in a directory (the leak signature).
static int count_temps(const fs::path& dir) {
  int c = 0;
  for (const auto& e : fs::directory_iterator(dir))
    if (e.path().filename().string().find(".tmp.") != std::string::npos) ++c;
  return c;
}

static bool sum_is(const std::string& path, double want) {
  return std::fabs(csqv::read_saved_sum(path) - want) < 1e-9;
}

int main() {
  const fs::path dir =
      fs::temp_directory_path() / ("csqv_save_test_" + std::to_string(csqv::unique_tag()));
  fs::remove_all(dir);
  fs::create_directories(dir);
  const std::string path = (dir / "champ.txt").string();

  const std::vector<double> x{0.25, 0.75}, y{0.5, 0.5}, r{0.1, 0.1};

  // 1. First save creates the file and leaves no temp.
  CHECK(csqv::save_champion(path, 2, 0, x, y, r, 1.000), "first save should write");
  CHECK(fs::exists(path), "file should exist after first save");
  CHECK(sum_is(path, 1.000), "on-disk sum should be 1.000");
  CHECK(count_temps(dir) == 0, "no temp after first save");

  // 2. A BETTER save must REPLACE the existing file. This is the Windows-critical path:
  //    plain std::rename fails when the target exists, freezing the file and leaking a temp.
  CHECK(csqv::save_champion(path, 2, 0, x, y, r, 2.000), "better save should write");
  CHECK(sum_is(path, 2.000), "on-disk sum should advance to 2.000");
  CHECK(count_temps(dir) == 0, "no temp after replacing save");

  // 3. A worse-or-equal save is a monotonic no-op, and never leaves a temp.
  CHECK(!csqv::save_champion(path, 2, 0, x, y, r, 2.000), "equal save should be a no-op");
  CHECK(!csqv::save_champion(path, 2, 0, x, y, r, 1.500), "worse save should be a no-op");
  CHECK(sum_is(path, 2.000), "on-disk sum stays 2.000 after no-ops");
  CHECK(count_temps(dir) == 0, "no temp after no-op saves");

  // 4. atomic_replace into a NON-EXISTENT directory must fail AND remove the temp (no leak).
  //    Locks the failure-cleanup path directly on any OS.
  const std::string tmp = (dir / "orphan.tmp.12345").string();
  { std::ofstream f(tmp); f << "junk\n"; }
  CHECK(!csqv::atomic_replace(tmp, (dir / "no_such_dir" / "target.txt").string()),
        "replace into a missing directory should fail");
  CHECK(!fs::exists(tmp), "failed replace must remove the temp (no leak)");

  // 5. A save whose temp cannot be written (parent dir missing) is a clean no-op: it
  //    returns false, creates no file, and leaves no temp. It must never replace a good
  //    champion with a partial/failed write.
  const std::string bad = (dir / "no_such_dir" / "champ.txt").string();
  CHECK(!csqv::save_champion(bad, 2, 0, x, y, r, 9.999), "save into a missing dir is a no-op");
  CHECK(!fs::exists(bad), "no file created on a failed write");
  CHECK(sum_is(path, 2.000), "the existing good champion is untouched by a failed save");
  CHECK(count_temps(dir) == 0, "no temp after a failed save");

  fs::remove_all(dir);
  if (g_failures == 0) std::printf("csqv_save_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}

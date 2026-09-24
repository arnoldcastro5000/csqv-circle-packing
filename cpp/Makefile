# Fallback build without CMake. On Windows with MinGW-w64:  mingw32-make -C cpp
# On Linux/macOS:  make -C cpp
# Produces cpp/build/csqv_worker and cpp/build/csqv_check.
CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops -std=c++17 -Icsqv
LDFLAGS ?= -pthread
# CSQV_ARCH is the -march target. The default, native, tunes the build for this host's CPU
# (+4-5% worker restarts on a Zen 3 host). A native binary can crash with an illegal
# instruction on an older CPU, so to copy a binary to another host, build it with a named
# target (make CSQV_ARCH=x86-64-v3) or with no -march (make CSQV_ARCH=).
# -ffp-contract=off is mandatory with any -march: a target with FMA lets the compiler fuse
# a * b + c, which changes the sealed checksums and the startup radii hash
# (csqv/fp_contract_test.cpp fails in such a build). arm64 always has FMA. The override
# appends both flags also to a CXXFLAGS from the environment or the command line.
CSQV_ARCH ?= native
override CXXFLAGS += $(if $(CSQV_ARCH),-march=$(CSQV_ARCH)) -ffp-contract=off

BUILD := build

# The compile command of the last native build. If CXX or CXXFLAGS (CSQV_ARCH included)
# change, the stamp is phony for this run: make rewrites it, and each native binary depends
# on it, so each one rebuilds with the new flags. GNU make functions only, no shell.
FLAGS_STAMP := $(BUILD)/.cxxflags
FLAGS_NOW := $(strip $(CXX) $(CXXFLAGS))
ifneq ($(strip $(file <$(FLAGS_STAMP))),$(FLAGS_NOW))
.PHONY: $(FLAGS_STAMP)
endif
$(FLAGS_STAMP): | $(BUILD)
	$(file >$@,$(FLAGS_NOW))

.PHONY: all clean test bench
all: $(BUILD)/csqv_worker $(BUILD)/csqv_check $(BUILD)/csqv_squeeze_check

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/csqv_worker: csqv/worker.cpp csqv/*.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/worker.cpp -o $@ $(LDFLAGS)

$(BUILD)/csqv_check: csqv/check.cpp csqv/*.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/check.cpp -o $@

$(BUILD)/csqv_squeeze_check: csqv/squeeze_check.cpp csqv/*.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/squeeze_check.cpp -o $@

$(BUILD)/csqv_save_test: csqv/save_test.cpp csqv/io.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/save_test.cpp -o $@

$(BUILD)/csqv_report_test: csqv/report_test.cpp csqv/report.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/report_test.cpp -o $@

$(BUILD)/csqv_center_polish_test: csqv/center_polish_test.cpp csqv/polish.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/center_polish_test.cpp -o $@

$(BUILD)/csqv_penalty_prefilter_test: csqv/penalty_prefilter_test.cpp csqv/polish.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/penalty_prefilter_test.cpp -o $@

$(BUILD)/csqv_fingerprint_test: csqv/fingerprint_test.cpp csqv/fingerprint.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/fingerprint_test.cpp -o $@

$(BUILD)/csqv_spread_test: csqv/spread_test.cpp csqv/spread.hpp csqv/search.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/spread_test.cpp -o $@

$(BUILD)/csqv_general_lp_test: csqv/general_lp_test.cpp csqv/general_lp.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/general_lp_test.cpp -o $@

$(BUILD)/csqv_lp_test: csqv/lp_test.cpp csqv/lp.hpp csqv/search.hpp csqv/benchmark.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/lp_test.cpp -o $@

$(BUILD)/csqv_lp_tableau_test: csqv/lp_tableau_test.cpp csqv/lp.hpp csqv/search.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/lp_tableau_test.cpp -o $@ $(LDFLAGS)

$(BUILD)/csqv_fp_contract_test: csqv/fp_contract_test.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/fp_contract_test.cpp -o $@

$(BUILD)/csqv_jam_slp_test: csqv/jam_slp_test.cpp csqv/jam_slp.hpp csqv/general_lp.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/jam_slp_test.cpp -o $@

# Bench tools for bit-identical speed work (see bench/README.md). The full gate against a
# git ref is bench/gate.sh; these targets build the tools from this tree only.
BENCH_TOOLS := perf_driver lp_bench trajectory_hash
bench: $(addprefix $(BUILD)/bench/,$(BENCH_TOOLS)) $(BUILD)/bench/pricing_check

$(BUILD)/bench:
	mkdir -p $(BUILD)/bench

$(addprefix $(BUILD)/bench/,$(BENCH_TOOLS)): $(BUILD)/bench/%: bench/%.cpp csqv/*.hpp | $(BUILD)/bench
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD)/bench/pricing_check: bench/pricing_check.cpp csqv/*.hpp | $(BUILD)/bench
	$(CXX) $(CXXFLAGS) -DCSQV_LP_VERIFY_PRICING $< -o $@

# Build and run the unit tests (save: atomic+monotonic, no leak; report: --record flag;
# center_polish: monotone exact-LP center squeeze; fingerprint: basin hash; spread: seed spreading;
# lp: incremental reduced-cost pricing, and a forced mid-solve fallback, are bit-identical to
# the full recompute; lp_tableau: each thread reuses one tableau buffer, with no stale state;
# fp_contract: CXXFLAGS do not fuse a * b + c into an FMA).
test: $(BUILD)/csqv_save_test $(BUILD)/csqv_report_test $(BUILD)/csqv_center_polish_test \
      $(BUILD)/csqv_fingerprint_test $(BUILD)/csqv_spread_test $(BUILD)/csqv_general_lp_test \
      $(BUILD)/csqv_jam_slp_test $(BUILD)/csqv_penalty_prefilter_test $(BUILD)/csqv_lp_test \
      $(BUILD)/csqv_lp_tableau_test $(BUILD)/csqv_fp_contract_test
	$(BUILD)/csqv_fp_contract_test
	$(BUILD)/csqv_save_test
	$(BUILD)/csqv_report_test
	$(BUILD)/csqv_center_polish_test
	$(BUILD)/csqv_penalty_prefilter_test
	$(BUILD)/csqv_fingerprint_test
	$(BUILD)/csqv_spread_test
	$(BUILD)/csqv_general_lp_test
	$(BUILD)/csqv_jam_slp_test
	$(BUILD)/csqv_lp_test
	$(BUILD)/csqv_lp_tableau_test

clean:
	rm -rf $(BUILD)

NATIVE_BINS := $(addprefix $(BUILD)/,csqv_worker csqv_check csqv_squeeze_check csqv_save_test \
  csqv_report_test csqv_center_polish_test csqv_penalty_prefilter_test csqv_fingerprint_test \
  csqv_spread_test csqv_general_lp_test csqv_lp_test csqv_lp_tableau_test csqv_fp_contract_test \
  csqv_jam_slp_test) \
  $(addprefix $(BUILD)/bench/,$(BENCH_TOOLS) pricing_check)
$(NATIVE_BINS): $(FLAGS_STAMP)

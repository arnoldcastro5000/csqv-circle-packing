# Fallback build without CMake. On Windows with MinGW-w64:  mingw32-make -C cpp
# On Linux/macOS:  make -C cpp
# Produces cpp/build/csqv_worker and cpp/build/csqv_check.
CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops -std=c++17 -Icsqv
LDFLAGS ?= -pthread

BUILD := build

.PHONY: all clean test
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

$(BUILD)/csqv_fingerprint_test: csqv/fingerprint_test.cpp csqv/fingerprint.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/fingerprint_test.cpp -o $@

$(BUILD)/csqv_spread_test: csqv/spread_test.cpp csqv/spread.hpp csqv/search.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/spread_test.cpp -o $@

$(BUILD)/csqv_general_lp_test: csqv/general_lp_test.cpp csqv/general_lp.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/general_lp_test.cpp -o $@

$(BUILD)/csqv_jam_slp_test: csqv/jam_slp_test.cpp csqv/jam_slp.hpp csqv/general_lp.hpp csqv/lp.hpp csqv/geometry.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) csqv/jam_slp_test.cpp -o $@

# Build and run the unit tests (save: atomic+monotonic, no leak; report: --record flag;
# center_polish: monotone exact-LP center squeeze; fingerprint: basin hash; spread: seed spreading).
test: $(BUILD)/csqv_save_test $(BUILD)/csqv_report_test $(BUILD)/csqv_center_polish_test \
      $(BUILD)/csqv_fingerprint_test $(BUILD)/csqv_spread_test $(BUILD)/csqv_general_lp_test \
      $(BUILD)/csqv_jam_slp_test
	$(BUILD)/csqv_save_test
	$(BUILD)/csqv_report_test
	$(BUILD)/csqv_center_polish_test
	$(BUILD)/csqv_fingerprint_test
	$(BUILD)/csqv_spread_test
	$(BUILD)/csqv_general_lp_test
	$(BUILD)/csqv_jam_slp_test

clean:
	rm -rf $(BUILD)

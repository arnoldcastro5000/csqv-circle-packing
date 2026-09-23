# CSQV native worker. Self-contained: no external dependencies, standard C++17 only.
# Builds two targets:
#   csqv_worker - the multi-threaded record search
#   csqv_check  - the parity harness (used by cpp/tools/run_parity.py)
#
# MinGW-w64 on Windows:
#   cmake -S cpp -B cpp/build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
#   cmake --build cpp/build -j
# Linux/macOS (default generator):
#   cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j
cmake_minimum_required(VERSION 3.15)
project(csqv_worker CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# -O3 plus -funroll-loops (measured ~+33% worker throughput, ticket 36; safe, no FP semantic
# change). Do NOT enable -ffast-math: the verifier relies on strict IEEE NaN/inf checks.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops")

find_package(Threads REQUIRED)

add_executable(csqv_worker csqv/worker.cpp)
target_include_directories(csqv_worker PRIVATE csqv)
target_link_libraries(csqv_worker PRIVATE Threads::Threads)

add_executable(csqv_check csqv/check.cpp)
target_include_directories(csqv_check PRIVATE csqv)

# squeeze_check: the C++-frame squeeze harness. Default = the first-order center_polish_analytic
# audit and analytic baseline (ticket 45); --jam = the SLP jammer + native .pck emit, the C++
# half of the acceptance gate (the acceptance harness re-scores it in the Python frame).
add_executable(csqv_squeeze_check csqv/squeeze_check.cpp)
target_include_directories(csqv_squeeze_check PRIVATE csqv)

# Unit tests. io.hpp: atomic + monotonic champion save, no temp leak. report.hpp:
# operator-supplied --record flag and advisory gap logic.
enable_testing()
add_executable(csqv_save_test csqv/save_test.cpp)
target_include_directories(csqv_save_test PRIVATE csqv)
add_test(NAME csqv_save_test COMMAND csqv_save_test)

add_executable(csqv_report_test csqv/report_test.cpp)
target_include_directories(csqv_report_test PRIVATE csqv)
add_test(NAME csqv_report_test COMMAND csqv_report_test)

# polish.hpp::center_polish: monotone exact-LP center squeeze on a champion.
add_executable(csqv_center_polish_test csqv/center_polish_test.cpp)
target_include_directories(csqv_center_polish_test PRIVATE csqv)
add_test(NAME csqv_center_polish_test COMMAND csqv_center_polish_test)

# polish.hpp::penalty_obj_grad: the squared-distance pre-filter is bit-exact vs the dense sweep.
add_executable(csqv_penalty_prefilter_test csqv/penalty_prefilter_test.cpp)
target_include_directories(csqv_penalty_prefilter_test PRIVATE csqv)
add_test(NAME csqv_penalty_prefilter_test COMMAND csqv_penalty_prefilter_test)

# fingerprint.hpp::basin_fingerprint: the WL contact-graph basin hash (stable, D4-invariant,
# discriminating). spread.hpp: the low-discrepancy seed spreading and defect kick. Both back the
# basin-dedup + seed-spread worker (ticket 42 / ADR 0002).
add_executable(csqv_fingerprint_test csqv/fingerprint_test.cpp)
target_include_directories(csqv_fingerprint_test PRIVATE csqv)
add_test(NAME csqv_fingerprint_test COMMAND csqv_fingerprint_test)

add_executable(csqv_spread_test csqv/spread_test.cpp)
target_include_directories(csqv_spread_test PRIVATE csqv)
add_test(NAME csqv_spread_test COMMAND csqv_spread_test)

# general_lp.hpp::general_lp: the general two-phase bounded-variable simplex (arbitrary rows,
# negative bounds, phase 1) that the jamming optimizer builds on (ticket 01).
add_executable(csqv_general_lp_test csqv/general_lp_test.cpp)
target_include_directories(csqv_general_lp_test PRIVATE csqv)
add_test(NAME csqv_general_lp_test COMMAND csqv_general_lp_test)

# jam_slp.hpp::jam_slp: the trust-region SLP jamming optimizer (Donev flex LP iterated to a
# jammed optimum), ported from problems/csqv/terminal_squeeze.py (ADR 0003).
add_executable(csqv_jam_slp_test csqv/jam_slp_test.cpp)
target_include_directories(csqv_jam_slp_test PRIVATE csqv)
add_test(NAME csqv_jam_slp_test COMMAND csqv_jam_slp_test)

# lp.hpp::detail::solve_reduced_lp: incremental reduced-cost pricing is bit-identical to the
# full recompute (radii and duals), and the half-integral guard never trips on real LPs.
add_executable(csqv_lp_test csqv/lp_test.cpp)
target_include_directories(csqv_lp_test PRIVATE csqv)
add_test(NAME csqv_lp_test COMMAND csqv_lp_test)

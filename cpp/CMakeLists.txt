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

# -O3 only. Do NOT enable -ffast-math: the verifier relies on strict IEEE NaN/inf checks.
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Threads REQUIRED)

add_executable(csqv_worker csqv/worker.cpp)
target_include_directories(csqv_worker PRIVATE csqv)
target_link_libraries(csqv_worker PRIVATE Threads::Threads)

add_executable(csqv_check csqv/check.cpp)
target_include_directories(csqv_check PRIVATE csqv)

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

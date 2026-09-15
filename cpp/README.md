# CSQV native worker (C++)

A native, multi-threaded C++ port of the scalable CSQV record search
(`results/run_scalable.py` + `results/scalable_polish.py`). It exists to run the search
at full core count on a host machine, past the vCPU cap of the sandbox the Python code
runs in. It targets **N=121** (the Specht-square record we beat) and works for any N.

The worker is **self-contained**: standard C++17 only, no external libraries, no network
at build time. It builds with a free compiler (MinGW-w64 on Windows).

## What it reproduces

One restart mirrors the Python chain exactly in structure:

1. `construct` - a varied seed (grid / jitter / hex / random)
2. `growpush` - LP-scaled overlap relaxation
3. `radii_lp` - the **exact** optimal radii for fixed centers (an LP)
4. `scalable_polish` - L-BFGS on a smooth overlap/wall penalty, rising penalty weight
5. `center_polish` - on a new best, an exact-LP center squeeze: a projected gradient
   ascent on the centers that maximizes the exact-LP radii sum (monotone). It recovers
   the last ~1e-5 of sum the penalty polish leaves unclaimed.
6. `verify_and_score` - the feasibility shrink and the score
7. perturb-the-best restarts for diversity; the best per seed persists atomically

The two trusted primitives are ported faithfully from `problems/csqv/`:

- **`radii_lp`** uses a compact bounded-variable primal simplex with lazy constraint
  generation. Only a few dozen pair constraints are ever tight (82 at N=121, 406 at
  N=484), so the active set stays tiny at any N. It is **strictly feasible** (overlap at
  machine precision) where scipy's HiGHS leaves ~1e-7 of overlap that must be shrunk
  away, so on degenerate near-tight packings it is slightly **more precise** than the
  Python solver (verified below).
- **`verify_and_score`** is a line-for-line port of the Python verifier.

## Verifying a packing

The C++ verifier scores the search, but a submission must be checked **independently** of
the search that produced it. This repository ships a standalone verifier that reads the
radii **as written** in a `.pck` (it does not re-derive them) and checks containment and
non-overlap at **zero tolerance**, in pure NumPy:

```sh
python3 tools/validate_pck.py discoveries/2026-09-14-15-csqv-120/csqv120.pck
```

Records move daily: re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` for the
LIVE best-known before any submission.

> Note: the original research monorepo also has a Python "source of truth" verifier
> (`verify_champion.py`) and a parity harness (`make_golden.py`, `run_parity.py`,
> `cross_validate.py`) that re-derive radii with the trusted `radii_lp`. Those depend on
> the monorepo's `problems/csqv/` package and are not included here; `tools/validate_pck.py`
> is the self-contained check for a finished `.pck`.

## Packomania `.pck` compliance (packomania.com/hints.html)

The worker writes `results-cpp/csqv<N>.pck` meeting the site format: line 1 the largest
radius (bare number), line 2 the author(s, comma-separated), then one `x y r` line per
circle **sorted by increasing radius**, in the **centered side-1 square** [-0.5, 0.5]^2
(container center at 0,0). Default precision is **15 dp** (full double precision for these
magnitudes; the site asks for as many decimals as possible); pass a trailing `pck_dp`
argument for a different value, e.g. `12`.

Packomania re-checks at zero tolerance, so `emit_pck` checks feasibility on the values AS
RE-PARSED FROM THE FILE and applies the smallest safety shrink that leaves a ROBUST
feasibility gap (overlap <= -1e-12, containment slack >= 1e-12), not merely <= 0. The
cushion sits well above the float noise floor (~1e-16 at 15 dp), so a re-check with a
different summation order cannot flip a residual positive. Cost to the sum is ~1e-10, far
below any record margin. This is why 15 dp is safe here where a naive fixed-precision emit
is not.

Submit by email to the maintainer, E. Specht. The current address (2026 CSQV page) is
`eckard.specht@ovgu.de`; the 2018 hints page lists `eckard.specht@physik.uni-magdeburg.de`.
Supply the coordinate file, not just the claimed sum.

## Build

MinGW-w64 (Windows), with CMake:

```sh
cmake -S cpp -B cpp/build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build -j
```

Without CMake (MinGW): `mingw32-make -C cpp`. Linux/macOS: `make -C cpp` or the CMake
lines with the default generator. Do **not** add `-ffast-math`: the verifier relies on
strict IEEE NaN/inf handling.

## Run

```sh
# csqv_worker <n> <budget_s> <base_seed> [threads] [out_dir] [author] [pck_dp] [--record <live_sum>]
cpp/build/csqv_worker 121 3600 0 12 results-cpp "Arnold Castro" 15 --record 5.797468812113
```

- `threads` independent workers use seeds `base_seed .. base_seed+threads-1`. Use your
  core count (e.g. 12 on a Ryzen 5 5600H).
- `--record <live_sum>` is OPTIONAL and advisory only. The search never uses it; it just
  labels the console with the gap to that value and prints `*** RECORD BEATEN ***` once the
  sum clears it. Records move daily, so the worker does NOT hardcode one: pass the live
  best-known from `https://www.packomania.com/csqv/txt/sumradii.txt`. Without it the worker
  omits the gap and never claims a record; `tools/validate_pck.py` remains the
  authority on whether a finished packing is feasible and beats the live record.
- Each worker writes `results-cpp/n<N>-s<seed>-best.txt` (centered frame) on every
  improvement, and the overall best to `results-cpp/n<N>-best.txt` plus a ready
  `results-cpp/csqv<N>.pck`.
- **Writes are atomic and MONOTONIC**: each save reads the sum already on disk and rewrites
  only if the new value is strictly higher, using a unique temp name. So a file never
  regresses and two writers cannot corrupt each other's temp, even if several worker
  PROCESSES share one `out_dir` (they cooperate: the global files converge to the true max).
  Prefer one process with many threads for efficiency, but concurrent processes are now safe.
- **Resume**: relaunch with the same `out_dir` and seeds. Each worker resumes from its saved
  champion, continues via perturb-best restarts, and immediately re-publishes the global
  `n<N>-best.txt` + `.pck` from the resumed champion (so they reflect it even with no new gain).
  To seed a fresh search from a known champion, copy it into `n<N>-s<seed>-best.txt` first.
- A cross of the live record prints `*** RECORD BEATEN ***`. Confirm the emitted `.pck`
  with `tools/validate_pck.py` before believing it.

## Startup accuracy/precision benchmark

Every worker start runs a fixed, deterministic suite of the exact primitives and writes it
to `results-cpp/benchmark-<YYYYMMDD-HHMMSS>.txt`. The cases are generated from the RAW
`mt19937_64` sequence (standard-specified, identical on every conforming compiler) mapped
to [0,1) with a fixed formula, so the SAME build produces the SAME cases on Linux and
Windows. The file records, per case, the centers and radii at 17 significant digits, the
sum/score/feasibility, and a bitwise FNV-1a hash of the radii. The startup line reports
`primitives_feasible=YES/NO`.

The benchmark file lets you compare two OS runs directly: the same build produces the
same cases on Linux and Windows, so their sums must agree to machine precision and their
radii hashes must be identical (bit-identical numerics, expected for the same compiler
family since `sqrt` and the basic operations are IEEE-correctly-rounded and the build uses
no `-ffast-math`). `check.cpp` (target `csqv_check`) reads unit-frame centers and prints
`radii_lp` + `verify_and_score`, for spot checks against any reference solver.

> The original research monorepo has a Python parity harness (`make_golden.py`,
> `run_parity.py`, `cross_validate.py`) that diffs these primitives against the trusted
> `problems/csqv/` oracle. It is not included here; this repository's self-contained check
> is `tools/validate_pck.py` on a finished `.pck`.

## Layout

```
cpp/
  csqv/
    geometry.hpp        frame, wall_slack, feasibility shrink, verify_and_score
    lp.hpp              radii_lp: bounded simplex + lazy constraint generation
    polish.hpp          penalty objective, projected L-BFGS (scalable_polish),
                        and center_polish (exact-LP center squeeze)
    search.hpp          construct + growpush
    benchmark.hpp       deterministic startup accuracy/precision self-benchmark
    io.hpp              atomic, monotonic champion save (no temp leak)
    report.hpp          the advisory --record flag
    worker.cpp          multi-threaded restart loop, resume, atomic save, .pck emit
    check.cpp           read centers -> radii_lp + verify_and_score
    *_test.cpp          unit tests (save, report, center_polish)
  CMakeLists.txt   Makefile   README.md
../tools/
    validate_pck.py     self-contained zero-tolerance .pck verifier (pure NumPy)
```

## Known limitation

At very large N (e.g. 484) `radii_lp` is slower than at N=121 (the dense inner simplex
grows with the active set). N=121, the requested target, is fast (~10 ms per `radii_lp`
call, faster than the Python HiGHS path). A sparse inner solve would lift the large-N
case; it is not needed for N=121.

# Bench tools: the bit-identical gate for speed work

A bit-identical speed change must give the same bytes for the same seed. These tools prove
that and measure the gain. `gate.sh` runs all of them against a git ref.

## Run the full gate

From the repo root:

```sh
cpp/bench/gate.sh <base-ref> [out_dir]
```

- Head is the `cpp/csqv` source in the working tree, uncommitted edits included.
- Base is `cpp/csqv` at `<base-ref>`, taken from a temporary git worktree. The script removes the
  worktree at exit.
- A third build, `fb`, is head with `-DCSQV_LP_FORCE_FALLBACK_AT=3`. In this build, the
  half-integral guard fails at the 3rd basis change of each LP solve, so the fallback path runs.
  The fallback build must also give the same bytes, and the worker must print the pricing
  WARNING line.
- The default run takes about 20 to 25 minutes on 4 vCPUs. For a quick check, set smaller
  values, for example `GATE_NS="30 90" GATE_SEEDS=2 GATE_BUDGET=3 GATE_SPEED=0`. The script
  header lists every variable.
- The exit status is 0 only if every gate holds. The last line is `GATE: PASS` or `GATE: FAIL`.
- Both sides build with the `cpp/Makefile` default flags (`-march=native -ffp-contract=off`). To
  gate a flag change, give each side its flags with `GATE_BASE_FLAGS` and `GATE_HEAD_FLAGS`.

| Gate | Tool | Pass condition |
| --- | --- | --- |
| Checksums | `perf_driver.cpp` | Seed 777: N=90 (20 iters) = 98.374190007874, N=143 (8 iters) = 49.910861299027, the same hex bits in base, head and fb. |
| G1 | the worker | The startup benchmark data lines (the radii hashes) are identical. |
| G2 | `trajectory_hash.cpp` | The hash of construct, growpush, the LP duals and the center polish is identical in base, head and fb. |
| G3 | the worker + `compare_runs.py` | The per-restart `best` events of each seed are identical up to the shorter run's restart count. |
| Pricing | `pricing_check.cpp` | Each incremental reduced cost equals the full recompute (`mismatches=0`), and the guard never fails (`fallbacks=0`). |
| G4 | `make test`, `make accept` | Both pass (the script skips `accept` if the Makefile has no such target). |

The speed section (not a gate) prints the `lp_bench` time and the 30 s worker restart count,
base vs head. It runs one job at a time, so the numbers are comparable.

## Build one tool

`make -C cpp bench` builds the tools from the working tree into `cpp/build/bench/`:

```sh
cpp/build/bench/perf_driver 90 20 777     # restarts/s + the sealed checksum
cpp/build/bench/lp_bench 121              # radii_lp time only
cpp/build/bench/trajectory_hash 121 20 5  # G2 hash: <n> <seeds> <polish steps>
cpp/build/bench/pricing_check 90 4        # in-loop pricing check: <n> <seeds>
```

## Profile the stages of a restart

The stage worker is the worker with the stage timers of `csqv/stage_timer.hpp`
(`-DCSQV_STAGE_TIMERS`). It takes the worker's arguments, follows the same search path for a
seed, and prints a stage table at exit. Each mark in the worker gives the time since the
previous mark of its thread to one stage. `radii_lp` adds its own time to a per-thread LP clock,
so the table also shows the `radii_lp` part of each stage. Without the macro, the marks expand
to nothing, and the production worker is the same binary as before.

```sh
make -C cpp stages
cpp/build/bench/csqv_worker_stages 121 60 0 1 /tmp/st121
```

The table has one row per stage. `loop%` is the part of the restart loop, `per restart` is the
time per restart, `radii_lp%` is the part of the stage in `radii_lp`, and `lp/rst` is the
`radii_lp` calls per restart. Startup, final polish and seal run once, so they show only
seconds. With several threads, the seconds are the sum over the threads.

## Rules

- Compile the tools only through `make bench` or `gate.sh`. Both give the headers with
  `-I<tree>/cpp/csqv` and nothing else. An old copy of a header next to a source file wins over
  `-I`. An old `lp.hpp` copy in `cpp/build/` once made a speed number measure the old code,
  while the checksums still matched.
- `pricing_check` needs `-DCSQV_LP_VERIFY_PRICING`. That macro adds a full recompute to each
  pricing step, so never define it in a worker build.
- `CSQV_LP_FORCE_FALLBACK_AT` is for tests and the `fb` build only. A worker built with it prints
  the pricing WARNING line at exit.

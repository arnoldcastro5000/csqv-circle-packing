# csqv-circle-packing

A solver and a set of verified record packings for the **CSQV** problem: pack `N`
circles of **variable radii** into a unit square with no overlap, maximizing the **sum of
the radii**. This is the "circles in a square, variable radii" (csqv) entry in Eckard Specht's
Packomania [record](https://www.packomania.com/csqv/csqv.html).

The solver is coded by Claude Code Opus 4.8 (the native C++ worker) and human-directed. It proves 
some targets are beatable and are offered as verified feasible packings. Each discovery folder 
carries its own `README.md` with the exact provenance.

The solver is a self-contained C++17 worker (no dependencies, no network at build time).
The `discoveries/` folder holds packings that beat, match, or extend the published
best-known, each independently verified feasible at zero tolerance.

## Records in this repository

Live best-known re-fetched 2026-09-16 (N=93, 94, 95, 96, 97, 98, 99, 118, 119) and from the
2026-09-15 page (N=120, 121, 141, 142, 143).

| N   | our packing (sum) | Packomania live best-known               | status |
|-----|-------------------|------------------------------------------|--------|
| 93  | 5.073447923447    | 5.071029091028 (Dutton [10], 2026-09-16) | **leads** live by +0.0477%; SEALED (jammed), held, submitted |
| 94  | 5.100695331105    | 5.098748716187 (Dutton [10], 2026-09-16) | **leads** live by +0.0382%; SEALED (jammed), held, submitted |
| 95  | 5.126321536373    | 5.125823937617 (Specht [2], 2026-09-16)  | **leads** live by +0.0097%; SEALED (jammed), held, submitted |
| 96  | 5.153337299999    | 5.152146722729 (Dutton [10], 2026-09-16) | **leads** live by +0.0231%; SEALED (jammed), held, submitted |
| 97  | 5.180768365834    | 5.180639279285 (Dutton [10], 2026-09-16) | **leads** live by +0.0025%; SEALED (jammed), held, submitted |
| 98  | 5.208554826179    | 5.206668546480 (Garg [18], 2026-09-16)   | **leads** live by +0.0362%; SEALED (jammed), held, submitted |
| 99  | 5.234274484443    | 5.233333680412 (Tasseff [13], 2026-09-16)| **leads** live by +0.0180%; SEALED (jammed), held, submitted |
| 118 | 5.724613195610    | 5.723940934671 (Tariq [20], 2026-09-16)  | **leads** live by +0.0117%; SEALED (jammed), held, submitted |
| 119 | 5.749508785837    | 5.749519103944 (Castro [21], 2026-09-16) | submitted + accepted 2026-09-16; site re-optimized higher (+1.0e-5) to a value our jamming optimizer independently reaches (5.749519103949) |
| 120 | 5.774548479515    | 5.774582323763 (Castro [21], 2026-09-15) | submitted + accepted; site re-optimized higher (+3.4e-5), so live now exceeds our recorded value |
| 121 | 5.799069987171    | 5.799103501951 (Castro [21], 2026-09-15) | submitted + accepted; site re-optimized higher (+3.4e-5), so live now exceeds our recorded value |
| 141 | 6.261881429857    | 6.261919395210 (Castro [21], 2026-09-15) | submitted + accepted; site re-optimized higher (+3.8e-5), so live now exceeds our recorded value |
| 142 | 6.289814405001    | 6.289851836935 (Castro [21], 2026-09-15) | submitted + accepted; site re-optimized higher (+3.7e-5), so live now exceeds our recorded value |
| 143 | 6.312963889026    | 6.312993919570 (Castro [21], 2026-09-15) | submitted + accepted; site re-optimized higher (+3.0e-5), so live now exceeds our recorded value |

Every N we submitted (119, 120, 121, 141, 142, 143) was accepted and then re-optimized
higher on the site (all credited to Castro [21]). N=141-143, previously untracked, are now 
published (introduced 14-Sep-2026); none is marked proven-optimal.

Records move daily. Re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` before
trusting any comparison above.

### The N=121 story

The submitted N=121 packing (sum 5.799069987171) was **accepted**. Packomania then
published it as **5.799103501951**, higher than the submitted value. This is not a counting
error: the published packing is the submitted arrangement (transposed) with the centers
nudged by ~1e-6, and running this repository's exact LP on the published centers reproduces
5.799103501951 to 1e-12. The maintainer re-optimized the accepted submission, standard
Packomania practice. That finding motivated the `center_polish` stage described below.

`discoveries/2026-09-13-20-csqv-121/` is the submitted packing; `2026-09-14-00-csqv-121/`
is the same centers with exactly-optimal (exact-LP) radii, 5.799074421369.

N=120 followed the same pattern: our submission (5.774548479515) was accepted and
then re-optimized higher on the site to 5.774582323763 (credited to Castro [21]).
So both N=120 and N=121 currently show a live value above the packing recorded here.

## The solver

The C++ worker (`cpp/`) runs, per restart: a varied construction, an LP-scaled overlap
relaxation, the **exact-LP radii** for fixed centers, an L-BFGS penalty polish, and then
a **`center_polish_analytic`** stage: a projected gradient ascent on the centers that
maximizes the exact-LP radii sum, using the LP dual as the **analytic** gradient (one LP solve
per step, monotone). It recovers the last ~1e-5 of sum the penalty polish leaves unclaimed, the
same gap the Packomania re-optimization closed on N=121. An optional `--spread` mode adds
cross-thread basin de-duplication and low-discrepancy seed spreading (ADR 0002). See
[`cpp/README.md`](cpp/README.md) for the full method, build, and run instructions.

For a finished champion there is an optional **terminal second-order squeeze**
(`problems/csqv/terminal_squeeze.py`, NumPy + SciPy): a full-space sparse NLP over (x, y, r)
plus a guarded contact-graph KKT-Newton, run once on the single best packing. It is not part of
the self-contained C++ worker. See `cpp/README.md` for the scope and limits.

Key idea: for **fixed centers**, the optimal radii are the solution of a linear program
(cheap and exact at any N). All the difficulty is in the **center arrangement**, so the
search explores center arrangements and the LP always assigns the best radii.

### Build

```sh
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j
# or, without CMake:
make -C cpp
```

### Run

```sh
# csqv_worker <n> <budget_s> <base_seed> [threads] [out_dir] [author] [pck_dp] [--record <live_sum>]
cpp/build/csqv_worker 121 3600 0 8 results-cpp "Your Name" 15 --record 5.797468812113
```

The worker writes the best centers and a ready-to-submit `.pck` into `out_dir`, resumes
from a saved champion on relaunch, and never regresses a saved file. `--record` is advisory
only (a console gap and banner); the search never uses it.

## Verifying a packing

`tools/validate_pck.py` reads the radii **as written** in a `.pck` (it does not re-derive
them) and checks containment and non-overlap at **zero tolerance**, in pure NumPy. It is an
independent second implementation, so it does not trust the solver that produced the file.

```sh
python3 tools/validate_pck.py discoveries/2026-09-14-15-csqv-120/csqv120.pck
```

## Submission format (Packomania)

A `.pck` file: line 1 the largest radius, line 2 the author, then one `x y r` line per
circle **sorted by increasing radius**, in the **centered side-1 square** [-0.5, 0.5]^2
(container center at the origin). Submit by email to the maintainer, E. Specht. Records
move daily, so re-fetch the live best-known at submission time. A submission is momentous
and outward-facing: verify rigorously first.

## Layout

```
cpp/           self-contained C++17 worker, tests, and its own README
tools/         validate_pck.py, the zero-tolerance .pck verifier
problems/csqv/ vendored trusted primitives + the terminal squeeze (NumPy + SciPy)
tests/         test_terminal_squeeze.py
discoveries/   verified record packings, each with a .pck and an honest README
docs/adr/      architecture decision records
docs/research/ the CSQV method, target-selection, and solver research notes
```

## License

[MIT](LICENSE). Record packings may also be freely reused; attribution is appreciated.

# csqv-circle-packing

A solver and a set of verified record packings for the **CSQV** problem: pack `N`
circles of **variable radii** into a unit square with no overlap, maximizing the **sum of
the radii**. This is the "circles in a square, variable radii" (csqv) entry in E. Specht's
[Packomania](https://www.packomania.com/) record book.

The solver is a self-contained C++17 worker (no dependencies, no network at build time).
The `discoveries/` folder holds packings that beat, match, or extend the published
best-known, each independently verified feasible at zero tolerance.

## Records in this repository

| N   | sum of radii     | Packomania best-known (2026-09-14) | status |
|-----|------------------|------------------------------------|--------|
| 120 | 5.774548479515   | 5.773179664810                     | **beats** the best-known by +0.0237% |
| 121 | 5.799069987171   | 5.797468812113 (at submission)     | **accepted**; Packomania re-optimized and published 5.799103501951 |
| 141 | 6.261881429857   | not published                      | feasible packing for an **untracked** N |
| 142 | 6.289814405001   | not published                      | feasible packing for an **untracked** N |
| 143 | 6.312963889026   | not published                      | feasible packing for an **untracked** N |

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

## Honest provenance

These packings were found by a **numerical search** (a hand-written cold multi-start, then
the native C++ worker), **not** by an LLM. They prove the targets are beatable and are
offered as verified feasible packings. Each discovery folder carries its own `README.md`
with the exact provenance.

## The solver

The C++ worker (`cpp/`) runs, per restart: a varied construction, an LP-scaled overlap
relaxation, the **exact-LP radii** for fixed centers, an L-BFGS penalty polish, and then
a **`center_polish`** stage: a projected gradient ascent on the centers that maximizes the
exact-LP radii sum (monotone). `center_polish` recovers the last ~1e-5 of sum the penalty
polish leaves unclaimed, the same gap the Packomania re-optimization closed on N=121. See
[`cpp/README.md`](cpp/README.md) for the full method, build, and run instructions.

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
discoveries/   verified record packings, each with a .pck and an honest README
docs/adr/      architecture decision records
docs/research/ the CSQV method, target-selection, and solver research notes
```

## License

[MIT](LICENSE). Record packings may also be freely reused; attribution is appreciated.

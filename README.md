# csqv-circle-packing

A solver and a set of verified record packings for the **CSQV** problem: pack `N`
circles of **variable radii** into a unit square with no overlap, maximizing the **sum of
the radii**. This is the "circles in a square, variable radii" (csqv) entry in Eckard Specht's
Packomania [record](https://www.packomania.com/csqv/csqv.html).

The solver is coded by Claude Code Opus 4.8 (the native C++ worker) and human-directed. It proves 
some targets are beatable and are submitted as verified feasible packings. Each discovery folder 
carries its own `README.md`.

The solver is a self-contained C++17 worker (no dependencies, no network at build time).
The `discoveries/` folder holds packings that beat, matched, or extended the published
best-known at the time they were found, each independently verified feasible at zero
tolerance.

## Records in this repository

Live best-known re-fetched 2026-09-23 for every N below.

| N   | our packing (sum) | Packomania live best-known   | status |
|-----|-------------------|------------------------------|--------|
| 84  | 4.816500878221    | 4.816500879756 (Castro [21]) | **still ours**: published and credited to us, the only record here still standing; the site re-optimized our submission by +1.5e-9 (float-scale), and our full seal independently reaches 4.816500879755 |
| 93  | 5.073447923447    | 5.076245164998 (Denoual)     | submitted; the live best-known now exceeds our packing |
| 94  | 5.100695331105    | 5.101975416147 (Denoual)     | submitted; the live best-known now exceeds our packing |
| 95  | 5.126321536373    | 5.129111105835 (Denoual)     | submitted; the live best-known now exceeds our packing |
| 96  | 5.153337299999    | 5.155289947655 (Heap [17])   | submitted; the live best-known now exceeds our packing |
| 97  | 5.180768365834    | 5.181609663002 (Heap [17])   | submitted; the live best-known now exceeds our packing |
| 98  | 5.208554826179    | 5.208578422017 (Heap [17])   | submitted; the live best-known now exceeds our packing |
| 99  | 5.234274484443    | 5.235812940861 (Heap [17])   | submitted; the live best-known now exceeds our packing |
| 118 | 5.724613195610    | 5.726609681279 (Heap [17])   | submitted; the live best-known now exceeds our packing |
| 119 | 5.749508785837    | 5.751784962150 (Heap [17])   | submitted + accepted, since surpassed on the live site |
| 120 | 5.774548479515    | 5.776687225070 (Heap [17])   | submitted + accepted, since surpassed on the live site |
| 121 | 5.799069987171    | 5.802127855949 (Heap [17])   | submitted + accepted, since surpassed on the live site |
| 141 | 6.261881429857    | 6.269795235354 (Heap [17])   | submitted + accepted, since surpassed on the live site |
| 142 | 6.289814405001    | 6.292381128997 (Heap [17])   | submitted + accepted, since surpassed on the live site |
| 143 | 6.312963889026    | 6.318041422739 (Denoual)     | submitted + accepted, since surpassed on the live site |

As of 2026-09-23, none of the packings listed above still leads. Other contributors have
surpassed every one of them on the live site: Wilfred Heap [17] now holds N=96, 97, 98, 99,
118, 119, 120, 121, 141, and 142; Jean-René Denoual holds N=93, 94, 95, and 143. The only
CSQV record still credited to us on Packomania is N=84 (4.816500879756, Arnold Castro [21]);
its packing is in `discoveries/2026-09-23-00-csqv-84/`. None of these packings is marked
proven-optimal.

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

### Requirements

- A C++17 compiler (g++ or clang; MinGW-w64 on Windows) and `make` or CMake. The worker has
  no other dependencies.
- Python 3 with NumPy for `tools/validate_pck.py`. SciPy is also necessary for the terminal
  squeeze (`problems/csqv/terminal_squeeze.py`), and pytest for `tests/`:
  `pip install numpy scipy pytest`.

### Build

```sh
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j
# or, without CMake:
make -C cpp
```

### Run

```sh
# csqv_worker <n> <budget_s> <base_seed> [threads] [out_dir] [author] [pck_dp] [--record <live_sum>]
cpp/build/csqv_worker 121 3600 0 8 results-cpp "Your Name" 15 --record 5.802127855949
```

The worker writes the best centers and a ready-to-submit `.pck` into `out_dir`, resumes
from a saved champion on relaunch, and never regresses a saved file. `--record` is advisory
only (a console gap and banner); the search never uses it. The example value is the live N=121
best-known on 2026-09-23. Records move, so take the current value from
`https://www.packomania.com/csqv/txt/sumradii.txt`.

## The search path

Each worker thread repeats a **restart** (one full attempt, from a starting arrangement to a
local optimum). A restart runs these stages in order.

1. **Pick the start.** A **warm restart** perturbs the thread's own best centers with small
   Gaussian noise. This is a **basin-hopping** step (basin hopping = perturb the current best,
   re-optimize, then keep the result only if it is better). It searches for **depth** (a better
   optimum near the current one). The search picks a warm restart most of the time once a thread
   has a best. A **cold restart** instead builds a fresh arrangement. In the default and
   `--count` modes a cold restart cycles through four constructions: grid, jitter, hex, and
   random. In `--spread` mode a **low-discrepancy sequence** (a number sequence whose points
   spread evenly and never cluster) picks the construction family and its rotation, aspect, and
   jitter, so each thread starts in a different region. Seven families are implemented (grid,
   jitter, hex, random, phyllotaxis, Archimedean spiral, and concentric rings); the build
   activates the first four by default.

2. **growpush.** This stage relaxes the overlap (it grows the circles and pushes each
   overlapping pair apart) until the arrangement is feasible and tight.

3. **The exact-LP radii (`radii_lp`).** For the **fixed** centers, this stage solves a **linear
   program** (an LP: it maximizes a linear objective under linear constraints). The LP gives each
   circle its largest radius with no overlap and full containment. The LP is exact and cheap at
   any `N`. So the search only moves the centers; the LP always assigns the best radii.

4. **The penalty polish (`scalable_polish`).** This stage moves the centers to raise the radii
   sum. It uses **L-BFGS** (a quasi-Newton optimizer: it approximates the curvature from recent
   gradients) on a **penalty** objective (it adds a soft cost for any overlap). It runs a
   **continuation ladder** (it solves an easy problem first, then repeats with a larger penalty
   weight each round), so the overlap falls to near zero.

5. **Bank the result.** If the polished sum beats the thread's best, the search saves it at once,
   so no improvement is ever lost.

6. **Basin bookkeeping (only in `--count` and `--spread`).** The search makes a **basin
   fingerprint** (a hash of the converged optimum; a **basin** is the set of starts that flow to
   one local optimum) and records it. The count of distinct fingerprints is the coverage metric.
   In `--spread` mode only, a cold seed that re-descends a seen basin triggers an **escape**: a
   **defect kick** (it removes the few smallest circles and reinserts them at the emptiest spots,
   which changes the **contact graph**, the graph of which circles touch) fires up to three times,
   to hop to an unseen basin.

7. **The final squeeze (`center_polish_analytic`).** At the budget end the thread runs one
   **projected gradient ascent** on the centers (projected = each step stays inside the box) that
   maximizes the exact-LP radii sum. It uses the **LP dual** (the shadow prices the LP returns) as
   the exact gradient, so one LP solve gives one exact step. This recovers the last ~1e-5 of sum
   the penalty polish leaves. A single finished champion may also get the optional terminal
   squeeze described above.

Two modes control the coverage machinery. The default mode is a plain **multi-start** (many
independent restarts). `--count` adds basin counting only. `--spread` adds the spread seeds, the
de-duplication, and the escape.

## Findings

- **The LP dominates the run time.** A profile of a record-band run shows the exact-LP solve is
  about 70% of all executed lines. Inside the LP, the **pricing** step is the hot spot (pricing =
  the sweep over the non-basic variables that picks the next one to enter the basis on each
  simplex iteration). The **pivot** update (it rewrites the tableau after each swap) is the next
  cost. The penalty polish is far smaller, and the basin and escape bookkeeping is negligible. So
  a faster LP is the fastest way to a better packing: more LP throughput means more restarts in
  the same budget, which means wider coverage.

- **Large `N` is throughput-bound.** At `N` above ~120 a single restart takes tens of seconds, so
  few restarts finish per hour. The warm depth walk then barely engages, and the LP speed sets the
  search rate.

- **The warm walk finds the champions.** Most record packings come from a warm basin-hopping
  refinement of a basin that a cold seed first found. The cold seed supplies the basin; the warm
  walk deepens it.

- **Cold spread seeds land in fresh basins.** The low-discrepancy spreading makes almost every
  cold seed reach a new basin at the record band, so the escape rarely fires there. The spreading,
  not the escape, gives `--spread` its coverage.

- **`--count` and `--spread` trade depth for coverage.** `--count` often reaches the single better
  basin at high `N`. `--spread` reaches more distinct basins. Each discovery README records which
  mode found that `N`.

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

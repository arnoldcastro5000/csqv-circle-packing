# Discovery: CSQV N=120 circle packing (beats best-known, exact-LP radii)

An improved packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 120. It beats the published
best-known.

IMPORTANT, honest provenance: this packing was found by an AI-authored numerical
search. It is a solver-tooling result that proves
the target is beatable.

## Provenance

| field | value |
|---|---|
| result files | `csqv120.pck` (15 dp, strictly feasible), `n120-best.txt` (the centers) |
| centers found by | the native C++ cold multi-start worker (`cpp/csqv/worker.cpp`, `search.hpp`), base_seed 0, thread seed 1, no LLM |
| radii + `.pck` by | the native C++ exact LP + rounding-aware emit in `cpp/` (`cpp/csqv/lp.hpp`, `worker.cpp::emit_pck`) |
| repo commit (code base) | 003734f |
| primitives used | exact radii-by-LP for fixed centers; independent verifier recomputes containment, overlap, and the sum |
| submitted to Packomania | NO. Held; a submission is the operator's call. |

## Result

| quantity | value |
|---|---|
| sum of radii (this packing, as re-parsed from the 15-dp `.pck`) | 5.774548479515367 |
| sum of radii (in-memory C++ champion) | 5.774548479635 |
| Packomania live best-known (fetched 2026-09-14) | 5.773179664810 |
| margin over the live record | +0.001368815 (about +0.0237%) |

The live best-known was confirmed twice on 2026-09-14: the master table
`csqv/txt/sumradii.txt` and the detail page `csqv/txt/csqv120.txt` both give
`sumradii = 5.773179664810`.

## Verification (independent, from the .pck)

Parsed the exact `.pck` and recomputed geometry from scratch (pure numpy, radii
read AS WRITTEN, no re-derivation), at ZERO tolerance:

- 120 circles, radii in [0.036, 0.063], sorted by increasing radius.
- Containment (centered side-1 square, frame center 0,0): max wall violation
  -1.0e-12 (inside, robust cushion).
- Non-overlap (all 7140 pairs): max pair overlap -2.0e-12 (no overlap, robust
  cushion).
- The first line (largest radius) equals the largest radius in the rows.
- STRICTLY FEASIBLE at zero tolerance: YES.

## Note on LP precision (exact-LP vs scipy)

The `.pck` uses the native C++ exact LP (`cpp/csqv/lp.hpp`), which is strictly
feasible to about 1e-17 on these centers and scores 5.774548. The Python (scipy
HiGHS) LP on the SAME centers leaves about 8.7e-8 of overlap on the raw radii and
so scores 5.774543 after the feasibility shrink. The gap 5.774548 vs 5.774543 is
LP PRECISION on the same centers, not a better arrangement. Both beat the live
record.

## Reproduce

The exact-LP radii and the strictly-feasible 15-dp `.pck` (from the native C++
worker, built per `cpp/README.md`):

```
csqv_worker 120 <budget_s> 0 <threads> <out_dir> "Arnold Castro" 15 --record 5.773179664810
```

The `--record <live_sum>` flag is optional and advisory (it labels the console
with the gap and prints `*** RECORD BEATEN ***`); pass the live best-known from
`sumradii.txt`. The search does not use it.

It resumes from `<out_dir>/n120-s<base_seed>-best.txt`, writes the champion
centers to `<out_dir>/n120-best.txt`, and emits `<out_dir>/csqv120.pck` at 15 dp.
Re-verify the round-trip at zero tolerance before trusting any value. Records move
daily; re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` at verification
time.

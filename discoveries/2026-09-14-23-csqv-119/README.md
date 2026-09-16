# Discovery: CSQV N=119 circle packing (beats best-known, exact-LP radii)

An improved packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 119. It beats the published
best-known.

IMPORTANT, honest provenance: this packing was found by an AI-authored numerical
search, NOT by the LLM-evolution loop. It is a solver-tooling result that proves
the target is beatable. Making it a true Crucible-Loop (LLM-evolved) discovery is
future work.

## Provenance

| field | value |
|---|---|
| result files | `csqv119.pck` (15 dp, strictly feasible), `n119-best.txt` (raw champion centers) |
| centers found by | the native C++ cold multi-start worker (`cpp/csqv/worker.cpp`, `search.hpp`), base_seed 0, thread seed 0, no LLM |
| radii + `.pck` by | the native C++ exact LP + rounding-aware emit in `cpp/` (`cpp/csqv/lp.hpp`, `worker.cpp::emit_pck`) |
| terminal squeeze | the exact-LP center squeeze via the analytic dual gradient (`cpp/csqv/polish.hpp::center_polish_analytic`); re-squeezed to the analytic optimum on 2026-09-16, so the residual headroom is now zero |
| repo commit (code base) | 537810f |
| primitives used | exact radii-by-LP for fixed centers; independent verifier recomputes containment, overlap, and the sum |
| submitted to Packomania | NO. Held; a submission is the operator's call. |

## Result

| quantity | value |
|---|---|
| sum of radii (this packing, as re-parsed from the 15-dp `.pck`) | 5.749508785837 |
| sum of radii (in-memory C++ champion, post re-squeeze) | 5.749508785956 |
| Packomania live best-known (fetched 2026-09-14) | 5.748642197922 |
| margin over the live record | +0.000866490 (about +0.0151%) |

The live best-known was confirmed twice on 2026-09-14: the master table
`csqv/txt/sumradii.txt` and the detail page `csqv/txt/csqv119.txt` both give
`sumradii = 5.748642197922`.

## Verification (independent, from the .pck)

Parsed the exact `.pck` and recomputed geometry from scratch (pure numpy, radii
read AS WRITTEN, no re-derivation), at ZERO tolerance:

- 119 circles, radii in [0.031043, 0.062312], sorted by increasing radius.
- Containment (centered side-1 square, frame center 0,0): max wall violation
  -1.0e-12 (inside, robust cushion).
- Non-overlap (all 7021 pairs): max pair overlap -2.0e-12 (no overlap, robust
  cushion).
- The first line (largest radius) equals the largest radius in the rows.
- STRICTLY FEASIBLE at zero tolerance: YES.

Reproduce with `python3 tools/validate_pck.py discoveries/2026-09-14-23-csqv-119/csqv119.pck`.

## Note on LP precision (exact-LP vs scipy)

The `.pck` uses the native C++ exact LP (`cpp/csqv/lp.hpp`), which is strictly
feasible to about 1e-17 on these centers and scores 5.749508. The Python (scipy
HiGHS) LP on the SAME centers leaves about 9.7e-8 of overlap on the raw radii and
so scores 5.749503 after the feasibility shrink. The gap 5.749508 vs 5.749503 is
LP PRECISION on the same centers, not a better arrangement. BOTH scores beat the
live record (5.749503 and 5.749508 both exceed 5.748642).

## Reproduce

The exact-LP radii and the strictly-feasible 15-dp `.pck` (from the native C++
worker, built per `cpp/README.md`):

```
csqv_worker 119 <budget_s> 0 <threads> <out_dir> "Arnold Castro" 15 --record 5.748642197922
```

The `--record <live_sum>` flag is optional and advisory (it labels the console
with the gap and prints `*** RECORD BEATEN ***`); pass the live best-known from
`sumradii.txt`. The search does not use it.

It resumes from `<out_dir>/n119-s<base_seed>-best.txt`, writes the champion
centers to `<out_dir>/n119-best.txt`, and emits `<out_dir>/csqv119.pck` at 15 dp.
Re-verify the round-trip at zero tolerance before trusting any value. Records move
daily; re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` at verification
time.

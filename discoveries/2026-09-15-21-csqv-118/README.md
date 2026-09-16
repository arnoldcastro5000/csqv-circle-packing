# Discovery: CSQV N=118 circle packing (beats best-known, exact-LP radii)

An improved packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 118. It beats the published
best-known.

IMPORTANT, honest provenance: this packing was found by a HAND-WRITTEN numerical
search, NOT by the LLM-evolution loop. It is a solver-tooling result that proves
the target is beatable. Making it a true Crucible-Loop (LLM-evolved) discovery is
future work.

## Provenance

| field | value |
|---|---|
| result files | `csqv118.pck` (15 dp, strictly feasible), `n118-best.txt` (the centers) |
| centers found by | the native C++ cold multi-start worker (`cpp/csqv/worker.cpp`, `search.hpp`), base_seed 2, thread seed 2, no LLM |
| radii + `.pck` by | the native C++ exact LP + rounding-aware emit in `cpp/` (`cpp/csqv/lp.hpp`, `worker.cpp::emit_pck`) |
| terminal squeeze | the exact-LP center squeeze via the analytic dual gradient (`cpp/csqv/polish.hpp::center_polish_analytic`); resume-squeezed to a strictly feasible optimum |
| repo commit (code base) | branch `feat/csqv-analytic-dual-gradient`, base 8c11af7 |
| primitives used | exact radii-by-LP for fixed centers; independent verifier recomputes containment, overlap, and the sum |
| submitted to Packomania | NO. Held; a submission is the operator's call. |

## Result

| quantity | value |
|---|---|
| sum of radii (this packing, as re-parsed from the 15-dp `.pck`) | 5.724587011021 |
| in-memory C++ champion | 5.724587011138 |
| Packomania live best-known (fetched 2026-09-16) | 5.723940934671 |
| margin over the live record | +0.000646076 (about +0.0113%) |

The live best-known was fetched 2026-09-16 from `csqv/csqv.html`: N=118 =
5.723940934671, credited to Zeeshan Tariq [20] (set 12-Sep-2026), unchanged from
the 2026-09-14 value except a last-digit rounding. N=118 is actively contested.

## Verification (independent, from the .pck)

`tools/validate_pck.py` re-parsed the exact `.pck` (radii read AS WRITTEN, no
re-derivation) at ZERO tolerance:

- 118 circles, radii in [0.031955, 0.061639], sorted by increasing radius.
- Containment (centered side-1 square, frame center 0,0): max wall violation
  -1.0e-12 (inside, robust cushion).
- Non-overlap (all 6903 pairs): max pair overlap -2.0e-12 (no overlap, robust
  cushion).
- The first line (largest radius) equals the largest radius in the rows.
- STRICTLY FEASIBLE at zero tolerance: YES.

Reproduce with `python3 tools/validate_pck.py discoveries/2026-09-15-21-csqv-118/csqv118.pck`.

## Reproduce

The exact-LP radii and the strictly-feasible 15-dp `.pck` (from the native C++
worker, built per `cpp/README.md`):

```
csqv_worker 118 <budget_s> 2 <threads> <out_dir> "Arnold Castro" 15 --record 5.723940934671
```

The `--record <live_sum>` flag is optional and advisory (it labels the console
with the gap and prints `*** RECORD BEATEN ***`); pass the live best-known from
`sumradii.txt`. The search does not use it.

It resumes from `<out_dir>/n118-s<base_seed>-best.txt`, writes the champion
centers to `<out_dir>/n118-best.txt`, and emits `<out_dir>/csqv118.pck` at 15 dp.
The worker squeezes the champion at clean exit and on resume; a run stopped early
must be relaunched (or resume-squeezed) before the `.pck` is trusted. Records move
daily; re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` at verification
time.

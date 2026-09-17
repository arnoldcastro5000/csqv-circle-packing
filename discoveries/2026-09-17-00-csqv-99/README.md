# Discovery: CSQV N=99 circle packing (beats best-known, exact-LP radii)

An improved packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 99. It beats the live best-known.

IMPORTANT, honest provenance: this packing was found by an AI-authored numerical
search. It is a solver-tooling result that proves the target is beatable.

## Provenance

| field | value |
|---|---|
| result files | `csqv99.pck` (15 dp, strictly feasible, SEALED at its jammed value), `n99-sealed.txt` (the jammed centers, the `.pck` source), `n99-best.txt` (the pre-seal worker centers) |
| centers found by | the native C++ cold multi-start worker (`cpp/csqv/worker.cpp`), base_seed 0, thread seed 0, no LLM |
| worker mode | `--count` (basin-dedup OFF). First seal from the high-N A/B sweep (N=99..110). Here `--count` was the ONLY mode to beat live: count 5.234266579 (beats) vs spread 5.231965147 (misses, +0.026%). Spread led distinct basins 37817 vs 34240, but its best basin fell short; count found the record-beating basin. |
| radii + `.pck` by | the exact-LP radii + rounding-aware feasibility shrink (the parse-back safety loop) |
| terminal squeeze | SEALED 2026-09-17: the first-order-squeezed champion was NOT jammed (Donev margin 9.54), so the jamming optimizer (`jam_slp`, a trust-region sequential-LP flex ascent) drove it to its same-basin jammed optimum (margin 0), then the `.pck` was emitted + zero-tol validated |
| primitives used | exact radii-by-LP for fixed centers; independent verifier recomputes containment, overlap, and the sum |
| submitted to Packomania | NO. Held; a submission is the operator's call. |

## Result

| quantity | value |
|---|---|
| sum of radii (SEALED, re-parsed from the 15-dp `.pck`) | 5.234274484442986 |
| sum before sealing (first-order squeezed) | 5.234262200774 |
| sealing gain | +1.23e-5 (the champion was un-jammed) |
| Packomania live best-known (fetched 2026-09-16) | 5.233333680412 |
| margin over the live record | +0.000940804031 (about +0.01798%) |

The live best-known was fetched 2026-09-16 from `csqv/csqv.html` (page updated
16-Sep-2026): N=99 = 5.233333680412, credited to Byron Tasseff [13] (private
communication, Aug 2026), not marked proven-optimal.

## Verification (independent, from the .pck)

`tools/validate_pck.py` re-parsed the exact `.pck` (radii read AS WRITTEN, no
re-derivation) at ZERO tolerance:

- 99 circles, radii in [0.036836, 0.069341], sorted by increasing radius.
- Containment (centered side-1 square, frame center 0,0): worst wall violation
  -1.0e-12 (inside, robust cushion).
- Non-overlap (all 4851 pairs): worst pair overlap -2.0e-12 (no overlap, robust
  cushion).
- The first line (largest radius) equals the largest radius in the rows.
- STRICTLY FEASIBLE at zero tolerance: YES.

Reproduce with `python3 tools/validate_pck.py discoveries/2026-09-17-00-csqv-99/csqv99.pck`.

## Reproduce

The exact-LP radii and the strictly-feasible 15-dp `.pck` (from the native C++
worker, built per `cpp/README.md`):

```
csqv_worker 99 <budget_s> 0 <threads> <out_dir> "Arnold Castro" 15 --record 5.233333680412
```

The `--record <live_sum>` flag is optional and advisory (it labels the console
with the gap and prints `*** RECORD BEATEN ***`); pass the live best-known from
`sumradii.txt`. The search does not use it.

It resumes from `<out_dir>/n99-s<base_seed>-best.txt`, writes the champion centers
to `<out_dir>/n99-best.txt`, and emits `<out_dir>/csqv99.pck` at 15 dp. The worker
squeezes the champion at clean exit and on resume; a run stopped early must be
relaunched (or resume-squeezed) before the `.pck` is trusted. The champion was
then SEALED (jammed) before recording. Records move daily; re-fetch
`https://www.packomania.com/csqv/txt/sumradii.txt` at verification time.

# Discovery: CSQV N=143 circle packing (feasible; N not in the record book)

A feasible packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 143.

Packomania does NOT publish a best-known for N = 143. The master table
`csqv/txt/sumradii.txt` is sparse in this range (it jumps 140 to 144), and the
detail page `csqv/txt/csqv143.txt` returns 404 (fetched live 2026-09-14). So this
is NOT a record beat; it is a verified feasible packing for an N the record book
does not track. A first-entry submission is possible but is the operator's call.

IMPORTANT, honest provenance: this packing was found by a numerical search with the
native C++ worker.

## Provenance

| field | value |
|---|---|
| result files | `csqv143.pck` (15 dp, strictly feasible), `n143-best.txt` (the centers) |
| centers found by | the native C++ cold multi-start worker (`cpp/csqv/worker.cpp`, `search.hpp`), base_seed 0, no LLM |
| radii + `.pck` by | the native C++ exact LP + rounding-aware emit in `cpp/` (`cpp/csqv/lp.hpp`, `worker.cpp::emit_pck`) |
| repo commit (code base) | e297555 |
| primitives used | exact radii-by-LP for fixed centers; independent verifier recomputes containment, overlap, and the sum |
| Packomania best-known | none published for N = 143 |
| submitted to Packomania | NO. Held; a first-entry submission is the operator's call. |

## Result

| quantity | value |
|---|---|
| sum of radii (this packing, as re-parsed from the 15-dp `.pck`) | 6.312963889026136 |
| sum of radii (in-memory C++ champion) | 6.312963889169 |
| Packomania live best-known (fetched 2026-09-14) | none (N not tracked) |

## Verification (independent, from the .pck)

Parsed the exact `.pck` and recomputed geometry from scratch (pure numpy, radii
read AS WRITTEN, no re-derivation), at ZERO tolerance:

- 143 circles, radii in [0.029, 0.058], sorted by increasing radius.
- Containment (centered side-1 square, frame center 0,0): max wall violation
  -1.0e-12 (inside, robust cushion).
- Non-overlap (all 10153 pairs): max pair overlap -2.0e-12 (no overlap, robust
  cushion).
- The first line (largest radius) equals the largest radius in the rows.
- STRICTLY FEASIBLE at zero tolerance: YES.

Note on LP precision: the `.pck` uses the native C++ exact LP, strictly feasible to
about 1e-17. The Python (scipy HiGHS) LP on the same centers leaves about 1e-7 of
overlap on the raw radii and scores marginally lower after the feasibility shrink.
The gap is LP precision on the same centers, not a different arrangement.

## Reproduce

```
csqv_worker 143 <budget_s> 0 <threads> <out_dir> "Arnold Castro" 15
```

It resumes from `<out_dir>/n143-s0-best.txt`, writes the champion centers to
`<out_dir>/n143-best.txt`, and emits `<out_dir>/csqv143.pck` at 15 dp. Re-verify the
round-trip at zero tolerance before trusting any value. The record book moves daily;
re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` at verification time to
confirm N = 143 is still untracked.

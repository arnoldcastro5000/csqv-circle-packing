# Discovery: CSQV N=121 circle packing (record-beating)

An improved packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 121. It beats the published
best-known and was submitted to Packomania on 2026-09-14.

IMPORTANT, honest provenance: this packing was found by an AI-authored numerical
search. It is a solver-tooling result that proves the target is beatable.

## Provenance

| field | value |
|---|---|
| result files | `csqv121.pck` (submitted), plus the coordinates it contains |
| generating code | `run_scalable.py` + `scalable_polish.py` (this directory) |
| method | AI-authored cold multi-start (no LLM) |
| repo commit (code base) | 34a6633 |
| primitives used | `problems.csqv.problem.radii_lp` (exact radii for fixed centers), `problems.csqv.data`/verifier |
| submitted to | E. Specht, packomania.com, 2026-09-14 |

## Result

| quantity | value |
|---|---|
| sum of radii (this packing) | 5.799069987171 |
| Packomania best-known at submit time | 5.797468812113 |
| margin | +0.001601175 (about +0.028%) |

The search later reached a marginally higher champion (5.799070342829); the
SUBMITTED value is the slightly lower 5.799069987171 in `csqv121.pck`. Both beat
the record.

## Verification (independent, from the .pck)

Parsed the exact `.pck` and recomputed geometry from scratch:

- 121 circles, radii in [0.031, 0.063], sorted by increasing radius.
- Containment (centered side-1 square): min wall slack +3.27e-8 (inside).
- Non-overlap (all 7260 pairs): min pair gap +4.03e-13 (no overlap).
- Strictly feasible at 12 decimals and at 1e-12; the reference verifier agrees
  (max_out 0.0, max_overlap 0.0).

Use 12-decimal precision. A 15-decimal emit surfaced a sub-epsilon (-3.9e-16)
contact, so 12 dp is the correct, strictly-feasible precision (and matches
Packomania's own csqv files).

## Method

For FIXED centers, `radii_lp` returns the EXACT optimal radii (a linear program),
so all slack lives in the center ARRANGEMENT. N=121 is a Specht square
CONSTRUCTION (a pattern rule plus a light solver), so it is construction-soft and
beatable by a different-basin COLD search. The chain, per restart:

1. construction (grid / jittered grid / hexagonal / random),
2. grow-and-push relaxation (inflate radii, push overlapping centers apart, clip
   to walls, iterate) to spread centers cheaply,
3. `scalable_polish`: L-BFGS-B on a smooth penalty objective (maximize sum of
   radii minus quadratic overlap and wall penalties, penalty continuation over
   rising lambda), then `radii_lp` projects to exact feasibility,
4. perturb-the-best restarts for diversity; keep the best feasible packing.

Why `scalable_polish` and not the repo's dense SLSQP `polish`: the dense polish
builds a dense jacobian (~5.6 GB and >10 min per call at N=484; it chokes on the
~n^2/2 constraints), so it is infeasible at large N. The L-BFGS-B penalty method
is O(n) memory (~105 MB at N=121) and sub-second, with comparable quality. This
is what made a record search tractable.

## Reproduce

From the repo root, with scipy installed:

```
PYTHONPATH=. python3 -u discoveries/2026-09-13-20-csqv-121/run_scalable.py 121 <budget_s> <seed>
```

It writes the best packing to `n121-s<seed>-best.txt` (centered frame) and prints
`*** RECORD BEATEN ***` when the sum exceeds the live record. Verify the `.pck`
with `tools/validate_pck.py` at zero tolerance before trusting any value. Records
move daily; re-fetch
`https://www.packomania.com/csqv/txt/sumradii.txt` at verification time.

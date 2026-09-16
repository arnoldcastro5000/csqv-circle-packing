# Discovery: CSQV N=121 circle packing (improved, exact-LP radii)

An improved packing for the Packomania CSQV problem (variable-radii circles in a
unit square, maximize the sum of radii) at N = 121. It beats the published
best-known. This is our BEST N=121 record.

This packing has the SAME centers as the earlier discovery
`discoveries/2026-09-13-20-csqv-121/`, but tighter, exactly-optimal radii. The
earlier discovery scored the centers with the Python (scipy HiGHS) LP, which
leaves about 6.8e-8 of overlap and so scores 5.799070 after the feasibility
shrink. The native C++ exact LP (`cpp/csqv/lp.hpp`) is strictly feasible to about
1e-17 on the same centers, so it scores HIGHER (5.799074) and stays strictly
feasible. The gap 5.799074 vs 5.799070 is LP PRECISION on the same centers, not a
better arrangement.

IMPORTANT, honest provenance: this packing was found by an AI-authored numerical
search, NOT by the LLM-evolution loop. It is a solver-tooling result that proves
the target is beatable. Making it a true Crucible-Loop (LLM-evolved) discovery is
future work.

## Provenance

| field | value |
|---|---|
| result files | `csqv121.pck` (15 dp, strictly feasible), `n121-best.txt` (the centers) |
| centers found by | `discoveries/2026-09-13-20-csqv-121/run_scalable.py` (AI-authored cold multi-start, seed 0, no LLM) |
| radii + `.pck` by | the native C++ exact LP + rounding-aware emit in `cpp/` (`cpp/csqv/lp.hpp`, `worker.cpp::emit_pck`) |
| repo commit (code base) | see `git log` at upload time |
| primitives used | exact radii-by-LP for fixed centers; independent verifier recomputes containment, overlap, and the sum |
| submitted to Packomania | NO. The earlier 12-dp 5.799069987171 file was submitted 2026-09-14; this improved file is held pending that outcome. |

## Result

| quantity | value |
|---|---|
| sum of radii (this packing, as re-parsed from the 15-dp `.pck`) | 5.799074421368918 |
| sum of radii (in-memory C++ champion) | 5.799074421490 |
| Packomania live best-known (fetched 2026-09-14) | 5.797468812113 |
| margin over the live record | +0.001605609 (about +0.0277%) |
| vs the SUBMITTED 5.799069987171 | +4.43e-6 higher (same centers, tighter radii) |

## Verification (independent, from the .pck)

Parsed the exact `.pck` and recomputed geometry from scratch, at ZERO tolerance:

- 121 circles, radii in [0.031, 0.063], sorted by increasing radius.
- Containment (centered side-1 square, frame center 0,0): max wall violation
  -1.0e-12 (inside, robust cushion).
- Non-overlap (all 7260 pairs): max pair overlap -2.0e-12 (no overlap, robust
  cushion).
- The first line (largest radius) equals the largest radius in the rows.
- STRICTLY FEASIBLE at zero tolerance: YES.

The 15-dp emit is safe here because the emit checks feasibility on the values AS
RE-PARSED from the file and shrinks to a robust -1e-12 cushion. This supersedes
the earlier discovery's "use 12 dp only" caution, which was a fragile zero-margin
emit.

## Reproduce

The centers (from the repo root, with scipy installed):

```
PYTHONPATH=. python3 -u discoveries/2026-09-13-20-csqv-121/run_scalable.py 121 <budget_s> 0
```

The exact-LP radii and the strictly-feasible 15-dp `.pck` (from the native C++
worker, built per `cpp/README.md`):

```
csqv_worker 121 <budget_s> 200 <threads> <out_dir> "Arnold Castro" 15 --record 5.797468812113
```

The `--record <live_sum>` flag is optional and advisory (it labels the console
with the gap and prints `*** RECORD BEATEN ***`); pass the live best-known from
`sumradii.txt`. The search does not use it.

It resumes from `<out_dir>/n121-s200-best.txt`, writes the champion centers to
`<out_dir>/n121-best.txt`, and emits `<out_dir>/csqv121.pck` at 15 dp. Re-verify
the round-trip at zero tolerance before trusting any value. Records move daily;
re-fetch `https://www.packomania.com/csqv/txt/sumradii.txt` at verification time.

## Note on the two records in this repo

- `discoveries/2026-09-13-20-csqv-121/` = the SUBMITTED packing (5.799069987171,
  12 dp, Python scipy LP).
- `discoveries/2026-09-14-00-csqv-121/` (this dir) = the BEST packing
  (5.799074421369, 15 dp, C++ exact LP), same centers, held pending the
  submission outcome.

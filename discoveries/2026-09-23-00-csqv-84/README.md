# Discovery: CSQV N=84 circle packing (published record)

A packing for the Packomania CSQV problem (variable-radii circles in a unit
square, maximize the sum of radii) at N = 84. It was submitted, accepted, and
published, credited to Arnold Castro. As of 2026-09-23 it is the only CSQV
record in this repository still credited to us on the live site.

IMPORTANT, honest provenance: found by the AI-authored numerical C++ worker.

## Provenance

| field | value |
|---|---|
| result file | `csqv84.pck` (strictly feasible, submittable, SEALED at its jammed value) |
| source champion | `n84-champion-spread.txt` (raw worker champion, centered frame) |
| sealed champion | `n84-sealed.txt` (jammed centers, the `.pck` source) |
| origin run | `cpp/build/84-20260920-024921-spread/`, base seed 0, a host run in the N=89..80 descending sweep |
| worker mode | `--spread` (basin-dedup/seed-spread ON) |
| method | AI-authored cold multi-start C++ worker + the exact-LP center squeeze `center_polish_analytic`, then SEALED in-run by the worker's terminal jamming optimizer (`jam_slp`), the `.pck` emitted and zero-tol validated |
| submitted | YES, accepted and published (credited to Arnold Castro) |

## Result

| quantity | value |
|---|---|
| sum of radii (SEALED, `.pck` as written) | 4.816500878221171 |
| raw worker champion sum (pre-seal) | 4.816493273962 |
| sealed sum (jammed, worker frame) | 4.816500878305 |
| sealing gain (jammed vs pre-seal baseline) | +7.60e-06 (the champion was un-jammed) |
| emit cost (`.pck` as written vs jammed) | -8.4e-11 |
| Packomania live best-known (fetched 2026-09-23) | 4.816500879756 (Arnold Castro, not proven-optimal) |

The live record was re-fetched 2026-09-23: N=84 = 4.816500879756, credited to
Arnold Castro. The site re-optimized the accepted submission by about +1.5e-9
(float-scale). This repository's full seal chain (the jamming optimizer plus a
guarded contact-graph Newton step) independently reaches 4.816500879755, matching
the published value to 1e-12, which confirms the site's re-optimizer is the
same-basin jamming squeeze. Re-fetch again before trusting any comparison.

## Verification (independent, from the `.pck` at zero tolerance)

Parsed the exact `.pck` and recomputed geometry from scratch (centered side-1
square, center 0,0, half-side 0.5), radii AS WRITTEN, zero tolerance:

- 84 circles, all radii positive, sorted by increasing radius.
- max wall violation = -1.000e-12 (every circle inside the square).
- max overlap = -2.000e-12 (no pair overlaps).
- STRICTLY FEASIBLE as written: YES (`validate_pck.py` VALID at zero tolerance).

## Next

- N=84 is published (credited to Arnold Castro); do not resubmit.
- Records move daily; re-fetch the live best-known before any comparison.

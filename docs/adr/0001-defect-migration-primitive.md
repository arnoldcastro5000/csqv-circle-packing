# Add a third fixed primitive: defect-migration

The CSQV solver contract fixed exactly two trusted primitives on `ctx`, `radii_lp`
and `polish`, kept out of the caller and the scoring path. We add a third, `defect_move`, a
combinatorial move that removes the k weakest circles and reinserts them into the largest empty
holes, then resets radii and polishes once.

## Why

The N=60 and N=90 calibrations showed our search re-derives strong basins but never beats them.
An iteration-depth probe ruled out polish under-convergence (raising SLSQP iterations 10x gave
zero gain). A local optimizer only reaches the floor of the basin its start sits in, so a
stronger or faster polish cannot break the plateau. Defect-migration is the one available move
that changes the COMBINATORIAL contact graph, so it can reach a basin the center-only
perturbation walk never visits. It is the classical record-chasing move (Grosso; Addis,
Locatelli, Schoen) and appeared in the evolved solvers reported in arXiv:2609.05093.

## The trade-off

We chose a FIXED, trusted primitive over caller-side shell code. A primitive is ungameable
(out of the scoring path), unit-testable with no model access, and deterministic given
(packing, rng, k). The cost is a wider contract than the deliberate two-primitive minimum and
one more piece of trusted surface to maintain. We accepted this because robustness and testability
outweigh contract minimalism here, consistent with the repo's quality-first mandate. The caller
still decides WHEN and how often to call it.

## Status

Accepted. If a defect-move-enabled re-calibration at N=90 does not surface a beating basin, the
effort declares the proven machinery the deliverable and stops (option d); this primitive stays
as trusted solver surface regardless.

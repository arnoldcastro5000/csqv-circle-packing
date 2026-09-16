# Basin dedup gates the seeds only, never the depth

The basin-dedup worker (tickets 41, 42) shares a [[Basin fingerprint]] across the worker threads
so they cover unique basins. This ADR fixes WHERE the fingerprint acts: it gates the SEEDS (the
fresh cold constructions) and never the perturbation depth inside a basin.

## Why

The count of distinct optima for N<=100 is super-exponential, so "exhaust all basins" means
maximum unique coverage, not enumeration. The natural but wrong design is to hard-skip ANY search
that lands in a seen basin. Records need DEPTH: the Monotonic Basin Hopping perturbation walk and
`defect_move` dig a promising basin out to its floor, and this depth is what produces a record. If
the gate drops a perturbation neighbour the instant its fingerprint is seen, it kills the funnel
walk, coverage rises, and P(record basin) falls. The metaheuristic evidence (MLSL, repulsion
multistart, taboo, clearing, Sobol multistart) all steer WHERE a search starts and never cap
depth, while basin hopping (Wales and Doye; Addis, Locatelli, Schoen) shows the perturbation walk
is the record-producing mechanism. See `docs/research/phase-2-csqv-basin-dedup-execution.md`.

## The trade-off

We chose SEEDS-ONLY dedup plus protected depth over pure coverage dedup. The cost is that two
threads can still spend depth-effort in the same basin if both start near it despite the Sobol
seed spreading, so raw distinct-basin coverage is not maximal. We accepted this because records,
not coverage, are the goal, and depth is where records come from. Coverage is only a proxy.

The evidence for the shared-hash dedup itself lifting P(record basin)/budget over an
already-diverse multistart is thin for this objective, so the fingerprint is a cheap add-on on top
of the well-supported core (Sobol [[Seed spreading]] + protected depth), and the full build is
gated on an empirical A/B (ticket 42 acceptance).

## Status

Accepted. Design decision for the basin-dedup worker; the verifier and the fixed-primitive
contract (ADR 0001) are unchanged. The build is [ticket 43](../../.scratch/csqv/issues/43-build-basin-dedup-seed-spread-worker.md).
Revisit if the A/B shows coverage dedup trades away depth even under the seeds-only gate.

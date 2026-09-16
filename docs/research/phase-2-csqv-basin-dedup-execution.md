# Phase 2 CSQV: how to execute cross-thread basin deduplication (N<=100)

Research date: 2026-09-15. This note answers [ticket 41](../../.scratch/csqv/issues/41-cross-thread-basin-deduplication.md).
It feeds the GO/NO-GO decision and spec in [ticket 42](../../.scratch/csqv/issues/42-basin-worker-go-no-go-and-spec.md).
It does not resolve ticket 42.

Read first for context: the basin definition in [CONTEXT.md](../../CONTEXT.md) (terms Basin,
Defect-migration, Multi-start), the current worker
[cpp/csqv/worker.cpp](../../cpp/csqv/worker.cpp), and the record-methods note
[phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md).

CSQV packs N variable-radius circles in a unit square. The objective maximizes the sum of the
radii. For fixed centers the optimal radii solve an exact LP. A local search reaches only its
start basin's optimum. The record is the single deepest basin. The count of distinct optima for
N<=100 is super-exponential, so "exhaust all basins" means MAXIMUM UNIQUE COVERAGE, not
enumeration. This note judges every option by P(record basin) / wall-clock budget, not by raw
distinct-basin count.

Method. This note uses primary sources: peer-reviewed papers, arXiv, first-party docs, and the
local worker code. Each non-obvious claim carries an inline URL or DOI. Some publisher full texts
sit behind a paywall gate. Those claims are marked "gate" and are corroborated by the abstract,
the index record, and a corroborating open source, not by a verbatim full-text quote.

## How the current worker picks seeds (the problem)

`worker_loop` in [cpp/csqv/worker.cpp:195-263](../../cpp/csqv/worker.cpp) runs one independent
thread per seed. Each restart either perturbs the thread's own best (`bx, by`, probability 0.6,
Gaussian scale 0.006 to 0.03) or builds a fresh construction `kinds[restarts % 4]` (grid, jitter,
hex, random). No thread reads another thread's basin. So two threads that build the same
construction kind fall into the same dominant basin, and a warm perturbation collapses back into
the incumbent basin (the map records N=60 and N=90 warm walks that re-derive the incumbent). Three
threads give less than 3x distinct coverage. The dedup idea shares a faithful combinatorial basin
FINGERPRINT so a thread can skip a basin another thread already found.

## 1. Fingerprint: contact graph, not quantized radii, never raw floats

**Raw-float fingerprints fail. Confirmed.** Two descents to the SAME basin stop at coordinates
that differ by the solver tolerance (this worker converges contacts to about 1e-9 to 1e-12; the
verifier margins are overlap 1.9e-12 and boundary 2.8e-17, per the map). Any hash of the float
coordinates then differs on every run, so every run looks distinct and the dedup gate never fires.
CONTEXT.md already fixes the correct rule: "Two searches reach the SAME basin when their converged
contact graphs match, not when their float coordinates match" ([CONTEXT.md](../../CONTEXT.md), term
Basin). So the fingerprint must be COMBINATORIAL.

**Contact graph (recommended).** Read the converged packing. Add an edge (i, j) when
`dist(c_i, c_j) - (r_i + r_j) <= tol`. Add a wall label to node i when i touches a side within
`tol`. This is exactly the physical contact network that the jamming literature reads off a
packing; Donev et al. build and analyze the same contact network to classify a packing and to
flag rattlers (free circles with no contacts)
([J. Comput. Phys. 197 (2004) 139-166, DOI 10.1016/j.jcp.2003.11.022](https://doi.org/10.1016/j.jcp.2003.11.022);
open PDF [math.nyu.edu](https://math.nyu.edu/inmemoriam/donev/Packing/Jamming_LP.pdf)).

**Contact tolerance.** A true contact converges to the solver floor (about 1e-9 to 1e-12 here). A
near-miss (a gap that is not a contact) sits at a macroscopic distance, orders of magnitude
larger. So a single `tol` near 1e-6 separates the two cleanly: it sits in the wide gap between the
contact floor and the smallest real gap. The choice is not precision-critical because the two
populations are orders apart. Rattlers get degree 0 and no wall label; keep them as isolated nodes
so a packing with a free circle still fingerprints stably.

**Canonical hash and collision behaviour.** Color each node by (wall-contact flags, quantized
radius bucket). Compute a canonical form of the colored graph, then hash it to 64 or 128 bits.
- An exact canonical labeling (nauty/Traces, McKay and Piperno, "Practical graph isomorphism, II",
  [J. Symbolic Computation 60 (2014) 94-112, DOI 10.1016/j.jsc.2013.09.003](https://doi.org/10.1016/j.jsc.2013.09.003))
  gives CF(X) = CF(Y) exactly when X and Y are isomorphic. Then the only collisions come from hash
  truncation, and a 128-bit hash makes a birthday collision negligible over the ~1e5 to 1e6
  fingerprints of a long run.
- A Weisfeiler-Leman color-refinement HASH is cheaper and never false-splits (same basin always
  hashes the same). It can rarely FALSE-MERGE two non-isomorphic graphs that WL cannot separate
  ([WL tutorial, NSF PAR 10299993](https://par.nsf.gov/servlets/purl/10299993)). Contact graphs of
  planar packings are near-planar and the node colors (radius bucket, wall flags) break most ties,
  so a false merge is rare. A false merge only skips a genuinely new basin, which is the safe
  direction (it never discards a found record; the incumbent is already saved).
- The canonical form is invariant under the square's D4 symmetry (rotations and reflections)
  automatically, because a rotated or reflected packing has an isomorphic colored graph. So
  symmetric copies of one basin hash together, which is correct: they are the same basin.

**Quantized sorted-radii (reject as the standalone fingerprint).** Round each radius, sort, hash
the list. It is combinatorial and cheap, and it is rotation/reflection invariant. But it
OVER-MERGES: two genuinely different arrangements can share the same radii multiset, so distinct
basins collide and the gate wrongly skips a new basin. It is also UNSTABLE at a bucket boundary:
two runs of one basin can straddle a rounding edge and split. It captures the optimum's value
signature, not its arrangement. Use it only as a cheap necessary-not-sufficient PRE-FILTER (equal
radii vectors are a precondition for the same basin), then confirm with the contact graph.

**Recommendation (1).** Fingerprint = canonical hash of the tolerance-thresholded, node-colored
CONTACT GRAPH. Use quantized sorted-radii only as an optional cheap pre-screen.

## 2. Shared-set architecture: a single lock is enough; Bloom is the fixed-footprint option

**Throughput reality.** The map estimates about 6 threads x about 1 restart/s, so about 6
fingerprint operations per second aggregate. A hash insert or lookup costs well under a
microsecond. A single global lock is therefore held for a negligible fraction of a second per
second of wall clock. The "global lock throttles throughput" worry is real only at millions of ops
per second; at 6 ops/s it is a non-issue. So a single `std::mutex` around one shared
`std::unordered_set<uint64_t/uint128>` is correct and simplest. Sharding or per-thread local sets
with a periodic merge add complexity for a contention problem that does not exist at this rate.

**Memory at N=100 (confirm the map estimate).** A long run is on the order of 1e5 to 1e6 restarts
(6 threads x 1/s x hours). Store only the 8 to 16 byte hash, not the full graph.
- Exact set: about 1e6 entries x (16 byte key + ~40 byte node overhead in `std::unordered_set`) is
  roughly 50 MB worst case; an open-addressing flat set is roughly 24 MB. So the map's "tens of MB"
  estimate is confirmed as a safe upper bound for the exact set.
- Bloom filter: for 1e6 elements at a 1% false-positive rate you need about 9.6 bits per element,
  so about 1.2 MB, a fixed footprint set at build time. The false-positive formula is
  p ~= (1 - e^(-kn/m))^k with optimal k = (m/n) ln 2; 1% needs ~9.6 bits/element and 0.1% needs
  ~14.4 bits/element (Bloom, "Space/time trade-offs in hash coding with allowable errors",
  [CACM 13(7) 1970, 422-426, DOI 10.1145/362686.362692](https://doi.org/10.1145/362686.362692)).
- A Bloom false positive means the worker RARELY skips a genuinely new basin. That is the safe
  direction and it costs one construction, so a 1% rate is cheap insurance for a two-order-smaller,
  fixed footprint.

**Recommendation (2).** Single global mutex + one shared set. Use an exact hash set (tens of MB,
audit-friendly) unless a fixed memory ceiling is wanted, in which case a Bloom filter at 1% FP
(~1.2 MB) is a clean drop-in. Sharding is not needed at this restart rate.

## 3. The dedup gate: gate SEEDS ONLY, protect depth

The literature splits cleanly along the coverage-vs-depth line the ticket names.

**Dedup the SEEDS (coverage side).** Multi-Level Single Linkage (MLSL) is the canonical
seed-dedup method. It clusters sampled points into regions of attraction and starts ONE local
search per region, and it provably finds all relevant local minima with the smallest possible
number of local searches (Rinnooy Kan and Timmer, "Stochastic global optimization methods, Part I:
Clustering methods" and "Part II: Multi level methods", Math. Programming 39 (1987) 27-56 and
57-78, [DOI 10.1007/BF02592070](https://doi.org/10.1007/BF02592070) open PDF, and
[DOI 10.1007/BF02592071](https://doi.org/10.1007/BF02592071)). This is exactly "skip a fresh
construction that would re-descend a seen basin." Crucially MLSL never caps DEPTH: each started
local search runs fully to convergence. The repulsion multistart method places each new start away
from already-found basins and reports fewer local searches than plain multistart to find the global
optimum ("The repulsion algorithm, a new multistart method for global optimization", Structural
Optimization, [DOI 10.1007/BF01197028](https://doi.org/10.1007/BF01197028); gate, corroborated by
the abstract and index). Taboo/tabu search for continuous minima keeps a memory of visited regions
to avoid re-entrapment, again a seed/steering mechanism, not a depth cap (Cvijovic and Klinowski,
"Taboo search: an approach to the multiple minima problem",
[Science 267 (1995) 664-666, DOI 10.1126/science.267.5198.664](https://doi.org/10.1126/science.267.5198.664)).
Niching by clearing keeps only the best individual per niche and spreads the rest, a
selection/seed diversity mechanism that still protects the depth of each niche winner (Petrowski,
"A clearing procedure as a niching method for genetic algorithms", ICEC'96 798-803,
[DOI 10.1109/ICEC.1996.542703](https://doi.org/10.1109/ICEC.1996.542703)). Low-discrepancy
(Sobol) multistart is the deterministic version of the same idea: it spreads the START points
evenly to avoid clustered, redundant descents, without touching depth (Bratley and Fox,
"Algorithm 659: Implementing Sobol's quasirandom sequence generator",
[ACM TOMS 14(1) 1988, DOI 10.1145/42288.214372](https://doi.org/10.1145/42288.214372)).

**Never dedup the DEPTH (record side).** Basin hopping perturbs the current local minimum and
re-minimizes, exploiting a funnel that funnels toward the global minimum (Wales and Doye, "Global
optimization by basin-hopping...", [J. Phys. Chem. A 101 (1997) 5111-5116, DOI 10.1021/jp970984n](https://doi.org/10.1021/jp970984n);
open [ORA copy](https://ora.ox.ac.uk/objects/uuid:6a36a972-20a7-4173-9943-f90f45c37a2c)). The whole
value of Monotonic Basin Hopping is the perturbation WALK inside and around a promising basin. If
the gate drops a perturbation neighbour the instant its fingerprint is "seen," it kills the funnel
walk that produces the record. Coverage rises and P(record basin) falls, the exact failure the
ticket warns of.

**Recommendation (3).** Gate SEEDS ONLY. Fingerprint each converged optimum and record it. Before
a thread commits to a fresh cold CONSTRUCTION, skip and redraw if that construction's basin is
already recorded (steer WHERE you start). Always fully refine (center_polish, defect_move) and
perturb-walk any basin a thread is actually in, and do not gate the perturbation neighbours. This
is MLSL/repulsion logic on the seeds plus MBH depth inside the basin.

## 4. Core evidence: honest reading for P(record basin)/budget

**What the literature supports strongly.** On a funneling landscape, perturbation-based search
(MBH/PBH) beats plain independent multistart by a wide margin. Disk packing in a square is
conjectured to have exactly this funneling landscape, and the MBH/PBH approach improved putative
optima for n<=130 in as many as 32 instances (Addis, Locatelli, Schoen, "Disk Packing in a Square:
A New Global Optimization Approach",
[optimization-online PDF](http://www.optimization-online.org/DB_FILE/2005/10/1229.pdf)). Population
Basin Hopping maintains a dissimilarity measure between population members to keep coverage wide
and is very efficient on hard packing problems (Grosso, Locatelli, Schoen, "A population-based
approach for hard global optimization problems based on dissimilarity measures", Math. Programming
110 (2007) 373-404, [DOI 10.1007/s10107-006-0006-3](https://doi.org/10.1007/s10107-006-0006-3);
open [optimization-online record](https://optimization-online.org/2005/02/1056/); dissimilarity
detail in Comput. Optim. Appl., [DOI 10.1007/s10589-008-9194-5](https://doi.org/10.1007/s10589-008-9194-5),
open [PDF](http://www.optimization-online.org/DB_FILE/2007/11/1841.pdf)). MLSL proves that
seed-dedup (start one search per region) reaches all minima with the fewest local searches, which
is strictly higher distinct-coverage per unit budget than plain multistart. Leary shows the funnel
structure is what lets perturbation search win over random restart (Leary, "Global optimization on
funneling landscapes", J. Global Optimization 18 (2000) 367-383,
[DOI 10.1023/A:1026500301312](https://doi.org/10.1023/A:1026500301312); gate, corroborated by
abstract and index).

**Where the literature is thin (be honest).** None of these papers measure the destination metric
directly: P(record basin) per wall-clock second for a CROSS-THREAD SHARED-FINGERPRINT dedup. The
MBH/PBH gains come from (a) perturbation on a funnel and (b) a population dissimilarity measure,
not from a hash set that forbids revisiting. PBH keeps DISTINCT members by a soft dissimilarity
distance, not by an exact combinatorial hash that hard-skips a basin. So the exact-hash dedup
mechanism's marginal value ON TOP of an already-diverse independent multistart is unmeasured for
this objective. Coverage is only a proxy: the record lives in ONE basin, so more distinct basins
help only if the record basin is among the newly covered ones AND each covered basin is still dug
deep enough. If the gate steals budget from depth, P(record) can fall. The evidence therefore
supports the MECHANISM SPLIT (seed-level diversity plus protected perturbation depth) more than it
supports a strict shared-fingerprint dedup as the primary lever. This is precisely the risk that
an empirical A/B at a known-hard N<=100 must settle before a full build.

## 5. Escape-on-hit: perturb with a combinatorial kick, then fall back to a fresh construction

When a fresh construction rolls into a seen basin, the funnel evidence favors a HOP from the seen
optimum over an immediate fully random restart: on a funneling landscape neighbouring basins are
correlated and lead downhill toward the global, so basin hopping reaches a good new basin faster
than random restart (Wales and Doye [DOI 10.1021/jp970984n](https://doi.org/10.1021/jp970984n);
Leary [DOI 10.1023/A:1026500301312](https://doi.org/10.1023/A:1026500301312)). But this worker's
own evidence qualifies "perturb": a SMALL center perturbation collapses back to the same basin
(the map records N=60 and N=90 warm walks re-deriving the incumbent). To hop to an ADJACENT UNSEEN
basin you must change the contact graph, which is exactly what `defect_move` does (remove the k
weakest circles, reinsert into the largest holes; it "changes the combinatorial contact graph, so
it can reach a basin the perturbation walk never visits", [CONTEXT.md](../../CONTEXT.md), term
Defect-migration). PBH already carries this pattern as a large-jump escape when the walk stalls.

**Recommendation (5).** On a hit, first HOP: apply a large kick (`defect_move` or a large-scale
perturbation) from the seen optimum, re-minimize, and check the new fingerprint. If the hop keeps
landing in seen basins after a small bounded number of tries, fall back to a fresh (Sobol-spread or
cold) construction. This is basin-hopping-with-restart, the standard PBH escape, and it reuses the
worker's existing `defect_move` primitive.

## Synthesis for ticket 42

- Fingerprint: canonical hash of the node-colored contact graph at `tol ~ 1e-6`; radii-vector as an
  optional pre-screen; never raw floats.
- Shared set: one global mutex + one shared set; exact set (tens of MB, confirmed) or a 1% Bloom
  filter (~1.2 MB fixed) for a memory ceiling. Contention is not a real problem at ~6 ops/s.
- Gate SEEDS ONLY (MLSL/repulsion/Sobol logic); protect perturbation and defect-move DEPTH (MBH).
- Escape-on-hit: hop with `defect_move` first, then fall back to a fresh construction.
- Evidence: strong for funnel-aware perturbation and for seed-level diversity beating plain
  multistart; THIN for an exact shared-hash dedup lifting P(record)/budget specifically. Coverage
  is a proxy only.

**One-line GO/NO-GO input (ticket 42 decides).** The evidence points to a QUALIFIED GO: build the
SEEDS-ONLY dedup worker (contact-graph fingerprint, single-lock shared set, escape by
`defect_move`, depth protected), but GATE the full build on an empirical P(record basin)/budget A/B
against the current independent-parallel worker at a known-hard N<=100, because the marginal gain
of exact shared-hash dedup over already-diverse multistart is unproven for this objective.

## Sources

- Rinnooy Kan, Timmer. Stochastic global optimization methods, Part I: Clustering methods. Math.
  Programming 39 (1987) 27-56. [DOI 10.1007/BF02592070](https://doi.org/10.1007/BF02592070) (open PDF).
- Rinnooy Kan, Timmer. Part II: Multi level methods. Math. Programming 39 (1987) 57-78.
  [DOI 10.1007/BF02592071](https://doi.org/10.1007/BF02592071).
- Wales, Doye. Global optimization by basin-hopping... J. Phys. Chem. A 101 (1997) 5111-5116.
  [DOI 10.1021/jp970984n](https://doi.org/10.1021/jp970984n);
  [open ORA copy](https://ora.ox.ac.uk/objects/uuid:6a36a972-20a7-4173-9943-f90f45c37a2c).
- Addis, Locatelli, Schoen. Disk Packing in a Square: A New Global Optimization Approach.
  [optimization-online PDF](http://www.optimization-online.org/DB_FILE/2005/10/1229.pdf).
- Grosso, Locatelli, Schoen. A population-based approach for hard global optimization problems based
  on dissimilarity measures. Math. Programming 110 (2007) 373-404.
  [DOI 10.1007/s10107-006-0006-3](https://doi.org/10.1007/s10107-006-0006-3) (gate);
  [open optimization-online record](https://optimization-online.org/2005/02/1056/).
- Dissimilarity measures for population-based global optimization algorithms. Comput. Optim. Appl.
  [DOI 10.1007/s10589-008-9194-5](https://doi.org/10.1007/s10589-008-9194-5);
  [open PDF](http://www.optimization-online.org/DB_FILE/2007/11/1841.pdf).
- Leary. Global optimization on funneling landscapes. J. Global Optimization 18 (2000) 367-383.
  [DOI 10.1023/A:1026500301312](https://doi.org/10.1023/A:1026500301312) (gate; abstract + index).
- Cvijovic, Klinowski. Taboo search: an approach to the multiple minima problem.
  [Science 267 (1995) 664-666, DOI 10.1126/science.267.5198.664](https://doi.org/10.1126/science.267.5198.664).
- "The repulsion algorithm, a new multistart method for global optimization". Structural
  Optimization. [DOI 10.1007/BF01197028](https://doi.org/10.1007/BF01197028) (gate; abstract + index).
- Petrowski. A clearing procedure as a niching method for genetic algorithms. ICEC'96 798-803.
  [DOI 10.1109/ICEC.1996.542703](https://doi.org/10.1109/ICEC.1996.542703).
- Bratley, Fox. Algorithm 659: Implementing Sobol's quasirandom sequence generator. ACM TOMS 14(1)
  (1988). [DOI 10.1145/42288.214372](https://doi.org/10.1145/42288.214372).
- Bloom. Space/time trade-offs in hash coding with allowable errors. CACM 13(7) (1970) 422-426.
  [DOI 10.1145/362686.362692](https://doi.org/10.1145/362686.362692).
- McKay, Piperno. Practical graph isomorphism, II. J. Symbolic Computation 60 (2014) 94-112.
  [DOI 10.1016/j.jsc.2013.09.003](https://doi.org/10.1016/j.jsc.2013.09.003).
- Weisfeiler-Leman test tutorial. [NSF PAR 10299993](https://par.nsf.gov/servlets/purl/10299993).
- Donev et al. A linear programming algorithm to test for jamming in hard-sphere packings. J.
  Comput. Phys. 197 (2004) 139-166. [DOI 10.1016/j.jcp.2003.11.022](https://doi.org/10.1016/j.jcp.2003.11.022);
  [open PDF](https://math.nyu.edu/inmemoriam/donev/Packing/Jamming_LP.pdf).

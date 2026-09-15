# Phase 2 CSQV: a solver-improvement catalog for beating a hard record

Research date: 2026-09-13. This note is ADDITIVE to the ranked levers note
([docs/research/phase-2-csqv-polish-and-solver-levers.md](phase-2-csqv-polish-and-solver-levers.md)).
That note ranked six STRATEGIC levers and put KKT-Newton polish fifth. This note goes DEEPER and
solver-specific. It enumerates every concrete way to raise our CSQV solver quality, with a primary
citation, an expected ceiling against an N=90-class target, an implementation shape inside our
fixed-primitive architecture, a throughput effect, and a contract or ADR risk. It ends with a
ranked build shortlist.

CSQV means variable-radius circles in the unit square. The objective maximizes the sum of the
radii. The circles must not overlap and must stay inside the square. The band is N = 28 to 114.
The hard case is an N=90-class target: genuinely hard, not proven optimal, where our current solver
re-derives a strong incumbent basin but never beats it.

Method. This note uses primary sources only: the sibling Discovery Loop paper (full text), the
SciPy, IPOPT, NLopt, and HiGHS solver docs, the Grosso and Addis packing papers (both PDFs read in
full locally with pymupdf), the Lubachevsky and Graham billiard papers, and our own code
(`problems/csqv/problem.py`, `problems/csqv/record.py`). Each quantified claim carries an inline
URL. The prose uses short active sentences and no em dash.

---

## Headline verdict: the top two levers are spent, and no LOCAL-SOLVER lever alone clears the bar

Context from the coordinator. The two highest levers are now exhausted. Target selection is spent:
calibration showed no anchored soft target is beatable. Global diversity is spent: population basin
hopping, dissimilar pools, and large-jump escape are built and tried, and they re-derive strong
basins but never beat them ([record.py `RECORD_SEED_SOURCE`](../../problems/csqv/record.py)). The
only unspent levers are the local-solver ones.

Answer, stated plainly. NO purely local-solver lever, on its own, has enough ceiling to beat a hard
non-proven N=90-class record. The reason is structural. A local optimizer descends from its start
to the floor of the basin that contains that start. It cannot reach a basin the global walk never
visits. So a stronger, more exact, or faster LOCAL solve can only (1) reach the current basin floor
more precisely, which does not beat a record the search already re-derives, or (2) reach that floor
FASTER, which raises the ceiling only INDIRECTLY by feeding the already-built diversity walk more
and truer basin samples. Effect (1) has ZERO ceiling for this target. Effect (2) has a real but
SMALL and unquantified ceiling.

This matches the strongest primary evidence. In the Discovery Loop all ten records were "First
Improved" at "Iter 0 (seed)", the plain "Penalty L-BFGS-B + LP radii" solver; the KKT-Newton exact
polish arrived at iter 8, inside the diminishing-returns phase where cost per unit improvement was
$1,138 versus $8.70 early, a 130x jump, and no ablation credits it with any record
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). The classical record method
contains NO KKT, Newton, or active-set exact polish anywhere in its text (verified by full read of
the PDF); it reaches putative optima with a strong NLP local solver inside a good basin-hopping walk
([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

One caveat sharpens effect (2). Our polish is TRUNCATED at N=90 (Section 6.1): the iteration cap is
`max(30, min(200, int(6000/n)))` = 66 SLSQP iterations at n=90, below even SciPy's SLSQP default of
100 ([problem.py](../../problems/csqv/problem.py);
[SciPy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)). So at N=90
the search may compare partially-polished points, not true basin floors. A faster second-order
polish would both converge each basin fully AND fit more restarts in the 150s budget (the search is
polish-bound at 0.12 to 0.19 restarts/s). That is the ONE local-solver family worth building, and
its value is throughput, not precision.

One move genuinely CHANGES basins and may be UNSPENT: defect migration (Section 3.3). It removes the
weakest circles and reinserts them into holes, changing the COMBINATORIAL contact structure, not
just the continuous center positions. Our built diversity perturbs centers only; it never changes
the contact graph combinatorially ([record.py](../../problems/csqv/record.py)). So the classical and
Discovery Loop defect-migration move is arguably not yet tried. It is a shell move, not a
local-solver swap, but it is the honest answer to "what can break a basin". Its ceiling is still
modest: the classical authors themselves "could not understand which part of our methods lead us to
win" ([Addis, Locatelli, Schoen 1343.pdf](https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf)).

Honest bottom line on the ceiling. Build order: a faster second-order polish (throughput enabler)
plus a defect-migration move (the only basin-changing move left), then a cheap calibration. If a
3x to 5x throughput gain and a working defect move do not surface a beating basin at N=90, then no
remaining lever clears the bar. In that case declare the proven, ungameable, record-re-deriving
machinery the DELIVERABLE and stop building. Do not chase a hard non-proven N=90 record with a
precision polish; its own precedent paper gives it a $1,138-per-unit price and no record credit.

---

## The basin test: how to judge every lever below

Every lever falls into exactly one of three effects. Judge its ceiling by which one it is.

- PRECISION inside a basin. It sharpens the floor of the basin the search already found. Ceiling
  against a re-derived hard record: ZERO. It cannot open a new basin. Examples: KKT-Newton exact
  polish, a tighter `ftol`, exact-Hessian convergence used only to add digits.
- THROUGHPUT and convergence reliability. It reaches the same basin floor faster or more reliably.
  Ceiling: INDIRECT and small. It raises the reachable ceiling only by giving the already-built
  diversity walk more and truer samples per second. Examples: sparse Jacobian, analytic Hessian,
  warm starts, a cheaper LP.
- BASIN CHANGE. It moves the search to a different basin the continuous walk does not reach.
  Ceiling: the only family that can beat a re-derived record, but still modest here. Examples:
  defect migration (combinatorial), a genuinely different construction that jams into a new basin,
  a barrier-crossing acceptance rule.

The current plateau is a basin problem, so only the third effect can beat the record directly, and
the second can only help it. The first cannot help at all.

---

## 1. Local optimizers for hard packing

The CSQV local solve has 3N variables (x, y, r per circle), a LINEAR objective (sum of r), 4N linear
wall constraints, and O(N^2) nonlinear pairwise non-overlap constraints. The Jacobian and constraint
count grow like N^2 ([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091)). Because the objective is
linear and the containment constraints are linear, ALL nonlinearity sits in the pairwise distance
constraints, and each has a simple sparse 4x4 Hessian block in (x_i, y_i, x_j, y_j). So the Lagrangian
Hessian is exact, analytic, and sparse. This fact drives the comparison below.

The standard record-chasing local solver is a strong SQP or interior-point NLP, not a light method.
Grosso et al. state the problem has "objective and constraint functions continuously differentiable
infinitely many times. Therefore, any local search method for this kind of problems can be employed.
However, our past experience ... suggests that SNOPT [9] is particularly well suited for these
problems" ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)). SNOPT
is a large-scale SPARSE SQP (Gill, Murray, Saunders, SIAM J. Optim. 2002). The classical record
engine is therefore sparse and second-order-aware, not dense.

The candidate engines, by class:

- SLSQP (our current polish). Sequential Least Squares Programming, the Kraft 1988 dense SQP; it
  handles bounds plus equality and inequality constraints and returns KKT multipliers; default
  `maxiter` is 100; it does NOT accept a Hessian
  ([SciPy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html);
  [SciPy minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)).
  It is DENSE, so it scales badly as the O(N^2) constraints grow. This is the throughput bottleneck
  at N=90.
- L-BFGS-B penalty (the Discovery Loop seed engine, and the engine that set all 10 records). Limited
  memory BFGS for BOUND-constrained problems only; it enforces non-overlap only through a penalty,
  so it needs the LP-radii step to price feasibility (Byrd, Lu, Nocedal 1995; Zhu, Byrd, Nocedal
  1997) ([SciPy L-BFGS-B](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html)).
  Cheap and scalable. Its penalty landscape has DIFFERENT basins of attraction than SLSQP, so an A/B
  is worth running: it was the ACTUAL record engine, and it may reach different floors than our SLSQP.
- trust-constr (in SciPy, no new dependency). A trust-region method that "switches between" a
  Byrd-Omojokun equality-constrained SQP and a "trust-region interior point method"; it is "the most
  appropriate for large-scale problems"; it accepts a SPARSE Jacobian (`sparse_jacobian`) and a
  Hessian in every form (exact callable, sparse, LinearOperator, or BFGS/SR1); default `maxiter`
  is 1000 (Byrd, Hribar, Nocedal 1999; Lalee, Nocedal, Plantenga 1998)
  ([SciPy trust-constr](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html);
  [SciPy NonlinearConstraint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html)).
  This is the single most attractive local-solver build: it is in-SciPy, sparse, and it can take our
  exact sparse Lagrangian Hessian for true second-order (near-quadratic) convergence.
- IPOPT (interior-point NLP, C++ dependency). "An interior point line search filter method that aims
  to find a local solution" for twice-differentiable functions; it solves "sparse, symmetric,
  indefinite linear systems" via a third-party sparse solver (MA27, MA57, MA97, MUMPS); it offers
  `hessian_approximation` = exact or limited-memory, and `warm_start_init_point` to start "from a
  previous optimization of a related problem"
  ([IPOPT overview](https://coin-or.github.io/Ipopt/);
  [IPOPT options](https://coin-or.github.io/Ipopt/OPTIONS.html); Wächter and Biegler,
  Math. Programming 106(1):25-57, 2006). Stronger and faster than trust-constr at higher N, but a
  heavy dependency (C++ plus HSL or MUMPS).
- NLopt (C library, lighter than IPOPT). "Only MMA and SLSQP support arbitrary nonlinear inequality
  constraints, and only SLSQP supports nonlinear equality constraints; the rest support
  bound-constrained or unconstrained problems only". NLopt recommends CCSAQ (`LD_CCSAQ`) over MMA,
  and it supports preconditioning "including a user-supplied Hessian approximation in the local
  model"; AUGLAG wraps any bound solver to add general constraints
  ([NLopt algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/)). CCSAQ is a
  gradient-based method-of-moving-asymptotes variant that natively handles the non-overlap
  inequalities and can be more robust than dense SLSQP.

(a) What it is: a menu of NLP engines for the fixed polish primitive, with primary docs above.
(b) Ceiling at N=90: THROUGHPUT and convergence reliability, not basin change. A second-order sparse
engine reaches each basin floor faster and more reliably, which lets the diversity walk sample more
basins per 150s. Real but indirect and modest. The L-BFGS-B penalty A/B has a MILD basin-diversity
side effect (a different landscape), so it is the one engine swap with a small direct-ceiling story.
(c) Implementation shape: swap or add a FIXED primitive alongside `polish`, kept out of the genome,
deterministic given centers, exactly like the current `polish` and `radii_lp`
([problem.py](../../problems/csqv/problem.py)).
(d) Throughput: trust-constr with a sparse Jacobian and exact Hessian should beat dense SLSQP at
N=90 both per iteration and in iterations to converge (Section 6). IPOPT is faster still at higher N.
(e) Contract or ADR risk: trust-constr adds NO dependency; IPOPT adds a heavy C++ plus HSL/MUMPS
dependency; NLopt adds a lighter C dependency. Benchmark trust-constr first before any C++ dependency.

---

## 2. Exact contact-graph / KKT-Newton polish: precision, not a basin lever

(a) What it is. Identify the active contact set (which walls and which pairwise constraints are tight
at the local optimum), then take Newton steps on the KKT stationarity system of that FIXED active
set. On a correct active set this gives quadratic or superlinear convergence to the exact tangency
configuration. This is the Discovery Loop iter-8 innovation, "replacing the general-purpose SLSQP
solver with a specialized Newton solver for the contact-graph KKT system"
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

(b) Ceiling at N=90: PRECISION inside a basin, so effectively ZERO against a re-derived record. It
fixes the active set from the current configuration, so the Newton step stays inside the current
basin and converges to that basin's exact floor. It squeezes the last digits of the feasibility
shrink; it cannot open a new basin. The primary data agrees: iter 8 sits in the diminishing phase at
$1,138 per unit improvement, its Table 2 total is flat (iter 7 to iter 8 both 59.97), and no
ablation credits it ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)). The classical record
paper contains no such method at all ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

Its only non-precision value is THROUGHPUT: a Newton step that converges in a handful of iterations
can be much cheaper than 66 truncated SLSQP iterations at N=90 (Section 6.1), which frees budget for
more restarts. So build it, if at all, as a fast FINAL step, not as a record lever.

(c) Implementation shape: a FIXED primitive (say `polish_kkt`) or a final stage inside `polish`, out
of the genome, deterministic given centers and the detected active set.
(d) Throughput: potentially high if it replaces many SLSQP iterations with few Newton steps; that is
its real case.
(e) Contract or ADR risk: active-set FRAGILITY. A near-tangent pair that is not truly active, or a
degenerate or over-determined contact graph, gives a singular or ill-conditioned KKT system and a
wrong or unstable step. It needs a guard: detect a singular system or a step that raises
infeasibility, and fall back to SLSQP or trust-constr. Contact detection needs a tolerance the paper
never specifies ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)).

---

## 3. Global and perturbation refinements beyond what we built

Our built shell does continuous center perturbation, an adaptive per-member jump size, a rare large
jump, and a dissimilar pool ([record.py `RECORD_SEED_SOURCE`](../../problems/csqv/record.py)). The
levers below go beyond that.

### 3.1 Monotonic basin hopping tuning (the delta rule)

The perturbation step size delta is decisive. Grosso et al.: choose it "in such a way that the
structure of the current local minimizer is not completely disrupted ... the method should move
between different but 'close' local minimizers". "If delta is too small, the starting point will be
very likely in the basin of attraction of the current local minimizer"; "if delta is too large, the
method becomes basically equivalent to a Multistart method". Their experiments make delta = 0.8 (in
a scaled container) "the most robust choice among the tested ones"
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)). Our seed uses a
Gaussian sigma of 0.02 to 0.2 on unit coordinates; the delta rule says the sweet spot is problem
scaled, so a scale sweep is cheap and in-genome. Ceiling: this is diversity tuning, part of the spent
global lever; low marginal ceiling. Effect: BASIN CHANGE tuning, but of an already-tried mechanism.

### 3.2 Funnel walks and acceptance rules (simulated annealing, threshold accepting)

The packing landscape is a FUNNEL: local optima "follow some sort of ordered structure ... known
under the name of 'funnel landscape'" ([Addis 1343.pdf](https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf)).
Monotonic acceptance descends the funnel but cannot cross a barrier UP. Simulated annealing or
threshold accepting can cross a barrier to a neighboring funnel branch, which is a BASIN CHANGE
mechanism our monotonic shell lacks. Ceiling: modest; it widens the reachable set beyond monotone
descent, but the funnel structure means the incumbent branch is usually the deepest. Implementation:
in-genome acceptance rule; no primitive change; contract-safe.

### 3.3 Defect migration (remove weak circles, re-equilibrate, reinsert): the one basin-changing move that may be UNSPENT

(a) What it is. Remove the weakest (smallest or most-constrained) circles, let the rest
re-equilibrate, then reinsert the removed circles into the holes that opened. This changes the
COMBINATORIAL contact structure, not just the continuous positions. Primary precedent is strong and
DOUBLE. Classical: Grosso et al. for unequal circles, "removing some circles, followed by insertion
of the missing circles, easily solves these instances", via "reduction of the search space during
the first phase where some circles are removed"
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)). Addis et al. built
a re-insertion-into-holes move that "gave an impressive improvement over our results, and many new
records could be found", and they observed that "improvements were consequences of the phenomena of
backtracking and survival"
([Addis 1343.pdf](https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf)). Discovery
Loop rediscovered it at iter 12: "removing weak circles, letting the packing re-equilibrate, then
reinserting into new holes ... a metaheuristic move that changes the combinatorial structure of the
packing" ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)).

(b) Ceiling at N=90: this is the ONLY remaining move with a true BASIN-CHANGE mechanism, and our
built diversity is continuous-only, so it may NOT be spent
([record.py](../../problems/csqv/record.py)). It is the honest best hope. Ceiling is still modest:
Addis et al. "could not understand which part of our methods lead us to win the contest", so no
primary source isolates its record contribution
([Addis 1343.pdf](https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf)).
(c) Implementation shape: IN-GENOME shell move using the fixed primitives (drop rows, call
`radii_lp` on the reduced set, `polish`, reinsert into the largest holes, `polish` again). No new
primitive, no verifier change.
(d) Throughput: each defect move costs one or two extra polish and LP calls; affordable given the
150s budget IF the polish is cheaper (Section 6).
(e) Contract or ADR risk: none to the contract; it is a pure genome move. It is exactly the kind of
combinatorial move the genome is meant to evolve.

### 3.4 Shrink-and-expand and billiard (Lubachevsky-Graham BLD) construction

(a) What it is. A billiard simulation: disks "move along straight lines, colliding with each other
and the region walls ... all the time maintaining a condition of no overlap", and "the disks are
uniformly allowed to gradually increase in size, until no significant growth can occur"
([Lubachevsky and Graham, DCG 18:179-194, 1997, arXiv:math/0406098](https://arxiv.org/abs/math/0406098);
event-driven engine in [Lubachevsky, J. Comput. Phys. 94(2):255-283, 1991,
arXiv:cond-mat/0503627](https://arxiv.org/abs/cond-mat/0503627)). Nurmela and Ostergard use a
repulsive-energy relaxation to reach dense equal-circle configurations (Discrete & Computational
Geometry 18(1):111-120, 1997, DOI 10.1007/PL00009306; primary text behind a publisher gate, so this
attribution is bibliographic only).

(b) Ceiling at N=90: as a CSQV SEED GENERATOR it jams into different configurations than the
grid-plus-LP seed, so it can seed a NEW basin (a basin-change mechanism at the construction stage).
But both methods target EQUAL radii, so they serve CSQV only as one initializer among many, and the
CSQV record configs are variable-radius and C1. Ceiling: low to modest, as a diversity seed.
(c) Implementation shape: an initializer inside the genome (produce centers), then `radii_lp` and
`polish` as usual. Optionally a fixed `bld_seed` primitive if a fast event-driven core is wanted.
(d) Throughput: a Python billiard sim is not cheap; use it to seed a few restarts, not every one.
(e) Contract or ADR risk: none if it only produces centers; the primitives and verifier are untouched.

---

## 4. Construction and initialization

(a) What it is. Better seeds than square-grid-plus-LP. Options from the precedent: hexagonal or
affine lattice templates "sized to hold exactly N circles" and an affine template bank (Discovery
Loop iters 4 and 7, mostly REJECTED), row-template cold starts (iter 14, rejected), golden-angle
spiral (ShinkaEvolve), greedy contact insertion, and the incumbent warm start we already use
([arXiv:2609.05093](https://arxiv.org/html/2609.05093);
[ShinkaEvolve arXiv:2509.19349](https://arxiv.org/html/2509.19349v1);
[record.py](../../problems/csqv/record.py)).

(b) Ceiling at N=90: a seed matters only if it seeds a DIFFERENT basin. That overlaps with the spent
global-diversity lever, and cold restarts at N=90 fail badly (49% gap), which says most fresh seeds
land in poor basins. Symmetry-aware seeds are a trap here: the record targets are C1, so a hard
symmetry constraint caps the sum below the record; use a symmetric layout only as ONE seed, then
relax to C1 in the polish
([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt);
[docs/research/phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md)).
Ceiling: low, except the incumbent warm start we already exploit.
(c) Implementation shape: in-genome initializers.
(d) Throughput: cheap.
(e) Contract or ADR risk: none.

---

## 5. The radii step: is the LP the bottleneck?

(a) What it is. `radii_lp` solves, for FIXED centers, the exact largest radii by a HiGHS linear
program: maximize sum of r subject to containment and r_i + r_j <= dist(i, j)
([problem.py](../../problems/csqv/problem.py);
[LP for radii, arXiv:1607.02184](https://arxiv.org/pdf/1607.02184)).

(b) Ceiling for QUALITY: ZERO. Given fixed centers the LP is EXACT and optimal, so no quality lives
here; all quality is in center placement. The radii step cannot beat a record. Its only lever is
THROUGHPUT.
(c)/(d) Throughput shape. The LP is re-solved from scratch on every restart and every polish seed,
and it builds a DENSE `A_ub` of shape (N(N-1)/2, N); at N=90 that is about 4005 by 90
([problem.py](../../problems/csqv/problem.py)). Two throughput wins: (1) build `A_ub` SPARSE (each
row has two nonzeros), which SciPy `linprog(method="highs")` accepts as a sparse matrix; (2) HOT
START the dual simplex from the previous basis, since between restarts only bounds and a few RHS
values change. HiGHS supports this natively: "Only the simplex solver can be hot-started", via
`setBasis` and `getBasis`, and it "reduce[s] the solution time by using data obtained on previous
runs" ([HiGHS further features](https://ergo-code.github.io/HiGHS/dev/guide/further/index.html);
[HiGHS overview, Huangfu and Hall 2018](https://ergo-code.github.io/HiGHS/dev/index.html)). BUT
`scipy.optimize.linprog` does NOT expose the HiGHS basis: its only warm-start parameter is `x0`,
"used only by the 'revised simplex' method", not by HiGHS
([SciPy linprog](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)).
To hot-start you must call the native `highspy` interface directly.
(e) Contract or ADR risk: sparse `A_ub` is a pure internal speedup, zero risk. A `highspy` dependency
is light and keeps the LP exact, so it does not touch the ungameable contract. Weigh it only if
profiling shows the LP is a real share of the per-restart cost; the polish is the larger cost.

---

## 6. Throughput and engineering: fit more restarts in 150s

The search is polish-bound at 0.12 to 0.19 restarts/s at N=90, so throughput directly multiplies the
diversity walk's basin samples. This is the section with the most concrete, code-level wins.

### 6.1 The polish is TRUNCATED at N=90 (a real code finding)

`_polish_scored` caps SLSQP at `eff_maxiter = max(30, min(maxiter, int(6000 / n)))`, which is 66
iterations at n=90 with `POLISH_MAXITER = 200` ([problem.py](../../problems/csqv/problem.py)). SciPy
SLSQP's own default is 100 ([SciPy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)).
At n=90 the solve has 270 variables and 4 * 90 + 4005 = 4365 constraints, so 66 dense SLSQP
iterations likely UNDER-converge. Consequence: the "basins" the search compares at N=90 may be
partially-polished points, not true floors, which BOTH leaks precision AND wastes the comparison.
This is a mixed precision-and-throughput symptom, and it is the strongest reason a better local
solve could matter at N=90 specifically. Fix: a sparse second-order solve converges in far fewer
iterations, so the cap stops binding.

### 6.2 The polish Jacobian is DENSE and rebuilt every iteration

`constraints_jac` builds a dense array of shape (4N + N(N-1)/2, 3N) on every SLSQP iteration; at n=90
that is about 4365 by 270, near 1.18 million entries, of which only about 26 thousand are nonzero
(each overlap row has 6 nonzeros, each wall row 2) ([problem.py](../../problems/csqv/problem.py)).
SLSQP cannot use a sparse Jacobian, but trust-constr and IPOPT can
([SciPy trust-constr](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html);
[SciPy NonlinearConstraint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html)).
A sparse Jacobian cuts the per-iteration linear algebra by roughly 45x at N=90.

### 6.3 Analytic sparse Hessian for true second-order steps

The objective is linear and containment is linear, so the Lagrangian Hessian is exactly the sum over
pairwise constraints of multiplier times a sparse 4x4 block in (x_i, y_i, x_j, y_j). It is analytic
and cheap to assemble. trust-constr accepts it as a callable returning a sparse matrix or
LinearOperator; IPOPT accepts `hessian_approximation = exact`
([SciPy trust-constr](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html);
[IPOPT options](https://coin-or.github.io/Ipopt/OPTIONS.html)). True second-order steps converge
near-quadratically, so each basin floor is reached in a handful of iterations, not 66-plus. This is
the throughput lever with the clearest ceiling, because it also removes the truncation of 6.1.

### 6.4 Warm-starting the local solve

MBH perturbs the incumbent, so each restart starts NEAR a known optimum. IPOPT can exploit this with
`warm_start_init_point` from "a previous optimization of a related problem"
([IPOPT options](https://coin-or.github.io/Ipopt/OPTIONS.html)). SLSQP and trust-constr have no
explicit warm-start, but the small perturbation already gives a good initial point, so the main win
is the LP hot start (Section 5) and the second-order convergence (6.3).

(a)-(e) for Section 6 as a whole. (a) engineering speedups with primary docs above. (b) ceiling:
THROUGHPUT, indirect, but the largest single-lever gain available, because at N=90 the polish is
both truncated and dense. (c) implementation: replace or augment the fixed `polish` primitive with a
sparse second-order trust-constr solve; make `radii_lp` sparse; optionally `highspy` hot start. All
FIXED, out of the genome. (d) throughput: plausibly 3x to 5x more restarts per 150s at N=90, and
truer floors; unverified until benchmarked. (e) contract or ADR risk: trust-constr and sparse `A_ub`
add no dependency and keep the primitives exact-or-deterministic; IPOPT and highspy add dependencies.
Benchmark the in-SciPy path first.

---

## 7. Second-order information

Covered in 6.3 and in Section 1. The one point to add: second-order info is a THROUGHPUT and
convergence lever, not a basin lever. A Newton or interior-point step from a given start lands in the
SAME basin as a first-order step from that start; it just gets there in fewer iterations and to a
truer floor. So exploiting the exact sparse Lagrangian Hessian (via trust-constr or IPOPT) raises the
ceiling only by enabling more and truer samples for the diversity walk. It cannot, by itself, beat a
re-derived record. The KKT-Newton contact-graph solve (Section 2) is the extreme second-order case
and shares this limit, plus active-set fragility.

---

## 8. What Discovery Loop's evolved solvers actually used

Separate what drove records from what was marginal, from the primary text
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

WHAT DROVE THE RECORDS: the SEED SOLVER and the TARGETS. Table 1 credits every one of the 10 records
(N = 101, 102, 103, 105, 106, 107, 108, 109, 111, 114) to "Iter 0 (seed)", the "Penalty L-BFGS-B + LP
radii" solver, a "simple multi-start penalty L-BFGS-B solver with LP-optimal radii". The gains were
2.43% to 5.39%. So a first-generation seed beat those best-knowns, which only happens when the
targets are soft.

WHAT WAS MARGINAL: the late additions. The productive phase was iters 0 to 5 at $8.70 per unit; the
diminishing phase was iters 6 to 14 at $1,138 per unit, a 130x jump. KKT-Newton (iter 8),
defect-migration (iter 12), and the island model (iter 5) all sit at or past this transition, and
Table 2's total is essentially flat across them (59.96 at iter 5 to 59.98 at iter 12). NO ablation
isolates any of them; early stopping after iter 9 "would have ... sacrificing only 0.006 on total sum
(0.01%)". The island model is "independent basin-hopping chains on separate CPU cores with periodic
migration of elite solutions". The independent verifier "shares no code with the solver, preventing
the LLM from gaming the evaluation", and it confirms containment and non-overlap "with zero
tolerance", recomputes the sum, and "applies a strict feasibility shrink". The run used a 120-second
per-target timeout, 6 parallel workers, an Intel Core i7-13700KF with 32 GB RAM, about 8 hours, and
$27.72 ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)).

Lesson for us. The precedent that broke records did it with a plain seed on soft targets, not with a
clever polish. Since our soft targets are spent, the precedent gives NO evidence that any late solver
add will beat a hard target. It weakly supports defect migration and the island model as accepted
building blocks, but with no measured record contribution.

---

## 9. Ranked build shortlist

Impact is the expected effect on BEATING a hard non-proven N=90-class record, given the top two
levers are spent. Effect is the basin-test class from Section 2.

| Rank | Improvement | Effect class | Ceiling | Layer | Contract-safe |
|------|-------------|--------------|---------|-------|---------------|
| 1 | Sparse second-order polish (trust-constr, exact sparse Hessian; fixes the N=90 truncation) | Throughput + convergence | Indirect, modest, but largest single lever | Fixed primitive | Yes |
| 2 | Defect-migration move (remove weak circles, reinsert into holes) | Basin change | The only basin-changing lever left; modest | Genome shell | Yes |
| 3 | L-BFGS-B penalty engine A/B (the actual Discovery Loop record engine) | Throughput + mild basin diversity | Low to modest; cheap to test | Fixed primitive | Yes |
| 4 | Barrier-crossing acceptance (simulated annealing or threshold accepting) | Basin change | Low; the funnel favors the incumbent branch | Genome shell | Yes |
| 5 | Sparse LP plus highspy hot start | Throughput | Indirect, small (LP is not the main cost) | Fixed primitive | Yes |
| 6 | KKT-Newton exact polish (as a FAST final step only, guarded) | Precision (+ throughput if few steps) | Zero for beating; some throughput | Fixed primitive | Yes |

Reasoning for the top build:

1. Build the sparse second-order polish FIRST. It is the largest single lever, it is in-SciPy with no
   new dependency, and it directly fixes the N=90 truncation (66 SLSQP iterations, dense Jacobian)
   that may be poisoning basin comparison. Its ceiling for beating the record is indirect and
   modest, but every other lever depends on the search actually reaching true basin floors and
   getting more samples per 150s
   ([SciPy trust-constr](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html);
   [problem.py](../../problems/csqv/problem.py)).
2. Build the defect-migration move SECOND. It is the ONLY remaining move with a true basin-change
   mechanism, our built diversity is continuous-only, and it has double primary precedent
   ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf);
   [Addis 1343.pdf](https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf);
   [arXiv:2609.05093](https://arxiv.org/html/2609.05093)).
3. Then CALIBRATE, do not keep building. Run the faster polish plus defect migration at N=90 and ask
   one question: does a better basin surface within the budget? If a 3x to 5x throughput gain and a
   working defect move do not beat the incumbent, no remaining lever clears the bar. Declare the
   proven, ungameable, record-re-deriving machinery the DELIVERABLE and stop.

THE SINGLE HIGHEST-IMPACT IMPROVEMENT: the sparse second-order polish (rank 1). It is the only lever
that both fixes a concrete N=90 defect (the truncated, dense SLSQP solve) and multiplies the sample
count of the already-built diversity walk, and it costs no new dependency. Its ceiling for beating
the record is honest and modest; it is the best available, not a guaranteed win.

---

## Sources

- Discovery Loop, full text (Table 1 all records at "Iter 0 (seed)"; seed "Penalty L-BFGS-B + LP
  radii"; iter 5 island model; iter 8 KKT-Newton, cost 3.22, total flat 59.97; iter 12 defect
  migration "removing weak circles ... reinserting into new holes"; productive $8.70 vs diminishing
  $1,138 per unit, 130x; verifier zero-tolerance plus feasibility shrink; 120 s, 6 workers,
  i7-13700KF, ~8 h, $27.72): https://arxiv.org/abs/2609.05093 and https://arxiv.org/html/2609.05093
- Grosso, Jamali, Locatelli, Schoen, "Solving the problem of packing equal and unequal circles in a
  circular container" (PDF read in full locally: SNOPT "particularly well suited"; "any local search
  method ... can be employed"; "not completely disrupted"; delta too small or too large; MBH
  "dramatically more efficient than Multistart"; Multistart "clearly inferior"; PBH dissimilarity
  cut dcut; delta = 0.8 most robust; "in 19 cases we could obtain an improvement" for equal circles
  vs Specht; defect move "removing some circles, followed by insertion of the missing circles"; NO
  KKT/Newton/active-set method present):
  https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf and
  https://optimization-online.org/2008/06/1999/
- Addis, Locatelli, Schoen, "Efficiently packing unequal disks in a circle" (PDF read in full
  locally: funnel-landscape conjecture; public contest "155 teams who submitted a total of 27 490
  tentative solutions"; re-insertion into holes "gave an impressive improvement ... many new records
  could be found"; "improvements were consequences of the phenomena of backtracking and survival";
  "it is quite difficult for us to understand which part of our methods lead us to win the contest"):
  https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf and
  https://optimization-online.org/2006/03/1343/ (ORL 36(1):37-42, 2008)
- SciPy SLSQP (Kraft 1988; default maxiter 100; no Hessian):
  https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html
- SciPy trust-constr (Byrd-Omojokun SQP plus trust-region interior point; sparse Jacobian; exact,
  sparse, LinearOperator, BFGS, or SR1 Hessian; maxiter 1000; Byrd-Hribar-Nocedal 1999,
  Lalee-Nocedal-Plantenga 1998):
  https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html
- SciPy L-BFGS-B (Byrd, Lu, Nocedal 1995; Zhu, Byrd, Nocedal 1997; bounds only):
  https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html
- SciPy NonlinearConstraint and LinearConstraint (sparse Jacobian; Hessian(x, v) with multipliers):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html and
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.LinearConstraint.html
- SciPy minimize (constraints only for COBYLA, COBYQA, SLSQP, trust-constr; hess only for the
  Newton and trust methods; trust-constr "most appropriate for large-scale problems"):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
- IPOPT (interior-point line-search filter; twice-differentiable; sparse symmetric indefinite
  systems via MA27, MA57, MA97, MUMPS; hessian_approximation exact or limited-memory;
  warm_start_init_point; Wächter and Biegler, Math. Programming 106(1):25-57, 2006):
  https://coin-or.github.io/Ipopt/ and https://coin-or.github.io/Ipopt/OPTIONS.html
- NLopt algorithms ("only MMA and SLSQP support arbitrary nonlinear inequality constraints, and only
  SLSQP supports nonlinear equality constraints"; CCSAQ recommended, supports a user Hessian
  approximation; AUGLAG wraps any bound solver):
  https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
- HiGHS (revised dual simplex the default; native hot start via setBasis/getBasis, simplex only;
  Huangfu and Hall, Math. Prog. Comp. 10(1):119-142, 2018):
  https://ergo-code.github.io/HiGHS/dev/index.html and
  https://ergo-code.github.io/HiGHS/dev/guide/further/index.html
- SciPy linprog (no HiGHS basis warm start; x0 only for revised simplex):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html
- Lubachevsky and Graham, "Curved Hexagonal Packings of Equal Disks in a Circle" (DCG 18:179-194,
  1997; the "billiards" grow-and-collide-until-jammed generator): https://arxiv.org/abs/math/0406098
- Lubachevsky, "How to Simulate Billiards and Similar Systems" (J. Comput. Phys. 94(2):255-283,
  1991; the event-driven collision engine): https://arxiv.org/abs/cond-mat/0503627
- Nurmela and Ostergard, "Packing up to 50 Equal Circles in a Square" (DCG 18(1):111-120, 1997;
  repulsive-energy relaxation; primary text behind a publisher gate, bibliographic only):
  https://doi.org/10.1007/PL00009306
- LP for radii given fixed centers, arXiv:1607.02184: https://arxiv.org/pdf/1607.02184
- Nonconvexity and O(N^2) constraint scaling, arXiv:2404.03091: https://arxiv.org/pdf/2404.03091
- ShinkaEvolve (golden-angle spiral init, hybrid gradient plus annealing), arXiv:2509.19349:
  https://arxiv.org/html/2509.19349v1
- Our code: problems/csqv/problem.py (polish eff_maxiter truncation, dense constraints_jac,
  radii_lp dense A_ub) and problems/csqv/record.py (RECORD_SEED_SOURCE continuous-perturbation shell)
- Companion notes: docs/research/phase-2-csqv-polish-and-solver-levers.md,
  docs/research/phase-2-csqv-record-methods-and-target-selection.md,
  docs/research/phase-2-csqv-method-and-precedent.md

## Open questions

- CORRECTION to a companion figure. The companion note states the Addis, Locatelli, Schoen disk
  paper "improve[d] over previously known putative optima in the range n <= 130 in as many as 32
  instances". A full local read of BOTH primary PDFs does NOT find that figure. The verified figures
  are Grosso's "in 19 cases we could obtain an improvement" for equal circles versus Specht
  ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)) and the Addis
  disk paper's contest win over 155 teams and 27,490 solutions, with no "32 / n<=130" count present
  ([1343.pdf](https://optimization-online.org/wp-content/uploads/2006/03/1343.pdf)). The "32 / n<=130"
  claim should be treated as UNVERIFIED and likely belongs to a different Addis-Locatelli-Schoen
  paper (possibly the JOGO disk-in-square work); track it down before reusing it.
- The N=90 throughput and truer-floor gain of a sparse second-order polish versus the current
  truncated dense SLSQP is unmeasured. Benchmark trust-constr with an exact sparse Hessian against
  SLSQP at N=30 and N=90: restarts per second, and best sum over a fixed 150 s.
- Whether defect migration surfaces a better basin at N=90 is untested; our built diversity never
  tried a combinatorial contact-graph move. A local A/B (continuous walk versus walk plus defect
  migration, same budget) would size it.
- Whether the N=90 plateau is partly an under-convergence artifact (Section 6.1) or purely basin
  scarcity is unknown. Re-run the plateau with an untruncated, converged polish before concluding
  the record is out of reach.
- IPOPT and highspy add dependencies. Confirm the in-SciPy trust-constr path is insufficient before
  taking either dependency.
- Nurmela-Ostergard primary text is behind a publisher gate; the repulsive-energy description here
  is bibliographic. Fetch the PDF through an authenticated source if the billiard or energy seed is
  built.

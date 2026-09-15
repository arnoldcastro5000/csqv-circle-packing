# Phase 2 CSQV: polish and solver levers, ranked

Research date: 2026-09-13. This note answers one question. Is a KKT-Newton exact polish on the
contact graph the BEST lever to raise our CSQV sum of radii, or are there stronger levers? It
ranks the candidate levers by expected impact on the achievable sum, with primary evidence, and
it states where KKT-Newton actually lands. Read the method and target notes first for the
pipeline, the soft-target procedure, and the MBH background
([docs/research/phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md);
[docs/research/phase-2-csqv-resume-and-warm-restart.md](phase-2-csqv-resume-and-warm-restart.md)).

CSQV means variable-radius circles in a unit square. The objective maximizes the sum of the
radii. The circles must not overlap and must stay inside the square. The relevant band is
N = 28 to 114.

This note uses primary sources: the sibling Discovery Loop paper (arXiv:2609.05093), the
Grosso, Jamali, Locatelli, Schoen circle-in-circle paper, the Addis, Locatelli, Schoen disk
paper, the scipy SLSQP reference, the IPOPT
project docs, and the Packomania data. Each quantified claim carries an inline URL. The prose
uses short active sentences and no em dash.

---

## Summary: the ranked verdict

The plateau at N=60 is a BASIN-DIVERSITY problem, not a polish-PRECISION problem. Every
candidate re-derives the same strong local optimum. A more exact polish inside that same basin
cannot beat the incumbent. So the highest-impact levers change WHICH basin the search reaches,
or change the TARGET so the basin the search reaches is already better than the record. A more
exact polish only sharpens the last digits of a basin the search already found.

Ranked by expected impact on the achievable sum of radii, strongest first:

1. TARGET SELECTION. Attack a soft N (an old, not-proven Cantrell holdout the 2026 wave
   skipped) instead of a hard strong optimum like N=60. Highest impact. The single decisive
   lever ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).
2. GLOBAL / PERTURBATION STRATEGY. Population Basin Hopping with a dissimilarity cut, adaptive
   and occasional large jumps, and explicit diversity injection. This is the lever that
   directly attacks the observed N=60 plateau. Monotonic Basin Hopping is "dramatically more
   efficient than Multistart" on circle packing
   ([Grosso, Jamali, Locatelli, Schoen 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)),
   and the population variant improved putative optima "in the range n <= 130 in as many as 32
   instances" on unequal disks
   ([Addis, Locatelli, Schoen, ORL 2008](https://optimization-online.org/2006/03/1343/)).
3. BETTER LOCAL NLP SOLVER as the MBH engine (a strong SQP or interior-point NLP, SNOPT or
   IPOPT class). Real but secondary. Grosso et al. report the local method is effective, not
   uniquely decisive: "any local search method ... can be employed", and past experience only
   "suggests that SNOPT ... is particularly well suited"
   ([Grosso, Jamali, Locatelli, Schoen 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
4. MORE / BETTER LLM SEARCH. Moderate. The LLM only shapes the search shell. At a plateau the
   NUMERIC layer binds, not the LLM layer ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).
5. KKT-NEWTON EXACT POLISH. Low to moderate. A precision refinement inside an already-found
   basin. It sharpens the final feasibility squeeze; it does not change the basin reached, so
   it cannot break a plateau. Cheap, safe, but not the bottleneck
   ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).
6. MORE COMPUTE. Lowest at a plateau. More restarts of a search that re-collapses into one
   basin returns the same optimum. Compute helps only after levers 1, 2 raise the ceiling
   ([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

KKT-Newton polish ranks FIFTH of six. It is not the best lever. The Discovery Loop evidence
argues strongly against it as the deciding factor (Section 1). It is a good, safe final-squeeze
primitive, but it targets the wrong bottleneck for our N=60 plateau.

---

## 1. The Discovery Loop evidence: the polish did NOT drive the records

The sibling paper is the direct precedent and the strongest evidence on this exact question
([arXiv:2609.05093 abstract](https://arxiv.org/abs/2609.05093);
[full text](https://arxiv.org/html/2609.05093)).

The decisive fact. In Table 1 every one of the 10 records was FIRST improved at "Iter 0
(seed)". The seed solver is "Penalty L-BFGS-B + LP radii"
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). The seed did the record
work. The KKT-Newton exact polish is a LATER discovery: "KKT-Newton exact polish (iter 8):
replacing the general-purpose SLSQP solver with a specialized Newton solver for the
contact-graph KKT system" ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).
Iteration 8 is well past iteration 0. The records were already set before KKT-Newton existed in
the run.

No ablation supports KKT-Newton. The paper lists the algorithmic innovations but gives NO
controlled comparison of solvers. It runs no experiment that isolates the effect of KKT-Newton,
the island model, or defect migration on the record sums
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). The cost analysis reports
that the late iterations (roughly iter 6 to 14) show diminishing returns
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). KKT-Newton (iter 8) sits
inside that diminishing-returns zone. So the primary evidence says the SEED SOLVER plus the
choice of TARGETS carried the records, and the exact polish was a marginal late add with no
proof it moved a record.

Read critically: the claim "KKT-Newton polish is the best lever" is unsupported by its own
source paper. The paper's own data ranks it below the seed and the targets.

---

## 2. Lever 1, target selection: the dominating lever

The evidence in Section 1 is also the evidence for target selection. Records fell at the seed
solver, so the binding variable was not solver quality but WHICH N the loop attacked. The
Discovery Loop targeted N in the 101 to 114 band and adjacent points, an area where soft
best-knowns existed ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). A
first-generation seed solver beat those best-knowns. That only happens when the target is soft.

Our own target note already builds the soft-target procedure from Packomania metadata: select N
still held by David W. Cantrell (reference [1], 2011/12), drop the proven-optimal bold set, keep
C1 targets the 2026 wave skipped, and confirm the neighbors changed hands but the target did
not ([author.txt](https://www.packomania.com/csqv/txt/author.txt);
[csqv.html](https://www.packomania.com/csqv/csqv.html);
[docs/research/phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md)).
The candidate soft N in our band are 28, 29, 30, 31, 37, 38, 41, 43, 45, 90, among others.

The N=60 lesson. N=60 shows a plateau because it is likely a strong or near-optimal local
optimum: every candidate re-derives it. That is a signal to change target, not to grind harder
on the same one. A soft N gives more headroom per unit of solver effort than a hard N gives per
unit of polish precision.

Impact: HIGHEST. Cost in our architecture: NONE structural. Target selection is a run input,
not a genome or primitive change. It does not touch the fixed-primitive or verifier contract.
Risk: the frontier moves weekly, so re-pull `author.txt` at attempt time
([author.txt](https://www.packomania.com/csqv/txt/author.txt)).

---

## 3. Lever 2, global and perturbation strategy: the plateau lever

This lever directly attacks the observed failure mode. The N=60 plateau is basin trapping: warm
restarts re-collapse into one basin. The cure is diversity in the SEARCH, which is exactly what
Monotonic Basin Hopping (MBH) and Population Basin Hopping (PBH) provide.

Two primary papers give the mechanics, and this note quotes both directly after reading the raw
PDF and the abstract. The mechanics paper is Grosso, Jamali, Locatelli, Schoen, "Solving the
problem of packing equal and unequal circles in a circular container"
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)). The
soft-optima result is Addis, Locatelli, Schoen, "Efficiently packing unequal disks in a circle"
([Optimization Online 1343](https://optimization-online.org/2006/03/1343/)). Note both use a
CIRCLE container, not a square, but the MBH/PBH method is container-general.

MBH mechanics, verbatim from the raw PDF. The perturbation "A good rule is to choose it in such
a way that the structure of the current local minimizer is not completely disrupted. The basic
idea is that the method should move between different but 'close' local minimizers." The move is
"a uniform random perturbation of each coordinate of each circle center within some interval
[-delta, delta]." The step size is decisive: "If delta is too small, the starting point will be
very likely in the basin of attraction of the current local minimizer", which is exactly our
plateau; if delta is too large "the method becomes basically equivalent to a Multistart method"
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)). PBH replaces
the single path with a population: it keeps N minimizers and, for each new candidate, competes
it with the most similar population member unless the dissimilarity exceeds a cut `dcut`, in
which case it competes with the worst member. This keeps members in different basins
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

The multistart comparison, verbatim, is the load-bearing evidence for this lever. "the rapid
increase in the number of local minimizers suggests that Multistart can not be an efficient
method for this problem", while MBH is "dramatically more efficient than Multistart." In the
experiments "the results obtained with Multistart are usually quite poor, and are clearly
inferior with respect to those obtained with MBH." Even when the authors gave Multistart twice
the local searches of MBH, Multistart reached the best-known "only in few cases, and only in a
single case (namely n = 61) the best known solution is reached quite regularly", and the MBH
cost per local search is "almost four times lower"
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)). On putative
optima, the same paper "in 19 cases we could obtain an improvement" over Specht's best-known for
equal circles in a circle ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

The landscape argument and the strongest soft-optima number come from the DISK paper, not the
1999.pdf. The 1999.pdf body does NOT state a funneling-landscape conjecture; it cites the funnel
idea only through a reference (Leary 2000). The disk paper states the conjecture and the number:
its stochastic search is "based on the conjecture that the problem possesses a so-called
funneling landscape", and it "could improve over previously known putative optima in the range
n <= 130 in as many as 32 instances", after a public contest with "155 teams ... 27,490
tentative solutions" ([Addis, Locatelli, Schoen, ORL 2008](https://optimization-online.org/2006/03/1343/)).
The parent basin-hopping work (Wales and Doye) established the funnel picture for molecular
clusters ([Wales and Doye 1997, DOI 10.1021/jp970984n](https://pubs.acs.org/doi/10.1021/jp970984n)).
So the lever that matters is the QUALITY OF THE WALK BETWEEN BASINS, not the sharpness of the
polish inside one basin. Adaptive and occasional large jumps, a diverse PBH population, and
island reset and migration are the concrete diversity mechanisms
([docs/research/phase-2-csqv-resume-and-warm-restart.md](phase-2-csqv-resume-and-warm-restart.md)).

Impact: HIGH. It is the lever that fixes the reported N=60 plateau. Cost in our architecture:
the search shell is exactly what the genome already evolves, so PBH dissimilarity, adaptive
jumps, and diversity injection are IN-GENOME changes. No new fixed primitive, no verifier
change. It does not break the ungameable contract; the verifier still recomputes the sum from
coordinates. Risk: none to the contract. The only risk is tuning (delta, `dcut`, jump
schedule), which the LLM can search.

---

## 4. Lever 3, better local NLP solver: real but secondary

The local solver is the engine each MBH step calls to descend a basin. A stronger engine finds
a deeper basin per step and needs fewer restarts.

The standard record-chasing local solver in the packing literature is a strong sparse NLP
solver, not a general-purpose light method. Grosso et al. state the problem "can be viewed as a
non-convex one with objective and constraint functions continuously differentiable", so "any
local search method ... can be employed. However, our past experience ... suggests that SNOPT
... is particularly well suited for these problems." SNOPT is "An SQP Algorithm for Large-Scale
Constrained Optimization" (Gill, Murray, Saunders, SIAM J. Optim. 2002), so the effective
classical choice is a large-scale sparse SQP, not a dense light solver
([Grosso, Jamali, Locatelli, Schoen 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
The candidate solvers, by class:

- SLSQP (our current polish). Sequential Least Squares Programming, an SQP method for
  constrained nonlinear optimization. It handles equality and inequality constraints through
  KKT multipliers and checks optimality on the gradient of the Lagrangian and the constraint
  violation ([scipy SLSQP reference](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)).
  It is a dense quasi-Newton SQP, so it scales poorly as the O(N^2) non-overlap constraint set
  grows past mid N.
- L-BFGS-B penalty (the Discovery Loop seed engine). A limited-memory quasi-Newton method on a
  penalty reformulation. It scales to more variables cheaply but only enforces constraints
  through the penalty, so it needs the LP-radii step to price feasibility
  ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).
- IPOPT (interior-point NLP). "A software package for large-scale nonlinear optimization" that
  finds "(local) solutions" of problems whose functions "can be nonlinear and nonconvex, but
  should be twice continuously differentiable"
  ([IPOPT project README](https://raw.githubusercontent.com/coin-or/Ipopt/stable/3.14/README.md);
  [coin-or.github.io/Ipopt](https://coin-or.github.io/Ipopt/)). It uses second derivatives and
  sparse linear algebra, so it scales to the large, sparse packing NLP far better than dense
  SLSQP.

Impact: MEDIUM. A stronger, sparse, second-order NLP engine deepens each basin and raises the
per-restart yield, which matters most at higher N in our band where dense SLSQP slows. But
Grosso frames the local method as effective, not uniquely decisive; the global walk carries the
result. Cost in our architecture: swapping or adding the local NLP is a FIXED-PRIMITIVE change
(like `polish`), kept out of the genome, exact-or-deterministic given centers. It does not break
the verifier contract. Risk: IPOPT adds a heavy C++ dependency; a pure-Python plus scipy box may
prefer a tuned SLSQP or a penalty L-BFGS-B first. Test SLSQP versus penalty L-BFGS-B versus an
interior-point NLP as the MBH engine before adding a heavy dependency.

---

## 5. Lever 4, more and better LLM search: moderate, and it binds the SHELL not the NLP

The LLM evolves only the search shell (basin hopping, dissimilar-basin pool, jump sizes, escape,
warm start). It does not touch the fixed numeric primitives. So the LLM layer can only improve
the DIVERSITY and SCHEDULING of the search, which is lever 2 expressed through the genome.

The Discovery Loop ran a single provider (Claude Fable 5.1 via the Anthropic CLI), no crossover
between providers, and still broke 10 records, with the records landing at the SEED
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). This says the LLM layer was
not the binding constraint there; the seed and the targets were. More iterations past the seed
showed diminishing returns ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

Impact: MODERATE. A better model or higher reasoning effort helps the LLM invent stronger
diversity mechanics (lever 2), so its value is largely lever 2 by another path. At a hard
plateau the NUMERIC layer binds, so more LLM effort on the shell hits a ceiling set by the
target and the NLP. Cost in our architecture: none to the contract; the genome still never sees
the best-known value, and the verifier still recomputes the sum. Risk: spending LLM budget on
shell variety while the target stays hard wastes it. Pair a stronger model with a soft target
and a diversity mandate.

---

## 6. Lever 5, KKT-Newton exact polish: where it actually lands

What it is. Identify the active contact set (which wall and pairwise constraints are tight at
the local optimum), then take Newton steps on the KKT stationarity conditions of that fixed
active set. On a correct active set this gives quadratic or superlinear convergence to the exact
tangency configuration, sharper than a general SQP that keeps re-checking the full constraint
set. This is the iter-8 innovation in the Discovery Loop
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

What it buys. Better FINAL PRECISION inside a basin the search already found. It squeezes the
last digits of feasibility and sum, so it improves the feasibility-shrink margin at submission.
It is a refinement of the polish step, at the same layer as SLSQP.

What it does NOT buy. It does not change WHICH basin the search reaches. The active contact set
is fixed from the current configuration, so the Newton step stays inside the current basin and
converges to that basin's exact optimum. Our N=60 failure is that every candidate reaches the
SAME basin. A more exact polish of that same basin returns the same optimum to more digits, not
a better optimum. So KKT-Newton targets the wrong bottleneck for the plateau. This matches the
paper: the records were set at the seed, before iter 8, with no ablation crediting the polish
(Section 1). It also matches the classical record method: the Grosso, Jamali, Locatelli, Schoen
paper contains NO KKT, Newton, or active-set exact polish anywhere in its text (verified by full
read of the raw PDF). The classical record chasers reached putative optima with a strong NLP
local solver inside a good MBH/PBH walk, never with an exact contact-graph Newton step
([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

Impact: LOW to MODERATE. Real value only as a final-squeeze primitive, and only after the
global walk finds a basin better than the record. It cannot open a new basin. Cost in our
architecture: LOW and SAFE. It fits the fixed-primitive contract exactly, like `radii_lp` and
`polish`. It is deterministic given the centers and the active set, so it stays OUT of the
genome, and the independent verifier still recomputes and shrinks the sum. It does not break the
ungameable contract. Risk: active-set fragility. A wrong or ill-conditioned active set (a
near-tangent that is not truly active, a degenerate contact graph) gives a wrong or unstable
Newton step, so it needs a guard that falls back to SLSQP when the KKT system is singular or the
step raises infeasibility. Adding it is defensible as a cheap final-precision primitive, but do
NOT expect it to move the plateau.

---

## 7. Lever 6, more compute: lowest at a plateau

More restarts, more cores, or a longer per-candidate budget multiply the number of times the
search runs. At a plateau each run re-collapses into the same basin, so more runs return the
same optimum. Grosso et al. warn that plain multistart is weak because "the rapid increase in
the number of local minimizers suggests that Multistart can not be an efficient method for this
problem", and that MBH wins by seeding each search near the current minimizer, not by more
independent random restarts
([Grosso, Jamali, Locatelli, Schoen 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

Impact: LOWEST until the ceiling rises. Compute pays off only AFTER lever 1 (a softer target)
or lever 2 (a more diverse walk) raises the reachable ceiling. Then more restarts sample more
distinct basins. Cost in our architecture: none to the contract, only wall-clock and cores.
Risk: spending compute on a hard target with a collapsing walk is near-zero yield.

---

## 8. Ranking table

Impact is the expected effect on the achievable sum of radii for our N=28 to 114 push, given the
observed N=60 plateau. Contract-safe means it keeps the fixed-primitive plus independent-verifier
contract intact.

| Rank | Lever | Impact | Layer | Contract-safe | Primary evidence |
|------|-------|--------|-------|---------------|------------------|
| 1 | Target selection (soft N) | Highest | Run input | Yes | Records set at seed, not by polish ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)) |
| 2 | Global / PBH + diversity injection | High | Genome (shell) | Yes | MBH "dramatically more efficient than Multistart" ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)); PBH improved 32 optima for n<=130 ([disk paper](https://optimization-online.org/2006/03/1343/)) |
| 3 | Stronger local NLP (SQP or interior-point) | Medium | Fixed primitive | Yes | SNOPT (large-scale SQP) "particularly well suited" ([1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)); IPOPT scales to large sparse NLP ([IPOPT](https://raw.githubusercontent.com/coin-or/Ipopt/stable/3.14/README.md)) |
| 4 | More / better LLM search | Moderate | Genome (shell) | Yes | Single-provider run broke 10 records at seed ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)) |
| 5 | KKT-Newton exact polish | Low to moderate | Fixed primitive | Yes | Iter-8 add, diminishing-returns zone, no ablation ([arXiv:2609.05093](https://arxiv.org/html/2609.05093)) |
| 6 | More compute | Lowest at plateau | Infra | Yes | Multistart weak; the MBH walk is the lever, not more restarts ([Grosso et al.](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)) |

---

## 9. Recommendation

Do not build KKT-Newton polish first. It is contract-safe and cheap, and it is worth adding
LATER as a final-precision squeeze, but it targets the wrong bottleneck for the N=60 plateau and
its own source paper gives no evidence it moved a record.

Order of work by expected payoff:

1. Change the target. Move off N=60 to a soft Cantrell holdout in the 28 to 45 band, chosen by
   the soft-target procedure, re-pulled at attempt time
   ([author.txt](https://www.packomania.com/csqv/txt/author.txt);
   [docs/research/phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md)).
2. Add diversity to the walk. Make the genome run a true PBH with a dissimilarity cut, adaptive
   and occasional large jumps, island reset and migration, so it stops re-collapsing into one
   basin ([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf);
   [docs/research/phase-2-csqv-resume-and-warm-restart.md](phase-2-csqv-resume-and-warm-restart.md)).
3. Test the local engine. Benchmark SLSQP versus penalty L-BFGS-B versus an interior-point NLP
   as the fixed MBH engine at N=30 and N=90 before adding a heavy dependency
   ([scipy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html);
   [IPOPT](https://raw.githubusercontent.com/coin-or/Ipopt/stable/3.14/README.md)).
4. Then, and only then, add KKT-Newton as a final-precision primitive with an SLSQP fallback
   guard, to sharpen the feasibility margin of a basin the walk already beat the record with
   ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

---

## Sources

- Sibling Discovery Loop paper (Table 1 "Iter 0 (seed)", iter-8 KKT-Newton, seed solver,
  targets, no ablation, diminishing returns, independent verifier):
  https://arxiv.org/abs/2609.05093 and https://arxiv.org/html/2609.05093
- Grosso, Jamali, Locatelli, Schoen, "Solving the problem of packing equal and unequal circles
  in a circular container" (raw PDF read in full: SNOPT as the well-suited large-scale SQP local
  solver, the perturbation "not completely disrupted" delta rule, "if delta is too small ...
  basin of attraction of the current local minimizer", MBH "dramatically more efficient than
  Multistart", Multistart "clearly inferior", "in 19 cases we could obtain an improvement", and
  NO KKT/Newton/active-set method present):
  https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf and
  https://optimization-online.org/2008/06/1999/ (published as Springer JOGO,
  https://link.springer.com/article/10.1007/s10898-009-9458-3)
- Addis, Locatelli, Schoen, "Efficiently packing unequal disks in a circle" (Operations Research
  Letters 36(1), 37-42, 2008: the funneling-landscape conjecture, "improve over previously known
  putative optima in the range n <= 130 in as many as 32 instances", public contest with 155
  teams and 27,490 tentative solutions): https://optimization-online.org/2006/03/1343/ and
  https://www.sciencedirect.com/science/article/abs/pii/S0167637707000314
- SNOPT (Gill, Murray, Saunders, "SNOPT: An SQP Algorithm for Large-Scale Constrained
  Optimization", SIAM J. Optim. 12, 979-1006, 2002), cited as reference [9] in the 1999.pdf.
- Wales and Doye 1997, basin hopping and the funneling landscape:
  https://pubs.acs.org/doi/10.1021/jp970984n
- scipy SLSQP reference (SQP with KKT multipliers, equality and inequality constraints,
  optimality on Lagrangian gradient and constraint violation):
  https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html
- IPOPT project (interior-point optimizer for large-scale nonlinear programming, local
  solutions, twice-differentiable, second derivatives):
  https://raw.githubusercontent.com/coin-or/Ipopt/stable/3.14/README.md and
  https://coin-or.github.io/Ipopt/
- Packomania CSQV data for target selection: https://www.packomania.com/csqv/csqv.html and
  https://www.packomania.com/csqv/txt/author.txt
- Companion notes: docs/research/phase-2-csqv-record-methods-and-target-selection.md and
  docs/research/phase-2-csqv-resume-and-warm-restart.md

## Open questions

- No primary source gives a head-to-head sum-of-radii ablation of KKT-Newton versus SLSQP on the
  SAME configuration. Our own A/B at a fixed contact graph (final sum and feasibility margin,
  KKT-Newton versus SLSQP) would size the final-squeeze gain precisely.
- The exact per-N record gain from a softer target versus a stronger walk is unquantified. A
  local run that fixes the walk and sweeps the target, then fixes the target and sweeps the
  walk, would separate lever 1 from lever 2 for our box.
- IPOPT versus tuned SLSQP versus penalty L-BFGS-B as the MBH engine is untested on our
  hardware. Benchmark at N=30 and N=90 before adding a C++ dependency.
- Resolved this session: the raw 1999.pdf was read in full. It CONFIRMS the SNOPT local solver,
  the delta rule, and MBH beating Multistart. It CORRECTS the attribution of the "32 improved
  optima for n <= 130" figure and the funneling-landscape conjecture, which belong to the Addis,
  Locatelli, Schoen disk paper (ORL 2008), NOT the 1999.pdf. The 1999.pdf reports "19 cases"
  improved for equal circles in a circle and does not state a funneling conjecture in its body.

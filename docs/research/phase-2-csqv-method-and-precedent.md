# Phase 2 CSQV: method and precedent for an LLM heuristic-discovery loop

This file answers one question. How do you attack CSQV with an LLM-as-mutation
heuristic-discovery loop, and what compute does it need per N? CSQV is the problem
of maximizing the sum of radii of N variable-sized non-overlapping circles inside
the unit square.

Method. This file uses primary sources only. It uses the original papers, their
author code repositories, and the official benchmark registry. It does not use
secondary write-ups. Each quantified or attributed claim carries an inline
citation with a URL. The prose uses short active sentences and no em dash.

Access notes. The `deepmind.google` AlphaEvolve blog does not state the circle
result. It returned only the kissing-number example
([DeepMind blog](https://deepmind.google/discover/blog/alphaevolve-a-gemini-powered-coding-agent-for-designing-advanced-algorithms/)).
The `arxiv.org/abs/2506.13131` abstract page also omits the circle result. So this
file cites the AlphaEvolve circle numbers by identity and quotes primary-adjacent
sources that reproduce them: the Packomania CSQV registry, the OpenEvolve author
repository, and the arXiv 2601.05943 primary paper. The Hugging Face OpenEvolve
blog returned HTTP 403. This file cites the OpenEvolve GitHub repository instead.

---

## Headline finding

CSQV is a benchmark with a live registry. E. Specht maintains the "best known
packings of variable-sized circles in a square with maximized sum of radii,"
complete up to N = 100
([packomania.com/csqv](https://www.packomania.com/csqv/csqv.html)). The registry
lists N = 26 sum of radii 2.635983084919 and N = 32 sum of radii 2.939572771205
([packomania.com/csqv](https://www.packomania.com/csqv/csqv.html)).

An LLM heuristic-discovery loop already broke real CSQV records. The Discovery Loop
system broke 10 Packomania CSQV records for N in 101-114 for a total LLM cost of
$27.72 ([arXiv:2609.05093](https://arxiv.org/abs/2609.05093);
[HTML](https://arxiv.org/html/2609.05093)). This is the closest primary precedent
for the target task.

The genome that wins on CSQV is a full search program, not a bare construction
rule. The winning programs pair a global-search shell (basin hopping or island
parallelism) with two fixed exact primitives: a linear program for radii given
centers, and an SLSQP or L-BFGS-B local polish
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). AlphaEvolve and
ShinkaEvolve evolve a similar hybrid, but they keep the refinement fixed and evolve
the initialization plus the local-search orchestration
([arXiv:2601.05943](https://arxiv.org/html/2601.05943v2);
[ShinkaEvolve arXiv:2509.19349](https://arxiv.org/html/2509.19349v1)).

---

## 1. State-of-the-art methods for CSQV and hard circle packing

CSQV is a nonconvex problem with many local optima. The non-overlap constraints
make it nonconvex ([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091);
[arXiv:1601.07277](https://arxiv.org/pdf/1601.07277)). The main solver families:

Nonlinear programming (direct QCQP). You state the problem with variables x_i,
y_i, r_i, a linear objective sum r_i, boundary box constraints, and O(N^2)
pairwise non-overlap constraints. General nonlinear solvers attack it directly. The
arXiv 2601.05943 study solved CSQV with the commercial solver FICO Xpress and the
open-source solver SCIP through "straightforward nonlinear programming
formulations" ([arXiv:2601.05943](https://arxiv.org/html/2601.05943v2)). The FICO
Xpress examples repository ships a CSQV example named "Maximize the sum of radii of
N non-overlapping circles in a square"
([examples.xpress.fico.com](https://examples.xpress.fico.com/example.pl?id=circlepacking)).
Direct NLP finds local optima. It needs multi-start or global search to reach
records ([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091)).

Radii-by-LP subproblem. Given fixed circle centers, the optimal radii follow from a
linear program. Maximizing the sum of radii of disjoint balls with fixed centers
"can be expressed as a linear program"
([arXiv:1607.02184](https://arxiv.org/pdf/1607.02184)). This is the key CSQV
primitive. It decouples the hard geometry (where to put centers) from the easy
sizing (how big each circle can grow). The Discovery Loop seed solver uses exactly
this: "Penalty L-BFGS-B + LP radii"
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

Basin hopping. This alternates a random perturbation of the incumbent with a local
NLP solve, then accepts by a monotonic or threshold rule. Grosso and coauthors
developed monotonic basin hopping (MBH) and population basin hopping (PBH) for
equal-circle packing, with the NLP solver SNOPT for the local step
([Szabo survey PDF](https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf)). Basin
hopping is the standard record-chasing method for hard packing. The Discovery Loop
evolved solvers used "basin hopping with perturbation of the incumbent, followed by
exact SLSQP optimization"
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

Physics-inspired relaxation (billiard and force models). Lubachevsky and Graham
proposed a billiard simulation. It treats each circle as a rigid billiard and grows
the circles under collision forces. Nurmela and Ostergard used a molecular
repulsion model ([Szabo survey PDF](https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf)).
These methods dominate the classic equal-circle-in-square records. They apply to
CSQV as a starting-configuration generator, but they were designed for equal radii.

Pattern search and action-space methods. Later record methods use structured local
search. Examples are the action-space based global optimization algorithm for
packing circles in a square
([ScienceDirect S0305054815000027](https://www.sciencedirect.com/science/article/abs/pii/S0305054815000027)),
threshold accepting, greedy vacancy search, and TAMSASS
([Szabo survey PDF](https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf)).

Evolved search programs (the LLM route). See sections 2 to 4. An LLM writes and
mutates the whole search program. An exact verifier scores each candidate. This
route now holds the small-N CSQV records and 10 large-N CSQV records
([arXiv:2506.13131](https://arxiv.org/abs/2506.13131);
[arXiv:2609.05093](https://arxiv.org/abs/2609.05093)).

---

## 2. AlphaEvolve on circle packing

Problem and result. AlphaEvolve attacked CSQV at N = 26. It reached sum of radii
2.63586276 and improved the prior best of 2.634 from Friedman, 2012
([arXiv:2506.13131](https://arxiv.org/abs/2506.13131), corroborated by the CSQV
registry [packomania.com/csqv](https://www.packomania.com/csqv/csqv.html) and by
[arXiv:2601.05943](https://arxiv.org/html/2601.05943v2)). The AlphaEvolve paper
also reports the rectangle-of-perimeter-4 relaxation. The arXiv 2601.05943 study
reads AlphaEvolve as reaching N = 26 sum 2.63930 and N = 27 sum 2.69015 in that
relaxed variant ([arXiv:2601.05943](https://arxiv.org/html/2601.05943v2)).

The genome. AlphaEvolve evolved "a program representing a search heuristic," not a
static list of coordinates ([arXiv:2601.05943](https://arxiv.org/html/2601.05943v2)).
The evolved program constructs a configuration and then refines it under a fixed
time budget on the order of 1000 seconds
([arXiv:2506.13131 HTML](https://arxiv.org/html/2506.13131v1)). The author-adjacent
OpenEvolve reproduction confirms the structure. It is a two-phase program: a
construction phase that sets positions and radii, then a refinement phase that
calls `scipy.optimize.minimize` with the SLSQP method against a negative-sum-of-radii
objective and non-overlap plus boundary constraints
([OpenEvolve GitHub](https://github.com/algorithmicsuperintelligence/openevolve/tree/main/examples/circle_packing)).
So the genome is a construction plus refinement program. It is not a bare
construction heuristic and not a static configuration. The refinement solver is a
fixed library call.

Compute. The AlphaEvolve paper does not give a per-task compute figure for the
circle result in the accessible sections
([arXiv:2506.13131 HTML](https://arxiv.org/html/2506.13131v1)). AlphaEvolve in
general uses a Gemini-driven distributed pipeline of many LLM samples over many
iterations ([arXiv:2506.13131](https://arxiv.org/abs/2506.13131)). The OpenEvolve
reproduction reached sum 2.634292402141039, which is 99.97% of the AlphaEvolve
value, on commodity hardware
([OpenEvolve GitHub](https://github.com/algorithmicsuperintelligence/openevolve/tree/main/examples/circle_packing)).

---

## 3. The sibling project: Discovery Loop (arXiv:2609.05093)

This is the direct CSQV precedent. The paper is "LLM-Guided Program Evolution for
Circle Packing: Breaking 10 Packomania Records for $28"
([arXiv:2609.05093](https://arxiv.org/abs/2609.05093);
[HTML](https://arxiv.org/html/2609.05093)).

Variant. It targets CSQV exactly: "maximize the sum of radii of N variable-radius
circles in the unit square"
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

The genome. The LLM evolves a complete replacement solver each iteration, not a
patch. The seed solver is "Penalty L-BFGS-B + LP radii," a multi-start penalty
method with L-BFGS-B and an LP for optimal radii given positions
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). The evolved solvers
grew richer structure: basin hopping with SLSQP polish, island-model parallelism
with independent basin-hopping chains, a KKT-Newton exact polish on the contact
graph, and defect-migration moves that remove weak circles and re-equilibrate
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

Method. It uses a single provider, Claude Fable 5.1 through the Anthropic CLI. It
uses no crossover. Each iteration yields one candidate replacement solver. An
independent zero-tolerance checker verifies each candidate: circles inside the unit
square, no overlaps, an independent recomputation of the sum, and a feasibility
shrink. It ran 15 iterations, index 0 to 14
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

Note. The task brief describes the method as "paired providers plus crossover." The
primary paper does not support that. The primary paper reports a single provider
and no crossover ([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). See
open questions.

Results and scale. It broke 10 Packomania CSQV records for N in the range 101-114,
with gains of 2.4% to 5.4%. Packomania accepted the results after independent
zero-tolerance verification. The total LLM cost was $27.72, or $13.95 with early
stopping ([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

Compute. The full run took "approximately 8 hours overnight" on an Intel Core
i7-13700KF with 32 GB RAM under Windows 11. The solver is about 400 lines of
Python. Evaluations ran in parallel with 6 workers, and each target had a 120
second timeout ([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

---

## 4. Other LLM-driven and evolutionary circle-packing results

ShinkaEvolve (arXiv:2509.19349). It attacked CSQV at N = 26 and reached sum of
radii 2.635983283. It beat AlphaEvolve. It used only 150 program evaluations, a
strong sample-efficiency claim. Its best program combined golden-angle spiral
initializations, a hybrid gradient plus annealing local search, and explicit
escape mechanisms for local optima
([ShinkaEvolve arXiv:2509.19349 HTML](https://arxiv.org/html/2509.19349v1),
corroborated by the CSQV registry
[packomania.com/csqv](https://www.packomania.com/csqv/csqv.html)). So the genome is
again a construction plus refinement program, with the initialization and the
local-search schedule as the evolved parts.

OpenEvolve (author repository). This is the open reimplementation of the AlphaEvolve
circle example. Its main example reached sum 2.634292402141039 at N = 26, 99.97% of
the AlphaEvolve value
([OpenEvolve GitHub](https://github.com/algorithmicsuperintelligence/openevolve/tree/main/examples/circle_packing)).
A reproduction issue reports a further result of 2.635977
([OpenEvolve issue #156](https://github.com/algorithmicsuperintelligence/openevolve/issues/156)).

Global optimization revisited (arXiv:2601.05943). This primary paper checks the
LLM records against classic global NLP solvers. It solved CSQV and the rectangle
relaxation with FICO Xpress and SCIP, and it reports one improving square solution
at N = 32 sum 2.93957 versus a prior 2.93794
([arXiv:2601.05943](https://arxiv.org/html/2601.05943v2)). The lesson is that a
plain NLP solver on a good formulation is competitive at small N.

FunSearch, EoH, and ReEvo. These systems do not target circle packing in their
primary papers. FunSearch evolves a scoring function for cap set and online bin
packing ([FunSearch is cited in the companion file
docs/research/phase-2-xl-llm-evolution-precedent.md]). EoH and ReEvo evolve
construction or guidance functions for routing and bin packing
([ReEvo arXiv:2402.01145](https://arxiv.org/abs/2402.01145)). They matter here as
the source of the genome pattern, an evolved function inside a fixed backbone, not
as circle-packing results.

---

## 5. The genome question for CSQV

The candidate genomes are:
(a) a construction or placement heuristic that places circles one at a time and
sets each radius;
(b) a local-move or perturbation rule inside a fixed refinement loop;
(c) a full search or optimization program;
(d) a symbolic starting configuration refined by a fixed solver.

What the precedent supports. The strongest CSQV precedent evolves option (c), a full
search program, but with fixed exact sub-primitives inside it. Discovery Loop
evolved the whole solver and broke 10 records
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). AlphaEvolve and
ShinkaEvolve sit between (a) and (c). They evolve the construction and the
local-search orchestration, and they keep a fixed refinement solver
([arXiv:2601.05943](https://arxiv.org/html/2601.05943v2);
[ShinkaEvolve arXiv:2509.19349 HTML](https://arxiv.org/html/2509.19349v1)).

Why option (a) alone is weak. CSQV needs a numeric refinement to reach a record. A
pure one-shot placement leaves the sum below the local optimum. Every record system
runs an optimizer after construction, either SLSQP or L-BFGS-B or a KKT polish
([OpenEvolve GitHub](https://github.com/algorithmicsuperintelligence/openevolve/tree/main/examples/circle_packing);
[arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)).

Why option (d) alone is weak. A fixed solver on one symbolic start finds one local
optimum. CSQV has many local optima
([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091)). Records need many restarts
and escape moves, which is search, not a single solve.

The two fixed primitives to keep out of the genome. First, the radii-by-LP step,
which is exact given centers
([arXiv:1607.02184](https://arxiv.org/pdf/1607.02184);
[arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). Second, the SLSQP or
L-BFGS-B local polish. The LLM should evolve the parts above these: the
initialization pattern, the perturbation and defect-migration moves, the restart
and acceptance schedule, and the island structure
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093);
[ShinkaEvolve arXiv:2509.19349 HTML](https://arxiv.org/html/2509.19349v1)).

---

## 6. Compute scaling with N

Cost of one evaluation. One evaluation places N circles and refines them. The
refinement is an NLP with 3N decision variables, x_i, y_i, r_i, and O(N^2) pairwise
non-overlap constraints ([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091)). The
constraint count and the Jacobian both grow like N^2. A dense SLSQP or interior
step then grows worse than N^2 per iteration. Multi-start and basin hopping
multiply this by the number of restarts. The radii-by-LP step is cheap. It is
linear in the constraints for fixed centers
([arXiv:1607.02184](https://arxiv.org/pdf/1607.02184)).

What is tractable on modest hardware. Discovery Loop ran CSQV up to N = 114 in pure
Python with scipy, at a 120 second per-target timeout, with 6 parallel workers, on
an Intel Core i7-13700KF with 32 GB RAM
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). This shows that
N up to roughly the low hundreds is tractable per evaluation in pure Python plus
scipy on a desktop. On a smaller 4-core, 7 GB machine in pure Python plus numpy,
the same N range stays feasible, but with fewer parallel workers and a longer
wall-clock. The main limits are the 120 second style per-evaluation budget and the
number of restarts you can afford, not memory. The O(N^2) constraint matrix at
N = 100 is about 10^4 pairs, which fits easily in 7 GB.

What needs heavy compute. Two things push past a small machine. First, very large N
in the hundreds to thousands, where the O(N^2) constraints and the dense linear
algebra dominate and call for a C++ or GPU solver. Second, a wide evolutionary
search with many islands and many restarts per candidate, where the total is
evaluations times restarts times per-solve cost. Discovery Loop kept this cheap by
capping each target at 120 seconds and 15 iterations
([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). AlphaEvolve used a
fixed budget near 1000 seconds per program and a large distributed LLM pipeline
([arXiv:2506.13131 HTML](https://arxiv.org/html/2506.13131v1)).

---

## Sources

- CSQV registry, best known variable-radius square packings, complete to N = 100,
  E. Specht: https://www.packomania.com/csqv/csqv.html
- AlphaEvolve, Novikov et al., arXiv:2506.13131:
  https://arxiv.org/abs/2506.13131 and https://arxiv.org/html/2506.13131v1
- AlphaEvolve DeepMind blog:
  https://deepmind.google/discover/blog/alphaevolve-a-gemini-powered-coding-agent-for-designing-advanced-algorithms/
- Discovery Loop, arXiv:2609.05093:
  https://arxiv.org/abs/2609.05093 and https://arxiv.org/html/2609.05093
- ShinkaEvolve, Lange et al., arXiv:2509.19349:
  https://arxiv.org/html/2509.19349v1
- Global optimization for combinatorial geometry revisited, arXiv:2601.05943:
  https://arxiv.org/html/2601.05943v2
- OpenEvolve circle packing example, author repository:
  https://github.com/algorithmicsuperintelligence/openevolve/tree/main/examples/circle_packing
- OpenEvolve reproduction issue #156:
  https://github.com/algorithmicsuperintelligence/openevolve/issues/156
- Maximizing the sum of radii of disjoint balls, LP for fixed centers,
  arXiv:1607.02184: https://arxiv.org/pdf/1607.02184
- Circle packing convexification techniques, arXiv:2404.03091:
  https://arxiv.org/pdf/2404.03091
- Heuristic convex-over-nonconvex system, arXiv:1601.07277:
  https://arxiv.org/pdf/1601.07277
- FICO Xpress CSQV example:
  https://examples.xpress.fico.com/example.pl?id=circlepacking
- Szabo and Csendes survey, circle packing into the square, classic methods:
  https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf
- Action-space based global optimization for packing circles in a square,
  ScienceDirect: https://www.sciencedirect.com/science/article/abs/pii/S0305054815000027
- ReEvo, Ye et al., arXiv:2402.01145: https://arxiv.org/abs/2402.01145

## Open questions

- The task brief says Discovery Loop uses "paired providers plus crossover." The
  primary paper reports a single provider (Claude Fable 5.1) and no crossover
  ([arXiv:2609.05093 HTML](https://arxiv.org/html/2609.05093)). Confirm which
  version the brief refers to. A later revision may add pairing.
- The exact AlphaEvolve per-task compute for the N = 26 circle result is not in the
  accessible sections of arXiv:2506.13131. The Colab notebook may hold it.
- The precise CSQV registry credit for the N = 26 and N = 32 values needs a check
  against the per-N detail pages, to see if AlphaEvolve or ShinkaEvolve holds the
  current record entry ([packomania.com/csqv](https://www.packomania.com/csqv/csqv.html)).
- No primary source gives a pure-numpy, no-scipy runtime curve for CSQV. All record
  systems use scipy or an NLP solver. A local benchmark is needed to size the
  4-core, 7 GB target precisely.

## Recommended genome and tractable N range

Evolve a full search program, option (c), with two fixed exact primitives kept out
of the genome: a linear program that sets radii given centers, and an SLSQP or
L-BFGS-B local polish. Let the LLM evolve the parts that matter for escaping local
optima: the initialization pattern (for example a golden-angle spiral), the
perturbation and defect-migration moves, the restart and acceptance schedule, and
an island structure. This mirrors the only system that broke real CSQV records
cheaply, Discovery Loop, and it matches the construction-plus-refinement pattern of
AlphaEvolve and ShinkaEvolve. Score every candidate with an independent
zero-tolerance verifier that checks containment, non-overlap, and the recomputed
sum. On a 4-core, 7 GB machine in pure Python plus numpy and scipy, target N from
about 26 to the low hundreds, with a per-evaluation timeout near 120 seconds and a
small number of parallel workers. This range covers the open-record frontier that
Discovery Loop reached (N in 101-114) at a total LLM cost near $28. Push past N in
the many hundreds only with a C++ or GPU solver, because the O(N^2) constraints and
dense linear algebra then dominate the per-evaluation cost.

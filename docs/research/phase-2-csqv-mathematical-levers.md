# Phase 2 CSQV: novel mathematical levers for a sum-of-radii record solver

Research date: 2026-09-15. This note adds NEW mathematical structure. It does not repeat the
solver-improvement catalog, the polish-and-solver-levers ranking, the method-and-precedent note,
or the target-selection note. Those already cover the KKT-Newton polish, basin hopping, defect
migration, target selection, and the NLP engine menu. Read them first:
[phase-2-csqv-solver-improvement-catalog.md](phase-2-csqv-solver-improvement-catalog.md),
[phase-2-csqv-polish-and-solver-levers.md](phase-2-csqv-polish-and-solver-levers.md),
[phase-2-csqv-method-and-precedent.md](phase-2-csqv-method-and-precedent.md),
[phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md).

CSQV means variable-radius circles in a unit square. The objective maximizes the sum of the radii.
The circles must not overlap and must stay inside the square. The centers are the free variables. For
fixed centers the optimal radii solve an exact linear program (LP): maximize sum r_i subject to
r_i + r_j <= dist(c_i, c_j) for every pair, r_i <= distance to each wall, and r_i >= 0.

Method. This note uses primary sources only: peer-reviewed papers, arXiv, first-party solver and
library docs, and standard optimization and computational-geometry references. Each non-obvious claim
carries an inline URL or DOI. Some primary full texts sit behind publisher gates. Those claims are
marked as confirmed by metadata and index, not by a verbatim full-text quote. The prose uses short
active sentences and no em dash.

The current pipeline (context). Construct centers, growpush relaxation, `radii_lp` (a bounded-variable
primal simplex with lazy pairwise constraints), `scalable_polish` (L-BFGS penalty continuation on
centers), `center_polish` (projected gradient ascent on centers that maximizes the exact-LP sum, with
a CENTRAL FINITE-DIFFERENCE gradient that costs 4n LP solves per gradient), `defect_move`, then the
exact verifier. The code is [cpp/csqv/lp.hpp](../../cpp/csqv/lp.hpp) and
[cpp/csqv/polish.hpp](../../cpp/csqv/polish.hpp). The `center_polish` comment itself states the gap
this note attacks: "analytic duals are not exposed by the lazy LP".

---

## A. Analytic LP-value gradient via duals (envelope / Danskin). Verdict: ESTABLISHED.

The math. `f(C)` is the optimal value of a parametric LP. The parameter is the center vector `C`. The
center enters only the right-hand sides: the pairwise RHS `d_ij(C) = ||c_i - c_j||` and the wall RHS
`w_i(C)`. Two standard facts give the gradient in closed form.

1. LP sensitivity. The optimal value of an LP is convex and piecewise linear in the right-hand side
   `b`. On each linear piece the gradient with respect to `b` equals the optimal DUAL vector `y*` (the
   shadow prices). This is the standard sensitivity-analysis result. Source: Bertsimas and Tsitsiklis,
   "Introduction to Linear Optimization," Athena Scientific 1997, Chapter 5, which "develops a
   characterization of dual optimal solutions as subgradients of a suitably defined optimal cost
   function" ([publisher preface](http://athenasc.com/linoptpreface.html); ISBN 9781886529199).

2. Envelope / Danskin theorem. For a value function `f(x) = max over z of phi(x, z)` with a smooth
   `phi`, the directional derivative is `f'(x; h) = max over the argmax set of grad_x phi`. When the
   maximizer is unique, `f` is differentiable and `grad f(x) = grad_x phi(x, z*)`. Source: Danskin,
   "The Theory of Max-Min," Springer 1967 ([DOI 10.1007/978-3-642-46092-0](https://doi.org/10.1007/978-3-642-46092-0));
   Bertsekas, "Nonlinear Programming," Athena Scientific.

Applied here, the chain rule gives an explicit gradient. Let `y_ij >= 0` be the optimal dual on the
active pairwise constraint `(i, j)` and `mu_i >= 0` the optimal dual on an active wall of circle `i`.
Then

  grad_{c_i} f(C) = sum over active neighbors j of y_ij * (c_i - c_j) / ||c_i - c_j||
                    + (wall unit-normal terms scaled by mu_i).

Every term is analytic. One LP solve plus a sparse assembly over the tight contacts replaces the 4n LP
solves of the central-difference gradient. The unit vector `(c_i - c_j)/||c_i - c_j||` is the same
contact direction already computed in the penalty gradient
([polish.hpp](../../cpp/csqv/polish.hpp), `penalty_obj_grad`).

Degeneracy is the one real caveat, and it is honest. A jammed CSQV optimum is highly degenerate: a
circle often touches three or more neighbors, and record contact counts run near 3N
([contacts.txt](https://www.packomania.com/csqv/txt/contacts.txt)). At a degenerate primal optimum the
DUAL is non-unique, so the value function has a kink and only ONE-SIDED directional derivatives exist.
The two-sided gradient equals a unique `y*` exactly when the primal is non-degenerate. Source:
"Induced optimal partition invariancy in linear optimization," which states "the primal and dual
optimal solution entries are the derivatives of the optimal value function with respect to the
objective and right-hand side vector entries, respectively," and that at break points one computes
left and right derivatives ([arXiv:2008.02305](https://arxiv.org/pdf/2008.02305)); and
Ghaffari-Hadigheh et al., EJOR 2007
([sensitivity analysis](https://www.sciencedirect.com/science/article/abs/pii/S0377221705002808)).
The practical fix is standard: any dual-optimal vertex gives a valid SUBGRADIENT, and projected
subgradient ascent still climbs. Danskin gives the exact directional derivative as a max over the dual
face, so a line search that reads the true ascent direction stays correct.

How practitioners compute this. Differentiable optimization differentiates the KKT (or conic)
optimality conditions implicitly. OptNet does it for a QP layer "using techniques from sensitivity
analysis, bilevel optimization, and implicit differentiation," at "virtually no additional cost on top
of the solve" (Amos and Kolter, ICML 2017,
[arXiv:1703.00443](https://arxiv.org/abs/1703.00443)). cvxpylayers differentiates through a cone
program via an affine-solver-affine map (Agrawal, Amos, Barratt, Boyd, Diamond, Kolter, NeurIPS 2019,
[arXiv:1910.12430](https://arxiv.org/abs/1910.12430)). Both note that a pure LP is polyhedral, so its
solution map is piecewise constant and its derivative is zero or undefined at kinks; the common remedy
is a small quadratic regularizer that smooths the LP into a strongly convex QP with a well-defined
gradient. For our VALUE gradient we do not need the solution derivative. We need the RHS gradient, which
is the dual, and which the simplex already computes.

Our own solver already computes what is needed. `solve_reduced_lp`
([lp.hpp](../../cpp/csqv/lp.hpp)) maintains a basis and prices reduced costs, so the dual vector is
`c_B^T B^{-1}` at optimum, available at the final basis for near-zero extra cost. The build is to return
it, not to add a solver.

ROI: HIGH, and DISTINCT from the existing catalog. The catalog optimized the NLP polish; it never
addressed the `center_polish` gradient. Replacing 4n LP solves per gradient with 1 solve plus a sparse
contact assembly is roughly a 4n-fold cut in the dominant cost of the terminal exact squeeze. At N=90
that is about 360 LP solves collapsed to 1 per gradient step. This does not open a new basin (the
existing catalog is right that precision alone cannot beat a re-derived record), but it converts the
exact center polish from an expensive final touch into a cheap step usable inside the search loop.

---

## B. Sparse contact-graph Newton as a terminal exact refiner. Verdict: ESTABLISHED (standard practice).

The existing notes already cover the KKT-Newton polish and rank it as precision-only. The NEW structure
here is the rigidity and tangency view, and what record refiners actually do to reach 15-plus digits.

The math. A jammed, locally maximal packing is a RIGID contact framework. Its tight contacts and wall
tangencies form a SQUARE system of algebraic (tangency and stationarity) equations in the centers and
radii. Newton's method on that square system converges quadratically to machine precision, once the
active set is correct. The rigidity and jamming basis is standard: Connelly rigidity of packings, and
the linear-programming jamming test of Donev, Torquato, Stillinger, and Connelly, "A linear programming
algorithm to test for jamming in hard-sphere packings," J. Comput. Phys. 197 (2004) 139-166
([DOI 10.1016/j.jcp.2003.10.047](https://doi.org/10.1016/j.jcp.2003.10.047)). That test linearizes the
impenetrability constraints about the packing and solves an LP on the contact network to classify local,
collective, and strict jamming. The Torquato and Stillinger review confirms "rigorous and efficient
linear-programming algorithms have been devised to assess whether a particular sphere packing is
locally, collectively, or strictly jammed" (Rev. Mod. Phys. 82 (2010) 2633,
[arXiv:1008.2982](https://arxiv.org/abs/1008.2982)).

What actual record refiners use. Nurmela and Ostergard reach high-precision equal-circle packings by
relaxing a repulsive-energy model, then stiffening the exponent stepwise until the potential acts as a
hard contact, which is equivalent to solving the nonlinear tangency system (Discrete and Computational
Geometry 18 (1997) 111-120, [DOI 10.1007/PL00009306](https://doi.org/10.1007/PL00009306)). Graham and
Lubachevsky grow rigid billiards until jammed ("Repeated patterns of dense packings of equal disks in a
square," Electron. J. Combin. 3(1) (1996) #R16,
[PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r16/pdf)). A modern record
refinement paper states record coordinates are refined "up to numerical machine precision" and warns of
the "false contact" problem, where two disks are extremely close but not truly tangent; it uses a
variance-minimizing refinement to about 1e-20 to avoid committing the wrong contact set (Amore et al.,
"Circle packing in regular polygons," [arXiv:2212.12287](https://arxiv.org/abs/2212.12287)). Certified,
computer-assisted OPTIMA come from interval branch-and-bound with interval Newton, which proved the
optimal equal-circle packings for N=28, 29, 30 (Markot and Csendes, "A new verified optimization
technique for the packing circles in a unit square problems," SIAM J. Optimization 16(1) (2005)
193-219).

ROI: LOW to MODERATE, matching the existing verdict. It is a precision and throughput lever, not a basin
lever. The genuinely new and cheap add here is a rigidity guard: use the Donev LP jamming test on the
detected contact graph BEFORE a Newton step, to reject a false-contact active set and to detect a
non-rigid (under-jammed) configuration that a Newton step would destabilize. That guard removes the
active-set fragility the catalog flagged, at the cost of one small LP. The "false contact" tolerance is
the same tolerance the catalog noted the sibling paper never specifies.

---

## C. Homotopy / numerical continuation across contact make-break events. Verdict: SPECULATIVE.

The math. Numerical continuation tracks a solution path of a parametrized system `H(x, lambda) = 0` by
predictor-corrector steps, and handles turning points and bifurcations (Allgower and Georg,
"Introduction to Numerical Continuation Methods," SIAM Classics 45, 2003,
[DOI 10.1137/1.9780898719154](https://doi.org/10.1137/1.9780898719154)). The packing analog is the
Lubachevsky-Stillinger inflation algorithm, where the radius is a continuation parameter: particles
start as points and grow uniformly while an event-driven molecular dynamics processes collisions, until
the system jams (Lubachevsky and Stillinger, "Geometric properties of random disk packings," J. Stat.
Phys. 60 (1990) 561-583, [DOI 10.1007/BF01025983](https://doi.org/10.1007/BF01025983)). The jamming
literature does treat packing generation as a parameter sweep in the growth rate, and it treats contact
make and break as non-smooth, discontinuous events that reorganize the contact network (Torquato and
Stillinger 2010, [arXiv:1008.2982](https://arxiv.org/abs/1008.2982)).

Why SPECULATIVE. The pieces are each sourced, but no primary source frames CSQV basin-to-basin motion as
"piecewise-smooth continuation with event detection at contact make and break," and none applies
predictor-corrector path following to the SUM-OF-RADII objective. Inflation targets equal or grown
radii, not a maximized sum of variable radii. The map from an inflation path to a better CSQV basin is a
synthesis, not a result. A continuation across a contact event is a non-smooth event that needs explicit
detection and restart of the corrector, which is exactly where the theory is hardest.

ROI: LOW to MODERATE, and higher RISK than the other levers. It is the one lever with a genuine
basin-change mechanism beyond defect migration, because a continuation path can walk a configuration
through a contact-topology change smoothly rather than by a random jump. But it is unproven for this
objective, and the existing defect-migration move is the cheaper, precedented basin-change lever. Treat
continuation as a research probe, not a near-term build.

---

## D. Laguerre / power-diagram-guided insertion. Verdict: ESTABLISHED geometry, with a KEY correction.

The correction first. For inserting a disk that must not overlap existing VARIABLE-radius disks, the
governing structure is NOT the power (Laguerre) diagram. The power diagram uses the squared power
distance `||x - c_i||^2 - r_i^2`, which does not equal the disjointness gap for variable radii
(Aurenhammer, "Power diagrams: properties, algorithms and applications," SIAM J. Comput. 16(1) (1987)
78-96, [DOI 10.1137/0216006](https://doi.org/10.1137/0216006);
[PDF](https://www.cs.jhu.edu/~misha/Spring16/Aurenhammer87.pdf)). The correct structure is the
ADDITIVELY WEIGHTED Voronoi diagram, also called the Apollonius diagram, whose distance is the linear
offset `||x - c_i|| - r_i`. CGAL states it directly: "delta(x, P_i) = ||x - c_i|| - w_i" and "the
Apollonius diagram, also called the Additively weighted Voronoi diagram," and when the weights are
disk radii the diagram "can be viewed as the Voronoi diagram of the set of circles"
([CGAL 2D Apollonius Graphs](https://doc.cgal.org/latest/Apollonius_graph_2/index.html)).

The math. The largest empty circle among point sites is centered at a Voronoi vertex if it is strictly
interior, or on a Voronoi edge or the hull boundary otherwise (Preparata and Shamos, "Computational
Geometry: An Introduction," Springer 1985; Chew and Drysdale, "Finding largest empty circles with
location constraints," Dartmouth TR 1986, [digital commons](https://digitalcommons.dartmouth.edu/cs_tr/29/)).
The weighted generalization is exact: the largest disk insertable among existing disks of radii `r_i` is
centered at a VERTEX of the additively weighted (Apollonius) diagram, and its radius is the maximum over
the diagram of the additively weighted clearance `||x - c_i|| - r_i` (also bounded by the wall distance).
CGAL computes the Apollonius graph EXACTLY when the kernel uses an exact number type
([CGAL reference](https://doc.cgal.org/latest/Apollonius_graph_2/group__PkgApolloniusGraph2Ref.html)).

ROI: MODERATE, and DISTINCT from the existing catalog. The catalog's defect-migration move reinserts
weak circles into "the largest holes" by a heuristic. The Apollonius diagram replaces that heuristic
with the exact largest-insertable-disk site. Two concrete uses. First, principled N to N+1 continuation:
to grow a good N-packing to N+1, insert the new disk at the Apollonius vertex of maximum clearance,
which is the provably largest single insertion. Second, a targeted defect reinsertion: after removing a
weak circle, place it at the best Apollonius vertex rather than a heuristic hole. This is a construction
and basin-seeding lever, not a precision lever, so it can seed a NEW basin. Cost caveat: an exact
Apollonius build is a C++ or CGAL dependency; a numerical additively weighted nearest-site scan over a
grid is a cheap first approximation that needs no dependency.

---

## E. Certified upper bounds for stopping. Verdict: SPECULATIVE, and honestly WEAK for this objective.

The honest finding. No published bound specifically upper-bounds the maximum SUM OF RADII of
non-overlapping variable disks with FREE centers in a fixed region. The famous packing bounds bound the
wrong quantity. The Cohn-Elkies linear-programming bound bounds the DENSITY of EQUAL spheres (Cohn and
Elkies, "New upper bounds on sphere packings I," Annals of Math. 157 (2003) 689-714,
[arXiv:math/0110009](https://arxiv.org/abs/math/0110009)). The Bachoc-Vallentin semidefinite bound
bounds KISSING numbers of EQUAL spheres (J. Amer. Math. Soc. 21 (2008) 909-924,
[arXiv:math/0608426](https://arxiv.org/abs/math/0608426)). Neither gives a sum-of-radii bound, and the
analogy to variable radii in a square is loose.

What exists is elementary or generic. The area bound is elementary: disjoint disks in a region of area
`A` give `pi * sum r_i^2 <= A`, so by Cauchy-Schwarz `sum r_i <= sqrt(N * A / pi)`. This is a
first-principles bound, and it is loose, since equality needs equal radii and full area coverage, which
disks cannot achieve. The exact-value companion, for FIXED centers, is Eppstein's LP (see Section F),
which is an exact solve, not a free-center bound
([arXiv:1607.02184](https://arxiv.org/abs/1607.02184)).

The one rigorous route to a stopping THEOREM is spatial branch-and-bound. The maximum sum of radii with
free centers is a nonconvex QCQP. A global solver builds McCormick convex underestimators or an RLT
relaxation and returns a valid DUAL BOUND, so a nonzero gap at termination is a rigorous optimality
certificate (McCormick, "Computability of global solutions to factorable nonconvex programs: Part I,"
Math. Programming 10 (1976) 147-175, [DOI 10.1007/BF01580665](https://doi.org/10.1007/BF01580665);
Sherali and Adams, "A Reformulation-Linearization Technique," Kluwer 1999). A recent paper solved CSQV
exactly this way with FICO Xpress and SCIP, "we consider the problem of packing n circles inside a
rectangle such that the sum of their radii is maximized" (Berthold, Kamp, Mexi, Pokutta, Polik, "Global
Optimization for Combinatorial Geometry Problems Revisited in the Era of LLMs,"
[arXiv:2601.05943](https://arxiv.org/abs/2601.05943)). Honest caveat: that paper explicitly does NOT
report the dual bounds or gaps ("we will not comment on the dual bounds and how they could be
improved"), so the certified-gap-as-theorem claim is sound in PRINCIPLE but is not yet evidenced for
CSQV, and spatial B&B is tractable only at small N.

ROI: LOW for beating a record, because a certified gap does not raise the achievable sum. But it is the
ONLY route that turns stopping into a theorem, and it has a clear niche: at small N (near the proven
band, N=28-33 in the equal-circle analogs) a spatial B&B run could CERTIFY that a soft Cantrell holdout
is or is not already optimal, which reframes a target from "hard" to "provably closed" and saves the
search from grinding a solved instance.

---

## F. Combinatorial structure of the radii LP. Verdict: ESTABLISHED, and the STANDOUT result.

The radii LP is NOT a generic LP. It has exactly two variables per pairwise constraint, both with a
plus sign, plus box bounds. This is a well-studied combinatorial class, and for the sum-of-radii
objective the exact combinatorial solver is already published.

The specific result. Eppstein, "Maximizing the Sum of Radii of Disjoint Balls or Disks," solves exactly
our fixed-center problem: for fixed centers, maximize the sum of radii of disjoint balls. He shows it is
an LP whose DUAL is a minimum-weight cycle cover or matching structure, and that in the plane it solves
in about O(n^{3/2}) time by exploiting that special structure, far faster than a general simplex
([arXiv:1607.02184](https://arxiv.org/abs/1607.02184)). This is the direct, primary, sum-of-radii
citation. It gives both a faster radii solve and the dual multipliers as a by-product, which are exactly
the multipliers Section A needs.

The general structure that explains it. Constraints with two variables each define the IP2 or 2VPI
class. The LP relaxation of a 2-variable-per-inequality system is solvable in the complexity of a
MIN-CUT or MAX-FLOW, and it yields HALF-INTEGRAL optima (Hochbaum, Megiddo, Naor, Tamir, "Tight bounds
and 2-approximation algorithms for integer programs with two variables per inequality," Math.
Programming 62 (1993) 69-83, [DOI 10.1007/BF01585160](https://doi.org/10.1007/BF01585160)). The
constraint `r_i + r_j <= 1` is precisely the fractional vertex-packing (stable-set) polytope, whose LP
relaxation is half-integral and computable by BIPARTITE matching or max-flow on the bipartite double
cover (Nemhauser and Trotter, "Vertex packings: structural properties and algorithms," Math. Programming
8 (1975) 232-248, [DOI 10.1007/BF01580444](https://doi.org/10.1007/BF01580444)); the fractional matching
polytope is half-integral (Balinski 1965, Management Science 12(3):253-313; Schrijver, "Combinatorial
Optimization," Springer 2003). The plus-plus sign is the non-monotone case; the bipartite double cover
converts each sum-constraint `r_i + r_j <= d` into difference form on the doubled variable set, which is
why the solve reduces to a min-cut. This also connects to the classic result that DIFFERENCE constraints
`r_i - r_j <= d` reduce to shortest paths by Bellman-Ford (Cormen, Leiserson, Rivest, Stein,
"Introduction to Algorithms," Section 24.4).

Duals for free. The reduction is itself a duality statement. Max-flow min-cut is LP strong duality, so
the combinatorial solve returns the optimal dual (the min-weight cover, or the potentials) alongside the
primal. So a combinatorial radii solve exposes the duals that hypothesis A needs at no extra cost. This
is the tight coupling: F gives A its inputs for free.

ROI: HIGH, and DISTINCT from the existing catalog, which treated `radii_lp` as an already-exact
black box whose only lever was throughput via sparse matrices and a HiGHS hot start. The new structure
says the radii LP has a purpose-built O(n^{3/2}) solver with a known dual, so the replacement is not a
generic speedup but an exact combinatorial algorithm that also unlocks the analytic gradient of A. The
current lazy simplex ([lp.hpp](../../cpp/csqv/lp.hpp)) is already fast because it keeps the active set
tiny, so the raw speed win over the simplex may be modest at our N; the load-bearing win is the DUAL it
hands back for the gradient.

---

## Ranked build-next short list (rigor times ROI), distinct from the existing catalog

The existing catalog ranked a sparse second-order NLP polish, defect migration, an L-BFGS-B engine A/B,
barrier-crossing acceptance, a sparse-plus-hot-start LP, and a KKT-Newton polish. The list below adds
NEW items and does not repeat those.

| Rank | New lever | Hypotheses | Verdict | Effect class | Contract-safe |
|------|-----------|------------|---------|--------------|---------------|
| 1 | Analytic dual gradient for `center_polish` (return the simplex duals, assemble the contact-direction gradient; subgradient fallback at degeneracy) | A + F | Established | Throughput (about 4n-fold cut in the exact center-polish gradient cost) | Yes, deterministic given centers |
| 2 | Apollonius (additively weighted Voronoi) insertion for N to N+1 growth and defect reinsertion, replacing the "largest holes" heuristic | D | Established geometry | Basin seeding (construction) | Yes |
| 3 | Eppstein exact combinatorial radii solver, replacing the general simplex and exposing the dual for lever 1 | F | Established | Throughput plus duals for free | Yes, exact |
| 4 | Rigidity guard on the contact-graph Newton step (Donev LP jamming test to reject false contacts before Newton) | B | Established | Precision plus robustness | Yes |
| 5 | Small-N certified upper bound via spatial branch-and-bound to close or confirm a soft target | E | Speculative for records, rigorous as a theorem | Stopping certificate | Yes, offline |
| 6 | Continuation across contact make-break events as a basin-change probe | C | Speculative | Basin change (unproven) | Yes, research only |

Reasoning for the top build. Build lever 1 first. It is the single new lever with both high rigor and
high ROI, it attacks a concrete measured cost (4n LP solves per `center_polish` gradient), the duals it
needs are already computed inside `solve_reduced_lp`, and it is contract-safe because the verifier still
recomputes the sum. Build lever 3 alongside it, because the Eppstein solver hands the exact dual to
lever 1 and is the primary sum-of-radii algorithm for this exact LP. Build lever 2 next, because it is
the only NEW basin-seeding lever here and it upgrades the existing defect move from a heuristic to an
exact largest-insertion. Levers 4, 5, 6 are guards and probes, not record levers.

---

## Bottom line

The strongest new mathematical structure is the tight A-plus-F coupling. The fixed-center radii problem
is not a generic LP but a purpose-built sum-of-radii combinatorial LP (Eppstein), whose dual is a
min-weight matching or cover, and whose dual is exactly the analytic gradient of the LP value with
respect to the centers (envelope and Danskin, LP sensitivity). So the current `center_polish`, which
spends 4n LP solves per gradient on central differences, can compute the same gradient from ONE solve
plus a sparse contact assembly, with a principled subgradient fallback at the degenerate contacts that
CSQV optima always have. This is a real, established, contract-safe throughput lever that the existing
catalog missed, because the catalog treated the radii LP and the center-polish gradient as fixed. It
does not, by itself, beat a re-derived record, which stays a basin problem. The one new basin lever here
is Apollonius-guided insertion, which turns the existing defect move from a heuristic into an exact
largest-insertion. Certified upper bounds and contact-event continuation are honest but weak for beating
a record: no specific sum-of-radii bound exists, spatial branch-and-bound certifies a gap only at small
N, and continuation across contact events is an unproven synthesis for this objective.

---

## Sources

- Bertsimas and Tsitsiklis, "Introduction to Linear Optimization," Athena Scientific 1997, Ch. 5
  (dual optimal solutions as subgradients of the optimal cost function): http://athenasc.com/linoptpreface.html
- Danskin, "The Theory of Max-Min and its Application to Weapons Allocation Problems," Springer 1967:
  https://doi.org/10.1007/978-3-642-46092-0 ; Bertsekas, "Nonlinear Programming," Athena Scientific.
- LP value-function derivatives and degeneracy: "Induced optimal partition invariancy in linear
  optimization," https://arxiv.org/pdf/2008.02305 ; Ghaffari-Hadigheh et al., EJOR 2007,
  https://www.sciencedirect.com/science/article/abs/pii/S0377221705002808
- OptNet, Amos and Kolter, ICML 2017: https://arxiv.org/abs/1703.00443
- Differentiable Convex Optimization Layers (cvxpylayers), Agrawal et al., NeurIPS 2019:
  https://arxiv.org/abs/1910.12430
- Donev, Torquato, Stillinger, Connelly, "A linear programming algorithm to test for jamming in
  hard-sphere packings," J. Comput. Phys. 197 (2004) 139-166: https://doi.org/10.1016/j.jcp.2003.10.047
- Torquato and Stillinger, "Jammed hard-particle packings," Rev. Mod. Phys. 82 (2010) 2633:
  https://arxiv.org/abs/1008.2982
- Nurmela and Ostergard, "Packing up to 50 equal circles in a square," DCG 18 (1997) 111-120:
  https://doi.org/10.1007/PL00009306
- Graham and Lubachevsky, "Repeated patterns of dense packings of equal disks in a square," EJC 3(1)
  (1996) #R16: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r16/pdf
- Amore et al., "Circle packing in regular polygons" (refinement to machine precision, false-contact
  problem): https://arxiv.org/abs/2212.12287
- Markot and Csendes, "A new verified optimization technique for the packing circles in a unit square
  problems," SIAM J. Optimization 16(1) (2005) 193-219 (interval branch-and-bound, computer-assisted
  proof for N=28, 29, 30).
- Lubachevsky and Stillinger, "Geometric properties of random disk packings," J. Stat. Phys. 60 (1990)
  561-583: https://doi.org/10.1007/BF01025983
- Allgower and Georg, "Introduction to Numerical Continuation Methods," SIAM Classics 45, 2003:
  https://doi.org/10.1137/1.9780898719154
- Aurenhammer, "Power diagrams: properties, algorithms and applications," SIAM J. Comput. 16(1) (1987)
  78-96: https://doi.org/10.1137/0216006 ; PDF https://www.cs.jhu.edu/~misha/Spring16/Aurenhammer87.pdf
- CGAL 2D Apollonius Graphs (additively weighted Voronoi, exact predicates):
  https://doc.cgal.org/latest/Apollonius_graph_2/index.html and
  https://doc.cgal.org/latest/Apollonius_graph_2/group__PkgApolloniusGraph2Ref.html
- Preparata and Shamos, "Computational Geometry: An Introduction," Springer 1985 (largest empty circle);
  Chew and Drysdale, Dartmouth TR 1986: https://digitalcommons.dartmouth.edu/cs_tr/29/
- Cohn and Elkies, "New upper bounds on sphere packings I," Annals of Math. 157 (2003) 689-714:
  https://arxiv.org/abs/math/0110009
- Bachoc and Vallentin, "New upper bounds for kissing numbers from semidefinite programming," JAMS 21
  (2008) 909-924: https://arxiv.org/abs/math/0608426
- McCormick, "Computability of global solutions to factorable nonconvex programs: Part I," Math.
  Programming 10 (1976) 147-175: https://doi.org/10.1007/BF01580665 ; Sherali and Adams,
  "A Reformulation-Linearization Technique," Kluwer 1999.
- Berthold, Kamp, Mexi, Pokutta, Polik, "Global Optimization for Combinatorial Geometry Problems
  Revisited in the Era of LLMs" (CSQV via FICO Xpress and SCIP; dual bounds not reported):
  https://arxiv.org/abs/2601.05943
- Eppstein, "Maximizing the Sum of Radii of Disjoint Balls or Disks" (exact fixed-center LP, dual =
  min-weight cycle cover, O(n^{3/2}) in the plane): https://arxiv.org/abs/1607.02184
- Hochbaum, Megiddo, Naor, Tamir, "Tight bounds and 2-approximation algorithms for integer programs
  with two variables per inequality," Math. Programming 62 (1993) 69-83:
  https://doi.org/10.1007/BF01585160
- Nemhauser and Trotter, "Vertex packings: structural properties and algorithms," Math. Programming 8
  (1975) 232-248: https://doi.org/10.1007/BF01580444 ; Balinski 1965, Management Science 12(3):253-313;
  Schrijver, "Combinatorial Optimization," Springer 2003.
- Cormen, Leiserson, Rivest, Stein, "Introduction to Algorithms," Section 24.4 (difference constraints
  and shortest paths).
- Our code: cpp/csqv/lp.hpp (`solve_reduced_lp` computes the basis and reduced costs) and
  cpp/csqv/polish.hpp (`center_polish` central-difference gradient, 4n LP solves per gradient).

## Open questions and flags

- Bertsimas and Tsitsiklis: the chapter (5) and the subgradient characterization are confirmed by the
  publisher preface and corroborating sources; the exact theorem index was not pulled from full text.
- "Duals for free" from the combinatorial solve follows rigorously from max-flow min-cut LP duality; no
  single verbatim primary sentence asserts it for this exact LP.
- Hochbaum-Megiddo-Naor-Tamir and Nemhauser-Trotter full texts sit behind a Springer gate; DOIs and the
  half-integrality and min-cut statements are confirmed by metadata and multiple corroborating sources.
- Nurmela-Ostergard energy wording and Markot-Csendes "interval Newton" phrasing are paraphrases from
  index and abstract; the methods (energy relaxation, verified interval branch-and-bound) are confirmed.
- No primary source frames CSQV basin motion as piecewise-smooth continuation with contact-event
  detection (Section C); present it as an analogy, not a cited result.
- arXiv:2601.05943 does not report CSQV dual bounds, so the certified-gap-as-theorem claim (Section E)
  is sound in principle but not yet evidenced for CSQV.
- The Eppstein O(n^{3/2}) speed win over our already-lazy simplex is unmeasured at our N; the primary
  value is the dual it exposes for the Section A gradient. Benchmark before replacing the simplex.
- The degenerate-dual subgradient behavior of the Section A gradient at true CSQV optima (many 3-plus
  contact circles) is untested; validate the subgradient ascent against the current finite-difference
  center polish on a fixed champion at N=30 and N=90.

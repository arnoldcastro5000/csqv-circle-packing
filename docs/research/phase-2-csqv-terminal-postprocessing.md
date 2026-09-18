# Phase 2 CSQV: terminal post-processing to the exact local optimum

Research date: 2026-09-16. This note answers ONE question. What algorithm squeezes the SINGLE best
CSQV champion to its EXACT local optimum, as a ONE-TIME terminal operation where cost is not a
constraint? This is not a per-restart inner-loop primitive. The one-shot framing removes the
throughput objection that ranked the contact-graph Newton polish fifth of six for in-loop use.
Precision-only is exactly what a terminal squeeze wants.

The primitives live in
[cpp/csqv/lp.hpp](../../cpp/csqv/lp.hpp) (`radii_lp`, `solve_reduced_lp` with duals) and
[cpp/csqv/polish.hpp](../../cpp/csqv/polish.hpp) (`center_polish`, `center_gradient`,
`center_polish_analytic`).

CSQV means variable-radius circles in the unit square. The objective maximizes the sum of the radii.
The circles must not overlap and must stay inside the square. For FIXED centers the optimal radii are
an exact linear program (LP), so all remaining slack sits in the center arrangement.

Method. This note uses primary sources only: the scipy and IPOPT solver docs and their underlying
method papers, the Nurmela-Ostergard and Amore refinement papers, the Donev jamming paper, the
interval-arithmetic packing proofs, the arbitrary-precision library docs, and the Packomania registry.
Each non-obvious claim carries an inline URL or DOI. Some full texts sit behind publisher gates. Those
are marked as confirmed by an open-access reproduction or by index and abstract, not by a verbatim
primary quote. The prose uses short active sentences and no em dash.

---

## The problem, stated precisely

The current terminal squeeze is `center_polish_analytic`
([polish.hpp](../../cpp/csqv/polish.hpp)). It is a FIRST-ORDER projected-gradient ascent on the
exact-LP value over the centers. The gradient is the LP optimal dual, which is the analytic gradient
of the value function by the envelope (Danskin) theorem (Danskin, "The Theory of Max-Min," Springer
1967, [DOI 10.1007/978-3-642-46092-0](https://doi.org/10.1007/978-3-642-46092-0); LP sensitivity in
Bertsimas and Tsitsiklis, "Introduction to Linear Optimization," Athena Scientific 1997, Ch. 5). It
reaches a first-order stationary point. It STALLS about 1 to 2e-5 short of the exact optimum at
degenerate contacts, because the value function has a KINK there and the optimum sits AT the kink
where the two-sided derivative does not exist (LP value-function derivatives and break points,
[arXiv:2008.02305](https://arxiv.org/pdf/2008.02305)).

The evidence that the missing gap is a SAME-BASIN second-order squeeze, not a new basin, is direct and
strong. Packomania accepted our N=121 submission at sum 5.799069987171, then published
5.799103501951, a gain of +3.35e-5 ([csqv.html](https://www.packomania.com/csqv/csqv.html), ref [21]).
The published packing is OUR arrangement, transposed, with centers nudged about 1e-6, and our own
exact LP on THEIR centers reproduces 5.799103 to 1e-12. So the
summation was correct and the gap was an un-squeezed optimum, not a bug and not a better basin. Our
own `center_polish_analytic` reproduces this behavior: on N=122 it lifted 5.816908404257 to
5.816928309336, a gain of +1.99e-5, feasible at zero tolerance. So the terminal
post-processor needs a TRUE SECOND-ORDER solve to the exact local optimum of the CURRENT basin. It
does not need a new basin.

The objective structure makes the second-order solve tractable. The objective sum(r) is LINEAR. The
containment constraints are LINEAR. ALL nonlinearity sits in the pairwise distance constraints
r_i + r_j <= dist(c_i, c_j), and each has a simple sparse 4x4 Hessian block in (x_i, y_i, x_j, y_j). So
the Lagrangian Hessian is exact, analytic, and sparse (catalog Section 1). A jammed CSQV optimum is a
RIGID contact framework: the tight contacts and wall tangencies form a square algebraic system in the
centers and radii, on which Newton converges quadratically once the active set is correct.

---

## A. Full-space sparse second-order NLP (trust-constr / IPOPT, exact Hessian). Verdict: BEST one-shot.

The math. State the full NLP in the 3N variables (x_i, y_i, r_i): a linear objective sum(r), 4N
linear wall constraints, and O(N^2) nonlinear non-overlap constraints, most of them slack. A
second-order constrained solver takes Newton-type steps on the KKT system with the EXACT sparse
Lagrangian Hessian. Near a good champion it converges quadratically, and it IDENTIFIES the active set
itself, so it needs no contact-graph guess.

scipy `trust-constr` is the in-library one-shot. It "is the most versatile constrained minimization
algorithm implemented in SciPy and the most appropriate for large-scale problems"
([trust-constr](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html)). It
runs a "Byrd-Omojokun Trust-Region SQP method" for equality constraints and switches to a
"trust-region interior point method" when inequalities are present
([minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)). It
accepts a SPARSE Jacobian and an EXACT Lagrangian Hessian as a callable `hess(x, v)`, where "`v` is
ndarray with shape (m,) containing Lagrange multipliers"
([NonlinearConstraint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html)).
Its default `maxiter` is 1000. The underlying methods are Byrd, Hribar, Nocedal, "An interior point
algorithm for large-scale nonlinear programming," SIAM J. Optim. 9(4) (1999) 877-900
([DOI 10.1137/S1052623497325107](https://doi.org/10.1137/S1052623497325107)); Lalee, Nocedal,
Plantenga, "On the implementation of an algorithm for large-scale equality constrained optimization,"
SIAM J. Optim. 8(3) (1998) 682-706
([DOI 10.1137/S1052623493262993](https://doi.org/10.1137/S1052623493262993)); Conn, Gould, Toint,
"Trust Region Methods," SIAM 2000
([DOI 10.1137/1.9780898719857](https://doi.org/10.1137/1.9780898719857)).

IPOPT is the stronger and more degeneracy-tolerant option at higher N. It "implements an interior
point line search filter method that aims to find a local solution" for functions that "should be
twice continuously differentiable" ([IPOPT overview](https://coin-or.github.io/Ipopt/)). It takes the
exact Hessian (`hessian_approximation = exact`, the default) and a warm start
(`warm_start_init_point` "from a previous optimization of a related problem")
([IPOPT options](https://coin-or.github.io/Ipopt/OPTIONS.html)). The method is Wachter and Biegler,
"On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear
programming," Math. Programming 106(1) (2006) 25-57
([DOI 10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y)).

Why this ranks first for a one-shot. It is the cleanest match to our situation. It handles the
DEGENERACY that defeats the first-order center polish, because an interior-point method walks through
the interior toward the highly-active jammed optimum rather than tripping on a kink. It handles the
active set automatically, so it avoids the contact-set misidentification that guards option B. It uses
the exact sparse Hessian, so it is genuinely second-order near the optimum. Cost is not a constraint
here, so the O(N^2) constraint set and a dense-ish factorization are acceptable for a single champion.
`trust-constr` adds NO dependency; it runs once, offline, on the one champion.

Honest caveat. SLSQP is the WRONG tool for the terminal squeeze. It is Kraft's dense SQP (Kraft, "A
software package for sequential quadratic programming," DFVLR-FB 88-28, 1988), its default `maxiter` is
100, and it does NOT take a Hessian
([SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)). It builds its own
BFGS approximation, so it is only quasi-second-order and it scales badly. Use `trust-constr` or IPOPT,
not SLSQP, for the exact terminal solve.

---

## B. Contact-graph KKT / Newton on the exact tangency system. Verdict: DEEPEST precision, but fragile.

The math. Detect the active contact set (which pairs touch, which walls are tangent). Write the square
system of tangency equations (dist(c_i, c_j) = r_i + r_j for each contact, r_i = wall distance for each
wall tangency) plus the KKT stationarity conditions. Take Newton steps on that fixed system. On a
correct active set this converges quadratically to machine precision.

This is what the documented record refiners use. Nurmela and Ostergard reach high-precision
equal-circle packings in TWO stages: first they relax a soft repulsive energy and anneal it, then they
feed the result into an exact nonlinear system solved by Newton to sharpen the precision (Nurmela and
Ostergard, "Packing up to 50 equal circles in a square," Discrete Comput. Geom. 18 (1997) 111-120,
[DOI 10.1007/PL00009306](https://doi.org/10.1007/PL00009306); primary PDF is Springer-gated). The
energy stage uses E = sum over pairs of (lambda / r_ij^2)^s with a RISING exponent s that stiffens the
soft potential toward a hard contact, accepting a new equilibrium only if it improves the packing;
this is quoted verbatim from an open-access reproduction of the Nurmela-Ostergard method (Amore et
al., "Circle packing in arbitrary domains," [arXiv:2308.16523](https://arxiv.org/pdf/2308.16523)). The
terminal Newton stage is summarized as "the solution obtained by the optimization algorithm is then
used as a starting point to solve a system of nonlinear equations with the aim of improving the
accuracy of the solution" (Hifi and M'Hallah, "A Literature Review on Circle and Sphere Packing
Problems," Advances in Operations Research 2009,
[DOI 10.1155/2009/150624](https://doi.org/10.1155/2009/150624); this exact clause is secondary, so
treat it as corroborated, not verbatim primary). Graham and Lubachevsky grow rigid billiards until
jammed, a different route to the same rigid tangency configuration (Graham and Lubachevsky, "Repeated
patterns of dense packings of equal disks in a square," Electron. J. Combin. 3(1) (1996) #R16,
[PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r16/pdf)).

The fragility is real and documented. Amore et al. DELIBERATELY REJECT the contact-graph Newton route
for exactly this reason. They write that a suspected contact "is in reality a 'false contact', meaning
that two disks are extremely close without touching each other. In such case the corresponding
equation should be eliminated by the set. Unfortunately one cannot know a priori whether a contact is
'true' or 'false'" ([arXiv:2212.12287](https://arxiv.org/abs/2212.12287), Sec. 2.3). CSQV optima are
highly degenerate, with many circles at three or more contacts (contacts run near 3N,
[contacts.txt](https://www.packomania.com/csqv/txt/contacts.txt)), so the active set is both large and
ambiguous. A wrong active set gives a wrong or singular Newton system.

So build B only as the DEEPEST final step, after A has cleaned the active set, and guard it (Section D
and E). It reaches machine precision that A may leave one or two digits short, but it must not run on
an unverified contact graph.

---

## C. Newton on the REDUCED center problem (radii = LP, dual = gradient). Verdict: elegant, but the optimum is a kink.

The math. Eliminate the radii: f(C) = the exact-LP value as a function of the centers only. The
gradient is the LP dual (envelope theorem), which `center_gradient` already assembles
([polish.hpp](../../cpp/csqv/polish.hpp)). A reduced Newton step would need the HESSIAN of f, the
second derivative of the LP value function in the centers. Within one optimal basis f is smooth and its
Hessian comes from differentiating the active tangency directions; across a basis change f has a KINK.

Why this is only third. The CSQV optimum sits AT a kink, because the jammed optimum is where contacts
make and the basis is degenerate (the same kink that stalls the first-order polish). So the reduced
Hessian is defined on each side but not at the optimum itself, and a reduced Newton must carry an
active-set or basis strategy through the non-smooth point. That is the same active-set problem as
option B, minus the clean full-space structure that a general NLP solver already handles. The reduced
form is the right MENTAL model, and it is the cheap in-loop gradient we already ship, but for the
one-time exact squeeze the full-space solver A is simpler and more robust. Keep C as the analytic
first-order conditioner that feeds A, not as the exact second-order engine.

---

## D. Rigidity / jamming certification of the contact graph. Verdict: REQUIRED guard for B.

The math. Before committing a contact set to a Newton solve, certify that the contacts form a rigid
framework and reject false contacts. The linear-programming jamming test linearizes the impenetrability
constraints about the packing and solves an LP on the contact network to classify the packing as
locally, collectively, or strictly jammed (Donev, Torquato, Stillinger, Connelly, "A linear
programming algorithm to test for jamming in hard-sphere packings," J. Comput. Phys. 197 (2004)
139-166, [DOI 10.1016/j.jcp.2003.10.047](https://doi.org/10.1016/j.jcp.2003.10.047)). The review
confirms "rigorous and efficient linear-programming algorithms have been devised to assess whether a
particular sphere packing is locally, collectively, or strictly jammed" (Torquato and Stillinger,
"Jammed hard-particle packings," Rev. Mod. Phys. 82 (2010) 2633,
[arXiv:1008.2982](https://arxiv.org/abs/1008.2982)).

How it guards. Run the Donev LP test on the detected contact graph. If the graph is not rigid
(under-jammed) or a contact carries no load, the Newton system for B is singular or wrong, so drop the
false contacts and re-detect before the Newton step. This test is cheap (one small LP) and it removes
the active-set fragility the catalog flagged. The false-contact tolerance is the same tolerance Amore
et al. never fully specify; treat it as a swept parameter, not a fixed constant.

---

## E. Extended / arbitrary-precision refinement. Verdict: use only if the float floor binds.

The float noise floor is real. Our `radii_lp` is strictly feasible to about 1e-12 and the `.pck`
cushion sits at that floor. A double-precision Newton residual cannot go below about 1e-16, and
round-off near a degenerate optimum can bind above that. The precision ladder is available and each
step is a primary tool:

- long double and `__float128`. Plain C++ `long double` (80-bit on x86, about 1e-19) and GCC
  `__float128` (quad, about 1e-33) are drop-in for a Newton refinement with no library. Boost
  exposes the quad type as `boost::multiprecision::float128`, and Boost.Multiprecision "provides
  integer, rational, floating-point, complex and interval number types having more range and precision
  than the language's ordinary built-in types," with precision "arbitrarily large (memory-limited),
  fixed at compile time, or variable at runtime"
  ([Boost.Multiprecision README](https://github.com/boostorg/multiprecision);
  [docs](https://www.boost.org/doc/libs/release/libs/multiprecision/)).
- mpmath. "A Python library for arbitrary-precision floating-point arithmetic" with root-finding
  (`findroot`, arbitrary-precision Newton) and linear algebra, version 1.4.1
  ([PyPI mpmath](https://pypi.org/project/mpmath/); [docs](https://mpmath.org/doc/current/)). Useful to
  run the B Newton step at a chosen `mp.dps` below the double floor for the single champion.
- Interval arithmetic for CERTIFIED refinement. The interval Newton / Krawczyk operator certifies that
  a box "contains a unique true zero" of a nonlinear system: Krawczyk 1969 built the interval Newton
  method, and Moore 1977 showed it "can be used to certify the existence and uniqueness of a solution
  to a system of nonlinear equations" (Breiding, Rose, Timme, "Certifying zeros of polynomial systems
  using interval arithmetic," [arXiv:2011.05000](https://arxiv.org/pdf/2011.05000)). This is the tool
  that PROVED the equal-circle optima for N=28, 29, 30 (Markot and Csendes, "A new verified
  optimization technique for the packing circles in a unit square problems," SIAM J. Optim. 16(1)
  (2005) 193-219, [DOI 10.1137/S1052623403425617](https://doi.org/10.1137/S1052623403425617)) and for
  N=31, 32, 33 (Markot, "Improved interval methods for solving circle packing problems in the unit
  square," J. Global Optim. 2021,
  [DOI 10.1007/s10898-021-01086-z](https://doi.org/10.1007/s10898-021-01086-z)); both texts are gated,
  so the result claims are confirmed by index and abstract. Kearfott is the foundational interval
  Newton with generalized bisection ("Interval Newton/generalized bisection when there are
  singularities near roots," [DOI 10.1007/BF02283694](https://doi.org/10.1007/BF02283694)).

When to reach for E. Only after A and B leave a residual at the float floor and a submission needs one
or two more digits. For CSQV at N in the low hundreds, a `long double` or `__float128` Newton in B is
the cheap first rung; full interval certification is a proof tool for small N, not a record-refiner for
N=142.

---

## F. Variance-minimizing refinement (Amore). Verdict: robust derivative-free fallback.

Amore et al. avoid the contact-graph Newton and instead minimize the VARIANCE of the squared
contact distances. When the variance vanishes "the contacts between the disks are fully enforced
within the numerical precision," and they reach a contact-distance variance of about 5.6e-21 in plain
double precision with python and numba, with no extended arithmetic
([arXiv:2212.12287](https://arxiv.org/abs/2212.12287), Sec. 2.3; the 5.6e-21 is a variance of squared
distances near 0.17, not a coordinate residual). The method declares contacts by a soft threshold that
it progressively tightens, so it never commits to a wrong exact equation set.

Why it earns a slot. It is the documented answer to the false-contact problem that makes B fragile. It
is derivative-free and stochastic, so it is slower than a Newton step, but it is robust where the
contact graph is ambiguous. Use it as the FALLBACK when the Donev test (D) cannot certify the graph,
or as a cross-check that B converged to the right configuration.

---

## G. Physical / relaxation conditioner (Lubachevsky-Stillinger, repulsive energy). Verdict: not needed here.

The Lubachevsky-Stillinger billiard inflates particles under event-driven collisions until the system
jams (Lubachevsky and Stillinger, "Geometric properties of random disk packings," J. Stat. Phys. 60
(1990) 561-583, [DOI 10.1007/BF01025983](https://doi.org/10.1007/BF01025983)), and the repulsive-energy
anneal is the Nurmela-Ostergard first stage (Section B). These CONDITION a cold or messy start into a
near-jammed configuration before an exact solve. We already hold a near-optimal, first-order-stationary
champion, so the conditioning work is done. A physical relaxation would DISTURB the champion and risk a
different basin, which the terminal squeeze must not do. Keep the physical conditioner for the SEARCH
stage (cold seeds), not for the terminal post-processor. Our `center_polish_analytic` is the correct,
non-disturbing first-order conditioner for this pipeline.

---

## Ranked candidates for a one-time terminal squeeze (best first)

Judge each by whether it reaches the EXACT local optimum of the current basin, robustly, in one shot.
Cost is not a constraint.

| Rank | Algorithm | Reaches exact optimum | Robustness | Evidence |
|------|-----------|-----------------------|------------|----------|
| 1 | Full-space sparse second-order NLP (trust-constr or IPOPT, exact Lagrangian Hessian) | Yes, quadratic near the optimum; self-identifies the active set | HIGH; interior-point handles the degeneracy that stalls the first-order polish | [trust-constr](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html); [IPOPT](https://coin-or.github.io/Ipopt/); Byrd-Hribar-Nocedal 1999 [DOI](https://doi.org/10.1137/S1052623497325107); Wachter-Biegler 2006 [DOI](https://doi.org/10.1007/s10107-004-0559-y) |
| 2 | Contact-graph KKT / Newton on the exact tangency system (Nurmela-Ostergard terminal step) | Yes, to machine precision, IF the active set is correct | LOW without guards; false contacts and degeneracy break it | Nurmela-Ostergard [DOI](https://doi.org/10.1007/PL00009306); false-contact caution [arXiv:2212.12287](https://arxiv.org/abs/2212.12287) |
| 3 | Rigidity guard (Donev LP jamming test) around rank 2 | Enables rank 2 safely; rejects false contacts | Required companion to rank 2 | Donev et al. [DOI](https://doi.org/10.1016/j.jcp.2003.10.047) |
| 4 | Variance-minimizing refinement (Amore) | Yes to about 5.6e-21 contact variance; derivative-free | HIGH where the contact graph is ambiguous; slower | [arXiv:2212.12287](https://arxiv.org/abs/2212.12287) |
| 5 | Extended-precision Newton (long double, __float128, mpmath) on rank 2 or 4 | Pushes below the 1e-12 float floor | Applies only when the float floor binds | [mpmath](https://pypi.org/project/mpmath/); [Boost](https://github.com/boostorg/multiprecision); interval Newton [arXiv:2011.05000](https://arxiv.org/pdf/2011.05000) |
| 6 | Reduced-center Newton (LP value Hessian) | In principle, but the optimum is a kink | Needs the same active-set handling as rank 2, minus the full-space simplicity | Envelope/Danskin |
| 7 | Physical conditioner (Lubachevsky-Stillinger, repulsive anneal) | No; conditions a start, risks a different basin | Wrong stage for a terminal squeeze | Lubachevsky-Stillinger [DOI](https://doi.org/10.1007/BF01025983) |

Reasoning for the top pick. Build rank 1 first. It is the single tool that reaches the exact optimum
of the current basin AND handles the degeneracy that defeats the first-order polish AND identifies the
active set itself, so it avoids the contact-set fragility that makes rank 2 dangerous. It costs no new
dependency in scipy. Rank 2 with the rank 3 guard is the deepest-precision follow-up, worth adding only
if rank 1 leaves a residual and a submission needs the last digits.

---

## Recommended terminal pipeline

Run this ONCE on the final champion, after the search selects it, before submission.

1. First-order conditioner. Run `center_polish_analytic`
   ([polish.hpp](../../cpp/csqv/polish.hpp)) to reach first-order stationarity and settle the contact
   graph. This is the current baseline and the correct, non-disturbing conditioner.
2. Exact second-order squeeze. Run a full-space `trust-constr` (or IPOPT) solve on (x, y, r) with a
   SPARSE Jacobian and the EXACT sparse Lagrangian Hessian, warm-started from step 1. This closes the
   1 to 2e-5 gap by taking true second-order steps through the degenerate region.
3. Deepest polish, optional and guarded. Detect the contact graph, run the Donev LP jamming test to
   reject false contacts and confirm rigidity, then take contact-graph Newton steps on the certified
   tangency system, in `long double` or `__float128` if the double residual binds. Fall back to the
   variance-minimizing refinement if the graph fails certification.
4. Exact radii. Re-solve `radii_lp` on the final centers for the exact optimal radii
   ([lp.hpp](../../cpp/csqv/lp.hpp)).
5. Feasibility shrink and zero-tolerance validation. Apply the existing verifier and feasibility
   shrink, so the emitted `.pck` is strictly feasible at zero tolerance.

Failure modes and their guards.
- Contact-set misidentification (false contacts). Guard with the Donev rigidity test (step 3) and the
  variance cross-check (step 3 fallback). Never run the contact-graph Newton on an uncertified graph.
- Newton divergence. The full-space solver in step 2 is trust-region safeguarded, so it does not
  diverge; the step-3 Newton must carry a line search or trust region and fall back to step 2 or the
  variance refinement on a bad step.
- Post-step infeasibility. Steps 4 and 5 always restore strict feasibility: `radii_lp` gives feasible
  radii for any centers, and the shrink banks the cushion. So a step that raises a tiny overlap cannot
  reach the submission.

---

## Acceptance test

Measure the pipeline against the current first-order `center_polish_analytic` baseline, on our OWN
submitted packings for N=142 and N=143, which Packomania re-optimized after acceptance.

- Target N=142: the published best-known is 6.289851836935, ref [21] (Arnold Castro, us)
  ([csqv.html](https://www.packomania.com/csqv/csqv.html)).
- Target N=143: the published best-known is 6.312993919570, ref [21] (us)
  ([csqv.html](https://www.packomania.com/csqv/csqv.html)).

Procedure. Take our SUBMITTED (pre-re-optimization) N=142 and N=143 packings. Record the baseline sum
after `center_polish_analytic`. Run the full pipeline. PASS if the pipeline recovers the published
values to within the `.pck` float floor (about 1e-9), that is, it closes the FULL gap that Specht's
re-optimization closed. Report the fraction of the gap recovered, against the baseline that currently
recovers only 40 to 66 percent of it. A
partial pass (for example 90 percent gap closure) still beats the baseline and is worth shipping; a
full pass means we can reach the re-optimized optimum BEFORE submission and never leave the gap for
Specht to close. Validate the final `.pck` at zero tolerance as always.

---

## Where the primary evidence is thin

- Specht's actual CSQV re-optimization method is UNDOCUMENTED. The Packomania front page and the CSQV
  page describe the problem and a closed "running search" program named "csqv", but NO algorithm for
  how the maintainer re-optimizes a submission ([packomania.com](https://www.packomania.com/);
  [csqv.html](https://www.packomania.com/csqv/csqv.html)). Specht's only DOCUMENTED refinement methods
  are the Modified Billiard Simulation (Szabo and Specht 2005, in the Szabo-Csendes survey,
  [45survey.pdf](https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf)) and Pulsating Disk Shaking
  (Szabo, Markot, Csendes, Specht, Casado, Garcia, "New Approaches to Circle Packing in a Square: With
  Program Codes," Springer 2007, [DOI 10.1007/978-0-387-45676-8](https://doi.org/10.1007/978-0-387-45676-8),
  Sec. 12.2). Both are STOCHASTIC, PHYSICAL methods for EQUAL circles, not a contact-graph Newton or an
  NLP polish, and neither is stated to be the CSQV variable-radii re-optimizer. So "match Specht" is an
  EMPIRICAL target, not a method we can copy. We INFER a "deeper same-basin second-order squeeze" from
  the N=121 evidence (published = our transposed arrangement, centers nudged about 1e-6, our LP
  reproduces the published sum to 1e-12), not from any Specht disclosure. The PDS verbatim algorithm
  stayed unreachable (Google Books and the citeseerx mirror both blocked).
- The Nurmela-Ostergard terminal Newton clause is from a secondary review (Hifi and M'Hallah 2009); the
  1997 primary PDF is Springer-gated. The energy-anneal stage is confirmed verbatim from an
  open-access reproduction ([arXiv:2308.16523](https://arxiv.org/pdf/2308.16523)), not from the 1997
  primary.
- Amore et al. explicitly REJECT the contact-graph Newton route because of the false-contact problem.
  This is primary and load-bearing: it is direct evidence that rank 2 is fragile and needs the rank 3
  guard, not a safe default.
- The Markot-Csendes and Markot interval-arithmetic proofs are behind SIAM and Springer gates; the
  N=28-30 and N=31-33 optimality claims are confirmed by index and abstract, not by verbatim primary
  text. They prove EQUAL-circle optima, and no interval proof exists for the CSQV sum-of-radii
  objective, so E is a precision rung, not a CSQV stopping theorem.
- The exact head-to-head gap closure of the full pipeline versus `center_polish_analytic` on N=142 and
  N=143 is UNMEASURED. The acceptance test above is the experiment that settles it.

---

## Sources

- scipy trust-constr (Byrd-Omojokun SQP plus trust-region interior point; sparse Jacobian; exact,
  sparse, or LinearOperator Hessian via hess(x, v) with multipliers; maxiter 1000; "most appropriate
  for large-scale"): https://docs.scipy.org/doc/scipy/reference/optimize.minimize-trustconstr.html ;
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html ;
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html
- Byrd, Hribar, Nocedal, "An interior point algorithm for large-scale nonlinear programming," SIAM J.
  Optim. 9(4) (1999) 877-900: https://doi.org/10.1137/S1052623497325107
- Lalee, Nocedal, Plantenga, "On the implementation of an algorithm for large-scale equality
  constrained optimization," SIAM J. Optim. 8(3) (1998) 682-706:
  https://doi.org/10.1137/S1052623493262993
- Conn, Gould, Toint, "Trust Region Methods," SIAM 2000: https://doi.org/10.1137/1.9780898719857
- scipy SLSQP (Kraft 1988 DFVLR-FB 88-28; maxiter 100; no Hessian):
  https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html
- IPOPT (interior-point filter line search; twice-differentiable; hessian_approximation exact default;
  warm_start_init_point): https://coin-or.github.io/Ipopt/ and https://coin-or.github.io/Ipopt/OPTIONS.html
- Wachter and Biegler, "On the implementation of an interior-point filter line-search algorithm for
  large-scale nonlinear programming," Math. Programming 106(1) (2006) 25-57:
  https://doi.org/10.1007/s10107-004-0559-y
- Nurmela and Ostergard, "Packing up to 50 equal circles in a square," Discrete Comput. Geom. 18 (1997)
  111-120 (energy anneal then terminal Newton on the exact contact system; primary Springer-gated):
  https://doi.org/10.1007/PL00009306
- Amore et al., "Circle packing in arbitrary domains" (open-access reproduction of the Nurmela-Ostergard
  energy with rising exponent s): https://arxiv.org/pdf/2308.16523
- Hifi and M'Hallah, "A Literature Review on Circle and Sphere Packing Problems," Advances in
  Operations Research 2009 (terminal-Newton summary of Nurmela-Ostergard, secondary):
  https://doi.org/10.1155/2009/150624
- Amore et al., "Circle packing in regular polygons" (variance-minimizing refinement to about 5.6e-21;
  false-contact caution; deliberate rejection of contact-graph Newton): https://arxiv.org/abs/2212.12287
- Graham and Lubachevsky, "Repeated patterns of dense packings of equal disks in a square," Electron.
  J. Combin. 3(1) (1996) #R16: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r16/pdf
- Donev, Torquato, Stillinger, Connelly, "A linear programming algorithm to test for jamming in
  hard-sphere packings," J. Comput. Phys. 197 (2004) 139-166: https://doi.org/10.1016/j.jcp.2003.10.047
- Torquato and Stillinger, "Jammed hard-particle packings," Rev. Mod. Phys. 82 (2010) 2633:
  https://arxiv.org/abs/1008.2982
- Lubachevsky and Stillinger, "Geometric properties of random disk packings," J. Stat. Phys. 60 (1990)
  561-583: https://doi.org/10.1007/BF01025983
- Danskin, "The Theory of Max-Min," Springer 1967 (envelope theorem):
  https://doi.org/10.1007/978-3-642-46092-0 ; LP value-function derivatives and break points:
  https://arxiv.org/pdf/2008.02305
- mpmath (arbitrary-precision floating point, findroot, linear algebra, v1.4.1):
  https://pypi.org/project/mpmath/ ; https://mpmath.org/doc/current/
- Boost.Multiprecision (cpp_bin_float, cpp_dec_float, mpfr, float128; arbitrary/fixed/runtime
  precision): https://github.com/boostorg/multiprecision ;
  https://www.boost.org/doc/libs/release/libs/multiprecision/
- Breiding, Rose, Timme, "Certifying zeros of polynomial systems using interval arithmetic" (interval
  Newton / Krawczyk existence and uniqueness, Krawczyk 1969, Moore 1977): https://arxiv.org/pdf/2011.05000
- Markot and Csendes, "A new verified optimization technique for the packing circles in a unit square
  problems," SIAM J. Optim. 16(1) (2005) 193-219 (interval proof N=28, 29, 30; gated):
  https://doi.org/10.1137/S1052623403425617
- Markot, "Improved interval methods for solving circle packing problems in the unit square," J. Global
  Optim. (2021) (interval proof N=31, 32, 33; gated): https://doi.org/10.1007/s10898-021-01086-z
- Kearfott, "Interval Newton/generalized bisection when there are singularities near roots":
  https://doi.org/10.1007/BF02283694
- Packomania CSQV registry (N=121 published 5.799103501951; N=142 published 6.289851836935; N=143
  published 6.312993919570; N=141/142/143 ref [21] = Arnold Castro): https://www.packomania.com/csqv/csqv.html ;
  contacts near 3N: https://www.packomania.com/csqv/txt/contacts.txt ;
  Szabo-Csendes survey (Modified Billiard Simulation): https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf
- Our code: cpp/csqv/lp.hpp (`radii_lp`, `solve_reduced_lp` with duals), cpp/csqv/polish.hpp
  (`center_polish`, `center_gradient`, `center_polish_analytic`).

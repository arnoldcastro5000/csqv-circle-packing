# Phase 2 CSQV: record methods and soft-target selection

Research date: 2026-09-13. This note goes deep on METHOD and TARGET SELECTION for a
single-N CSQV record push. It does not repeat the landscape overview. Read the two
companion notes first:
[docs/research/phase-2-csqv-target-and-landscape.md](phase-2-csqv-target-and-landscape.md)
and [docs/research/phase-2-csqv-method-and-precedent.md](phase-2-csqv-method-and-precedent.md).

CSQV means variable-radius circles in a unit square. The objective maximizes the sum of
the radii. The circles must not overlap and must stay inside the square.

Method. This note uses primary sources only. It pulls the raw Packomania data files
(`author.txt`, `symmetry.txt`, `sumradii.txt`, `contacts.txt`), the CSQV main page and its
reference legend, the Hints page, sample coordinate files, the sibling paper full text
(arXiv:2609.05093), and the primary monotonic-basin-hopping packing literature. Each
quantified claim carries an inline URL citation. The prose uses short active sentences and
no em dash.

---

## 1. The single-N record pipeline

### 1a. What actually sets CSQV records now

Every CSQV record system runs the same two-layer shape. A global-search shell explores where
to put the centers. Two fixed exact primitives price and polish each candidate: a linear
program that sets the optimal radii given fixed centers, and a nonlinear local polish (SLSQP,
L-BFGS-B, or a KKT-Newton step on the contact graph). The LLM route and the classical route
agree on this shape ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093);
[Grosso, Addis, Locatelli, Schoen, optimization-online 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

The radii-by-LP primitive is the CSQV lever. Given fixed centers, the maximum feasible radii
solve a linear program: each radius is bounded by half the wall distance and by the pairwise
gaps to neighbors. This decouples the hard nonconvex geometry (center placement) from the
cheap sizing (radii). The sibling seed solver uses exactly this
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

### 1b. The sibling discovery-loop pipeline, iteration by iteration

The sibling paper is the direct precedent: "LLM-Guided Program Evolution for Circle Packing:
Breaking 10 Packomania Records for $28"
([arXiv:2609.05093 abstract](https://arxiv.org/abs/2609.05093);
[full text](https://arxiv.org/html/2609.05093)). It evolves a complete replacement solver
each iteration, no crossover between providers, single model (Claude Fable 5.1 via the
Anthropic CLI). The 15 iterations (index 0 to 14), verbatim labels
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)):

- Iter 0 (seed): "Penalty L-BFGS-B + LP radii".
- Iter 1: "Basin hopping + SLSQP contact polish".
- Iter 2: "Sparse penalty + adaptive operators + elite pool".
- Iter 3: failed code generation (no code).
- Iter 4: "Hexagonal lattice init + crossover".
- Iter 5: "Island-model parallel basin hopping".
- Iter 6: "Ruin-and-recreate LNS" (rejected).
- Iter 7: "Affine lattice template bank".
- Iter 8: "KKT-Newton exact polish".
- Iter 9: "Lattice-aware slip moves".
- Iter 10: "Surrogate-screened basin hopping" (rejected).
- Iter 11: "Sparse KKT polish" (rejected).
- Iter 12: "Defect-migration moves".
- Iter 13: "Formulation-space search" (rejected).
- Iter 14: "Row-template cold starts" (rejected).

Seeding strategies the sibling used: random init (seed), hexagonal-lattice templates "sized
to hold exactly N circles" (iter 4), an affine lattice template bank (iter 7), and
row-template cold starts (iter 14) ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

The winning building blocks are the accepted iterations: basin hopping with SLSQP contact
polish, an elite pool, island-model parallel basin hopping, a KKT-Newton exact polish on the
contact graph, and defect-migration moves that drop weak circles and re-equilibrate. The
lattice and surrogate iterations were mostly rejected
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

### 1c. The classical global-optimization pipeline (Addis, Locatelli, Schoen)

The classical record method for packing unequal circles is Monotonic Basin Hopping (MBH) and
its population variant Population Basin Hopping (PBH). Exact mechanics
([Grosso, Addis, Locatelli, Schoen, optimization-online 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)):

- Local search: they use the NLP solver SNOPT for the local step. Any local method works,
  but "past experience suggests that SNOPT is" the effective choice.
- Perturbation move: "a uniform random perturbation of each coordinate of each" circle by a
  step delta. The rule is to "choose it in such a way that the structure of the current local
  minimizer is not completely disrupted", so the method walks between "close" local minimizers.
  If delta is too small, the perturbed point falls back into the same basin.
- Acceptance: monotonic. MBH keeps the incumbent, perturbs it, runs a local search, and
  accepts only if the new local minimum is better.
- Landscape conjecture: the packing problem has a "funneling landscape", the same feature
  found in molecular-conformation problems. This is why perturb-and-polish beats plain
  multistart, whose local-minima count "increases quite quickly".
- PBH: keeps a population of N local minimizers. It perturbs each member, then a child either
  competes with the most-similar population member (dissimilarity `d` below a cut `dcut`) or
  with the worst member (if dissimilar to all). A dissimilarity measure `d` over circle
  distances from the barycenter, plus `dcut`, keeps the population diverse.
- Stop rule: stop if the best member does not change for a fixed number of iterations.
- Result strength: the approach "could improve over previously known putative optima in the
  range n <= 130 in as many as 32 instances".

Physics-inspired seeds (Lubachevsky-Graham billiard, Nurmela-Ostergard molecular repulsion)
generate good equal-circle starts and feed the polish, but they target equal radii, so they
serve CSQV only as an initializer ([Szabo, Csendes survey](https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf)).

### 1d. Pipeline synthesis for our per-N push

The robust single-N pipeline is MBH or PBH with a fixed LP-radii pricing step and a fixed
SLSQP polish, plus defect migration to escape stuck contact graphs. The initialization
pattern, the perturbation step and subset, the acceptance schedule, and any island structure
are the parts worth evolving. Keep the LP-radii and the local polish out of the genome; they
are exact given centers ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093);
[optimization-online 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

---

## 2. Symmetry exploitation

The finding is decisive and it points AGAINST a symmetry constraint for most target N.

`symmetry.txt` records the symmetry group of each best-known CSQV packing
([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt)). The pattern:

- N=1 to 12: high symmetry (D4, D2, D1) because the packings are near-trivial and forced.
- N=13 and up: almost every best-known record is C1, meaning NO symmetry. The record
  configuration is asymmetric.
- The only symmetric records above N=12 are scattered special N: D1 at 14, 27, 34, 38, 51,
  64, 100, 196, 256, 324, and the alternating D1/D2 series on the large grid-like N; D2 at
  21, 25, 36, 49, 81, 121, 225, 289, 361, 441, and the perfect-square series.

Read against the proven-optimal (bold) set, the symmetric entries are mostly the proven or
grid-construction cases. The whole live frontier (the dense N=13 to 124 band, all of N=101 to
125, and the between-square high N) is C1 ([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt);
[csqv.html reference legend](https://www.packomania.com/csqv/csqv.html)).

Implication for the genome. A symmetry constraint (D1 mirror, D2, D4) cuts the free variables
roughly in half (D1) to a quarter (D4), which shrinks the NLP and speeds each solve. But the
data says the true CSQV optima at record N are asymmetric, so a hard symmetry constraint
would cap the reachable sum below the record. Do NOT hard-constrain to a symmetry group for a
generic C1 target. Two safe uses of symmetry remain: use a symmetric configuration only as
one SEED among many (cheap warm start, then relax to full C1 in the polish), and consider a
symmetry constraint only if the specific target N already holds a symmetric record (rare
above N=12) or is a perfect-square grid case ([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt)).

---

## 3. Soft-target identification

### 3a. The metadata Packomania exposes

Per-N data files under `csqv/txt/` ([directory index](https://www.packomania.com/csqv/)):

- `author.txt`: the current record holder per N. This is the staleness lever
  ([author.txt](https://www.packomania.com/csqv/txt/author.txt)).
- `symmetry.txt`: the symmetry group per N (Section 2)
  ([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt)).
- `sumradii.txt`: the current best sum of radii per N, 12 decimals
  ([sumradii.txt](https://www.packomania.com/csqv/txt/sumradii.txt)).
- `contacts.txt`: the contact count per N. It grows about linearly, roughly 3N contacts,
  a sanity check on how jammed a record is
  ([contacts.txt](https://www.packomania.com/csqv/txt/contacts.txt)).
- The CSQV main page carries the numbered Reference legend that maps each contributor to a
  DATE, plus the bold (proven-optimal) marks and the History of updates
  ([csqv.html](https://www.packomania.com/csqv/csqv.html)).
- Caveat 1: `records.txt` currently returns a server-side parser error, not data
  ([records.txt](https://www.packomania.com/csqv/txt/records.txt)).
- Caveat 2: the per-N HTML page date (for example "03-Sep-2026" on the N=30 page) is a site
  regeneration date, not the record date. Do not read it as staleness
  ([csqv30.html](https://www.packomania.com/csqv/csqv30.html)). Use the Reference legend date
  and the author identity instead.

### 3b. The reference legend dates each contributor

The legend maps contributor to era ([csqv.html](https://www.packomania.com/csqv/csqv.html)):

- [1] David W. Cantrell, sci.math forum, 2011/12. This is the pre-LLM classical wave. "About
  2011: first complete presentation from N=1 to N=100 due to Cantrell [1]."
- [2] E. Specht, program csqv, 2011 to 2026. Grid-style constructions, mostly the high N.
- [3] E. Friedman, Circles in Squares site. [4] A. Novikov et al., AlphaEvolve, Jun 2025.
- [5] to [20]: the 2026 LLM and human wave. [8] Haowei Lin, mid-Jul 2026. [10] Everett
  Dutton, Jul/Aug 2026. [11] Jason Liang, Aug 2026. [13] Byron Tasseff, Aug 2026. [14] Wes
  Sander (the sibling discovery-loop), Sep 2026. [15] Yue Huang, Sep 2026. [18] Anant Garg,
  Sep 2026. [20] Zeeshan Tariq, Sep 2026.

Any entry still credited to [1] Cantrell is OLD (2011/12) and predates the entire 2026 wave.
Any entry credited to [1] AND not bold (not proven optimal) is a soft target: old author, not
proven, and skipped by the wave that swept the neighboring N.

### 3c. The proven-optimal (bold) lock

Bold sum-of-radii entries (proven optimal, unbeatable) from the main page
([csqv.html](https://www.packomania.com/csqv/csqv.html)):
N = 1 to 12, 14, 16, 18, 25, 36, 49, 64, 81, 100, 121. Note N=27 is NOT bold; Jason Liang
beat Cantrell there with a D1 configuration ([author.txt](https://www.packomania.com/csqv/txt/author.txt)).

### 3d. Concrete procedure to find a soft N

1. Pull `author.txt`. Select N where the holder is "David W. Cantrell" (reference [1], 2011/12).
2. Remove the bold N (proven optimal): 1 to 12, 14, 16, 18, 25, 36, 49, 64, 81, 100, 121.
3. Cross-check `symmetry.txt`. Prefer C1 targets for a generic push, or a symmetric target
   only if you will constrain to that group.
4. Confirm the neighbors were swept by the 2026 wave but the target was skipped. This shows
   the wave found the target hard or uninteresting, not that it is optimal.
5. Re-pull `author.txt` at attempt time. The frontier moves weekly, so a soft N can flip.
6. Size the compute against N (Section 4). Prefer low N for a pure-Python plus scipy box.

### 3e. Candidate soft N

Applying the procedure to the current `author.txt`
([author.txt](https://www.packomania.com/csqv/txt/author.txt)),
`symmetry.txt` ([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt)), and the
bold set ([csqv.html](https://www.packomania.com/csqv/csqv.html)):

Mid-band Cantrell holdouts, not proven, C1 symmetry, skipped by the 2026 wave that took every
neighbor: N = 28, 29, 30, 31, 37, 41, 43, 45, 60, 90. Add N = 38 (Cantrell, not bold, but D1
symmetry). These are the strongest soft targets. They are old (2011/12), single classical
contributor, not proven optimal, uncontested through 2026-09-12, and low N so cheap to solve.
The surrounding N (33, 34, 35, 39, 40, 42, 47, 48, 50 and up) all changed hands in the 2026
wave, which frames these holdouts as overlooked, not optimal
([author.txt](https://www.packomania.com/csqv/txt/author.txt);
[csqv.html History of updates](https://www.packomania.com/csqv/csqv.html)).

Small Cantrell holdouts (N = 13, 15, 17, 19, 20, 21, 22, 23) are also old and not bold, but
small N is heavily studied and likely near-optimal, so headroom is thin. Treat them as lower
priority.

Sparse high band (N > 250): the perfect-square N (289, 324, 361, 400, 441, 484, 529, 576,
625, and up) are held by [2] E. Specht as grid-style constructions. These are single-author,
not proven, and plausibly soft, since a naive grid leaves sum-of-radii headroom for variable
radii. But N is large, so O(N^2) constraints make each solve heavy for a pure-Python box (see
Section 4). The between-square high N (265, 313, 365, 421, 481, 545, 613, 685, 761, 925) are
recent Tariq [20] entries, so they are contested, not soft
([author.txt](https://www.packomania.com/csqv/txt/author.txt);
[symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt)).

Avoid, per the brief and the data: N = 101 to 114 (the sibling and Tariq/Huang see-saw; the
sibling [14] currently holds 101, 103, 105, 106, 107, 108, 109, 111, while Tariq [20] and
Huang [15] retook 102, 104, 110, 112, 113, 114), and the hot 100 to 250 race generally
([author.txt](https://www.packomania.com/csqv/txt/author.txt)).

Headroom read. `sumradii.txt` shows a smooth, near-linear growth of the best sum with N
(for example N=30 is 2.842668747462, N=90 is 4.996999276706)
([sumradii.txt](https://www.packomania.com/csqv/txt/sumradii.txt)). A target that sits below
the local trend of its neighbors, or that carries fewer contacts than neighbors of similar N
in `contacts.txt`, is a jamming-headroom signal worth a closer look
([contacts.txt](https://www.packomania.com/csqv/txt/contacts.txt)).

---

## 4. Compute reality for a single N

Cost of one evaluation. A CSQV local solve has 3N decision variables (x, y, r per circle) and
O(N^2) pairwise non-overlap constraints plus 4N wall constraints. The constraint count and
the Jacobian grow like N^2, and a dense SLSQP step grows worse than N^2 per iteration. The
LP-radii step is cheap and near-linear given fixed centers
([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091);
[arXiv:1607.02184](https://arxiv.org/pdf/1607.02184)).

Reference compute from the sibling. The full loop reached N up to 114 in pure Python with
scipy, at a 120-second per-target timeout, with 6 parallel workers, on an Intel Core
i7-13700KF with 32 GB RAM, in "approximately 8 hours overnight" for 15 iterations. The solver
is about 400 lines ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). That 8
hours covers the whole evolutionary loop across all targets and LLM iterations, not one N.

Reference compute from the classical route. MBH and PBH with SNOPT improved 32 putative
optima for n <= 130. The method runs many perturb-and-polish restarts per N; the count scales
with the funnel depth, not a fixed budget
([optimization-online 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

Where a 4-core, 7.8 GB pure-Python plus scipy box is competitive. The soft band N = 28 to 90
is well inside reach. At N=30 the constraint set is about 30^2 = 900 pairs; at N=90 it is
about 8100 pairs. Both fit in far under 7.8 GB (the N=114 case is about 1.3 x 10^4 pairs and
the sibling ran it on a desktop). A single SLSQP polish at N=30 is sub-second and at N=90 is a
few seconds, so a 4-worker MBH or PBH run can afford hundreds to thousands of restarts inside
a few CPU-hours. The 120-second-per-restart, 4-worker style budget is realistic here.

Wall-clock budget to propose per N:

- N=28 to 45: about 1 to 2 CPU-hours of MBH restarts, roughly 30 to 60 minutes wall on 4
  cores. Comfortable.
- N=46 to 90: about 2 to 6 CPU-hours, roughly 1 to 2 hours wall on 4 cores. Comfortable.
- N=91 to 130: feasible but slower per restart; budget 4 to 12 CPU-hours and expect fewer
  restarts inside the window. This is also the contested band, so weigh the crowding.
- N > 250 (Specht grid constructions): each solve is heavy (O(N^2) dense linear algebra) and
  pure Python is not competitive against a C++ or GPU solver. Do not target this band on the
  4-core box unless you exploit a symmetry or grid structure to shrink the solve
  ([arXiv:2404.03091](https://arxiv.org/pdf/2404.03091);
  [arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

Bottom line. The 4-core pure-Python plus scipy box is genuinely competitive in the N=28 to 90
soft band, which is exactly where the Cantrell holdouts sit. Memory is not the limit; the
per-restart solve time and the restart count are.

---

## 5. Submission and verification

Format for a credible submission to Eckard Specht ([hints.html](https://www.packomania.com/hints.html)):

- File extension `.pck`, named after the subdirectory prefix (CSQV files are `csqvN.pck`).
- Line 1: the radius of the LARGEST object. No letters, no equal sign.
- Line 2: the author name(s), comma separated for multiple authors.
- Following lines: `x y radius`, columns separated by white space or tabs.
- Sort all coordinate lines by INCREASING radius. The first line is the smallest circle, the
  last line is the largest.
- Provide "as many decimal places as possible".
- Rescale the configuration to the standard container: square of side 1, centered at the
  origin ([hints.html](https://www.packomania.com/hints.html)).

The published coordinate files confirm the exact layout. The N=114 file has a header, then
one line per circle with index, x, y, radius at 12 decimals, then a trailer with the sum
([csqv114.txt](https://www.packomania.com/csqv/txt/csqv114.txt)). Verbatim first data lines:

```
   1   0.016278166455  -0.467106148930   0.032893851070
   2   0.386933093095   0.385019052197   0.032923559086
   3   0.359907669527  -0.467004518472   0.032995481528
```

Trailer:

```
# -----------------------------------------------------
#                           sumradii =   5.628706112373
```

The published files are centered at the origin, so the container runs from -0.5 to 0.5 on
each axis, and the lines are sorted by increasing radius, matching the Hints rules
([csqv114.txt](https://www.packomania.com/csqv/txt/csqv114.txt);
[csqv30.txt](https://www.packomania.com/csqv/txt/csqv30.txt)).

Feasibility and precision. The Hints page states no numeric tolerance, but the credibility
check is implicit in the 12-decimal coordinates: any reader recomputes containment (every
circle inside the square), pairwise non-overlap (center distance >= sum of the two radii), and
the sum of radii, all from the coordinates ([hints.html](https://www.packomania.com/hints.html);
[csqv114.txt](https://www.packomania.com/csqv/txt/csqv114.txt)). The sibling verifier sets the
practical bar: it "(1) Confirms all circles are within [0,1]^2 (with zero tolerance) (2)
Confirms no pair of circles overlaps (with zero tolerance) (3) Recomputes the sum of radii
independently (4) Applies a strict feasibility shrink to eliminate numerical-precision
artifacts" ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). Adopt the same
zero-tolerance check plus a feasibility shrink before submission, and report the symmetry
group so it matches the `symmetry.txt` convention (C1, D1, D2, D4)
([symmetry.txt](https://www.packomania.com/csqv/txt/symmetry.txt)).

---

## Design questions this raises for our per-N push

- Target N: do we commit to a low mid-band Cantrell holdout (28, 29, 30, 31, 37, 41, 43, 45,
  60, 90, or the D1 case 38) where the author is old, the entry is not proven, the wave
  skipped it, and the compute is cheap? Or do we chase a Specht high-N grid construction (>
  250) that is soft but compute-heavy? The 4-core box argues strongly for the low mid-band.
- Symmetry constraint: keep it OFF for the C1 targets, since the true record configs are
  asymmetric. Use a symmetric configuration only as one seed, and relax to C1 in the polish.
  Reconsider a hard constraint only for a target whose current record is already symmetric.
- Seed strategy: which initializers do we include? At least random plus a hexagonal or grid
  template sized to N, plus a warm start from the current published record coordinates (start
  from the incumbent and perturb, MBH style).
- Budget and multi-start: how many MBH or PBH restarts per N, at what per-restart timeout, on
  4 workers? Propose 120 seconds per restart and a total wall-clock of 30 to 120 minutes for
  N in 28 to 90. Do we add defect-migration and a KKT contact polish, the two accepted
  sibling moves, or keep the solver minimal first?
- Selection metric: score on the exact verifier (zero-tolerance containment and non-overlap,
  independently recomputed sum, feasibility shrink). Do we also track contact count against
  `contacts.txt` as a jamming-progress signal?
- Submission workflow: who runs the final zero-tolerance check and the feasibility shrink, and
  do we emit the `.pck` file directly (line 1 largest radius, line 2 author, sorted by
  increasing radius, 12+ decimals, centered unit square) so a record is one step from Specht?
- Staleness guard: do we re-pull `author.txt` and the reference legend at attempt time, since
  the frontier moves weekly and a soft N can flip before we submit?

---

## Sources

- CSQV main page, reference legend, bold proven-optimal marks, History of updates:
  https://www.packomania.com/csqv/csqv.html
- CSQV data directory index: https://www.packomania.com/csqv/
- author.txt (current holder per N): https://www.packomania.com/csqv/txt/author.txt
- symmetry.txt (symmetry group per N): https://www.packomania.com/csqv/txt/symmetry.txt
- sumradii.txt (best sum per N): https://www.packomania.com/csqv/txt/sumradii.txt
- contacts.txt (contact count per N): https://www.packomania.com/csqv/txt/contacts.txt
- records.txt (currently returns a server-side parser error):
  https://www.packomania.com/csqv/txt/records.txt
- Hints and submission format: https://www.packomania.com/hints.html
- Sample coordinate files: https://www.packomania.com/csqv/txt/csqv114.txt and
  https://www.packomania.com/csqv/txt/csqv30.txt
- Per-N page (date caveat): https://www.packomania.com/csqv/csqv30.html
- Sibling discovery-loop paper: https://arxiv.org/abs/2609.05093 and
  https://arxiv.org/html/2609.05093
- Grosso, Addis, Locatelli, Schoen, packing equal and unequal circles, MBH and PBH:
  https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf and
  https://optimization-online.org/2008/06/1999/
- Addis, Locatelli, Schoen, Efficiently packing unequal disks in a circle:
  https://www.semanticscholar.org/paper/Efficiently-packing-unequal-disks-in-a-circle-Addis-Locatelli/2ce59b910ded2f50f1798cb0171c57e4f4c06f95
- Szabo, Csendes survey, classic packing methods (billiard, molecular, TAMSASS):
  https://www.inf.u-szeged.hu/~pszabo/Pub/45survey.pdf
- LP for radii given fixed centers, arXiv:1607.02184: https://arxiv.org/pdf/1607.02184
- Nonconvexity and constraint scaling, arXiv:2404.03091: https://arxiv.org/pdf/2404.03091

## Open questions

- The per-N record DATE is not in a clean data file. The Reference legend dates the
  contributor, and `author.txt` names the holder, but a precise per-N last-change date needs
  the main table cells re-read from raw HTML. `records.txt` would help but currently errors
  ([records.txt](https://www.packomania.com/csqv/txt/records.txt)).
- The Cantrell holdout list is a snapshot from `author.txt` on 2026-09-13. Re-pull before any
  attempt ([author.txt](https://www.packomania.com/csqv/txt/author.txt)).
- The exact restart count and per-N wall-clock for a real 4-core pure-Python run is an
  estimate from the constraint scaling and the sibling desktop numbers. A local benchmark at
  N=30 and N=90 would size it precisely.
</content>
</invoke>

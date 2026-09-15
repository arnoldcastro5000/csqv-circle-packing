# Phase-2 target study: CSQV (variable circles in a square, maximize sum of radii)

Research date: 2026-09-12. Author context: this note evaluates CSQV as a target for the
LLM-as-mutation discovery loop. It uses primary sources only: the Packomania site, its
Hints page, coordinate files, the CSQV update log, and the primary paper of the sibling
project.

Repo context read first: `.scratch/crucible-loop/map.md` "Out of scope" section flags
Packomania as an open frontier that is "now partly crowded by the sibling project
`ucsandman/discovery-loop`, arXiv:2609.05093". This note tests that flag against the live
site. The finding is stronger than the note: the CSQV frontier is now HOT and crowded, not
just partly crowded (see Section 4).

## 1. Objective, constraints, and verification

CSQV means "circles in a square, variable radii". The objective is to pack N circles of
unequal radii in a unit square and to MAXIMIZE the sum of the radii. The container is a
square of side 1. The circles must not overlap and must stay inside the square.
[csqv.html](https://www.packomania.com/csqv/csqv.html),
[packomania.com main page](https://www.packomania.com/).

The site reports the objective as a single "sum of radii" value per N. The site publishes
per-N data files: `txt/sumradii.txt` (sum of radii), `txt/density.txt`, `txt/contacts.txt`,
`txt/boundary.txt`, `txt/core.txt`, `txt/symmetry.txt`, `txt/author.txt`, and
`txt/records.txt`.
[csqv directory index](https://www.packomania.com/csqv/).

Verification uses the coordinate files. Each packing has an ASCII coordinate file. A
sampled file (N=114) has one line per circle with four columns: index, x, y, and individual
radius, at 12 decimal places. The coordinates are centered at the origin, so the container
is the square from -0.5 to 0.5 on each axis. The file ends with the sum of all radii.
[csqv114.txt coordinates](https://www.packomania.com/csqv/txt/csqv114.txt).
First data lines, verbatim:

```
   1   0.016278166455  -0.467106148930   0.032893851070
   2   0.386933093095   0.385019052197   0.032923559086
   3   0.359907669527  -0.467004518472   0.032995481528
```

The submission format is defined on the Hints page. Line 1 gives the radius of the largest
object. Line 2 gives the author name(s). The following lines give x, y, and the individual
radius. The rule is "Always provide as many decimal places as possible". The rule is also
"Sort all coordinate lines by increasing radii (only for unequal objects)". Objects are
rescaled to a standard container (square side 1, centered at the origin). Columns use white
space or tabs. Files use the `.pck` extension named after the subdirectory.
[hints.html](https://www.packomania.com/hints.html).

The Hints page does not spell out an automated non-overlap or containment check. The check
is implicit in the coordinate files: any reader can recompute containment (all circles
inside the square), pairwise non-overlap (center distance greater than or equal to the sum
of the two radii), and the sum of radii, all from the 12-decimal coordinates.
[hints.html](https://www.packomania.com/hints.html),
[csqv114.txt](https://www.packomania.com/csqv/txt/csqv114.txt).

Independent-verifier note for design reuse: the sibling project applies exactly this check.
Its verifier "(1) Confirms all circles are within [0,1]^2 (with zero tolerance) (2) Confirms
no pair of circles overlaps (with zero tolerance) (3) Recomputes the sum of radii
independently (4) Applies a strict feasibility shrink to eliminate numerical-precision
artifacts".
[arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093).

## 2. N coverage: locked vs open

The CSQV table covers N from 1 to 10000, with individual pages `csqv1.html` through
`csqv10000.html`. The site states the table is complete up to N=100 and sparser above.
[csqv.html](https://www.packomania.com/csqv/csqv.html),
[csqv directory index](https://www.packomania.com/csqv/).

Proven-optimal entries appear in bold in the sum-of-radii column. The site marks these as
proven optimal, so they cannot be beaten. The reported bold entries include the small cases
N=1 through N=12, plus N=27 and N=100, and selected perfect-square configurations such as
N=196, N=256, and N=324.
[csqv.html](https://www.packomania.com/csqv/csqv.html).

Density of proven-optimal vs best-known: proven-optimal entries are SPARSE. The lock covers
the small N (about N=1 to 12) plus scattered special N. The large majority of the table,
including nearly all of the dense N=13 to 100 band and everything above N=100, is
best-known only, so it is beatable. This matches the site framing that "for certain values
of N several distinct optimal configurations exist" and that most entries are current best
records, not proofs.
[csqv.html](https://www.packomania.com/csqv/csqv.html).

## 3. Staleness and winnable niches

The site keeps a "History of updates" section on the CSQV page. It shows heavy recent
activity, not staleness. The CSQV section reports last update 2026-09-09 on the main page,
and per-N pages carry a 12-Sep-2026 date.
[packomania.com main page](https://www.packomania.com/),
[csqv114.html](https://www.packomania.com/csqv/csqv114.html).

Update-log timeline (from the CSQV "History of updates"):

- 2026-07-23: "Revival of the old web page as of 2012. In the last time, new approaches
  driven by LLM yield better results".
- 2026-07-24: N=26 credited to Yiping Wang.
- 2026-07-27: N=26 improved from 2.635977394754 to 2.635983084918.
- 2026-08-01: Haowei Lin credited for N=26, 33, 34, 35, 39, 40; Everett Dutton for N=42, 47,
  48, 53 to 59, 61, 69, 73 to 75, 77, 78, 80, 83, 86 to 88, 91, 95, 96, 99.
- 2026-08-03: Jason Liang beat David Cantrell's N=27 record with a D1-symmetry configuration.
- 2026-08-04: Dutton submitted many candidates across N=48 to 98.
- 2026-08-08: New symmetric packings with full D4 symmetry for N=113, 145, 181, 221.
- 2026-08-11: Jason Liang earned his own color in the table with 21 new records.
- 2026-09-12: Zeeshan Tariq contributed extensive improvements across many N.

[csqv.html History of updates](https://www.packomania.com/csqv/csqv.html).

Most recently changed N (per the log): 102, 104, 110, 112, 114 to 120, 125, 135, 140, 144,
145, 150, 155, 169, 170, 221, 265, 313, 365, 421, 481, 545, 613, 685, 761, 925.
[csqv.html](https://www.packomania.com/csqv/csqv.html).

Reading for winnable niches: the actively improved band is N in the low hundreds (about 100
to 250) and a scattered set of higher N. These N are contested, so a record there is short
lived. The complete but less-touched band is the middle N (about 13 to 100) where the small
cases are locked but the mid cases are best-known only. The sparse high-N band (above about
250, up toward 10000) has fewer entries and fewer contributors, so it may hold softer
targets, but the site data thins out there. See Section 6 for the synthesis.
[csqv.html](https://www.packomania.com/csqv/csqv.html),
[csqv directory index](https://www.packomania.com/csqv/).

Caveat: the per-N last-change date and contributor reference number sit in the Reference
column and the `txt/author.txt` file. An automated fetch reads the main table and the update
log, but did not enumerate every per-N reference number. Treat the per-N staleness map as
partial (see Open questions).
[csqv.html](https://www.packomania.com/csqv/csqv.html).

## 4. Competitive landscape

The CSQV frontier is crowded and moving fast, driven mainly by LLM optimizers since July
2026. Specht himself notes on 2026-07-23 that "new approaches driven by LLM yield better
results".
[csqv.html](https://www.packomania.com/csqv/csqv.html).

Named record setters, July to September 2026: Yiping Wang, Haowei Lin, Everett Dutton,
Jason Liang, Zeeshan Tariq, and the classical contributor David Cantrell (whose N=27 record
was beaten). Jason Liang alone landed 21 new records in about a week (2026-08-11). Records
land at a rate of many per week across this window.
[csqv.html](https://www.packomania.com/csqv/csqv.html).

The sibling project is the discovery loop by Wes Sander, code at
`github.com/ucsandman/discovery-loop`, paper arXiv:2609.05093. It broke 10 CSQV records in
the range N=101 to 114, with gains of 2.4% to 5.4% over the prior best, at a total LLM cost
of 27.72 US dollars across 15 iterations. It used a single model, "Claude Fable 5.1
(Anthropic), accessed via the Claude CLI". Its method is program evolution from a seed
solver. The seed was "Penalty L-BFGS-B + LP radii". Later iterations added "Basin hopping +
SLSQP contact polish", "Hexagonal lattice init + crossover", and "Island-model parallel
basin hopping".
[arXiv:2609.05093 abstract](https://arxiv.org/abs/2609.05093),
[arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093).

Exact N and gains from the sibling paper (prior record -> their best, gain):

- N=101: 5.163845 -> 5.289154, +2.43%
- N=102: 5.055187 -> 5.318238, +5.20%
- N=103: 5.085509 -> 5.345481, +5.11%
- N=105: 5.125967 -> 5.401298, +5.37%
- N=106: 5.151736 -> 5.429079, +5.38%
- N=107: 5.180124 -> 5.453952, +5.29%
- N=108: 5.205806 -> 5.481819, +5.30%
- N=109: 5.231096 -> 5.507926, +5.29%
- N=111: 5.278427 -> 5.554909, +5.24%
- N=114: 5.336683 -> 5.624188, +5.39%

[arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093).

Method notes: the sibling paper cites AlphaEvolve as the foundational LLM-guided
program-evolution work. It does NOT mention Gurobi in the fetched text, and it does not name
other solver authors beyond crediting Specht for the database. The task prompt mentions
"AlphaEvolve, Gurobi-based, and individual LLM optimizers" as CSQV contributors; the primary
CSQV log confirms the individual optimizers by name (Wang, Lin, Dutton, Liang, Tariq) and
the LLM-driven wave, but an explicit AlphaEvolve or Gurobi credit inside the Packomania CSQV
log was not confirmed by the fetch (see Open questions).
[csqv.html](https://www.packomania.com/csqv/csqv.html),
[arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093).

Crowding verdict: HIGH in the N=100 to 250 band. Multiple LLM optimizers now push this band
weekly, and the sibling loop already holds 10 records in N=101 to 114 that our loop would
have to re-beat. This is the most crowded slice. The frontier is wide (N up to 10000), so
crowding drops at high N, but so does the published baseline density.
[csqv.html](https://www.packomania.com/csqv/csqv.html),
[arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093).

## 5. Submission and recognition

A new record is submitted to Eckard Specht (the maintainer) as a coordinate file. The file
is a `.pck` file. Line 1 is the largest radius. Line 2 is the author name(s). The remaining
lines are x, y, and individual radius, sorted by increasing radius, at maximum decimal
places, rescaled to the unit square centered at the origin. Coordinates are REQUIRED, since
the record is only credible from the coordinates.
[hints.html](https://www.packomania.com/hints.html),
[csqv114.txt](https://www.packomania.com/csqv/txt/csqv114.txt).

Specht updates the site actively. The CSQV section shows a dense update log across July to
September 2026, with the main page dated 2026-09-09 and per-N pages dated 12-Sep-2026.
Recognition is public: a contributor is named in the update log and the Reference column,
and a prolific contributor can earn a dedicated color in the table (Jason Liang, 2026-08-11).
[packomania.com main page](https://www.packomania.com/),
[csqv.html](https://www.packomania.com/csqv/csqv.html).

## Where the winnable niches are

The crowded slice is N=100 to about 250. Several LLM optimizers, including the sibling
discovery loop, already contest it, and the sibling holds N=101 to 114. Do not open there. A
record there is fragile and needs to out-run active competitors. The locked slice (small N
plus scattered proven-optimal N such as 27 and 100) is unwinnable by definition. The best
remaining niches are two. First, the mid band N=13 to about 100 that is complete but is not
the current LLM battleground, where a strong local-search-plus-restart method may still lift
best-known values that pre-date the 2026 LLM wave. Second, the sparse high band above about
250, where fewer contributors compete and some entries look softer, at the cost of a heavier
harness and thinner published baselines. The realistic read: CSQV is a live, credible record
book, but it is now a fast, crowded race dominated by LLM optimizers since July 2026, so any
target N needs a fresh staleness check at attempt time and a method stronger than a bare
construction heuristic (contact polish plus basin hopping plus multi-start, per the sibling
recipe).
[csqv.html](https://www.packomania.com/csqv/csqv.html),
[arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093).

## Sources

- Packomania CSQV page and update log:
  https://www.packomania.com/csqv/csqv.html
- Packomania main page (section dates, links):
  https://www.packomania.com/
- Packomania Hints / submission format:
  https://www.packomania.com/hints.html
- Packomania CSQV directory index (data files, per-N pages):
  https://www.packomania.com/csqv/
- Sampled coordinate file, N=114:
  https://www.packomania.com/csqv/txt/csqv114.txt
- Per-N page, N=114 (date, coordinate link):
  https://www.packomania.com/csqv/csqv114.html
- Sibling project paper, abstract:
  https://arxiv.org/abs/2609.05093
- Sibling project paper, full text (records table, models, verifier):
  https://arxiv.org/html/2609.05093
- Sibling project code (referenced in paper): github.com/ucsandman/discovery-loop

## Open questions

- Exact proven-optimal set: the fetch confirmed bold (proven) entries at N=1 to 12, 27, 100,
  and some perfect squares (196, 256, 324), but did not enumerate the full bold set. Read
  `txt/sumradii.txt` or the full table to list every locked N.
  https://www.packomania.com/csqv/
- Per-N staleness map: the Reference column and `txt/author.txt` hold the per-N contributor
  reference number and date. The fetch read the update log and main table, not every per-N
  reference. Pull `txt/author.txt` and `txt/records.txt` to build a precise last-change map.
  https://www.packomania.com/csqv/
- AlphaEvolve and Gurobi on CSQV: the task named these as contributors, but the fetched CSQV
  log and sibling paper did not confirm a named AlphaEvolve or Gurobi CSQV record. Confirm
  from the CSQV Reference legend (entries [1] to [20]).
  https://www.packomania.com/csqv/csqv.html
- Automated fetch was rendered by a summarizing model, not raw HTML. Verbatim table cells
  (per-N bold flags, reference numbers) should be re-read from raw source before any record
  attempt.

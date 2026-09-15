# Phase 2 CSQV: wider-band target scan (N=115 to 250 and N>250)

## Scope

This note scans packomania CSQV (circles in a square, variable radii, maximize
the sum of radii of N non-overlapping circles in the unit square) for SOFT,
un-contested, warm-startable target N in a band wider than our first scan.
Source table and update log:
https://www.packomania.com/csqv/csqv.html

Prior scan covered N in {13, 20, 26, 32, 40, 50, 60, 75, 90}. It found SOFT
(pre-2026 Cantrell records) at 13, 20, 60, 90 and FRESH 2026-wave records at
26, 32, 40, 50, 75. The sibling "discovery loop" holds N=101 to 114 and is out
of scope here.

## Method

We read the Results table and the dated "History of updates" log at
https://www.packomania.com/csqv/csqv.html . We judged freshness from the
reference tag and the dated log, not from per-N page footers (those show a
page-regeneration date, not a record date). We spot-checked coordinate files at
https://www.packomania.com/csqv/txt/csqv<N>.txt to confirm warm-start data.

Reference tags used below:
- [1] David W. Cantrell, sci.math 2011/12 (old, pre-wave).
- [2] E. Specht, program csqv, 2011 to 2026 (site owner; old OR fresh 2026).
- [8]..[20] are 2026-wave contributors (Lin, Dutton, Liang, Abiyev, Tasseff,
  Sander, Huang, Tallbacka, Heap, Garg, Zhuang, Tariq), see the References
  block at https://www.packomania.com/csqv/csqv.html .

## Headline finding: the wider band is NOT soft by date

The 2026 wave swept the whole wider band. There are ZERO Cantrell-era ([1])
untouched records above N=100. Every N that exists in the Results table from
115 to 250, and every N above 250, was touched on or after 2026-07-23.

Two sweeps account for it, both from the update log at
https://www.packomania.com/csqv/csqv.html :
- 2026-09-09, Specht [2] generated all perfect-square packings N=n^2 for
  n=11..100 by generalizing Cantrell's square pattern. He describes it as a
  construction where he lets "the solver do the rest".
- 2026-09-12, Zeeshan Tariq [20] delivered a large bundle including N=221, 265,
  313, 365, 421, 481, 545, 613, 685, 761, 925 (plus 144, 169 and the 115..170
  range).

The specific high-N spot-checks requested (365, 421, 481, 545, 613, 685, 761,
925) are all Tariq [20] FRESH. N=250 and N=300 have no table entry at all.

## Status of every relevant table entry

Band 1 (N=115 to 250, excluding the 101 to 114 sibling range):

| N | sum of radii | tag | freshness (log date, who) | coords |
|---|---|---|---|---|
| 115 | 5.652195250495 | [20] | FRESH 12-Sep Tariq, 09-Sep Heap | yes |
| 116 | 5.676465244072 | [20] | FRESH Tariq, Heap, Garg | yes |
| 117 | 5.699606654188 | [20] | FRESH Tariq, Heap | yes |
| 118 | 5.723940934671 | [20] | FRESH Tariq, Heap | yes |
| 119 | 5.748642197923 | [20] | FRESH Tariq, Heap | yes |
| 120 | 5.773179664812 | [20] | FRESH Tariq, Abiyev, Zhuang, Heap | yes |
| 121 (11^2) | 5.797190396246 | [2] | FRESH 09-Sep Specht square, un-reclaimed | yes (verified) |
| 122 | 5.822533741922 | [18] | FRESH 10-Sep Garg | yes |
| 123 | 5.846120697887 | [18] | FRESH Garg | yes |
| 124 | 5.870322310717 | [18] | FRESH 10-Sep Garg | yes |
| 125 | 5.901645773238 | [20] | FRESH Tariq, Abiyev, Heap | yes |
| 130 | 6.020600850901 | [15] | FRESH 03-Sep Huang, Abiyev | yes |
| 135 | 6.132840789573 | [20] | FRESH Tariq, Heap | yes |
| 140 | 6.245907944668 | [20] | FRESH Tariq, Huang | yes |
| 144 (12^2) | 6.334929299516 | [20] | FRESH 12-Sep Tariq (reclaimed) | yes |
| 145 | 6.362566039232 | [20] | FRESH Tariq, Garg, Heap | yes |
| 150 | 6.475472914692 | [20] | FRESH Tariq, Heap | yes |
| 155 | 6.577414557048 | [20] | FRESH Tariq, Heap | yes |
| 160 | 6.686377484769 | [15] | FRESH Huang | yes |
| 169 (13^2) | 6.875964055290 | [20] | FRESH 12-Sep Tariq (reclaimed) | yes |
| 170 | 6.899613322825 | [20] | FRESH Tariq, Heap | yes |
| 181 | 7.116795383460 | [15] | FRESH Huang | yes |
| 196 (14^2) | 7.404833875833 | [2] | FRESH 09-Sep Specht square, un-reclaimed | yes (verified) |
| 221 | 7.851345119145 | [20] | FRESH 12-Sep Tariq | yes |
| 225 (15^2) | 7.940727647130 | [2] | FRESH 09-Sep Specht square, un-reclaimed | yes (verified) |

Note: N between the listed values (for example 126 to 129, 131 to 134, 250)
have no table entry. There is no record to beat and no coordinate file to
warm-start.

Band 2 (N>250), pattern of the entire high band:
- Every perfect square (256, 289, 324, 361, 400, 441, 484, 529, 576, 625, 676,
  729, 784, 841, 900, 961, 1024, ... up to 10000) is a Specht [2] construction
  dated 2026-09-09. FRESH by date, weak by construction, but un-reclaimed.
- Every non-square listed above 250 (265, 313, 365, 421, 481, 545, 613, 685,
  761, 925) is Tariq [20], 2026-09-12. FRESH and actively contested.
- 841 (29^2) also appeared in Specht's 08-Aug D4 "far from being records" set,
  then stayed [2].

## Ranked shortlist (softest warm-startable N in the wider band)

All picks are Specht [2] perfect-square constructions dated 2026-09-09, the only
tier in the wider band that is single-touch, un-reclaimed, and carries an
explicit "construction, not record" signal. Ranking favors the in-band range
first, then smaller N for compute feasibility.

1. N=121 (11^2). Sum 5.797190396246. Coords verified at
   https://www.packomania.com/csqv/txt/csqv121.txt . Smallest, most feasible.
2. N=196 (14^2). Sum 7.404833875833. Coords verified at
   https://www.packomania.com/csqv/txt/csqv196.txt .
3. N=225 (15^2). Sum 7.940727647130. Coords verified at
   https://www.packomania.com/csqv/txt/csqv225.txt .
4. N=256 (16^2). Sum 8.476802859291. First above 250. Coords verified at
   https://www.packomania.com/csqv/txt/csqv256.txt .
5. N=289 (17^2). Sum 9.012780594064. Coords expected at
   https://www.packomania.com/csqv/txt/csqv289.txt (not individually verified).

Exclude N=144 and N=169: they are squares but Tariq reclaimed them on
2026-09-12, so they are contested.

## Caveat: soft is not beatable, and these are not even soft by date

These squares are FRESH 2026 constructions, not stale pre-2026 records. Their
appeal is a construction-quality signal (Specht built them by a pattern rule and
a light solver polish, not by heavy optimization), not a staleness signal. Small
or highly symmetric packings can be near-optimal, so this list is input to a
compute-feasibility check and a per-N beatability calibration, not a final
target set. Two active contributors (Specht [2] and Tariq [20]) are still moving
records across the band, so any candidate can be reclaimed before we act.

## Sources

- CSQV Results table, update log, and References:
  https://www.packomania.com/csqv/csqv.html
- Coordinate files (warm-start), verified:
  https://www.packomania.com/csqv/txt/csqv121.txt ,
  https://www.packomania.com/csqv/txt/csqv196.txt ,
  https://www.packomania.com/csqv/txt/csqv225.txt ,
  https://www.packomania.com/csqv/txt/csqv256.txt

## Open questions

- Do the un-reclaimed square constructions (196, 225, 256, 289 and larger) hold
  measurable slack? Run a per-N polish on the warm start and compare to the
  published sum before committing compute.
- Where does the compute ceiling fall? N=121 to 289 is far larger than our prior
  targets. Confirm the solver scales before picking a square target.
- The wider band has no pre-2026 soft targets. Should we redirect to the
  remaining SOFT sub-100 Cantrell entries (13, 20, 60, 90) instead, where the
  record is genuinely old and untouched?
- Is the square construction pattern uniformly weak, or only weak for some n?
  Calibrate a few squares before assuming the whole family is beatable.

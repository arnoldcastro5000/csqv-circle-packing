# Phase 2 CSQV: resume and warm restart

Research date: 2026-09-13. This note answers whether and how a second-wave record run should
build on a prior session, so the Opus wave CONTINUES the search instead of re-rolling the same
fixed start. It backs ticket 11
([.scratch/csqv/issues/11-resume-or-seed-a-record-run-from-a-prior-session.md](../../.scratch/csqv/issues/11-resume-or-seed-a-record-run-from-a-prior-session.md)).
Read the method note first for the pipeline and the sibling precedent
([docs/research/phase-2-csqv-record-methods-and-target-selection.md](phase-2-csqv-record-methods-and-target-selection.md)).

This note uses primary sources: the framework papers (FunSearch, AlphaEvolve, EoH, ReEvo), the
official code (google-deepmind/funsearch, algorithmicsuperintelligence/openevolve), the
basin-hopping literature (Wales and Doye; Grosso, Locatelli, Schoen), the novelty-search
paper, and the Gurobi warm-start manual. Each claim carries an inline URL. The prose uses short
active sentences and no em dash.

## Summary

Checkpoint and resume is an ENGINEERING feature, not a research contribution. FunSearch and
AlphaEvolve describe the run state (islands, programs, scores) but the reference FunSearch code
holds it in memory with no save or load. OpenEvolve, the widely used open AlphaEvolve, adds a
concrete disk checkpoint and a `--checkpoint` resume that reloads the whole program database and
the iteration counter. The unit these frameworks carry forward is the PROGRAM population, not
the concrete artifact. The inner local search (MBH, PBH, basin hopping) is the layer that
carries the ARTIFACT: it seeds each step from the current best CONFIGURATION and recomputes the
score. Warm restart from a prior best is standard and legitimate in both layers, because an
independent verifier recomputes the objective from the coordinates. The crux finding: a plain
warm restart from a single plateaued point re-collapses into the same optimum. The lever that
beats a plateau is DIVERSITY, delivered by island reset, migration, perturbation, and novelty
pressure, not the resume mechanism by itself.

---

## 1. Resume mechanics in LLM program-evolution frameworks

### FunSearch: islands and program database, but no persistence in the reference code

FunSearch splits the population into islands. It stores correct programs in a program database
and samples k programs from an island to build each prompt
([Nature s41586-023-06924-6](https://www.nature.com/articles/s41586-023-06924-6)). The
reference code confirms the run state: `ProgramsDatabase` holds `self._islands`, per-island
clusters keyed by score signature, and `self._best_program_per_island`, and it draws with
`np.random` ([programs_database.py](https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py)).
The same file has NO save, load, checkpoint, or serialize method. The pipeline is
single-threaded and in-memory, so a run cannot resume or seed from a prior run's database out of
the box ([programs_database.py](https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py);
[github.com/google-deepmind/funsearch](https://github.com/google-deepmind/funsearch)).

### AlphaEvolve: an evolutionary database, described but not released

AlphaEvolve stores solutions "with evaluation results (scores and program outputs) attached to
them" and manages them with an algorithm "inspired by a combination of the MAP elites
algorithm and island-based population models". The stated goal is "to optimally resurface
previously explored ideas in future generations" and to balance exploitation of the best
programs against diversity for exploration
([arXiv:2506.13131](https://arxiv.org/abs/2506.13131);
[ar5iv full text](https://ar5iv.labs.arxiv.org/html/2506.13131)). The white paper gives no
checkpoint schema, no RNG detail, and no cross-run seeding protocol. AlphaEvolve is not open, so
its persistence is unverifiable from a primary source.

### OpenEvolve: the concrete, verifiable resume implementation

OpenEvolve is the open AlphaEvolve. It is the strongest primary source for a real resume path.
The default config saves a checkpoint every 10 iterations and seeds every component from one
random seed:

```
checkpoint_interval: 10   # Save checkpoints every N iterations
random_seed: 42           # Random seed for reproducibility (null = random, 42 = default)
```

([configs/default_config.yaml](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/configs/default_config.yaml)).

The CLI exposes `--checkpoint`, "Path to checkpoint directory to resume from". On resume it
calls `openevolve.database.load(args.checkpoint)` and reports "Checkpoint loaded successfully
(iteration {database.last_iteration})"
([cli.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/cli.py)).
The `save` method writes each program to disk and a single `metadata.json` that persists the
full database state ([database.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py)):

```
metadata = {
    "island_feature_maps": ...,
    "islands": ...,
    "archive": ...,
    "best_program_id": ...,
    "island_best_programs": ...,
    "last_iteration": ...,
    "current_island": ...,
    "island_generations": ...,
    "last_migration_generation": ...,
    "feature_stats": ...,
}
```

Two facts matter for us. First, the checkpoint IS the whole program database. A resume points at
any prior checkpoint directory, so a run CAN seed from a DIFFERENT run's database: you pass its
checkpoint path ([cli.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/cli.py)).
Second, `metadata.json` does NOT store the RNG stream. The database only re-seeds `random` from
`config.random_seed` at construction
([database.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py)).
So a resume restores the POPULATION and the COUNTERS, but not the exact random stream. It is a
warm restart, not a bit-exact continuation. This is the key pitfall for Section 6.

### EoH and ReEvo: initialize from scratch, no cross-run resume as a feature

EoH builds the initial population by prompting the LLM to produce heuristics "entirely from
scratch" through one initialization operator plus five variation operators
([arXiv:2401.02051](https://arxiv.org/abs/2401.02051)). ReEvo adds reflective evolution with
"verbal gradients" over an evolutionary population, and reports sample efficiency, not a resume
protocol ([arXiv:2402.01145](https://arxiv.org/abs/2402.01145)). Neither paper treats
checkpoint or cross-run seeding as a contribution. The pattern holds across the lineage: resume
is an implementation convenience, and only OpenEvolve documents it as a first-class feature.

Bottom line for question 1. Checkpoint and restart is not standard in the papers. When a system
does persist, it persists the PROGRAM DATABASE (islands, programs, scores, elite archive,
iteration and migration counters). RNG state is usually NOT persisted. Seeding one run from a
different run's database is possible where the resume takes a database path (OpenEvolve), but no
paper studies it.

---

## 2. Program vs solution transfer

The frameworks answer this by their architecture. The evolving unit is the PROGRAM, and the
database carries programs with their scores, not concrete artifacts. FunSearch stores programs
in islands and re-runs them to score them
([programs_database.py](https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py)).
AlphaEvolve stores "solutions...with evaluation results (scores and program outputs) attached",
where the solution is the algorithm and the output is what running it produces
([ar5iv:2506.13131](https://ar5iv.labs.arxiv.org/html/2506.13131)). The artifact is
reproducible from the program, so carrying the program carries the capability, and the
population carries the search diversity too.

The local-search layer answers the opposite way. MBH, PBH, and basin hopping carry the ARTIFACT:
the concrete configuration is the state that seeds the next step (Section 3). This is the right
level for a warm start, because a good configuration is an expensive discovery that the program
cannot cheaply re-derive.

No primary source gives a head-to-head ablation of "carry the program" versus "carry the
artifact"; the two live at different layers and are complementary. The defensible synthesis for
CSQV: carry the PROGRAM and its population to keep the evolved SEARCH capability and its
diversity, and, separately, hand the inner solver the prior best PACKING as a warm-start seed so
each program does not re-discover it. Carrying the artifact alone, without the population, throws
away the search machinery the wave learned. Carrying the program alone, without a warm-start
packing, forces every candidate to re-find a jamming that the incumbent already encodes.

---

## 3. Warm restart in MBH, PBH, and basin hopping

Warm restart from a prior best configuration is the DEFINITION of basin hopping, not an add-on.

Wales and Doye 1997 transform the energy landscape so that "any point in configuration space" is
associated with "the local minimum obtained by geometry optimization started from that point",
which turns the surface into interpenetrating basins. The method steps between minima, and it
found the lowest known structures for all Lennard-Jones clusters up to 110 atoms, including some
never found before in unbiased searches
([Wales and Doye, J. Phys. Chem. A 1997, DOI 10.1021/jp970984n](https://pubs.acs.org/doi/10.1021/jp970984n);
[abstract, ADS](https://ui.adsabs.harvard.edu/abs/1997JPCA..101.5111W/abstract);
[abstract, Oxford ORA](https://ora.ox.ac.uk/objects/uuid:6a36a972-20a7-4173-9943-f90f45c37a2c)).
Each iteration perturbs the CURRENT minimized coordinates, minimizes again, and accepts by a
Metropolis rule on the energy. The current best configuration is the seed of the next step,
which is exactly a warm restart.

The scipy reference implementation states the same loop as a clean primary. `basinhopping`
runs a "random perturbation of the coordinates", then a "local minimization", then "accept or
reject" the new minimum on the "minimized function value" by a Metropolis criterion with
temperature T: a step that lowers the value is always accepted, a step that raises it is
accepted with a probability set by T. The `x0` argument is the starting coordinates and the
`seed` argument fixes the RNG, and the docs credit the method to "David Wales and Jonathan
Doye" 1997 ([scipy.optimize.basinhopping](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html)).
So a warm restart in this family is literally: pass the prior best coordinates as `x0`.

Grosso, Addis, Locatelli, and Schoen apply the same idea to circle packing. Monotonic Basin
Hopping (MBH) keeps the incumbent, perturbs "each coordinate of each" circle by a step delta,
runs a local search, and accepts only if the new local minimum is better. The perturbation must
"not completely disrupt" the current minimizer, so the walk moves between close minima; too small
a delta falls back into the same basin. Population Basin Hopping (PBH) keeps a population of
minimizers and uses a dissimilarity measure with a cut `dcut` to keep members diverse. The
approach improved putative optima for n <= 130 in as many as 32 instances
([Grosso, Locatelli, Schoen, optimization-online 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

When does continuation beat a fresh restart? The landscape decides. Both sources report a
"funneling" landscape, where good minima cluster in a funnel and a perturb-and-polish walk
descends it faster than independent multistart, whose local-minima count grows quickly
([Wales and Doye 1997](https://pubs.acs.org/doi/10.1021/jp970984n);
[Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
So warm restart wins when the incumbent already sits in the right funnel and headroom remains
near it. It does not win when the incumbent is a hard local optimum in the wrong funnel; then the
walk needs a large jump or a fresh start to reach a different funnel. This is the same tension
Section 4 resolves.

---

## 4. Basin-trapping and plateau escape (the crux)

The plateau is real and predictable. A warm restart that perturbs a single plateaued point with
a small step re-minimizes into the SAME basin. Grosso et al. state it directly: if the
perturbation delta is too small, the perturbed point "falls back into the same basin"
([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
A same-start re-roll and a naive warm restart therefore converge to the same optimum. The
mechanisms that break this are all forms of DIVERSITY, not the resume itself.

- Larger and adaptive jumps. MBH tunes delta so the move escapes the current basin without
  destroying the structure. A plateaued search needs a bigger or occasional very large jump to
  leave the funnel ([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
- Population dissimilarity. PBH holds a diverse population with a `dcut` cut, so members occupy
  different basins and the search does not collapse to one point
  ([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
- Acceptance temperature. A Metropolis accept with a temperature T lets the walk climb uphill
  and leave a basin. T raises the probability of accepting a worse minimum, so a higher T loosens
  the walk off a plateau; monotonic MBH sets T to zero and only descends
  ([scipy.optimize.basinhopping](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html);
  [Wales and Doye 1997](https://pubs.acs.org/doi/10.1021/jp970984n)).
- Island reset (diversity injection). FunSearch periodically discards the weaker half of the
  islands and reseeds them from survivors. The reference method is `reset_islands`, which "Resets
  the weaker half of islands", sorts by score with noise to break ties, and seeds each reset
  island from a "founder" program of a surviving island
  ([programs_database.py](https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py)).
  This is explicit diversity injection to prevent premature convergence.
- Island migration. OpenEvolve migrates top programs between islands on a schedule
  (`migration_interval` default 10, `migration_rate` default 0.1) to spread good material without
  merging the populations ([database.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py)).
- Novelty pressure. Lehman and Stanley show that ambitious objectives are deceptive: "most
  ambitious objectives do not illuminate a path to themselves", and the objective gradient leads
  to dead-end local optima. Rewarding behavioral novelty alone escapes this deception
  ([Lehman and Stanley 2011, Evolutionary Computation 19(2)](https://dl.acm.org/doi/abs/10.1162/evco_a_00025);
  [PDF](https://www.cs.swarthmore.edu/~meeden/DevelopmentalRobotics/lehmanNoveltySearch11.pdf)).
  For us this argues for a diversity or novelty term when the sum-of-radii objective plateaus.

Does resuming a plateaued search help? Only if the resume also injects diversity. A pure resume
restores the plateaued state, which is exactly the trap. The real lever is diversity injection:
island reset, migration, larger jumps, and population dissimilarity. Resume is the vehicle;
diversity is the fuel.

---

## 5. Objective-value hiding under warm start

The precedent for giving an optimizer the prior best COORDINATES while withholding the target
VALUE is standard warm-start practice, and it is not gaming when an independent check recomputes
the score.

Basin hopping and MBH seed from the prior best CONFIGURATION and recompute the energy or the
sum of radii themselves. The value is never an input; it is an output of scoring the
configuration ([Wales and Doye 1997](https://pubs.acs.org/doi/10.1021/jp970984n);
[Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).

Mixed-integer solvers make the same split explicit. Gurobi accepts a MIP start, "a starting
solution with values for continuous and discrete variables", and it recommends providing values
only for the integer variables and letting the solver compute the rest. The user supplies the
variable values, and the solver computes and validates the objective
([Gurobi MIP starts](https://support.gurobi.com/hc/en-us/articles/360043834831-How-do-I-use-MIP-starts);
[Gurobi warm-start manual](https://docs.gurobi.com/projects/optimizer/en/current/features/warmstart.html)).
A warm start that carries variables but not the optimum is the norm.

The legitimacy argument for CSQV is the independent verifier. The sibling discovery loop scores
every candidate with a verifier that "(1) Confirms all circles are within [0,1]^2 (with zero
tolerance) (2) Confirms no pair of circles overlaps (with zero tolerance) (3) Recomputes the sum
of radii independently (4) Applies a strict feasibility shrink"
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). Because the verifier
recomputes the sum from the coordinates, a warm start that hands over coordinates cannot inflate
the reported score, and a program that never sees the target VALUE cannot hardcode or fit to it.
This is the standard defense against benchmark gaming: keep the target value out of the solver
and let an independent recompute settle the score. A carried PACKING is a legitimate MBH
continuation. A carried PROGRAM is safe only if the genome never receives the best-known value.

---

## 6. Resumable-run checkpoint schema

The only fully specified primary schema is OpenEvolve's `metadata.json` plus per-program files
(Section 1). It persists the population and the archive (`islands`, `island_feature_maps`,
`archive`, `best_program_id`, `island_best_programs`, `feature_stats`) and the counters
(`last_iteration`, `current_island`, `island_generations`, `last_migration_generation`), and it
stores each program with its metrics and optional prompts
([database.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py)).
That set maps cleanly onto a general checkpoint schema for an evolutionary plus local-search run:

- Population and islands: every program with its source, score, and island id.
- Elite archive: the MAP-Elites or elite pool that guards diversity.
- Best-so-far: the global best and the per-island best, program and packing both.
- Counters: iteration, generation, current island, last migration.
- History: the prompt or few-shot context per program, if reflection is used.
- RNG state: the random stream, to make a resume deterministic.

Pitfalls, from the primary evidence:

- Nondeterminism from missing RNG. OpenEvolve saves only `random_seed` in config, not the live
  RNG stream, so a resume continues with a re-seeded generator, not the original stream
  ([database.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py);
  [default_config.yaml](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/configs/default_config.yaml)).
  A bit-exact resume must persist the RNG state of every stochastic component. The LLM sampler is
  nondeterministic anyway, so treat a resume as a WARM restart, not a replay.
- Stale scores. A score is only meaningful under the exact scoring code and budget. If the
  verifier, the timeout, or the solver changes between sessions, carried scores are stale and can
  mislead selection. Re-score the carried elite under the current verifier before you trust it.
  The sibling verifier already recomputes the sum, so make a re-score on resume mandatory
  ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).
- Feasibility drift. Carry the packing at full precision and re-run the feasibility shrink on
  load, because a serialized configuration can pick up numeric artifacts
  ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)).

---

## Recommendation for ticket 11

WHAT to carry. Carry BOTH layers, at their right level.
- Carry the PROGRAM population, not one champion. Reseed the elite pool from the prior session's
  programs so the Opus wave inherits the SEARCH capability and its diversity, the same unit
  FunSearch and AlphaEvolve carry
  ([programs_database.py](https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py);
  [ar5iv:2506.13131](https://ar5iv.labs.arxiv.org/html/2506.13131)).
- Carry the best PACKING as an inner-solver warm start via `ctx.incumbent`, MBH style, so no
  candidate re-discovers the incumbent jamming
  ([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)).
- Carry HISTORY (prompt context) only if the operator uses reflection. It is optional and lower
  payoff than the population and the packing.

CLI shape. Provide BOTH shapes, because they serve different intents.
- `--resume <session-dir>` continues one run: it reloads that session's population, elite pool,
  counters, and incumbent. This mirrors OpenEvolve `--checkpoint`
  ([cli.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/cli.py)).
- `--seed-from <program.py>` and `--warm-from <packing>` compose a fresh run from named assets,
  which covers the Sonnet to Opus handoff and any cross-run seed. This matches OpenEvolve's
  ability to resume from any checkpoint path, that is, from a different run
  ([cli.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/cli.py)).
- Persist a resumed run in a NEW session dir that records its parent, so the lineage stays
  auditable and no prior run is overwritten.

Cross-model safety (Sonnet to Opus). Nothing carried is model-specific. The population is Python
source, the packing is coordinates, the scores are numbers, and the history is text. None encodes
the generating model. Re-score the carried elite under the current verifier on load, because a
score is only valid under its scoring code
([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093)). Keep the best-known VALUE out
of the genome; a carried packing is a legitimate warm start because the verifier recomputes the
sum ([arXiv:2609.05093 full text](https://arxiv.org/html/2609.05093);
[Gurobi MIP starts](https://support.gurobi.com/hc/en-us/articles/360043834831-How-do-I-use-MIP-starts)).

Diversity injection: needed. State it plainly. A resume alone will NOT beat Sonnet's plateau. A
warm restart from a plateaued point re-collapses into the same basin
([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)),
and an ambitious objective is deceptive near its optimum
([Lehman and Stanley 2011](https://dl.acm.org/doi/abs/10.1162/evco_a_00025)). The Opus wave must
add diversity: island reset that reseeds weak islands from survivors
([programs_database.py](https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py)),
migration between islands
([database.py](https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py)),
larger or adaptive perturbation jumps and a diverse PBH population
([Grosso et al. 1999.pdf](https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf)),
and, if the plateau persists, a novelty or diversity term alongside the sum-of-radii objective
([Lehman and Stanley 2011](https://www.cs.swarthmore.edu/~meeden/DevelopmentalRobotics/lehmanNoveltySearch11.pdf)).
Build resume for continuity, but treat diversity injection as the actual mechanism that beats the
plateau.

---

## Sources

- FunSearch, Nature: https://www.nature.com/articles/s41586-023-06924-6
- FunSearch reference code, ProgramsDatabase and reset_islands:
  https://raw.githubusercontent.com/google-deepmind/funsearch/main/implementation/programs_database.py
  and https://github.com/google-deepmind/funsearch
- AlphaEvolve: https://arxiv.org/abs/2506.13131 and https://ar5iv.labs.arxiv.org/html/2506.13131
- OpenEvolve default config (checkpoint_interval, random_seed):
  https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/configs/default_config.yaml
- OpenEvolve CLI (--checkpoint resume, database.load, last_iteration):
  https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/cli.py
- OpenEvolve database (metadata.json schema, migration, RNG seeding):
  https://raw.githubusercontent.com/algorithmicsuperintelligence/openevolve/main/openevolve/database.py
- EoH (Evolution of Heuristics): https://arxiv.org/abs/2401.02051
- ReEvo (Reflective Evolution): https://arxiv.org/abs/2402.01145
- Wales and Doye 1997, basin hopping for Lennard-Jones clusters:
  https://pubs.acs.org/doi/10.1021/jp970984n and
  https://ui.adsabs.harvard.edu/abs/1997JPCA..101.5111W/abstract and
  https://ora.ox.ac.uk/objects/uuid:6a36a972-20a7-4173-9943-f90f45c37a2c
- scipy basinhopping reference (perturb, minimize, Metropolis with temperature T; x0, seed):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html
- Grosso, Locatelli, Schoen, MBH and PBH for circle packing:
  https://optimization-online.org/wp-content/uploads/2008/06/1999.pdf
- Lehman and Stanley 2011, novelty search:
  https://dl.acm.org/doi/abs/10.1162/evco_a_00025 and
  https://www.cs.swarthmore.edu/~meeden/DevelopmentalRobotics/lehmanNoveltySearch11.pdf
- Gurobi MIP starts (warm start with variable values):
  https://support.gurobi.com/hc/en-us/articles/360043834831-How-do-I-use-MIP-starts and
  https://docs.gurobi.com/projects/optimizer/en/current/features/warmstart.html
- Sibling discovery loop, independent verifier: https://arxiv.org/abs/2609.05093 and
  https://arxiv.org/html/2609.05093

## Open questions

- The exact FunSearch `reset_islands` period and OpenEvolve migration cadence are code defaults,
  not tuned for CSQV. A local sweep of reset period, migration rate, and jump size at the target
  N would size the diversity knobs.
- No primary source ablates "carry program" versus "carry artifact" head to head. Our own A/B on
  the target N (resume with population only, packing only, both) would settle it for CSQV.
- AlphaEvolve persistence is undescribed in the white paper. OpenEvolve is the closest verifiable
  proxy, and its RNG is not fully checkpointed.

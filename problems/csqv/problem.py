"""Variable-radii circle packing (CSQV): the Phase-2b problem module for the loop.

This module plugs into `loop.py` exactly like the Phase-1 bin-packing and Phase-2
CVRP modules: it exposes `train_instances`, `test_instances`, `evaluate_source`,
`SEED_SOLVER_SOURCE`, `SYSTEM_INSTRUCTION`, and `build_prompt`. The orchestrator
stays problem-agnostic and imports this module by name.

The vocabulary (search program, packing, sum of radii, feasibility shrink,
radii-LP, polish, selection score, train N, held-out N, naive baseline) is defined
in the project glossary CONTEXT.md. The locked decisions live in the wayfinder
tickets under `.scratch/csqv/`:
- genome, verifier, and fixed-primitive contract (ticket 01),
- fitness, baseline, N-band, and staleness protocol (ticket 02),
- data, best-known values, and acceptance tolerance (ticket 03).

Genome contract (the Discovery-Loop framing, wider than the Phase-1/2 selection
genome): the LLM writes exactly one function `pack(n, ctx)` that returns an
(n, 3) array of [x, y, r] rows in the unit square [0, 1] x [0, 1]. The genome is a
whole search program: it initializes, perturbs, restarts, and refines. It may call
two FIXED trusted primitives exposed on `ctx`, but never redefines them:
- `ctx.radii_lp(centers)`: the exact largest radii for fixed centers, by a linear
  program;
- `ctx.polish(packing)`: a local numeric refinement (SLSQP) of centers and radii.
The genome also gets `ctx.rng` (a per-N-seeded generator), `ctx.deadline`, and
`ctx.time_left()`. It never sees the best-known values, so it cannot game the
objective.

Scoring maps onto the loop's generic `excess` (lower is better). For CSQV the
verified objective is the exact sum of radii; the selection quantity is the gap to
best-known, `(bks - sum) / bks`, averaged across the N in the set. A packing that
matches best-known scores 0; one that beats it scores negative (a record). The
independent verifier recomputes containment and non-overlap, applies a FEASIBILITY
SHRINK, and scores the shrunk sum, so a near-feasible candidate scores on an honest
valid value and a non-finite or out-of-bounds candidate is INVALID.

Security note: `evaluate_source` execs solver source in-process, so the loop calls
it only inside the sandbox subprocess with a hard timeout. The seed solver and the
tests are trusted, so they may call it directly. scipy is a hard dependency of this
module (install: `pip install --user --break-system-packages scipy`).
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import numpy.typing as npt
from scipy.optimize import linprog, minimize

from . import data

Array = npt.NDArray[np.float64]

# The unit square the genome and verifier use. The stored anchors live in the
# centered frame; data.py converts. TOL is the acceptance tolerance (data.py).
SIDE = data.SIDE
TOL = data.TOL

# --- N band and split (ticket 02) ------------------------------------------
# The anchored mid band, validation first: every target N has ground-truth
# coordinates, so the gap to best-known is exactly checkable. Train N drive
# selection; the disjoint held-out N are interleaved and reported only, never
# selected on (the overfitting guard; CSQV has one packing per N).
TRAIN_NS = [20, 40, 60, 90]
TEST_NS = [32, 50, 75]

# One soft time budget for the whole candidate, split across its N inside
# evaluate_source (ticket 02: one total budget, not a fixed per-N budget). Each N
# gets an equal share of the budget still left, so a cheap early N leaves more for
# an expensive later one. The subprocess timeout in the loop is the HARD backstop
# and must exceed this budget plus a margin; a CSQV run passes a large --timeout
# (the loop default of a few seconds is for the fast Phase-1/2 problems, not CSQV).
TOTAL_BUDGET_S = 90.0
# The SLSQP iteration cap for the polish primitive, scaled down for large N inside
# polish so a single solve stays bounded.
POLISH_MAXITER = 200

# The rng seed is a fixed function of N, so evaluation is reproducible and the
# genome cannot special-case a known seed to leak the answer (it never sees the
# best-known anyway).
BASE_SEED = 2026


class SolverError(Exception):
    """A malformed genome. It makes the candidate INVALID, never a loop crash."""


class SolverInterrupted(BaseException):
    """Raised INTO the genome (from a SIGTERM handler) to stop it on a deadline while
    keeping its best feasible work. It subclasses BaseException, not Exception, so a
    genome's own `except Exception` cannot swallow the stop; `run_pack` catches it and
    returns the best feasible packing seen so far instead of losing the whole run."""


@dataclass(frozen=True)
class Instance:
    """One CSQV instance: a circle count N and its published best-known sum."""

    n: int
    bks: float


@dataclass
class Result:
    """The outcome of scoring one candidate on one instance set.

    `excess` is the mean gap-to-best-known across the N (lower is better; negative
    means the candidate beat best-known). `worst_gap` is the single worst per-N gap,
    reported alongside so a candidate that sacrifices one N is visible; it rides only
    in-process (the sandbox serializes only valid/excess/runtime/error). An INVALID
    result is the worst possible score.
    """

    valid: bool
    excess: float | None = None
    runtime: float = 0.0
    error: str | None = None
    worst_gap: float | None = None


# --- Instance data (ticket 02, 03) -----------------------------------------
def _instances(ns: list[int]) -> list[Instance]:
    best_known = data.load_best_known()
    out: list[Instance] = []
    for n in ns:
        if n not in best_known:
            raise ValueError(f"no best-known sum for N={n}")
        out.append(Instance(n=n, bks=best_known[n]))
    return out


def train_instances() -> list[Instance]:
    """The set the search selects on."""
    return _instances(TRAIN_NS)


def test_instances() -> list[Instance]:
    """The held-out set, reported and watched for the plateau. Never selected on."""
    return _instances(TEST_NS)


# --- Fixed primitive: radii-LP ---------------------------------------------
def radii_lp(centers: Array) -> Array:
    """Return the largest radii that keep the packing valid for FIXED centers.

    This is a linear program: maximize the sum of radii subject to containment
    (each circle inside the unit square) and pairwise non-overlap
    (r_i + r_j <= dist(i, j)). It is trusted harness code, never part of the genome.
    A center outside the square yields a zero containment bound, so its radius is 0.
    On the rare solver failure it falls back to the containment bounds scaled to
    feasibility, so the returned radii never overlap.
    """
    centers = np.asarray(centers, dtype=np.float64)
    if centers.ndim != 2 or centers.shape[1] != 2:
        raise SolverError(f"radii_lp expects (n, 2) centers, got {centers.shape}")
    n = centers.shape[0]
    x = centers[:, 0]
    y = centers[:, 1]
    # Containment upper bound per circle: the distance to the nearest wall.
    cbound = np.maximum(data.wall_slack(x, y, SIDE), 0.0)
    if n < 2:
        return cbound.astype(np.float64)

    i, j, dist = data.pairwise_distances(x, y)

    # One inequality row per pair: +1 at column i and column j, bounded by dist.
    rows = np.repeat(np.arange(i.size), 2)
    cols = np.empty(i.size * 2, dtype=np.int64)
    cols[0::2] = i
    cols[1::2] = j
    vals = np.ones(i.size * 2, dtype=np.float64)
    a_ub = np.zeros((i.size, n), dtype=np.float64)
    a_ub[rows, cols] = vals

    bounds = [(0.0, float(c)) for c in cbound]
    res = linprog(
        c=-np.ones(n),
        A_ub=a_ub,
        b_ub=dist,
        bounds=bounds,
        method="highs",
    )
    if not res.success or res.x is None:
        # Fallback: shrink the containment bounds to a feasible, non-overlapping
        # set rather than return the overlapping bounds themselves.
        _, s = _feasible_scale(np.clip(x, 0.0, SIDE), np.clip(y, 0.0, SIDE), cbound)
        return (cbound * s).astype(np.float64)
    return np.asarray(res.x, dtype=np.float64)


# --- Fixed primitive: polish -----------------------------------------------
def _pairs(n: int) -> tuple[Array, Array]:
    i, j = np.triu_indices(n, k=1)
    return i.astype(np.int64), j.astype(np.int64)


def polish(packing: Array, maxiter: int = POLISH_MAXITER, deadline: float | None = None) -> Array:
    """Locally refine centers and radii toward a higher sum, keeping feasibility.

    An SLSQP refinement over all [x, y, r] with analytic gradients: box bounds keep
    centers in the square and radii non-negative; the nonlinear constraints enforce
    containment and non-overlap. It is trusted harness code, never part of the
    genome. It never worsens the score: it returns whichever of the input and the
    refined packing has the larger feasibility-shrunk sum, so a failed or unhelpful
    refinement is a no-op. A non-finite or wrong-shaped input is returned unchanged.
    SLSQP cannot be interrupted mid-solve, so the iteration cap is scaled down for
    large N to keep a single call bounded, and the deadline is checked before the
    call. The genome should test ctx.time_left() before it calls polish.
    """
    return _polish_scored(packing, maxiter, deadline)[0]


def _polish_scored(
    packing: Array, maxiter: int = POLISH_MAXITER, deadline: float | None = None
) -> tuple[Array, float]:
    """`polish` plus the returned packing's feasibility-shrunk sum, so callers that need
    the score (Ctx.polish, to track the best feasible packing) do not verify it again.
    A no-op return carries `-inf` so the caller skips it without a second verify."""
    packing = np.asarray(packing, dtype=np.float64)
    n = packing.shape[0]
    if packing.shape != (n, 3) or not np.isfinite(packing).all() or n == 0:
        return packing, float("-inf")
    if deadline is not None and time.perf_counter() >= deadline:
        _, in_sum = verify_and_score(packing)
        return packing, in_sum

    # Bound a single solve: SLSQP cost grows fast in n, so cap iterations by n.
    eff_maxiter = max(30, min(maxiter, int(6000 / n)))

    i, j = _pairs(n)

    def neg_sum(z: Array) -> float:
        return -float(z[2::3].sum())

    def neg_sum_grad(z: Array) -> Array:
        g = np.zeros_like(z)
        g[2::3] = -1.0
        return g

    def constraints(z: Array) -> Array:
        p = z.reshape(n, 3)
        x, y, r = p[:, 0], p[:, 1], p[:, 2]
        contain = np.concatenate([x - r, SIDE - x - r, y - r, SIDE - y - r])
        if i.size == 0:
            return contain
        dx = x[i] - x[j]
        dy = y[i] - y[j]
        dist = np.sqrt(dx * dx + dy * dy)
        overlap = dist - r[i] - r[j]
        return np.concatenate([contain, overlap])

    def constraints_jac(z: Array) -> Array:
        p = z.reshape(n, 3)
        x, y, r = p[:, 0], p[:, 1], p[:, 2]
        m = 4 * n + i.size
        jac = np.zeros((m, 3 * n), dtype=np.float64)
        idx = np.arange(n)
        # x - r >= 0
        jac[idx, 3 * idx] = 1.0
        jac[idx, 3 * idx + 2] = -1.0
        # SIDE - x - r >= 0
        jac[n + idx, 3 * idx] = -1.0
        jac[n + idx, 3 * idx + 2] = -1.0
        # y - r >= 0
        jac[2 * n + idx, 3 * idx + 1] = 1.0
        jac[2 * n + idx, 3 * idx + 2] = -1.0
        # SIDE - y - r >= 0
        jac[3 * n + idx, 3 * idx + 1] = -1.0
        jac[3 * n + idx, 3 * idx + 2] = -1.0
        if i.size:
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dist = np.sqrt(dx * dx + dy * dy)
            dist = np.where(dist > 0, dist, 1.0)  # guard; coincident pairs are rejected later
            base = 4 * n + np.arange(i.size)
            jac[base, 3 * i] = dx / dist
            jac[base, 3 * j] = -dx / dist
            jac[base, 3 * i + 1] = dy / dist
            jac[base, 3 * j + 1] = -dy / dist
            jac[base, 3 * i + 2] = -1.0
            jac[base, 3 * j + 2] = -1.0
        return jac

    bounds = [(0.0, SIDE), (0.0, SIDE), (0.0, None)] * n
    try:
        res = minimize(
            neg_sum,
            packing.flatten(),
            jac=neg_sum_grad,
            method="SLSQP",
            bounds=bounds,
            constraints=[{"type": "ineq", "fun": constraints, "jac": constraints_jac}],
            options={"maxiter": eff_maxiter, "ftol": 1e-10},
        )
    except Exception:  # noqa: BLE001  a numeric failure is a no-op, never a crash
        _, in_sum = verify_and_score(packing)
        return packing, in_sum

    candidate = np.asarray(res.x, dtype=np.float64).reshape(n, 3)
    if not np.isfinite(candidate).all():
        _, in_sum = verify_and_score(packing)
        return packing, in_sum
    _, cand_sum = verify_and_score(candidate)
    _, in_sum = verify_and_score(packing)
    if cand_sum > in_sum:
        return candidate, cand_sum
    return packing, in_sum


# --- Fixed primitive: defect-migration -------------------------------------
def _defect_move_scored(
    packing: Array,
    rng: np.random.Generator,
    k: int = 1,
    deadline: float | None = None,
    polish_maxiter: int = POLISH_MAXITER,
) -> tuple[Array, float]:
    """`defect_move` plus the returned packing's feasibility-shrunk sum, so Ctx.defect_move
    tracks the best packing without a second verify. A no-op or invalid input carries its
    own score (`-inf` when invalid), so the caller skips it cleanly.

    One move: remove the k weakest circles (smallest radius), reinsert each into the largest
    empty hole (a jittered coarse grid probed for the biggest admissible radius, greedily so
    later insertions see earlier ones), reset ALL radii by the exact LP, then one polish.
    Deterministic given (packing, rng state, k). Unlike a center perturbation it changes the
    combinatorial contact graph, so it can reach a basin the perturbation walk never visits.

    NOT monotonic: unlike `polish` (which never worsens the score), this removes and reinserts
    circles, so the refit-and-polish result CAN score below the input. Callers must compare and
    keep the better packing; `Ctx.defect_move` does this through the best-feasible tracker, so a
    genome that funnels its output there never regresses. A genome that treats the return value
    as unconditionally better than its input can go backward.
    """
    packing = np.asarray(packing, dtype=np.float64)
    n = packing.shape[0]
    # n <= 1 has no meaningful defect to migrate; a wrong shape or non-finite is a no-op.
    if packing.shape != (n, 3) or not np.isfinite(packing).all() or n <= 1:
        _, s = verify_and_score(packing)
        return packing, s
    k = int(max(1, min(k, n - 1)))  # keep at least one circle in place

    r = packing[:, 2]
    weak = np.argsort(r)[:k]  # the k smallest radii are the defects
    keep = np.ones(n, dtype=bool)
    keep[weak] = False
    kept = packing[keep, :2].copy()  # (n - k, 2) centers that stay

    # Candidate holes: a jittered coarse grid, so repeated calls explore (deterministic
    # given rng). The grid density grows with n so larger packings still find gaps.
    g = max(6, int(np.ceil(np.sqrt(n))) * 3)
    axis = np.linspace(0.02 * SIDE, 0.98 * SIDE, g)
    gx, gy = np.meshgrid(axis, axis)
    cand = np.column_stack([gx.ravel(), gy.ravel()])
    cand = np.clip(cand + rng.normal(0.0, 0.5 * SIDE / g, size=cand.shape), 0.0, SIDE)

    placed = kept.copy()
    placed_r = packing[keep, 2].copy()  # the kept circles' existing radii (free and accurate)
    new_centers = np.empty((k, 2), dtype=np.float64)
    for m in range(k):
        # Admissible radius at each candidate: the largest circle that fits inside the
        # square and touches no already-placed circle.
        adm = np.maximum(data.wall_slack(cand[:, 0], cand[:, 1], SIDE), 0.0)
        if placed.shape[0] > 0:
            dx = cand[:, 0:1] - placed[:, 0][None, :]
            dy = cand[:, 1:2] - placed[:, 1][None, :]
            radial = np.sqrt(dx * dx + dy * dy) - placed_r[None, :]
            adm = np.minimum(adm, np.min(radial, axis=1))
        best = int(np.argmax(adm))
        new_centers[m] = cand[best]
        placed = np.vstack([placed, cand[best][None, :]])
        placed_r = np.append(placed_r, max(0.0, float(adm[best])))

    centers = np.vstack([kept, new_centers])
    refit = np.column_stack([centers, radii_lp(centers)])
    return _polish_scored(refit, maxiter=polish_maxiter, deadline=deadline)


def defect_move(
    packing: Array,
    rng: np.random.Generator,
    k: int = 1,
    deadline: float | None = None,
    polish_maxiter: int = POLISH_MAXITER,
) -> Array:
    """Remove the k weakest circles and reinsert them into the largest empty holes, then
    refit radii and polish. A fixed trusted primitive (ADR 0001), never part of the genome.
    See `_defect_move_scored` for the contract. Returns the input unchanged when it cannot
    improve the geometry (n <= 1, wrong shape, non-finite)."""
    return _defect_move_scored(packing, rng, k, deadline, polish_maxiter)[0]


# --- Independent verifier with the feasibility shrink -----------------------
def _feasible_scale(x: Array, y: Array, r: Array, tol: float = TOL) -> tuple[bool, float]:
    """Return (ok, s): the largest uniform factor s in [0, 1] that makes the packing
    strictly valid, and whether the geometry is correctable at all.

    Centers must already sit inside the square (clip first). `ok` is False only for
    two coincident centers, which no uniform shrink can separate. Otherwise s scales
    every radius so no circle leaves the square and no two overlap. The containment
    and non-overlap geometry comes from data.py, so there is one reference.
    """
    slack = np.maximum(data.wall_slack(x, y, SIDE), 0.0)
    s = 1.0
    pos = r > 0
    if pos.any():
        s = min(s, float((slack[pos] / r[pos]).min()))
    if x.shape[0] >= 2:
        i, j, dist = data.pairwise_distances(x, y)
        if (dist < tol).any():
            return False, 0.0  # coincident centers, uncorrectable
        pair_sum = r[i] + r[j]
        m = pair_sum > 0
        if m.any():
            s = min(s, float((dist[m] / pair_sum[m]).min()))
    return True, max(0.0, min(1.0, s))


def _validate_and_shrink(
    packing: Array, tol: float = TOL
) -> tuple[bool, float, Array, Array, Array]:
    """Core verifier: return (ok, s, x, y, r).

    `ok` is False (with the other fields undefined) when the packing is INVALID: the
    wrong shape, non-finite, a negative radius, a center outside the square by more
    than `tol`, or two coincident centers (uncorrectable by a shrink). Otherwise `x`
    and `y` are the centers clipped into the square, `r` is the raw radii clipped to
    non-negative (NOT yet scaled), and `s` in [0, 1] is the largest uniform factor
    that makes every circle fit and no two overlap. Both `verify_and_score` and
    `canonical_feasible` build on this, so there is one reference and no drift.
    """
    packing = np.asarray(packing, dtype=np.float64)
    n = packing.shape[0]
    empty = np.empty(0, dtype=np.float64)
    if packing.shape != (n, 3):
        return False, 0.0, empty, empty, empty
    if not np.isfinite(packing).all():
        return False, 0.0, empty, empty, empty

    x = packing[:, 0].copy()
    y = packing[:, 1].copy()
    r = packing[:, 2]
    if (r < -tol).any():
        return False, 0.0, empty, empty, empty
    r = np.clip(r, 0.0, None)

    # Centers must sit inside the square (a center outside cannot be shrunk valid).
    if (x < -tol).any() or (x > SIDE + tol).any() or (y < -tol).any() or (y > SIDE + tol).any():
        return False, 0.0, empty, empty, empty
    x = np.clip(x, 0.0, SIDE)
    y = np.clip(y, 0.0, SIDE)

    ok, s = _feasible_scale(x, y, r, tol)
    if not ok:
        return False, 0.0, empty, empty, empty
    return True, s, x, y, r


def verify_and_score(packing: Array, tol: float = TOL) -> tuple[bool, float]:
    """Recompute feasibility and score a unit-square packing on the shrunk sum.

    Trust nothing from the genome. A packing is INVALID (returns (False, -inf)) when
    it is the wrong shape, non-finite, has a negative radius, has a center outside
    the square by more than `tol`, or has two coincident centers (uncorrectable by a
    shrink). Otherwise the verifier finds the largest uniform factor s in [0, 1] that
    makes every circle fit inside the square and no two overlap, then scores the
    shrunk sum s * sum(r). A packing already feasible scores its raw sum (s = 1); a
    slightly-violating one loses only a little (the numerical-artifact shrink); a
    badly-violating one is penalized hard, down to a valid zero (ticket 02: no
    INVALID cliff at s = 0). The score is always on a strictly-valid packing, so the
    objective cannot be gamed.
    """
    ok, s, _x, _y, r = _validate_and_shrink(packing, tol)
    if not ok:
        return False, float("-inf")
    return True, s * float(r.sum())


def canonical_feasible(packing: Array, tol: float = TOL) -> Array | None:
    """Return the strictly-valid packing the verifier scores, or None if INVALID.

    The centers are clipped into the square and every radius is scaled by the uniform
    feasibility factor s, so the result at most TOUCHES (no overlap, inside the
    square) and its sum of radii equals `verify_and_score`'s score. Use this to EMIT a
    packing (for example a .pck submission): the raw genome packing may hold tiny
    violations that a zero-tolerance external re-check would reject.
    """
    ok, s, x, y, r = _validate_and_shrink(packing, tol)
    if not ok:
        return None
    return np.column_stack([x, y, s * r])


# --- Context and solver loading --------------------------------------------
@dataclass
class Ctx:
    """The per-N context the genome receives. It exposes only N, a seeded rng, the
    two fixed primitives, a soft deadline, and (record mode only) the incumbent
    published packing to warm-start from; never the best-known sum. Seeing the
    incumbent coordinates is a legitimate MBH warm start, not gaming: the verifier
    still recomputes and shrinks, so a candidate must genuinely improve the sum."""

    n: int
    rng: np.random.Generator
    deadline: float
    _polish_maxiter: int = POLISH_MAXITER
    incumbent: Array | None = None
    # The best feasible packing (and its shrunk sum) seen through the primitives. It
    # lets run_pack preserve real work when the genome returns something worse, raises,
    # or is interrupted on a deadline. Private: the genome never reads it.
    _best_packing: Array | None = field(default=None, repr=False)
    _best_sum: float = field(default=float("-inf"), repr=False)

    def time_left(self) -> float:
        """Seconds left before the soft per-N deadline (never negative)."""
        return max(0.0, self.deadline - time.perf_counter())

    def _note_scored(self, packing: Array, s: float) -> None:
        """Remember `packing` (with its already-known feasibility-shrunk sum `s`) if it
        is the best seen so far, as a copy so a later genome mutation cannot corrupt the
        preserved best. `-inf` (an infeasible or no-op result) is never kept."""
        if s > self._best_sum:
            self._best_sum = s
            self._best_packing = np.array(packing, dtype=np.float64, copy=True)

    def _note(self, packing: Array) -> None:
        """Verify `packing` and remember it if it is the best feasible one seen. Used for
        the genome's own return, whose score is not yet known."""
        _, s = verify_and_score(packing)
        self._note_scored(packing, s)

    def radii_lp(self, centers: Array) -> Array:
        """The largest radii for fixed centers (the exact LP primitive)."""
        return radii_lp(centers)

    def polish(self, packing: Array) -> Array:
        """A local numeric refinement of centers and radii (the SLSQP primitive). Every
        polished result funnels through the best-feasible tracker (reusing the score
        polish already computed, so no extra verify), so the best packing is never
        lost even if the genome is later interrupted or returns something worse."""
        result, s = _polish_scored(packing, maxiter=self._polish_maxiter, deadline=self.deadline)
        self._note_scored(result, s)
        return result

    def defect_move(self, packing: Array, k: int = 1) -> Array:
        """Remove the k weakest circles and reinsert them into the largest empty holes,
        then refit and polish (the defect-migration primitive, ADR 0001). Like polish, the
        result funnels through the best-feasible tracker, so a strong packing is never lost
        even if the genome is later interrupted or returns something worse."""
        result, s = _defect_move_scored(
            packing, self.rng, k, deadline=self.deadline, polish_maxiter=self._polish_maxiter
        )
        self._note_scored(result, s)
        return result


def load_solver(source: str) -> Callable[[int, Ctx], Array]:
    """Compile solver source and return its `pack` callable.

    Raises on any source that defines no callable named `pack`.
    """
    namespace: dict = {"np": np, "numpy": np, "math": math}
    exec(source, namespace)  # noqa: S102  the loop sandboxes this in a subprocess
    pack = namespace.get("pack")
    if not callable(pack):
        raise ValueError("source defines no callable named 'pack'")
    return pack


def run_pack(
    pack: Callable[[int, Ctx], Array],
    n: int,
    deadline: float,
    polish_maxiter: int = POLISH_MAXITER,
    incumbent: Array | None = None,
    seed: int | None = None,
) -> Array:
    """Build the context, run the genome once for N, and return its packing.

    `deadline` is an absolute perf_counter time; the genome sees it through
    ctx.time_left(). `incumbent` (record mode) is a read-only warm-start packing.
    `seed` overrides the default (BASE_SEED + n) rng seed, so independent repeats
    of the SAME candidate explore different basins instead of one fixed stream.

    Never loses feasible work: the returned packing is the BEST feasible one seen,
    whether that is the genome's own return or a better packing it produced through
    ctx.polish. A SolverInterrupted (a deadline SIGTERM) or any genome fault returns
    that best feasible packing instead of failing; only a genome that produced NO
    feasible packing raises SolverError, which the caller turns into an INVALID result.
    """
    ro_incumbent: Array | None = None
    if incumbent is not None:
        ro_incumbent = np.array(incumbent, dtype=np.float64)
        ro_incumbent.flags.writeable = False  # the genome cannot mutate the incumbent
    ctx = Ctx(
        n=n,
        rng=np.random.default_rng(BASE_SEED + n if seed is None else seed),
        deadline=deadline,
        _polish_maxiter=polish_maxiter,
        incumbent=ro_incumbent,
    )
    try:
        raw = pack(n, ctx)
    except SolverInterrupted:
        if ctx._best_packing is not None:
            return ctx._best_packing  # deadline stop: keep the best feasible work
        raise SolverError("interrupted before any feasible packing") from None
    except Exception as error:  # noqa: BLE001  any genome fault is INVALID
        if ctx._best_packing is not None:
            return ctx._best_packing  # keep the best feasible work despite the fault
        raise SolverError(f"pack raised: {error}") from error
    try:
        arr: Array | None = np.asarray(raw, dtype=np.float64)
    except (ValueError, TypeError):
        arr = None  # a non-array return cannot be scored; it just is not the best
    if arr is not None:
        ctx._note(arr)  # count the genome's own return
    if ctx._best_packing is not None:
        return ctx._best_packing  # the best feasible packing seen, never lost
    if arr is None:
        raise SolverError("pack returned a non-array and produced no feasible packing")
    return arr  # no feasible packing; let the verifier reject it


# --- Scoring ---------------------------------------------------------------
def gap_on(sum_radii: float, bks: float) -> float:
    """The fraction under best-known for one N (lower is better; negative beats it)."""
    return (bks - sum_radii) / bks


def evaluate_source(
    source: str,
    instances: list[Instance],
    total_budget: float = TOTAL_BUDGET_S,
    polish_maxiter: int = POLISH_MAXITER,
) -> Result:
    """Compile, run, verify, and score a candidate on an instance set.

    `total_budget` is one soft wall-clock budget for the whole candidate; each N
    gets an equal share of the budget still left (ticket 02). Returns a valid Result
    with the mean gap-to-best-known across the N and the worst per-N gap, or an
    INVALID Result on any error or verifier rejection. This never raises for a bad
    genome; it returns INVALID instead.
    """
    start = time.perf_counter()
    if not instances:
        return Result(False, error="no instances provided")
    try:
        pack = load_solver(source)
    except Exception as error:  # noqa: BLE001  any load failure is INVALID
        return Result(False, error=f"load failed: {error}", runtime=time.perf_counter() - start)

    overall_deadline = start + total_budget
    gaps: list[float] = []
    try:
        for idx, inst in enumerate(instances):
            # Split the budget still left equally across the N not yet run, so a
            # cheap early N leaves more time for an expensive later one.
            remaining = overall_deadline - time.perf_counter()
            share = remaining / (len(instances) - idx)
            deadline = time.perf_counter() + max(0.0, share)
            packing = run_pack(pack, inst.n, deadline, polish_maxiter)
            if packing.shape != (inst.n, 3):
                return Result(
                    False,
                    error=f"pack returned shape {packing.shape}, expected {(inst.n, 3)}",
                    runtime=time.perf_counter() - start,
                )
            feasible, shrunk_sum = verify_and_score(packing)
            if not feasible:
                return Result(
                    False,
                    error=f"verifier rejected the packing for N={inst.n}",
                    runtime=time.perf_counter() - start,
                )
            gaps.append(gap_on(shrunk_sum, inst.bks))
    except SolverError as error:
        return Result(False, error=f"solver error: {error}", runtime=time.perf_counter() - start)
    except Exception as error:  # noqa: BLE001  any runtime fault is INVALID
        return Result(False, error=f"runtime error: {error}", runtime=time.perf_counter() - start)

    return Result(
        True,
        excess=float(np.mean(gaps)),
        runtime=time.perf_counter() - start,
        worst_gap=float(np.max(gaps)),
    )


# --- Seed solver (iteration 0): the naive baseline (ticket 02) --------------
# Square-grid centers, radii by the exact LP, then one polish. Deterministic and
# repeatable. It is the floor the evolved program must beat, and the evolution seed.
SEED_SOLVER_SOURCE = '''\
import math

import numpy as np


def pack(n, ctx):
    """Naive baseline: square-grid centers, exact radii by LP, then one polish.

    Place n centers on the smallest square grid that holds them, inset by half a
    cell so every center sits inside the unit square. The harness primitives set the
    largest feasible radii and refine locally. The genome only orchestrates; it
    never sets a radius by hand or scores.
    """
    k = int(math.ceil(math.sqrt(n)))
    step = 1.0 / k
    centers = np.empty((n, 2), dtype=float)
    filled = 0
    for row in range(k):
        for col in range(k):
            if filled == n:
                break
            centers[filled, 0] = (col + 0.5) * step
            centers[filled, 1] = (row + 0.5) * step
            filled += 1
    radii = ctx.radii_lp(centers)
    packing = np.column_stack([centers, radii])
    return ctx.polish(packing)
'''


def baseline_sum(n: int, total_budget: float = TOTAL_BUDGET_S, polish_maxiter: int = POLISH_MAXITER) -> float:
    """Run the seed solver for one N and return its feasibility-shrunk sum of radii.

    Used to record the baseline the evolved program must beat (ticket 02).
    """
    pack = load_solver(SEED_SOLVER_SOURCE)
    packing = run_pack(pack, n, time.perf_counter() + total_budget, polish_maxiter)
    feasible, shrunk_sum = verify_and_score(packing)
    if not feasible:
        raise SolverError(f"seed solver produced an INVALID packing for N={n}")
    return shrunk_sum


# --- Prompt building -------------------------------------------------------
SYSTEM_INSTRUCTION = '''\
You write one Python function for a variable-radii circle-packing (CSQV) solver.

The problem: place n non-overlapping circles of unequal radii inside the unit \
square [0, 1] x [0, 1] to MAXIMIZE the sum of the radii.

Contract:
- Define exactly one function: pack(n, ctx).
- Return a numpy array of shape (n, 3). Each row is [x, y, r]: a circle center \
(x, y) and its radius r, in the unit square. Circles must stay inside the square \
and must not overlap.
- pack is a whole search program. Initialize centers, perturb, restart, and refine \
to push the sum of radii up. Use ctx to search harder within the time budget.
- ctx.radii_lp(centers) takes an (n, 2) array of centers and returns the largest \
radii that keep the packing valid for those centers. Use it; do not set radii by hand.
- ctx.polish(packing) takes an (n, 3) packing and returns a locally refined one \
with a higher or equal sum. Call it on your best candidates.
- ctx.rng is a numpy random generator for reproducible restarts. ctx.time_left() \
returns the seconds left in the per-N budget; stop searching before it hits zero.
- A verifier recomputes feasibility, shrinks any tiny violation, and scores the \
shrunk sum. Non-finite or out-of-bounds output scores worst.

Output exactly one fenced python code block and no prose.
'''


def build_prompt(
    champion_source: str,
    champion_excess: float,
    baseline_excess: float,
    history: list[dict] | None = None,
    peers: list[str] | None = None,
) -> str:
    """Assemble the per-iteration prompt: champion source, peers, scoreboard, history.

    The fixed contract lives in SYSTEM_INSTRUCTION. This carries the changing
    content the search learns from. `peers` are other strong programs from the top-K
    population, shown so the model can combine ideas. Runtime rides in the history so
    the model learns to write fast programs; runtime never enters the objective. Gap
    is the mean fraction under best-known; lower is better, and a negative gap beats
    the record.
    """
    lines = [
        "You improve a CSQV circle-packing search program by evolution.",
        "Lower gap is better. Gap is the mean fraction under the best-known sum of "
        "radii across the target N. A negative gap beats the record.",
        "",
        f"Naive baseline gap (held-out, square grid): {baseline_excess:.4f}",
        f"Current champion gap (train): {champion_excess:.4f}",
        "",
        "Current champion source:",
        "```python",
        champion_source.strip(),
        "```",
        "",
    ]
    if peers:
        lines.append("Other strong programs so far (ideas to combine or improve on):")
        for peer in peers:
            lines.extend(["```python", peer.strip(), "```"])
        lines.append("")
    if history:
        lines.append("Recent attempts (most recent last):")
        for attempt in history:
            iteration = attempt.get("iteration", "?")
            runtime = attempt.get("runtime", 0.0)
            note = attempt.get("note", "")
            if attempt.get("valid"):
                gap = attempt.get("excess")
                gap_str = f"{gap:.4f}" if isinstance(gap, (int, float)) else "n/a"
                line = f"- iter {iteration}: gap {gap_str}, runtime {runtime:.3f}s"
            else:
                line = f"- iter {iteration}: INVALID ({attempt.get('error', '')}), runtime {runtime:.3f}s"
            lines.append(line + (f". {note}" if note else ""))
        lines.append("")
    lines.append(
        "Write a better pack function. Output exactly one fenced python block, no prose."
    )
    return "\n".join(lines)

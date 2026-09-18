"""Variable-radii circle packing (CSQV): the trusted primitives.

This module exposes the trusted numeric primitives for placing n non-overlapping
circles of unequal radii inside the unit square [0, 1] x [0, 1] to maximize the sum
of the radii. A caller supplies circle centers or a whole packing; this module
computes radii, refines a packing, and scores it. The primitives are:

- `radii_lp(centers)`: the exact largest radii for FIXED centers, by a linear
  program. It maximizes the sum of radii subject to containment (each circle inside
  the square) and pairwise non-overlap.
- `polish(packing)`: a local numeric refinement (SLSQP) of centers and radii that
  never worsens the feasibility-shrunk sum.
- `defect_move(packing, rng)`: a local move that removes the weakest circles and
  reinserts them into the largest empty holes, then refits radii and polishes.
- `verify_and_score(packing)` / `canonical_feasible(packing)`: an independent
  verifier that recomputes containment and non-overlap, applies a feasibility
  shrink, and scores (or emits) the strictly-valid packing.

The frame is the unit square [0, 1] x [0, 1]. TOL is the acceptance tolerance
(data.py). scipy is a hard dependency of this module (install:
`pip install --user --break-system-packages scipy`).
"""

from __future__ import annotations

import time

import numpy as np
import numpy.typing as npt
from scipy.optimize import linprog, minimize

from . import data

Array = npt.NDArray[np.float64]

# The unit square the caller and verifier use. The stored anchors live in the
# centered frame; data.py converts. TOL is the acceptance tolerance (data.py).
SIDE = data.SIDE
TOL = data.TOL

# The SLSQP iteration cap for the polish primitive, scaled down for large N inside
# polish so a single solve stays bounded.
POLISH_MAXITER = 200


class SolverError(Exception):
    """A malformed packing. It makes the candidate INVALID, never a crash."""


class SolverInterrupted(BaseException):
    """Raised INTO a running solver (from a SIGTERM handler) to stop it on a deadline
    while keeping its best feasible work. It subclasses BaseException, not Exception, so
    a solver's own `except Exception` cannot swallow the stop."""


# --- Fixed primitive: radii-LP ---------------------------------------------
def radii_lp(centers: Array) -> Array:
    """Return the largest radii that keep the packing valid for FIXED centers.

    This is a linear program: maximize the sum of radii subject to containment
    (each circle inside the unit square) and pairwise non-overlap
    (r_i + r_j <= dist(i, j)). It is trusted primitive code. A center outside the
    square yields a zero containment bound, so its radius is 0. On the rare solver
    failure it falls back to the containment bounds scaled to feasibility, so the
    returned radii never overlap.
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
    containment and non-overlap. It is trusted primitive code. It never worsens the
    score: it returns whichever of the input and the refined packing has the larger
    feasibility-shrunk sum, so a failed or unhelpful refinement is a no-op. A
    non-finite or wrong-shaped input is returned unchanged. SLSQP cannot be
    interrupted mid-solve, so the iteration cap is scaled down for large N to keep a
    single call bounded, and the deadline is checked before the call. The caller
    should test the time left before it calls polish.
    """
    return _polish_scored(packing, maxiter, deadline)[0]


def _polish_scored(
    packing: Array, maxiter: int = POLISH_MAXITER, deadline: float | None = None
) -> tuple[Array, float]:
    """`polish` plus the returned packing's feasibility-shrunk sum, so a caller that needs
    the score does not verify it again. A no-op return carries `-inf` so the caller skips
    it without a second verify."""
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
    """`defect_move` plus the returned packing's feasibility-shrunk sum, so a caller
    tracks the best packing without a second verify. A no-op or invalid input carries its
    own score (`-inf` when invalid), so the caller skips it cleanly.

    One move: remove the k weakest circles (smallest radius), reinsert each into the largest
    empty hole (a jittered coarse grid probed for the biggest admissible radius, greedily so
    later insertions see earlier ones), reset ALL radii by the exact LP, then one polish.
    Deterministic given (packing, rng state, k). Unlike a center perturbation it changes the
    combinatorial contact graph, so it can reach a basin the perturbation walk never visits.

    NOT monotonic: unlike `polish` (which never worsens the score), this removes and reinserts
    circles, so the refit-and-polish result CAN score below the input. A caller must compare and
    keep the better packing; a caller that treats the return value as unconditionally better than
    its input can go backward.
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
    refit radii and polish. A fixed trusted primitive (ADR 0001). See `_defect_move_scored`
    for the contract. Returns the input unchanged when it cannot improve the geometry
    (n <= 1, wrong shape, non-finite)."""
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

    Trust nothing about how the packing was produced. A packing is INVALID (returns
    (False, -inf)) when it is the wrong shape, non-finite, has a negative radius, has
    a center outside the square by more than `tol`, or has two coincident centers
    (uncorrectable by a shrink). Otherwise the verifier finds the largest uniform
    factor s in [0, 1] that makes every circle fit inside the square and no two
    overlap, then scores the shrunk sum s * sum(r). A packing already feasible scores
    its raw sum (s = 1); a slightly-violating one loses only a little (the
    numerical-artifact shrink); a badly-violating one is penalized hard, down to a
    valid zero (no INVALID cliff at s = 0). The score is always on a strictly-valid
    packing, so the objective cannot be gamed.
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
    packing (for example a .pck submission): the raw packing may hold tiny violations
    that a zero-tolerance external re-check would reject.
    """
    ok, s, x, y, r = _validate_and_shrink(packing, tol)
    if not ok:
        return None
    return np.column_stack([x, y, s * r])

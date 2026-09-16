"""Terminal second-order NLP squeeze for a single CSQV champion (ticket 45).

This is a ONE-TIME terminal post-processor. It takes the best candidate the search
found and drives it to the EXACT local optimum of its current basin, the deeper
same-basin optimum that Packomania's own re-optimization reaches after a submission.
It is NOT an inner-loop primitive, so cost is not a constraint here.

Why a second-order solve. The first-order conditioner `center_polish_analytic`
(cpp/csqv/polish.hpp) ascends the exact-LP value over the centers by its LP dual
gradient. It reaches a first-order stationary point but STALLS about 1 to 2e-5 short
of the optimum, because the jammed optimum sits AT a kink of the LP value function
where the two-sided derivative does not exist. Closing that last gap needs true
second-order steps. See docs/research/phase-2-csqv-terminal-postprocessing.md
(section A, ranked first) for the full derivation and the primary sources.

Method. State the full non-linear program in the 3N variables (x_i, y_i, r_i):
a LINEAR objective sum(r), 4N LINEAR wall (containment) constraints, and the
non-linear pairwise non-overlap constraints dist(c_i, c_j) >= r_i + r_j. All the
non-linearity sits in the pair distances, and each pair has a simple sparse 4x4
Hessian block in (x_i, y_i, x_j, y_j), so the Lagrangian Hessian is exact, analytic,
and sparse. scipy `trust-constr` takes the sparse Jacobian and the exact Lagrangian
Hessian, walks the interior toward the highly-active jammed optimum (so the kink that
stalls the first-order polish does not trip it), and self-identifies the active set.

The verifier stays the source of truth. This module returns a unit-frame packing;
the caller re-solves the exact radii with `radii_lp`, applies the feasibility shrink
with `canonical_feasible`, and validates the emitted .pck at zero tolerance. So a
step that raises a tiny overlap can never reach a submission.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint, linprog, minimize, nnls
from scipy.sparse import coo_matrix, csr_matrix

from . import data, problem

Array = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]

SIDE = problem.SIDE

# Wall ids for the four containment constraints of a circle. Each wall residual is >= 0 when
# the circle is inside the square; a residual of 0 means the circle is TANGENT to that wall.
_WALL_LEFT, _WALL_RIGHT, _WALL_BOTTOM, _WALL_TOP = 0, 1, 2, 3


def _containment_constraint(n: int) -> LinearConstraint:
    """The 4N linear wall constraints, as a sparse LinearConstraint on z = [x,y,r] * n.

    Per circle i the four rows keep it inside the square [0, SIDE]^2:
      x_i - r_i >= 0,  x_i + r_i <= SIDE,  y_i - r_i >= 0,  y_i + r_i <= SIDE.
    Written as bounded linear forms A z in [lb, ub], so trust-constr treats them as the
    exact linear constraints they are (no Hessian, no Jacobian recompute).
    """
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    lb = np.empty(4 * n, dtype=np.float64)
    ub = np.empty(4 * n, dtype=np.float64)
    inf = np.inf
    for i in range(n):
        xi, yi, ri = 3 * i, 3 * i + 1, 3 * i + 2
        # x_i - r_i in [0, inf]
        rows += [4 * i, 4 * i]
        cols += [xi, ri]
        vals += [1.0, -1.0]
        lb[4 * i], ub[4 * i] = 0.0, inf
        # x_i + r_i in [-inf, SIDE]
        rows += [4 * i + 1, 4 * i + 1]
        cols += [xi, ri]
        vals += [1.0, 1.0]
        lb[4 * i + 1], ub[4 * i + 1] = -inf, SIDE
        # y_i - r_i in [0, inf]
        rows += [4 * i + 2, 4 * i + 2]
        cols += [yi, ri]
        vals += [1.0, -1.0]
        lb[4 * i + 2], ub[4 * i + 2] = 0.0, inf
        # y_i + r_i in [-inf, SIDE]
        rows += [4 * i + 3, 4 * i + 3]
        cols += [yi, ri]
        vals += [1.0, 1.0]
        lb[4 * i + 3], ub[4 * i + 3] = -inf, SIDE
    a = coo_matrix((vals, (rows, cols)), shape=(4 * n, 3 * n)).tocsr()
    return LinearConstraint(a, lb, ub)


def _overlap_constraint(n: int, i: IntArray, j: IntArray) -> NonlinearConstraint:
    """The pairwise non-overlap constraints dist(c_i, c_j) - r_i - r_j >= 0.

    Supplies the SPARSE Jacobian and the EXACT sparse Lagrangian Hessian. The Hessian
    of the constraint dot product v . c(z) is a sum of per-pair 4x4 blocks over the
    center variables (the radius part is linear, so it adds nothing). Each block is
    v_p * J^T (1/d)(I - uu^T) J with u the unit contact direction and J the map from
    (x_i, y_i, x_j, y_j) to (dx, dy); it is exactly the Hessian of the pair distance.
    """
    m = i.size
    # Constant Jacobian sparsity: 6 nonzeros per row (x_i, y_i, x_j, y_j, r_i, r_j).
    jrows = np.repeat(np.arange(m), 6)
    jcols = np.empty(6 * m, dtype=np.int64)
    jcols[0::6] = 3 * i      # x_i
    jcols[1::6] = 3 * i + 1  # y_i
    jcols[2::6] = 3 * j      # x_j
    jcols[3::6] = 3 * j + 1  # y_j
    jcols[4::6] = 3 * i + 2  # r_i
    jcols[5::6] = 3 * j + 2  # r_j

    # Constant Hessian sparsity: the 4x4 center block per pair, i.e. 16 entries.
    hidx = np.empty((m, 4), dtype=np.int64)
    hidx[:, 0] = 3 * i      # x_i
    hidx[:, 1] = 3 * i + 1  # y_i
    hidx[:, 2] = 3 * j      # x_j
    hidx[:, 3] = 3 * j + 1  # y_j
    hrows = np.repeat(hidx, 4, axis=1).ravel()          # each of 4 rows over 4 cols
    hcols = np.tile(hidx, (1, 4)).ravel()

    def fun(z: Array) -> Array:
        p = z.reshape(n, 3)
        x, y, r = p[:, 0], p[:, 1], p[:, 2]
        d = np.hypot(x[i] - x[j], y[i] - y[j])
        return d - r[i] - r[j]

    def jac(z: Array) -> csr_matrix:
        p = z.reshape(n, 3)
        x, y = p[:, 0], p[:, 1]
        dx = x[i] - x[j]
        dy = y[i] - y[j]
        d = np.hypot(dx, dy)
        d = np.where(d > 0.0, d, 1.0)  # guard; a valid champion has no coincident pair
        ux, uy = dx / d, dy / d
        vals = np.empty(6 * m, dtype=np.float64)
        vals[0::6] = ux
        vals[1::6] = uy
        vals[2::6] = -ux
        vals[3::6] = -uy
        vals[4::6] = -1.0
        vals[5::6] = -1.0
        return coo_matrix((vals, (jrows, jcols)), shape=(m, 3 * n)).tocsr()

    def hess(z: Array, v: Array) -> csr_matrix:
        p = z.reshape(n, 3)
        x, y = p[:, 0], p[:, 1]
        dx = x[i] - x[j]
        dy = y[i] - y[j]
        d = np.hypot(dx, dy)
        d = np.where(d > 0.0, d, 1.0)
        inv = 1.0 / d
        d3 = inv / (d * d)
        # M = (1/d) I - (1/d^3) [[dx^2, dx dy],[dx dy, dy^2]], the Hessian of d in (dx, dy).
        m00 = inv - dx * dx * d3
        m01 = -dx * dy * d3
        m11 = inv - dy * dy * d3
        # 4x4 center block = [[M, -M], [-M, M]], scaled by the multiplier v_p.
        blk = np.empty((m, 4, 4), dtype=np.float64)
        blk[:, 0, 0], blk[:, 0, 1], blk[:, 1, 0], blk[:, 1, 1] = m00, m01, m01, m11
        blk[:, 2, 2], blk[:, 2, 3], blk[:, 3, 2], blk[:, 3, 3] = m00, m01, m01, m11
        blk[:, 0, 2], blk[:, 0, 3], blk[:, 1, 2], blk[:, 1, 3] = -m00, -m01, -m01, -m11
        blk[:, 2, 0], blk[:, 3, 0], blk[:, 2, 1], blk[:, 3, 1] = -m00, -m01, -m01, -m11
        blk *= v[:, None, None]
        return coo_matrix((blk.ravel(), (hrows, hcols)), shape=(3 * n, 3 * n)).tocsr()

    return NonlinearConstraint(fun, 0.0, np.inf, jac=jac, hess=hess)


def _near_contact_pairs(centers: Array, radii: Array, margin: float) -> tuple[IntArray, IntArray]:
    """Upper-triangle pair indices (i, j) whose non-overlap constraint can bind.

    A pair with distance d and radii r_i, r_j is near contact when d - (r_i + r_j) < margin.
    A pair far outside that band cannot become tight under the tiny center motion of a
    same-basin terminal squeeze, so it is dropped from the NLP. Correctness does NOT depend
    on this prune: the final `radii_lp` and `canonical_feasible` re-check ALL pairs and
    shrink, so a dropped pair that somehow approaches contact is still caught before emit.
    The prune only bounds the constraint count (~3N near contacts, not O(N^2)), which keeps
    the second-order solve fast and local.
    """
    x, y = centers[:, 0], centers[:, 1]
    i, j, d = data.pairwise_distances(x, y)
    keep = (d - (radii[i] + radii[j])) < margin
    return i[keep].astype(np.int64), j[keep].astype(np.int64)


def terminal_squeeze(
    packing_unit: Array,
    maxiter: int = 1000,
    gtol: float = 1e-12,
    xtol: float = 1e-14,
    trust_box: float = 1e-2,
    verbose: int = 0,
) -> Array:
    """Squeeze a unit-frame champion to its exact same-basin local optimum.

    `packing_unit` is an (n, 3) array of [x, y, r] rows in the unit square. Returns the
    refined (n, 3) unit-frame packing. It is MONOTONE by construction at the call site:
    the caller re-solves radii with `radii_lp` and keeps the better of input and output,
    so a solver that fails to improve is a no-op. On any solver fault it returns the
    input unchanged.

    `trust_box` bounds each center to a small box of that half-width around the champion.
    That encodes the SAME-BASIN intent explicitly: the true terminal motion is tiny
    (~1e-5), so the box gives ample headroom yet stops the interior-point solve from
    wandering into a different arrangement. It also makes the near-contact prune SAFE: a
    pair farther apart than 2 * trust_box cannot reach contact inside the boxes, so it can
    be dropped from the NLP. The final `radii_lp` + `canonical_feasible` re-check every
    pair regardless, so the prune and the box cost no correctness and keep the solve fast.
    """
    packing_unit = np.asarray(packing_unit, dtype=np.float64)
    n = packing_unit.shape[0]
    if packing_unit.shape != (n, 3) or not np.isfinite(packing_unit).all() or n == 0:
        return packing_unit
    if n < 2:
        return packing_unit  # a single circle is already at its wall-slack optimum

    # Seed the radii from the exact LP so the start is strictly FEASIBLE and near-optimal,
    # whatever radii the caller passed. An infeasible start (stale radii from a moved
    # champion) makes the interior-point solve stall; the LP radii are the right warm start.
    centers = packing_unit[:, :2]
    r_lp = problem.radii_lp(centers)
    z0 = np.empty(3 * n, dtype=np.float64)
    z0[0::3] = centers[:, 0]
    z0[1::3] = centers[:, 1]
    z0[2::3] = r_lp

    # Include any pair that could reach contact inside the trust boxes (both centers move
    # up to trust_box toward each other, so 2 * trust_box of closing), plus a safety band.
    pair_margin = 2.0 * trust_box + 1e-3
    i, j = _near_contact_pairs(centers, r_lp, pair_margin)
    if i.size == 0:
        return packing_unit  # no active pairs to squeeze against

    grad = np.zeros(3 * n, dtype=np.float64)
    grad[2::3] = -1.0  # minimize -sum(r)
    zero_hess = csr_matrix((3 * n, 3 * n))

    def neg_sum(z: Array) -> float:
        return -float(z[2::3].sum())

    def neg_sum_jac(_z: Array) -> Array:
        return grad

    def neg_sum_hess(_z: Array) -> csr_matrix:
        return zero_hess

    # Per-center trust box: x, y within trust_box of the champion (clipped to the square);
    # radii free below (containment caps them from above). This keeps the solve same-basin.
    lb = np.empty(3 * n, dtype=np.float64)
    ub = np.empty(3 * n, dtype=np.float64)
    lb[0::3] = np.maximum(centers[:, 0] - trust_box, 0.0)
    ub[0::3] = np.minimum(centers[:, 0] + trust_box, SIDE)
    lb[1::3] = np.maximum(centers[:, 1] - trust_box, 0.0)
    ub[1::3] = np.minimum(centers[:, 1] + trust_box, SIDE)
    lb[2::3], ub[2::3] = 0.0, np.inf
    bounds = Bounds(lb, ub)
    constraints = [_containment_constraint(n), _overlap_constraint(n, i, j)]

    try:
        res = minimize(
            neg_sum,
            z0,
            method="trust-constr",
            jac=neg_sum_jac,
            hess=neg_sum_hess,
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": maxiter, "gtol": gtol, "xtol": xtol,
                     "initial_tr_radius": trust_box, "verbose": verbose},
        )
    except Exception:  # noqa: BLE001  a numeric failure is a no-op, never a crash
        return packing_unit

    out = np.asarray(res.x, dtype=np.float64).reshape(n, 3)
    if not np.isfinite(out).all():
        return packing_unit
    return out


def _exact_score(centers: Array) -> float:
    """The exact-LP feasibility-shrunk sum for fixed centers, the true CSQV objective."""
    r = problem.radii_lp(centers)
    _, s = problem.verify_and_score(np.column_stack([centers, r]))
    return s


def squeeze_champion(
    packing_unit: Array,
    maxiter: int = 1000,
    hops: int = 0,
    trust_box: float = 1e-2,
    escape_scales: tuple[float, ...] = (3e-6, 1e-5, 3e-5),
    escape_box: float = 1e-3,
    patience: int = 60,
    seed: int = 0,
    newton: bool = True,
    verbose: int = 0,
) -> tuple[Array, float]:
    """Run the terminal squeeze pipeline on a champion; return (canonical, score).

    Pipeline (research note section "Recommended terminal pipeline"):
      1. second-order NLP squeeze over (x, y, r) via `terminal_squeeze`;
      2. exact radii: re-solve `radii_lp` on the squeezed centers;
      3. feasibility shrink + verify: `canonical_feasible` gives the strictly-valid
         packing the verifier scores, so the result at most touches.

    A near-optimal champion sits AT the degenerate kink of the LP value function, which is
    an approximate stationary point of the smooth NLP, so the direct second-order step
    alone recovers only the first-order optimum. To cross the kink to the deeper same-basin
    optimum, `hops` > 0 runs a bounded MONOTONIC BASIN HOP: perturb the current best by a
    tiny `escape_scales` jitter, re-run the second-order squeeze in a small `escape_box`,
    and keep the result only if it scores higher. `seed` makes the hop deterministic;
    `patience` stops early after that many hops with no gain.

    `newton` runs the rank-2 guarded contact-graph Newton (`contact_newton`) as a final
    precision step. It sharpens a JAMMED packing to machine precision; on a champion that is
    not yet jammed (Donev margin above its guard) it is a no-op, so it only ever helps.

    The returned packing is the BEST seen, so the score never decreases (monotone: safe
    before any compare or submit).
    """
    packing_unit = np.asarray(packing_unit, dtype=np.float64)
    centers = packing_unit[:, :2].copy()

    best_centers = centers
    best_sum = _exact_score(centers)  # no-op floor: a non-improving squeeze cannot regress

    squeezed = terminal_squeeze(packing_unit, maxiter=maxiter, trust_box=trust_box, verbose=verbose)
    s_sq = _exact_score(squeezed[:, :2])
    if s_sq > best_sum:
        best_centers, best_sum = squeezed[:, :2].copy(), s_sq

    if hops > 0:
        rng = np.random.default_rng(seed)
        n = centers.shape[0]
        no_gain = 0
        for _ in range(hops):
            eps = float(rng.choice(escape_scales))
            trial = best_centers + rng.normal(0.0, eps, size=(n, 2))
            trial = np.clip(trial, 0.0, SIDE)
            # terminal_squeeze re-seeds radii from radii_lp, so a placeholder column suffices
            # (avoids a redundant LP solve on this hot path).
            packing = np.column_stack([trial, np.zeros(n)])
            hopped = terminal_squeeze(packing, maxiter=maxiter, trust_box=escape_box)
            s = _exact_score(hopped[:, :2])
            if s > best_sum + 1e-12:
                best_centers, best_sum = hopped[:, :2].copy(), s
                no_gain = 0
            else:
                no_gain += 1
                if no_gain >= patience:
                    break

    if newton:
        # contact_newton re-derives radii from the centers, so a placeholder column suffices.
        polished = contact_newton(np.column_stack([best_centers, np.zeros(best_centers.shape[0])]))
        s = _exact_score(polished[:, :2])
        if s > best_sum:
            best_centers, best_sum = polished[:, :2].copy(), s

    r = problem.radii_lp(best_centers)
    canon = problem.canonical_feasible(np.column_stack([best_centers, r]))
    assert canon is not None  # best_centers scored finite above, so it is always feasible
    return canon, best_sum


# --- Rank-2: guarded contact-graph KKT-Newton (research note section B) ------
#
# At a jammed CSQV optimum the tight contacts and wall tangencies form a rigid framework.
# The optimum is the KKT point of "maximize sum(r) s.t. active contacts hold with equality":
#   stationarity   grad(-sum r) - A^T nu = 0        (3N equations, nu >= 0 the multipliers)
#   primal         c_active(z) = 0                   (M equations, the tangency system)
# with unknowns z = [x, y, r] (3N) and nu (M). This is a SQUARE system on which Newton
# converges quadratically once the active set is correct. The KKT Jacobian is the saddle
#   J = [ H   -A^T ]   with  H = -sum_p nu_p Hess(dist_p)  (the exact sparse Lagrangian
#       [ A    0   ]        Hessian, reused from the second-order squeeze).
#
# The fragility is FALSE CONTACTS: a pair ~1e-9 apart that does not carry load. A wrong
# active set gives a singular or wrong Newton system. Three guards handle it:
#   1. drop non-load-bearing contacts by complementary slackness (nnls multiplier ~ 0);
#   2. the Donev LP jamming test certifies the pruned contact set is rigid before Newton;
#   3. the caller keeps the better of input and output, and the terminal radii_lp + shrink
#      + zero-tol validate is the hard feasibility backstop, so a bad step is only wasted
#      work, never a regression or an infeasible emit.


def _active_set(
    centers: Array, radii: Array, pair_tol: float, wall_tol: float
) -> tuple[IntArray, IntArray]:
    """Detect the active contacts and wall tangencies at the given packing.

    Returns (pairs, walls): `pairs` is a (P, 2) array of contacting circle index pairs
    (surface gap below `pair_tol`); `walls` is a (W, 2) array of [circle index, wall id]
    for each wall the circle is tangent to (residual below `wall_tol`).
    """
    x, y = centers[:, 0], centers[:, 1]
    n = centers.shape[0]
    i, j, d = data.pairwise_distances(x, y)
    keep = (d - (radii[i] + radii[j])) < pair_tol
    pairs = np.column_stack([i[keep], j[keep]]).astype(np.int64)

    residuals = np.stack([x - radii, SIDE - x - radii, y - radii, SIDE - y - radii], axis=1)
    wi, wid = np.nonzero(residuals < wall_tol)
    walls = np.column_stack([wi, wid]).astype(np.int64)
    return pairs, walls


def _active_jacobian(centers: Array, radii: Array, pairs: IntArray, walls: IntArray) -> csr_matrix:
    """The Jacobian A (M x 3N) of the active constraints, rows = pairs then walls.

    A pair row is the gradient of dist(i, j) - r_i - r_j; a wall row is the gradient of the
    tangent wall residual (x_i - r_i, SIDE - x_i - r_i, y_i - r_i, or SIDE - y_i - r_i).
    """
    n = centers.shape[0]
    x, y = centers[:, 0], centers[:, 1]
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    for row, (i, j) in enumerate(pairs):
        dx, dy = x[i] - x[j], y[i] - y[j]
        d = np.hypot(dx, dy)
        d = d if d > 0.0 else 1.0
        ux, uy = dx / d, dy / d
        rows += [row] * 6
        cols += [3 * i, 3 * i + 1, 3 * j, 3 * j + 1, 3 * i + 2, 3 * j + 2]
        vals += [ux, uy, -ux, -uy, -1.0, -1.0]

    p = pairs.shape[0]
    for w, (i, wid) in enumerate(walls):
        row = p + w
        if wid == _WALL_LEFT:
            rows += [row, row]; cols += [3 * i, 3 * i + 2]; vals += [1.0, -1.0]
        elif wid == _WALL_RIGHT:
            rows += [row, row]; cols += [3 * i, 3 * i + 2]; vals += [-1.0, -1.0]
        elif wid == _WALL_BOTTOM:
            rows += [row, row]; cols += [3 * i + 1, 3 * i + 2]; vals += [1.0, -1.0]
        else:  # _WALL_TOP
            rows += [row, row]; cols += [3 * i + 1, 3 * i + 2]; vals += [-1.0, -1.0]

    m = pairs.shape[0] + walls.shape[0]
    return coo_matrix((vals, (rows, cols)), shape=(m, 3 * n)).tocsr()


def _active_residual(centers: Array, radii: Array, pairs: IntArray, walls: IntArray) -> Array:
    """The active-constraint values c_a (M,), rows = pairs then walls. Zero at a jammed point."""
    x, y = centers[:, 0], centers[:, 1]
    out = np.empty(pairs.shape[0] + walls.shape[0], dtype=np.float64)
    if pairs.shape[0]:
        i, j = pairs[:, 0], pairs[:, 1]
        out[: pairs.shape[0]] = np.hypot(x[i] - x[j], y[i] - y[j]) - radii[i] - radii[j]
    for w, (i, wid) in enumerate(walls):
        if wid == _WALL_LEFT:
            out[pairs.shape[0] + w] = x[i] - radii[i]
        elif wid == _WALL_RIGHT:
            out[pairs.shape[0] + w] = SIDE - x[i] - radii[i]
        elif wid == _WALL_BOTTOM:
            out[pairs.shape[0] + w] = y[i] - radii[i]
        else:
            out[pairs.shape[0] + w] = SIDE - y[i] - radii[i]
    return out


def donev_jamming_margin(centers: Array, radii: Array, pairs: IntArray, walls: IntArray) -> float:
    """The Donev LP jamming margin: the best first-order objective gain from an unjamming flex.

    Linearize the active contacts about the packing and solve
        maximize  sum(delta_r)   s.t.  A . delta_z >= 0,  -1 <= delta_z <= 1,
    where A . delta_z >= 0 keeps every active contact non-overlapping to first order (Donev,
    Torquato, Stillinger, Connelly, J. Comput. Phys. 197 (2004) 139). A margin near 0 means
    no feasible flex increases the radii: the contact set is RIGID (collectively jammed), so
    the Newton target is a true jammed optimum. A clearly positive margin means an improving
    motion exists that the active set does not yet capture, so Newton on it is not trustworthy.
    Returns the margin (>= 0), or +inf if the LP fails (treated as not certified).
    """
    n = centers.shape[0]
    a = _active_jacobian(centers, radii, pairs, walls)
    c = np.zeros(3 * n, dtype=np.float64)
    c[2::3] = -1.0  # maximize sum(delta_r) == minimize -sum(delta_r)
    res = linprog(
        c=c,
        A_ub=(-a).toarray(),
        b_ub=np.zeros(a.shape[0]),
        bounds=[(-1.0, 1.0)] * (3 * n),
        method="highs",
    )
    if not res.success or res.x is None:
        return np.inf
    return float(res.x[2::3].sum())


def contact_newton(
    packing_unit: Array,
    max_iter: int = 30,
    pair_tol: float = 1e-7,
    wall_tol: float = 1e-7,
    load_tol: float = 1e-9,
    jam_tol: float = 1e-6,
    ftol: float = 1e-15,
) -> Array:
    """Sharpen a jammed champion to machine precision by KKT-Newton on its contact graph.

    Detects the active set, drops non-load-bearing (false) contacts by complementary
    slackness, certifies the pruned set is rigid with the Donev jamming test, then takes
    guarded Newton steps on the KKT system. Returns the refined (n, 3) unit-frame packing,
    or the input unchanged when it cannot certify or improve (a no-op is safe: the caller
    keeps the better of input and output).
    """
    packing_unit = np.asarray(packing_unit, dtype=np.float64)
    n = packing_unit.shape[0]
    if packing_unit.shape != (n, 3) or not np.isfinite(packing_unit).all() or n < 2:
        return packing_unit

    centers = packing_unit[:, :2].copy()
    radii = problem.radii_lp(centers)
    z = np.empty(3 * n, dtype=np.float64)
    z[0::3], z[1::3], z[2::3] = centers[:, 0], centers[:, 1], radii

    pairs, walls = _active_set(centers, radii, pair_tol, wall_tol)
    if pairs.shape[0] + walls.shape[0] == 0:
        return packing_unit

    # Stationarity is linear in the multipliers: A^T nu = grad(-sum r). Solve it non-negative
    # and drop contacts that carry no load (nu ~ 0) as false contacts (complementary slackness).
    a = _active_jacobian(centers, radii, pairs, walls)
    grad0 = np.zeros(3 * n, dtype=np.float64)
    grad0[2::3] = -1.0
    nu, _ = nnls(a.T.toarray(), grad0)
    load = nu > load_tol * max(1.0, nu.max())
    npairs = pairs.shape[0]
    pairs = pairs[load[:npairs]]
    walls = walls[load[npairs:]]
    nu = nu[load]
    if pairs.shape[0] + walls.shape[0] == 0:
        return packing_unit

    # Guard: certify the pruned contact set is rigid before trusting a Newton solve on it.
    if donev_jamming_margin(centers, radii, pairs, walls) > jam_tol:
        return packing_unit

    # Guarded Newton on F = [grad(-sum r) - A^T nu ; c_active] = 0.
    w = np.concatenate([z, nu])
    for _ in range(max_iter):
        zc = w[:3 * n]
        cur_centers = np.column_stack([zc[0::3], zc[1::3]])
        cur_radii = zc[2::3]
        nu = w[3 * n:]
        a = _active_jacobian(cur_centers, cur_radii, pairs, walls)
        resid = _active_residual(cur_centers, cur_radii, pairs, walls)
        f_stat = grad0 - a.T @ nu
        f = np.concatenate([f_stat, resid])
        if np.linalg.norm(f) < ftol:
            break

        # KKT Jacobian J = [[H, -A^T], [A, 0]] with H = -sum_p nu_p Hess(dist_p).
        h = _pair_hessian(cur_centers, pairs, nu[: pairs.shape[0]], n)
        m = a.shape[0]
        jac = np.zeros((3 * n + m, 3 * n + m), dtype=np.float64)
        jac[:3 * n, :3 * n] = -h
        jac[:3 * n, 3 * n:] = -a.T.toarray()
        jac[3 * n:, :3 * n] = a.toarray()
        try:
            step = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            return packing_unit  # singular KKT: bail, the caller keeps the input

        # Fraction-to-boundary so the multipliers stay non-negative, then a damped update.
        dnu = step[3 * n:]
        alpha = 1.0
        neg = dnu < 0.0
        if neg.any():
            alpha = min(1.0, 0.99 * float(np.min(-nu[neg] / dnu[neg])))
        w = w + alpha * step
        if not np.isfinite(w).all():
            return packing_unit

    out = np.column_stack([w[:3 * n][0::3], w[:3 * n][1::3], w[:3 * n][2::3]])
    if not np.isfinite(out).all():
        return packing_unit
    return out


def _pair_hessian(centers: Array, pairs: IntArray, nu_pairs: Array, n: int) -> Array:
    """H = sum_p nu_p Hess(dist_p) as a dense (3N x 3N) matrix (walls are linear, no Hessian).

    Returns +sum nu_p Hess(dist_p); the caller negates it for the Lagrangian Hessian block.
    """
    h = np.zeros((3 * n, 3 * n), dtype=np.float64)
    x, y = centers[:, 0], centers[:, 1]
    for p, (i, j) in enumerate(pairs):
        dx, dy = x[i] - x[j], y[i] - y[j]
        d = np.hypot(dx, dy)
        if d <= 0.0:
            continue
        inv = 1.0 / d
        d3 = inv / (d * d)
        m00 = inv - dx * dx * d3
        m01 = -dx * dy * d3
        m11 = inv - dy * dy * d3
        block = np.array([[m00, m01], [m01, m11]]) * nu_pairs[p]
        idx = [3 * i, 3 * i + 1]
        jdx = [3 * j, 3 * j + 1]
        for a_local, ai in enumerate(idx):
            for b_local, bj in enumerate(idx):
                h[ai, bj] += block[a_local, b_local]
            for b_local, bj in enumerate(jdx):
                h[ai, bj] -= block[a_local, b_local]
        for a_local, ai in enumerate(jdx):
            for b_local, bj in enumerate(idx):
                h[ai, bj] -= block[a_local, b_local]
            for b_local, bj in enumerate(jdx):
                h[ai, bj] += block[a_local, b_local]
    return h

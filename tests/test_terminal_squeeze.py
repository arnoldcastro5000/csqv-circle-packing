"""Tests for the terminal second-order NLP squeeze (ticket 45). No model access; scipy required.

These lock the mechanics of the one-time terminal post-processor: the exact sparse Jacobian
and Lagrangian Hessian of the non-overlap constraints, same-basin recovery of a known optimum,
monotonicity (it never regresses a champion), strict feasibility of the emitted packing, and
the small edge cases. The heavy N=142/143 acceptance run is a separate opt-in script
(cpp/tools/terminal_squeeze_accept.py), too slow for the unit suite.
"""

from __future__ import annotations

import numpy as np
import pytest

from problems.csqv import data, problem, terminal_squeeze as ts


def _exact_score(centers: np.ndarray) -> float:
    r = problem.radii_lp(centers)
    _, s = problem.verify_and_score(np.column_stack([centers, r]))
    return s


# --- exact derivatives of the non-overlap constraints -----------------------
def test_overlap_jacobian_matches_finite_difference() -> None:
    rng = np.random.default_rng(0)
    n = 5
    z = rng.uniform(0.2, 0.8, size=3 * n)
    z[2::3] = rng.uniform(0.02, 0.08, n)
    i, j = np.triu_indices(n, k=1)
    nc = ts._overlap_constraint(n, i.astype(np.int64), j.astype(np.int64))
    analytic = nc.jac(z).toarray()

    h = 1e-7
    f0 = nc.fun(z)
    numeric = np.zeros_like(analytic)
    for k in range(3 * n):
        zp = z.copy()
        zp[k] += h
        numeric[:, k] = (nc.fun(zp) - f0) / h
    assert np.abs(analytic - numeric).max() < 1e-5


def test_overlap_hessian_matches_finite_difference_and_is_symmetric() -> None:
    rng = np.random.default_rng(1)
    n = 5
    z = rng.uniform(0.2, 0.8, size=3 * n)
    z[2::3] = rng.uniform(0.02, 0.08, n)
    i, j = np.triu_indices(n, k=1)
    i, j = i.astype(np.int64), j.astype(np.int64)
    v = rng.uniform(0.1, 1.0, size=i.size)
    nc = ts._overlap_constraint(n, i, j)
    hess = nc.hess(z, v).toarray()
    assert np.abs(hess - hess.T).max() == 0.0  # assembled symmetric by construction

    # The Hessian of dot(fun, v) is the derivative of jac(z)^T v.
    def grad(zz: np.ndarray) -> np.ndarray:
        return ts._overlap_constraint(n, i, j).jac(zz).toarray().T @ v

    h = 1e-7
    g0 = grad(z)
    numeric = np.zeros((3 * n, 3 * n))
    for k in range(3 * n):
        zp = z.copy()
        zp[k] += h
        numeric[:, k] = (grad(zp) - g0) / h
    assert np.abs(hess - numeric).max() < 1e-4


# --- same-basin recovery ----------------------------------------------------
def test_recovers_a_perturbed_optimum() -> None:
    """A small perturbation off a stored optimum is squeezed back to nearly its sum."""
    unit = data.centered_to_unit(data.load_packing_centered(20))
    optimum = _exact_score(unit[:, :2])
    rng = np.random.default_rng(3)
    perturbed = unit.copy()
    perturbed[:, :2] += rng.normal(0.0, 1e-4, size=(20, 2))
    perturbed_sum = _exact_score(perturbed[:, :2])

    squeezed = ts.terminal_squeeze(perturbed, trust_box=1e-3)
    squeezed_sum = _exact_score(squeezed[:, :2])

    # It recovers the great majority of the gap the perturbation opened.
    recovered = (squeezed_sum - perturbed_sum) / (optimum - perturbed_sum)
    assert recovered > 0.9


# --- monotonicity and feasibility of the full pipeline ----------------------
def test_squeeze_champion_never_regresses() -> None:
    unit = data.centered_to_unit(data.load_packing_centered(13))
    base = _exact_score(unit[:, :2])
    _, squeezed = ts.squeeze_champion(unit, hops=0)
    assert squeezed >= base - 1e-12


def test_squeeze_champion_emits_strictly_feasible_packing() -> None:
    unit = data.centered_to_unit(data.load_packing_centered(20))
    rng = np.random.default_rng(5)
    unit[:, :2] += rng.normal(0.0, 1e-4, size=(20, 2))
    canon, _ = ts.squeeze_champion(unit, hops=3, seed=0)

    x, y, r = canon[:, 0], canon[:, 1], canon[:, 2]
    i, j, d = data.pairwise_distances(x, y)
    max_overlap = float(((r[i] + r[j]) - d).max())
    max_wall = float((r - np.maximum(data.wall_slack(x, y), 0.0)).max())
    assert max_overlap <= 0.0  # strictly non-overlapping at zero tolerance
    assert max_wall <= 0.0     # strictly contained at zero tolerance


def test_basin_hop_is_deterministic_given_a_seed() -> None:
    unit = data.centered_to_unit(data.load_packing_centered(13))
    _, a = ts.squeeze_champion(unit, hops=5, seed=7)
    _, b = ts.squeeze_champion(unit, hops=5, seed=7)
    assert a == b


# --- rank-2: guarded contact-graph Newton -----------------------------------
def test_active_jacobian_matches_finite_difference() -> None:
    unit = data.centered_to_unit(data.load_packing_centered(20))
    centers, radii = unit[:, :2], problem.radii_lp(unit[:, :2])
    pairs, walls = ts._active_set(centers, radii, 1e-7, 1e-7)
    a = ts._active_jacobian(centers, radii, pairs, walls).toarray()

    z = np.empty(3 * 20)
    z[0::3], z[1::3], z[2::3] = centers[:, 0], centers[:, 1], radii

    def resid(zz: np.ndarray) -> np.ndarray:
        c = np.column_stack([zz[0::3], zz[1::3]])
        return ts._active_residual(c, zz[2::3], pairs, walls)

    h = 1e-7
    r0 = resid(z)
    numeric = np.zeros_like(a)
    for k in range(3 * 20):
        zp = z.copy()
        zp[k] += h
        numeric[:, k] = (resid(zp) - r0) / h
    assert np.abs(a - numeric).max() < 1e-5


def test_donev_margin_is_zero_at_a_jammed_optimum() -> None:
    unit = data.centered_to_unit(data.load_packing_centered(32))
    centers, radii = unit[:, :2], problem.radii_lp(unit[:, :2])
    pairs, walls = ts._active_set(centers, radii, 1e-7, 1e-7)
    margin = ts.donev_jamming_margin(centers, radii, pairs, walls)
    assert margin < 1e-9  # a jammed packing admits no unjamming flex


def test_donev_margin_is_positive_when_unjammed() -> None:
    unit = data.centered_to_unit(data.load_packing_centered(32))
    rng = np.random.default_rng(4)
    centers = unit[:, :2] + rng.normal(0.0, 3e-5, size=(32, 2))
    radii = problem.radii_lp(centers)
    pairs, walls = ts._active_set(centers, radii, 1e-7, 1e-7)
    margin = ts.donev_jamming_margin(centers, radii, pairs, walls)
    assert margin > 1e-3  # an unjammed packing has an improving flex


def test_contact_newton_sharpens_a_near_jammed_packing() -> None:
    """A jammed optimum nudged by 1e-9 keeps its contact graph; Newton restores tangency."""
    unit = data.centered_to_unit(data.load_packing_centered(32))
    rng = np.random.default_rng(9)
    centers = unit[:, :2] + rng.normal(0.0, 1e-9, size=(32, 2))
    base = _exact_score(centers)

    polished = ts.contact_newton(np.column_stack([centers, problem.radii_lp(centers)]))
    pc, pr = polished[:, :2], problem.radii_lp(polished[:, :2])
    pairs, walls = ts._active_set(pc, pr, 1e-7, 1e-7)
    residual = np.linalg.norm(ts._active_residual(pc, pr, pairs, walls))
    assert residual < 1e-13                    # sharpened to near machine precision
    assert _exact_score(pc) >= base - 1e-12    # never regresses


def test_contact_newton_noops_on_an_unjammed_input() -> None:
    """The Donev guard refuses a non-jammed packing, so Newton returns it unchanged."""
    unit = data.centered_to_unit(data.load_packing_centered(32))
    rng = np.random.default_rng(4)
    centers = unit[:, :2] + rng.normal(0.0, 3e-5, size=(32, 2))
    packing = np.column_stack([centers, problem.radii_lp(centers)])
    out = ts.contact_newton(packing)
    assert np.array_equal(out, packing)


# --- edge cases -------------------------------------------------------------
def test_single_circle_is_a_noop() -> None:
    unit = np.array([[0.5, 0.5, 0.5]], dtype=float)
    out = ts.terminal_squeeze(unit)
    assert np.array_equal(out, unit)


def test_non_finite_input_is_returned_unchanged() -> None:
    bad = np.array([[0.5, 0.5, np.nan], [0.3, 0.3, 0.1]], dtype=float)
    out = ts.terminal_squeeze(bad)
    assert out is bad or np.array_equal(out[np.isfinite(out)], bad[np.isfinite(bad)])

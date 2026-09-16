"""CSQV data layer: Packomania best-known values, coordinate anchors, and the
reference cross-check that locks the acceptance convention.

This module anchors the Phase-2b CSQV problem (wayfinder ticket "Anchor the CSQV
data, best-known values, and acceptance tolerance"). It does NOT build the problem
module; that is a separate ticket. It provides the parts a later verifier and
scorer reuse:

  - the best-known sum-of-radii values (Packomania, complete to N=100),
  - a loader for the anchor coordinate files,
  - the frame convention (the centered side-1 square the files use, versus the unit
    square the genome and verifier use),
  - a reference verifier that recomputes containment, non-overlap, and the sum.

Convention (pinned from the Packomania coordinate files, fetched 2026-09-12; see
data/README.md):

  - Storage frame: centered origin, the square [-0.5, 0.5] x [-0.5, 0.5], side 1. A
    boundary circle satisfies abs(x) + r = 0.5 exactly (for example csqv26 circle
    4: 0.415360499304 + 0.084639500696 = 0.5).
  - Genome and verifier frame (CONTEXT.md, ticket "Lock the CSQV genome, verifier,
    and fixed-primitive contract"): the unit square [0, 1] x [0, 1]. Convert by
    adding 0.5 to x and y; the radius does not change.
  - Values and coordinates carry 12 decimals.
  - Acceptance tolerance TOL: a packing is feasible when its worst boundary
    violation and its worst pairwise overlap are both at most TOL. The stored
    best-known packings are feasible far inside TOL (the test reports the real
    margin).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import numpy.typing as npt

Array = npt.NDArray[np.float64]

# The container is a square of side 1. Storage centers it on the origin.
SIDE = 1.0

# The acceptance tolerance for containment and non-overlap. The Packomania files
# carry 12 decimals, so a touching contact or a boundary circle lands within a few
# rounding units of exact. TOL sits well above that rounding floor and well below
# any real overlap, so it accepts every published packing and rejects a real
# violation. The test asserts the measured margin stays far inside TOL.
TOL = 1e-9

_DATA_DIR = Path(__file__).parent / "data"
_COORDS_DIR = _DATA_DIR / "coords"


@dataclass(frozen=True)
class VerifyResult:
    """The outcome of the reference verifier on one packing."""

    n: int
    feasible: bool
    max_out: float  # worst containment violation (0 when every circle is inside)
    max_overlap: float  # worst pairwise overlap (0 when no two circles overlap)
    sum_radii: float


@lru_cache(maxsize=1)
def _best_known() -> dict[int, float]:
    values: dict[int, float] = {}
    for line in (_DATA_DIR / "sumradii.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            values[int(parts[0])] = float(parts[1])
    return values


def load_best_known() -> dict[int, float]:
    """Return the Packomania best-known sum of radii per N, from sumradii.txt.

    The file is complete for N = 1 to 100 and sparse above. Keys are the circle
    count N; values are the published sum of radii. A fresh dict is returned on each
    call, so a caller may mutate it without poisoning the shared cache.
    """
    return dict(_best_known())


def available_coords() -> list[int]:
    """Return the N for which an anchor coordinate file is stored, sorted.

    A file whose name after `csqv` is not a plain integer is ignored, so a stray
    file does not break callers or test collection.
    """
    ns: list[int] = []
    for p in _COORDS_DIR.glob("csqv*.txt"):
        suffix = p.stem.removeprefix("csqv")
        if suffix.isdigit():
            ns.append(int(suffix))
    return sorted(ns)


def load_packing_centered(n: int) -> Array:
    """Load the stored packing for N in the centered frame as an (N, 3) array.

    Columns are [x, y, r] in the square [-0.5, 0.5] x [-0.5, 0.5]. The row indices
    must be exactly 1..N in order, so a dropped, duplicated, or reordered row is a
    hard error rather than a silent wrong packing.
    """
    indices: list[int] = []
    rows: list[tuple[float, float, float]] = []
    for line in (_COORDS_DIR / f"csqv{n}.txt").read_text().splitlines():
        if line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) == 4 and parts[0].isdigit():
            indices.append(int(parts[0]))
            rows.append((float(parts[1]), float(parts[2]), float(parts[3])))
    if indices != list(range(1, n + 1)):
        raise ValueError(f"csqv{n}.txt row indices are not 1..{n} in order")
    packing = np.array(rows, dtype=np.float64)
    if packing.shape != (n, 3):
        raise ValueError(f"csqv{n}.txt parsed to {packing.shape}, expected ({n}, 3)")
    return packing


def _shift_centers(packing: Array, delta: float) -> Array:
    """Return a copy with the x and y centers shifted by `delta`; radii unchanged."""
    out = packing.copy()
    out[:, 0] += delta
    out[:, 1] += delta
    return out


def centered_to_unit(packing: Array) -> Array:
    """Convert a centered-frame packing to the unit square [0, 1] x [0, 1]."""
    return _shift_centers(packing, 0.5)


def unit_to_centered(packing: Array) -> Array:
    """Convert a unit-square packing back to the centered frame."""
    return _shift_centers(packing, -0.5)


def sum_radii(packing: Array) -> float:
    """Return the sum of the circle radii (the CSQV objective, frame-invariant)."""
    return float(packing[:, 2].sum())


def wall_slack(x: Array, y: Array, side: float = SIDE) -> Array:
    """Return each center's distance to the nearest wall of the square [0, side]^2.

    This is the containment bound: a circle with radius r at (x, y) fits inside the
    square exactly when r <= wall_slack. The build verifier and the radii-LP
    primitive both reuse this, so the containment geometry lives in one place.
    """
    return np.minimum.reduce([x, side - x, y, side - y])


def pairwise_distances(x: Array, y: Array) -> tuple[Array, Array, Array]:
    """Return the upper-triangle pair indices and their center distances.

    (i, j, dist) with i < j. Only the upper triangle is formed, so no full n x n
    matrix is materialized. The build verifier and scorer reuse this, so the
    non-overlap geometry lives in one place.
    """
    n = x.shape[0]
    i, j = np.triu_indices(n, k=1)
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    return i.astype(np.int64), j.astype(np.int64), np.sqrt(dx * dx + dy * dy)


def verify(packing_unit: Array, side: float = SIDE, tol: float = TOL) -> VerifyResult:
    """Recompute containment, pairwise overlap, and the sum for a unit-square packing.

    The packing is feasible when its worst containment violation and its worst
    pairwise overlap are both at most `tol`. This is the reference check the later
    verifier reuses; it trusts nothing about how the packing was produced. A
    non-finite coordinate or radius is never feasible; an evolved genome that emits
    NaN or inf is rejected, not silently accepted.
    """
    n = int(packing_unit.shape[0])
    if not np.isfinite(packing_unit).all():
        return VerifyResult(
            n=n,
            feasible=False,
            max_out=float("inf"),
            max_overlap=float("inf"),
            sum_radii=float("nan"),
        )

    x = packing_unit[:, 0]
    y = packing_unit[:, 1]
    r = packing_unit[:, 2]

    # Containment: every circle must fit inside its distance to the nearest wall.
    containment = r - wall_slack(x, y, side)
    max_out = max(0.0, float(containment.max())) if r.size else 0.0

    # Non-overlap: for every pair, the center distance must be at least r_i + r_j.
    max_overlap = 0.0
    if n >= 2:
        i, j, dist = pairwise_distances(x, y)
        overlap = (r[i] + r[j]) - dist
        max_overlap = max(0.0, float(overlap.max()))

    feasible = max_out <= tol and max_overlap <= tol
    return VerifyResult(
        n=n,
        feasible=feasible,
        max_out=max_out,
        max_overlap=max_overlap,
        sum_radii=sum_radii(packing_unit),
    )

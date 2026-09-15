"""A memory-light, scalable center+radius polish for large N (ticket 27/28 large-N path).

The dense SLSQP polish is infeasible at N=484 (~5.6 GB/call, >10 min, SLSQP chokes on ~118k
constraints). This replaces it with L-BFGS-B on a smooth penalty objective: maximize sum of radii
minus a quadratic penalty for overlap and wall violations. L-BFGS-B stores only a few correction
vectors (O(n) memory), and the objective/gradient are vectorized over pairs, so it scales to N=484
in low memory. After the optimize, radii_lp projects to an EXACTLY feasible packing and the verifier
scores it. Penalty continuation (rising lambda) tightens the constraints across rounds.
"""
import numpy as np
from scipy.optimize import minimize

from problems.csqv import problem

SIDE = 1.0


def _obj_grad(z, n, i, j, lam):
    p = z.reshape(n, 3)
    x, y, r = p[:, 0], p[:, 1], p[:, 2]
    g = np.zeros((n, 3))
    f = -r.sum()
    g[:, 2] -= 1.0
    # pairwise overlap penalty: relu(r_i + r_j - d)^2
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    d = np.sqrt(dx * dx + dy * dy) + 1e-12
    o = (r[i] + r[j]) - d
    m = o > 0
    if m.any():
        oi, jj, oo = i[m], j[m], o[m]
        f += lam * np.sum(oo * oo)
        w = 2.0 * lam * oo
        ux, uy = (dx[m] / d[m]), (dy[m] / d[m])
        np.add.at(g[:, 2], oi, w)
        np.add.at(g[:, 2], jj, w)
        np.add.at(g[:, 0], oi, w * (-ux))
        np.add.at(g[:, 0], jj, w * (ux))
        np.add.at(g[:, 1], oi, w * (-uy))
        np.add.at(g[:, 1], jj, w * (uy))
    # wall penalties: relu(r - x), relu(x + r - 1), relu(r - y), relu(y + r - 1)
    for viol, dvx, dvr in (
        (r - x, -1.0, 1.0),          # left
        (x + r - SIDE, 1.0, 1.0),    # right
        (r - y, 0.0, 1.0),           # bottom (dv on y handled below)
        (y + r - SIDE, 0.0, 1.0),    # top
    ):
        mm = viol > 0
        if not mm.any():
            continue
        f += lam * np.sum(viol[mm] ** 2)
        w = 2.0 * lam * viol[mm]
        g[mm, 2] += w * dvr
    # y-wall center gradients (kept separate for clarity)
    vb = r - y
    mb = vb > 0
    g[mb, 1] += 2.0 * lam * vb[mb] * (-1.0)
    vt = y + r - SIDE
    mt = vt > 0
    g[mt, 1] += 2.0 * lam * vt[mt] * (1.0)
    # x-wall center gradients
    vl = r - x
    ml = vl > 0
    g[ml, 0] += 2.0 * lam * vl[ml] * (-1.0)
    vr = x + r - SIDE
    mr = vr > 0
    g[mr, 0] += 2.0 * lam * vr[mr] * (1.0)
    return f, g.ravel()


def scalable_polish(packing, rounds=(1e2, 1e3, 1e4, 1e5), maxiter=200):
    """L-BFGS-B penalty optimize, then project via radii_lp. Returns (packing, feasible_sum)."""
    packing = np.asarray(packing, dtype=np.float64)
    n = packing.shape[0]
    i, j = np.triu_indices(n, k=1)
    z = packing.copy().ravel()
    bounds = [(0.0, SIDE), (0.0, SIDE), (0.0, None)] * n
    for lam in rounds:
        res = minimize(_obj_grad, z, args=(n, i, j, lam), jac=True, method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": maxiter, "ftol": 1e-12})
        z = res.x
    c = z.reshape(n, 3)[:, :2]
    r = problem.radii_lp(c)
    ok, s = problem.verify_and_score(np.column_stack([c, r]))
    return np.column_stack([c, r]), (s if ok else -1.0)


if __name__ == "__main__":
    import sys
    import time

    from problems.csqv import data
    LIVE = {121: 5.797442192252, 144: 6.334929299516, 484: 11.693330291415}
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 121
    rec = LIVE[n]
    rng = np.random.default_rng(0)
    s = int(np.ceil(np.sqrt(n)))
    xs = np.linspace(0.03, 0.97, s)
    gx, gy = np.meshgrid(xs, xs)
    c = np.column_stack([gx.ravel(), gy.ravel()])[:n]
    # grow-push warm-up
    i, j = np.triu_indices(n, k=1)
    for _ in range(10):
        r = problem.radii_lp(c) * 1.06
        for _ in range(15):
            dx = c[i, 0] - c[j, 0]; dy = c[i, 1] - c[j, 1]; d = np.sqrt(dx*dx+dy*dy)+1e-12
            ov = np.maximum(0.0, (r[i]+r[j])-d); ux, uy = dx/d, dy/d
            fx = np.zeros(n); fy = np.zeros(n)
            np.add.at(fx, i, .25*ov*ux); np.add.at(fy, i, .25*ov*uy)
            np.add.at(fx, j, -.25*ov*ux); np.add.at(fy, j, -.25*ov*uy)
            c = c + np.column_stack([fx, fy]); c[:, 0] = np.clip(c[:, 0], r, 1-r); c[:, 1] = np.clip(c[:, 1], r, 1-r)
    gp = problem.verify_and_score(np.column_stack([c, problem.radii_lp(c)]))[1]
    t = time.perf_counter()
    _, ps = scalable_polish(np.column_stack([c, problem.radii_lp(c)]))
    dt = time.perf_counter() - t
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"N={n}: growpush={gp:.4f} +scalable_polish={ps:.6f} record={rec:.6f} "
          f"gap={100*(rec-ps)/rec:+.4f}% time={dt:.0f}s peakRSS={peak:.0f}MB")

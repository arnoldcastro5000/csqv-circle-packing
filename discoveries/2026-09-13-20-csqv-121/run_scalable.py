"""Long parallel record search using the memory-light scalable polish (ticket 28/30 large-N path).

Per restart: varied construction -> grow-push -> scalable_polish (L-BFGS-B penalty) -> keep best,
with perturb-the-best restarts for diversity. Memory-light, so many workers run in parallel. Best
packing persists per-seed atomically, so a long unattended run never loses a record. A cross of the
live record prints *** RECORD BEATEN ***.

Usage: run_scalable.py <n> <budget_s> <seed>
"""
import os
import sys
import time
from pathlib import Path

import numpy as np

from problems.csqv import data, problem

# scalable_polish.py sits next to this file in the discovery archive.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scalable_polish import scalable_polish

LIVE_RECORD = {121: 5.797442192252, 144: 6.334929299516, 484: 11.693330291415}


def construct(n, kind, rng):
    s = int(np.ceil(np.sqrt(n)))
    if kind == "grid":
        m = 0.02 + 0.03 * rng.random()
        xs = np.linspace(m, 1 - m, s)
        gx, gy = np.meshgrid(xs, xs)
        c = np.column_stack([gx.ravel(), gy.ravel()])[:n]
    elif kind == "jitter":
        xs = np.linspace(0.03, 0.97, s)
        gx, gy = np.meshgrid(xs, xs)
        c = np.column_stack([gx.ravel(), gy.ravel()])[:n]
        c = c + rng.normal(0, 0.02, c.shape)
    elif kind == "hex":
        dx = 1.0 / s
        pts, y, row = [], dx * 0.6, 0
        while y < 1.0 and len(pts) < n * 2:
            off = dx / 2 if row % 2 else 0.0
            x = dx * 0.6 + off
            while x < 1.0:
                pts.append((x, y))
                x += dx
            y += dx * 0.866
            row += 1
        c = np.array(pts)[:n] if len(pts) >= n else construct(n, "grid", rng)
    else:
        c = rng.random((n, 2))
    return np.clip(c, 0.001, 0.999)


def growpush(c, n, i, j, rng, outer=12, inner=15):
    for _ in range(outer):
        r = problem.radii_lp(c) * (1.0 + 0.05 + 0.05 * rng.random())
        for _ in range(inner):
            dx = c[i, 0] - c[j, 0]; dy = c[i, 1] - c[j, 1]
            d = np.sqrt(dx * dx + dy * dy) + 1e-12
            ov = np.maximum(0.0, (r[i] + r[j]) - d); ux, uy = dx / d, dy / d
            fx = np.zeros(n); fy = np.zeros(n)
            np.add.at(fx, i, 0.25 * ov * ux); np.add.at(fy, i, 0.25 * ov * uy)
            np.add.at(fx, j, -0.25 * ov * ux); np.add.at(fy, j, -0.25 * ov * uy)
            c = c + np.column_stack([fx, fy])
            c[:, 0] = np.clip(c[:, 0], r, 1 - r); c[:, 1] = np.clip(c[:, 1], r, 1 - r)
    return c


def save_best(n, seed, packing_unit, s):
    out = Path(f"results/n{n}-s{seed}-best.txt")
    tmp = out.with_suffix(".txt.tmp")
    np.savetxt(tmp, data.unit_to_centered(packing_unit), fmt="%.12f",
               header=f"N={n} seed={seed} sum={s:.12f}")
    tmp.replace(out)


def main():
    n = int(sys.argv[1]); budget = float(sys.argv[2]); seed = int(sys.argv[3])
    rec = LIVE_RECORD[n]
    rng = np.random.default_rng(seed)
    i, j = np.triu_indices(n, k=1)
    kinds = ["grid", "jitter", "hex", "random"]
    best, bestc, restarts = -1.0, None, 0
    # Resume: if a saved champion exists for this (n, seed), warm-start from it so a relaunch
    # continues the search (via perturb-best) instead of starting cold. The best file is the
    # centered frame the worker wrote; convert to unit and re-score exactly.
    resume_path = Path(f"results/n{n}-s{seed}-best.txt")
    if resume_path.exists():
        try:
            saved = np.array([[float(v) for v in ln.split()]
                              for ln in resume_path.read_text().splitlines()
                              if not ln.lstrip().startswith("#") and len(ln.split()) == 3])
            if saved.shape == (n, 3):
                pk_u = data.centered_to_unit(saved)
                r = problem.radii_lp(pk_u[:, :2])
                ok, s = problem.verify_and_score(np.column_stack([pk_u[:, :2], r]))
                if ok and s > best:
                    best, bestc = s, pk_u[:, :2].copy()
                    print(f"  RESUMED from {resume_path.name}: best {s:.9f} "
                          f"gap {100*(rec-s)/rec:+.5f}%", flush=True)
        except Exception as e:  # noqa: BLE001  a bad/partial file must not block a fresh start
            print(f"  resume skipped ({e}); starting cold", flush=True)
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < budget:
        try:
            if bestc is not None and rng.random() < 0.6:
                scale = rng.choice([0.006, 0.012, 0.03])
                c = np.clip(bestc + rng.normal(0, scale, bestc.shape), 0.001, 0.999)
            else:
                c = construct(n, kinds[restarts % len(kinds)], rng)
            c = growpush(c, n, i, j, rng)
            pol, s = scalable_polish(np.column_stack([c, problem.radii_lp(c)]))
            if s > best:
                best, bestc = s, pol[:, :2].copy()
                save_best(n, seed, pol, s)
                el = time.perf_counter() - t0
                cross = "  *** RECORD BEATEN ***" if s > rec + 1e-6 else ""
                print(f"  [{el:6.0f}s r{restarts}] best {s:.9f} gap {100*(rec-s)/rec:+.5f}%{cross}",
                      flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  r{restarts} error: {e}", flush=True)
        restarts += 1
    print(f"\nN={n} seed={seed}: best={best:.9f} record={rec:.9f} restarts={restarts} "
          f"gap={100*(rec-best)/rec:+.5f}%", flush=True)


if __name__ == "__main__":
    main()

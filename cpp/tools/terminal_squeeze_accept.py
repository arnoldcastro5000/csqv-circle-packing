"""Acceptance measurement for the terminal second-order squeeze (ticket 45).

Runs the terminal pipeline (problems.csqv.terminal_squeeze.squeeze_champion) on our submitted
N=142 and N=143 champions and reports the fraction of the gap-to-published it closes. All sums
are measured in ONE frame, the trusted Python verifier (problem.verify_and_score), so the
percentages are consistent (the C++ exact LP scores ~5e-6 higher than Python's HiGHS on the same
centers, so the C++ center_polish_analytic baseline is NOT mixed in here; measure that separately
with cpp/csqv/squeeze_check, which reports it in the C++ frame).

The published values are Packomania's re-optimized best-known (ref [21], us). They are the
empirical bar: Specht's re-optimizer is undocumented. A partial pass that beats the first-order
baseline still ships (ticket 44/45).

Usage: PYTHONPATH=. python3 cpp/tools/terminal_squeeze_accept.py [hops]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from problems.csqv import data, terminal_squeeze  # noqa: E402

# N -> (submitted champion file, Packomania published re-optimized best-known).
CASES = {
    142: ("discoveries/2026-09-14-15-csqv-142/n142-best.txt", 6.289851836935),
    143: ("discoveries/2026-09-14-15-csqv-143/n143-best.txt", 6.312993919570),
}


def _load_unit(path: str) -> np.ndarray:
    rows = [[float(v) for v in ln.split()]
            for ln in open(path)
            if not ln.lstrip().startswith("#") and len(ln.split()) == 3]
    return data.centered_to_unit(np.array(rows, dtype=np.float64))


def main() -> int:
    hops = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    for n, (path, published) in CASES.items():
        unit = _load_unit(path)
        raw = terminal_squeeze._exact_score(unit[:, :2])
        gap = published - raw
        t = time.time()
        _, squeezed = terminal_squeeze.squeeze_champion(unit, hops=hops, seed=0)
        dt = time.time() - t
        print(f"N={n}")
        print(f"  radii_lp baseline        = {raw:.12f}")
        print(f"  terminal squeeze         = {squeezed:.12f}  ({100 * (squeezed - raw) / gap:.1f}% of gap)")
        print(f"  published (target)       = {published:.12f}")
        print(f"  remaining gap            = {published - squeezed:.3e}   ({dt:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

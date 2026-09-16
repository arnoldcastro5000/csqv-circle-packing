#!/usr/bin/env python3
"""Independent .pck validator. Reads radii AS WRITTEN (no re-derivation).

Centered frame: square [-0.5, 0.5]^2, center (0,0), side 1.
Checks the artifact exactly as a submitter/maintainer would re-parse it, at zero tolerance:
geometry (no overlap, inside the box) PLUS Packomania format compliance (radii sorted
increasing, line-1 max equals the actual max, filename N matches the row count, all radii
positive and finite). Returns a nonzero exit code on any invalid file, so it is safe as a
`validate && submit` guard. Pure numpy, no project imports: an ungameable second implementation.
"""
import re
import sys

import numpy as np

HALF = 0.5  # square is [-HALF, HALF]^2
DECLARED_MAX_TOL = 1e-12  # line-1 max radius vs actual, at 15 dp emit precision


def load_pck(path):
    with open(path) as fh:
        raw = fh.readlines()
    if len(raw) < 2:
        raise ValueError(f"{path}: not a .pck (need a max-radius line, an author line, then rows)")
    declared_max = float(raw[0].split()[0])
    author = raw[1].strip()
    rows = []
    for ln in raw[2:]:
        parts = ln.split()
        if len(parts) == 3:
            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
    return declared_max, author, np.array(rows, dtype=np.float64)


def validate(path):
    m = re.search(r"csqv(\d+)\.pck$", path)
    fileN = int(m.group(1)) if m else None
    declared_max, author, a = load_pck(path)
    n = a.shape[0]
    x, y, r = a[:, 0], a[:, 1], a[:, 2]

    res = {"path": path, "N": n, "fileN": fileN, "author": author}
    res["N_ok"] = (fileN is None) or (fileN == n)
    res["all_finite"] = bool(np.all(np.isfinite(a)))
    res["all_r_pos"] = bool(np.all(r > 0.0))

    # declared max radius vs actual
    res["max_r"] = float(r.max())
    res["declared_max"] = declared_max
    res["declared_max_diff"] = abs(declared_max - float(r.max()))

    # sorted increasing by radius (non-decreasing)
    dr = np.diff(r)
    res["sorted_incr"] = bool(np.all(dr >= -1e-15))
    res["worst_inversion"] = float(dr.min()) if n > 1 else 0.0

    # containment: worst wall violation (>0 means poke outside)
    wall_slack = HALF - np.maximum(np.abs(x), np.abs(y))
    wall_viol = r - wall_slack
    res["worst_wall"] = float(wall_viol.max())
    res["worst_wall_idx"] = int(np.argmax(wall_viol))

    # overlap: worst pairwise violation (>0 means overlap). N<2 has no pairs.
    i, j = np.triu_indices(n, k=1)
    if len(i):
        dx = x[i] - x[j]
        dy = y[i] - y[j]
        dist = np.sqrt(dx * dx + dy * dy)
        ovl = (r[i] + r[j]) - dist
        k = int(np.argmax(ovl))
        res["worst_overlap"] = float(ovl.max())
        res["worst_overlap_pair"] = (int(i[k]), int(j[k]))
    else:
        res["worst_overlap"] = float("-inf")
        res["worst_overlap_pair"] = (-1, -1)

    res["sum_radii"] = float(r.sum())
    res["feasible"] = (res["worst_wall"] <= 0.0) and (res["worst_overlap"] <= 0.0)

    # The submission gate: geometric feasibility PLUS Packomania format compliance. A file that
    # is non-overlapping but unsorted, mis-declares its max radius, has a wrong N, or carries a
    # non-positive/non-finite radius is still rejected by the maintainer, so it is not valid.
    res["valid"] = bool(
        res["feasible"] and res["all_finite"] and res["all_r_pos"]
        and res["sorted_incr"] and res["N_ok"] and res["declared_max_diff"] <= DECLARED_MAX_TOL
    )
    return res


def main(argv=None):
    """Validate each .pck path. Return 0 if all are valid, 1 otherwise (the submission gate)."""
    paths = sys.argv[1:] if argv is None else list(argv)
    all_valid = True
    for p in paths:
        r = validate(p)
        all_valid = all_valid and r["valid"]
        print("=" * 78)
        print(f"FILE      {r['path']}")
        print(f"author    {r['author']!r}")
        print(f"N         {r['N']}  (filename says {r['fileN']})  "
              f"match={'OK' if r['N_ok'] else 'FAIL'}")
        print(f"finite    {r['all_finite']}   radii>0 {r['all_r_pos']}")
        print(f"sorted    increasing-by-r={r['sorted_incr']}  "
              f"worst_step={r['worst_inversion']:.2e}")
        print(f"maxR      file-line1={r['declared_max']:.15f}  "
              f"actual-max={r['max_r']:.15f}  diff={r['declared_max_diff']:.2e}")
        print(f"WALL      worst violation = {r['worst_wall']:+.3e}  "
              f"(circle {r['worst_wall_idx']})  {'FEASIBLE' if r['worst_wall']<=0 else 'VIOLATED'}")
        print(f"OVERLAP   worst violation = {r['worst_overlap']:+.3e}  "
              f"pair {r['worst_overlap_pair']}  {'FEASIBLE' if r['worst_overlap']<=0 else 'VIOLATED'}")
        print(f"SUM       {r['sum_radii']:.15f}")
        print(f"VERDICT   {'*** VALID at zero tolerance ***' if r['valid'] else '!!! INVALID (see fields above) !!!'}")
    return 0 if all_valid else 1


if __name__ == "__main__":
    sys.exit(main())

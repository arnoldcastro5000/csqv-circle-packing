"""Compare worker runs of a head build (and a forced-fallback build) with a base build.

Gate G1: the startup benchmark data lines (the radii hashes) are identical.
Gate G3: the per-restart "best" events of a seed are identical up to the shorter run's
restart count. A time-budget run is never reproducible as a whole, so the faster build
only extends the base sequence.
A fallback run must also print the pricing WARNING line; a head run must not.

Usage: python3 compare_runs.py <runs_dir>
<runs_dir> holds log_<build>_<n>_<seed>.txt and o_<build>_<n>_<seed>/ from gate.sh, where
<build> is base, head or fb. Exit status 0 only if every gate holds.
"""
import collections
import glob
import os
import re
import sys

EVENT = re.compile(r"\[\s*\d+s s(\d+) r(\d+)\] best (\S+)")
WARNING = "WARNING radii-LP pricing fell back"


def events(log):
    with open(log) as f:
        return [(int(m.group(2)), m.group(3)) for m in map(EVENT.search, f) if m]


def restarts(log):
    with open(log) as f:
        m = re.search(r"restarts=(\d+)", f.read())
    return int(m.group(1)) if m else -1


def benchmark_lines(out_dir):
    files = glob.glob(os.path.join(out_dir, "benchmark-*.txt"))
    if len(files) != 1:
        return None
    with open(files[0]) as f:
        return [line for line in f if not line.startswith("#")]


def main(runs_dir):
    failures = 0
    compared = 0
    totals = collections.defaultdict(lambda: [0, 0, 0])  # n -> runs, base restarts, head restarts
    for build in ("head", "fb"):
        for log in sorted(glob.glob(os.path.join(runs_dir, f"log_{build}_*.txt"))):
            name = os.path.basename(log)[len(f"log_{build}_"):-len(".txt")]
            n = name.split("_")[0]
            base_log = os.path.join(runs_dir, f"log_base_{name}.txt")
            if not os.path.exists(base_log):
                print(f"MISSING base run for {log}")
                failures += 1
                continue
            compared += 1
            rb, rv = restarts(base_log), restarts(log)
            limit = min(rb, rv)
            b = [e for e in events(base_log) if e[0] <= limit]
            v = [e for e in events(log) if e[0] <= limit]
            if not b or b != v:
                print(f"G3 MISMATCH {log}: base {b[:3]} vs {v[:3]}")
                failures += 1
            bb = benchmark_lines(os.path.join(runs_dir, f"o_base_{name}"))
            bv = benchmark_lines(os.path.join(runs_dir, f"o_{build}_{name}"))
            if bb is None or bb != bv:
                print(f"G1 DIFF {log}")
                failures += 1
            with open(log) as f:
                warned = WARNING in f.read()
            if warned != (build == "fb"):
                print(f"WARNING line {'missing' if build == 'fb' else 'present'} in {log}")
                failures += 1
            if build == "head":
                t = totals[n]
                t[0] += 1
                t[1] += rb
                t[2] += rv
    for n, (runs, rb, rv) in sorted(totals.items(), key=lambda kv: int(kv[0])):
        ratio = rv / rb if rb > 0 else float("nan")
        print(f"N={n} runs={runs} restarts base={rb} head={rv} ratio={ratio:.2f}")
    print(f"runs compared: {compared}")
    ok = failures == 0 and compared > 0
    print("G1+G3", "IDENTICAL" if ok else f"FAILED ({failures} problems)")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))

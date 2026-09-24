#!/usr/bin/env bash
# Bit-identical gate for a speed change: compares the working tree (head) with a git ref (base).
#
# Usage: cpp/bench/gate.sh <base-ref> [out_dir]
#
# Head = the csqv sources in THIS working tree (uncommitted edits included). Base = the
# csqv sources at <base-ref>, from a temporary git worktree. The bench tools always come from
# this tree's cpp/bench/, and each build uses -I<tree>/cpp/csqv only, so no stale header can
# shadow the source. A third build, fb, is head with the half-integral guard forced to fail
# at the 3rd basis change of each LP solve (the fallback path).
#
# Gates (each must hold; the exit status is 0 only if all hold):
#   checksums  perf_driver seed 777 (N=90 20 iters, N=143 8 iters): base == head == fb ==
#              the sealed values.
#   G1 + G3    worker runs (1 thread): the startup benchmark lines are identical, and the
#              per-restart "best" events extend the base sequence exactly.
#   G2         trajectory_hash (radii, duals, center polish): base == head == fb.
#   pricing    pricing_check: each incremental reduced cost equals the full recompute, and
#              the entering pick from the pricing keys equals the pick of the scan.
#   G4         make test and make accept in this tree (accept only if the target exists).
# Speed (not a gate): lp_bench and a 30 s worker run per N, base vs head, one job at a time.
#
# Environment (defaults in brackets):
#   GATE_NS       worker and G2 sizes              [30 60 90 121 143]
#   GATE_SEEDS    worker seeds per N               [12]
#   GATE_BUDGET   worker budget per run, seconds   [10]
#   GATE_JOBS     parallel worker runs             [2]
#   GATE_G2_SEEDS trajectory_hash seeds per N      [20]
#   GATE_SPEED    1 = run the speed section        [1]
#   GATE_BASE_FLAGS  compiler flags for base       [the cpp/Makefile default, below]
#   GATE_HEAD_FLAGS  compiler flags for head + fb  [the cpp/Makefile default, below]
#   CXX           compiler                         [g++]
# The default flags are "-O3 -funroll-loops -std=c++17 -march=native -ffp-contract=off", as in
# cpp/Makefile. To gate a flag change, give the old flags to base, for example:
#   GATE_BASE_FLAGS="-O3 -funroll-loops -std=c++17" cpp/bench/gate.sh HEAD
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  sed -n '2,34p' "$0"
  exit 2
fi

BASE_REF=$1
HERE=$(cd "$(dirname "$0")" && pwd)
CPP=$(cd "$HERE/.." && pwd)
REPO=$(git -C "$CPP" rev-parse --show-toplevel)
OUT=${2:-/tmp/csqv-gate-$(date +%Y%m%d-%H%M%S)}
NS=${GATE_NS:-30 60 90 121 143}
SEEDS=${GATE_SEEDS:-12}
BUDGET=${GATE_BUDGET:-10}
JOBS=${GATE_JOBS:-2}
G2_SEEDS=${GATE_G2_SEEDS:-20}
SPEED=${GATE_SPEED:-1}
CXX=${CXX:-g++}
DEFAULT_FLAGS="-O3 -funroll-loops -std=c++17 -march=native -ffp-contract=off"
read -r -a BASE_FLAGS <<< "${GATE_BASE_FLAGS:-$DEFAULT_FLAGS}"
read -r -a HEAD_FLAGS <<< "${GATE_HEAD_FLAGS:-$DEFAULT_FLAGS}"
SEALED_90=98.374190007874
SEALED_143=49.910861299027

mkdir -p "$OUT"
BASE_TREE="$OUT/base-src"
git -C "$REPO" worktree add --detach "$BASE_TREE" "$BASE_REF" > /dev/null
trap 'git -C "$REPO" worktree remove --force "$BASE_TREE"' EXIT
echo "base: $BASE_REF = $(git -C "$BASE_TREE" rev-parse --short HEAD)"
echo "head: working tree at $(git -C "$REPO" rev-parse --short HEAD)$([[ -z $(git -C "$REPO" status --porcelain -- cpp) ]] || echo ' + local edits')"
echo "base flags: ${BASE_FLAGS[*]}"
echo "head flags: ${HEAD_FLAGS[*]}"
echo "out:  $OUT"

FAIL=0
fail() {
  echo "FAIL: $*"
  FAIL=1
}

# build <name> <base|head> <csqv_dir> <source> [extra flags...]
build() {
  local name=$1 side=$2 inc=$3 src=$4
  shift 4
  local flags=("${HEAD_FLAGS[@]}")
  [[ $side == base ]] && flags=("${BASE_FLAGS[@]}")
  "$CXX" "${flags[@]}" -I"$inc" "$@" "$src" -o "$OUT/bin/$name" -pthread
}

echo "== build"
mkdir -p "$OUT/bin"
pids=()
for side in base head fb; do
  case $side in
    base) inc="$BASE_TREE/cpp/csqv"; flags=base; extra=() ;;
    head) inc="$CPP/csqv"; flags=head; extra=() ;;
    fb) inc="$CPP/csqv"; flags=head; extra=(-DCSQV_LP_FORCE_FALLBACK_AT=3) ;;
  esac
  build "worker_$side" "$flags" "$inc" "$inc/worker.cpp" ${extra[@]+"${extra[@]}"} &
  pids+=($!)
  for tool in perf_driver lp_bench trajectory_hash; do
    build "${tool}_$side" "$flags" "$inc" "$HERE/$tool.cpp" ${extra[@]+"${extra[@]}"} &
    pids+=($!)
  done
done
build pricing_check head "$CPP/csqv" "$HERE/pricing_check.cpp" -DCSQV_LP_VERIFY_PRICING &
pids+=($!)
for pid in "${pids[@]}"; do
  wait "$pid" || { echo "FAIL: a build failed"; exit 1; }
done

echo "== checksums (perf_driver seed 777)"
for spec in "90 20 $SEALED_90" "143 8 $SEALED_143"; do
  read -r n iters sealed <<< "$spec"
  hexes=()
  for side in base head fb; do
    line=$("$OUT/bin/perf_driver_$side" "$n" "$iters" 777)
    echo "  $side $line"
    [[ $line == *"checksum=$sealed "* ]] || fail "N=$n $side checksum is not the sealed $sealed"
    hexes+=("${line##*hex=}")
  done
  [[ ${hexes[0]} == "${hexes[1]}" && ${hexes[1]} == "${hexes[2]}" ]] ||
    fail "N=$n checksum bits differ: ${hexes[*]}"
done

echo "== G2 (trajectory_hash, $G2_SEEDS seeds, 5 polish steps)"
for n in $NS; do
  for side in base head fb; do
    "$OUT/bin/trajectory_hash_$side" "$n" "$G2_SEEDS" 5 > "$OUT/g2_${side}_$n.txt" &
  done
  wait
  hashes=()
  for side in base head fb; do
    hashes+=("$(sed 's/ time=.*//' "$OUT/g2_${side}_$n.txt")")
  done
  echo "  ${hashes[0]}"
  [[ ${hashes[0]} == "${hashes[1]}" && ${hashes[1]} == "${hashes[2]}" ]] ||
    fail "G2 N=$n: base '${hashes[0]}' head '${hashes[1]}' fb '${hashes[2]}'"
done

echo "== pricing check (incremental d_j vs the full recompute; key pick vs the scan)"
for n in 30 90 143; do
  if ! "$OUT/bin/pricing_check" "$n" 4 | sed 's/^/  /'; then fail "pricing check N=$n"; fi
done

echo "== G1 + G3 (worker, ${BUDGET} s, 1 thread, $SEEDS seeds per N, $JOBS jobs)"
RUNS="$OUT/runs"
mkdir -p "$RUNS"
export OUT RUNS BUDGET
read -r -a ns <<< "$NS"
fb_ns=$(printf '%s\n' "${ns[0]}" "${ns[${#ns[@]} / 2]}" "${ns[${#ns[@]} - 1]}" | sort -nu)
{
  for n in $NS; do
    for ((s = 0; s < SEEDS; ++s)); do
      echo "base $n $s"
      echo "head $n $s"
    done
  done
  # The fallback build: the smallest, a middle and the largest N, up to 4 seeds each.
  for n in $fb_ns; do
    for ((s = 0; s < SEEDS && s < 4; ++s)); do echo "fb $n $s"; done
  done
} | xargs -P "$JOBS" -n 3 bash -c \
  '"$OUT/bin/worker_$0" "$1" "$BUDGET" "$2" 1 "$RUNS/o_$0_$1_$2" > "$RUNS/log_$0_$1_$2.txt" 2>&1'
if ! python3 "$HERE/compare_runs.py" "$RUNS" | sed 's/^/  /'; then fail "G1 + G3"; fi

echo "== G4 (make test, make accept)"
if make -C "$CPP" test > "$OUT/make_test.txt" 2>&1; then
  echo "  make test: PASS"
else
  fail "make test (see $OUT/make_test.txt)"
fi
if ! make -C "$CPP" -n accept > /dev/null 2>&1; then
  echo "  make accept: no such target (skipped)"
elif make -C "$CPP" accept > "$OUT/make_accept.txt" 2>&1; then
  echo "  make accept: PASS"
else
  fail "make accept (see $OUT/make_accept.txt)"
fi

if [[ $SPEED == 1 ]]; then
  echo "== speed (one job at a time; not a gate)"
  for n in 90 121 143; do
    tb=$("$OUT/bin/lp_bench_base" "$n")
    th=$("$OUT/bin/lp_bench_head" "$n")
    [[ ${tb##*acc=} == "${th##*acc=}" ]] || fail "lp_bench N=$n acc differs"
    sb=$(sed 's/.*time=\([^ ]*\).*/\1/' <<< "$tb")
    sh=$(sed 's/.*time=\([^ ]*\).*/\1/' <<< "$th")
    echo "  lp_bench N=$n base ${sb}s head ${sh}s speedup $(awk "BEGIN{printf \"%.2f\", $sb/$sh}")x"
  done
  for n in 90 121; do
    for side in base head; do
      "$OUT/bin/worker_$side" "$n" 30 0 1 "$OUT/speed_${side}_$n" > "$OUT/speed_${side}_$n.txt" 2>&1
    done
    rb=$(grep -o 'restarts=[0-9]*' "$OUT/speed_base_$n.txt" | head -1 | cut -d= -f2)
    rh=$(grep -o 'restarts=[0-9]*' "$OUT/speed_head_$n.txt" | head -1 | cut -d= -f2)
    echo "  worker 30 s N=$n restarts base $rb head $rh ($(awk "BEGIN{printf \"%.2f\", $rh/$rb}")x)"
  done
fi

if [[ $FAIL == 0 ]]; then
  echo "GATE: PASS"
else
  echo "GATE: FAIL"
fi
exit $FAIL

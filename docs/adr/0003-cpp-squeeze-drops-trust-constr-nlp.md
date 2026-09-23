# The native C++ squeeze drops the trust-constr NLP

The Python terminal squeeze (`problems/csqv/terminal_squeeze.py`) has three stages: a
`trust-constr` second-order NLP, the SLP jammer (`jam_slp`), and the KKT-Newton
(`contact_newton`). When we port the squeeze into the self-contained C++17 worker so it
seals its own record in-process, we port the EFFECTIVE pipeline (the existing
`center_polish_analytic`, then `jam_slp`, then an optional gated `contact_newton`) and we
DROP the `trust-constr` NLP.

## Why

The NLP does not earn a dependency-free port. On the real champions N=142 and N=143 the
SLP jammer drives the Donev jamming margin to ~0 and reaches the sealed value to the float
floor (100% of the gap), while the `trust-constr` NLP plateaus at ~72-75% and is
dominated. So the NLP is the weakest stage on exactly the inputs the squeeze runs on.

The cost is also lopsided. `jam_slp` needs a general bounded-variable LP, a bounded and
testable add. A dependency-free `trust-constr` equivalent is a sparse trust-region
interior-point NLP solver written from scratch, a large build, for a stage a cheaper stage
already beats. The worker's headline property is that it is self-contained with no
dependencies, so a heavyweight vendored solver is exactly what we must avoid.

## The trade-off

We chose behavioral parity on the OUTPUT (the sealed value) over a literal stage-for-stage
port. A future reader will see a stage missing and wonder why, so we record it here. The
cost is that the C++ path does not reproduce the NLP's intermediate iterates; we accept
this because the NLP result is strictly dominated by `jam_slp` on these inputs, so the
sealed value is unchanged.

`contact_newton` is kept but gated: it lands only if it adds measurable digits over
`jam_slp` alone in the C++ frame (`jam_slp` already reached the float floor on N=142/143,
so its marginal value may be nil). This keeps the critical path short.

## Status

Accepted. The Python pipeline is unchanged and stays as the reference and the
source-of-truth verifier for submissions. The build is the general bounded-variable
simplex, then the native `jam_slp` port. Revisit if a case appears where `jam_slp` leaves
a gap the NLP would have closed.

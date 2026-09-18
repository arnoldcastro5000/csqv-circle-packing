# CSQV data (Packomania), provenance and convention

This directory holds the CSQV problem data: variable-radii circles packed in
a square of side 1, maximizing the sum of radii. It holds the best-known values, a
set of anchor coordinate files, and the pinned acceptance convention. The loader,
the frame conversion, and the reference verifier live in `problems/csqv/data.py`.

## Source

- Registry: Packomania CSQV, "best known packings of variable-sized circles in a
  square with maximized sum of radii", maintained by E. Specht.
  Index: https://www.packomania.com/csqv/csqv.html
- Fetched: 2026-09-12.

## Files

- `sumradii.txt` : the best-known sum of radii per N. Complete for N = 1 to 100 and
  sparse above. Source: https://www.packomania.com/csqv/txt/sumradii.txt
- `author.txt` : the contributor reference per N. Kept for provenance and for the
  contested-versus-open landscape read. Source:
  https://www.packomania.com/csqv/txt/author.txt
- `coords/csqvN.txt` : the coordinate file for N circles, one anchor per N in
  {13, 20, 26, 32, 40, 50, 60, 75, 90, 100}, spanning the less-contested mid band.
  Source: https://www.packomania.com/csqv/txt/csqvN.txt

The stored anchor set is a sample for the cross-check, not a final target band.

## Coordinate format

Each `coords/csqvN.txt` file has comment lines that start with `#`, then N data
rows `index x y r`, then a footer `sumradii = <value>`. Values carry 12 decimals.

## Frame convention

- Storage frame: centered origin, the square [-0.5, 0.5] x [-0.5, 0.5], side 1. A
  boundary circle satisfies `abs(x) + r = 0.5` exactly (for example `csqv26`
  circle 4: `0.415360499304 + 0.084639500696 = 0.5`).
- Solver and verifier frame: the unit square [0, 1] x [0, 1]. Convert by adding
  0.5 to x and to y; the radius does not change.

## Acceptance tolerance

A packing is feasible when its worst containment violation and its worst pairwise
overlap are both at most `TOL = 1e-9` (see `problems/csqv/data.py`). This sits well
above the rounding floor of the 12-decimal files and well below any real overlap.

Measured across every stored anchor packing (see the test):

- worst boundary violation: 2.78e-17
- worst pairwise overlap: 1.89e-12
- worst gap between our recomputed sum and the published best-known: 1.00e-12

So every published best-known packing is feasible far inside `TOL`, and our sum of
radii reproduces the published value to the rounding floor. This locks the
convention the later verifier must match, the analog of the CVRP BKS cross-check.

## Submission format (for later, when a champion beats a live best-known)

Packomania accepts a `.pck` file: line 1 the largest radius, line 2 the author,
then `x y r` rows sorted by increasing radius at maximum decimals, coordinates in
the centered frame. Convert from the unit-square frame with
`problems/csqv/data.py:unit_to_centered` before writing.

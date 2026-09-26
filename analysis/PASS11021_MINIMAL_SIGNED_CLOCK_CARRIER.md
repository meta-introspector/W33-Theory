# Pass 11021 — The minimal exact signed clock carrier is 24-dimensional

Producer: `analysis/w33_pass11021_minimal_signed_clock_carrier.py`
Certificate: `data/w33_pass11021_minimal_signed_clock_carrier.json`
Regression: `tests/test_w33_pass11021_minimal_signed_clock_carrier.py`

Pass 11020 proved that the canonical signed E6 cubic has an exact hidden
GL2(3) monomial symmetry on its full 27-coordinate H27 carrier.

The remaining question was how large a signed-equivariant field is forced once
the coarse four-clock amplitudes are embedded into that exact action.

The answer is exact: the four clock indicators generate a 24-dimensional
invariant subspace, and the three-dimensional clock augmentation module already
generates the same 24 dimensions.
## Exact coordinate split

The GL2(3) support action on the 27 H27 coordinates has four orbits of sizes

    1 + 2 + 8 + 16 = 27.

The three coordinates with quotient address (a,b)=(0,0) form the central
three-dimensional invariant sector.

The other 24 coordinates form the noncentral invariant sector:

    27 = 24 noncentral + 3 central.

The clock fibres lie entirely in the noncentral sector, so their signed orbit
closure can never require the central three coordinates. The executable rank
calculation proves that they require all 24 noncentral coordinates.
## Every clock direction resolves as 2 + 4

Each of the four projective clock directions contains six H27 coordinates.
Against the two noncentral support orbits, every fibre has the same split:

    6 = 2 + 4.

Thus the four-amplitude clock quotient suppresses a uniform internal
two-layer structure in every direction.

The eight-point orbit contributes two coordinates per clock direction.
The sixteen-point orbit contributes four coordinates per direction.
## The eight-point orbit is a quadratic phase graph

In the current physical-Clifford H27 gauge, the eight-point support orbit has
exactly one phase lift over every nonzero quotient vector (a,b).

An exhaustive quadratic fit over F3 gives one and only one solution:

    c = a*b + a + b  (mod 3).

The other two phase lifts over each nonzero (a,b) make the complementary
sixteen-point orbit.

So the 24-dimensional carrier is not an anonymous collection of coordinates:
it is the union of one quadratic phase section and its two-sheet complement.
## Minimality

Let v_1,...,v_4 be the four six-point fibre indicator vectors.

Closing their signed GL2(3) orbit gives rank 24.

More strongly, take only the three augmentation generators

    v_1-v_4,  v_2-v_4,  v_3-v_4.

Their signed orbit span also has rank 24.

Therefore the exact signed completion of the three-dimensional clock order
parameter is already the full noncentral H27 sector.

The center remains separately invariant with rank 3.
## Character diagnostics

For the 24-dimensional signed sector, the character trace histogram is

    trace 0: 27 elements
    trace 2: 12 elements
    trace 6:  8 elements
    trace 24: 1 element.

The invariant subspace has dimension 2 and the commutant has dimension 19.

For comparison, the central three-dimensional sector has invariant dimension 2
and commutant dimension 5; the full 27-dimensional representation has invariant
dimension 4 and commutant dimension 34.

These are exact representation diagnostics, not physical mode counts.

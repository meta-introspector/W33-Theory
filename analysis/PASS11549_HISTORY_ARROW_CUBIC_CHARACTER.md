# Pass 11549 — The finite history arrow is a cubic Hamming character

**Status:** exact finite theorem.  The old global orientation cycle now has a
closed local formula.

Pass 11547 turns the 27 histories into ternary Hamming addresses.  A null step
then changes all three coordinates, so its displacement is
[
d=(d_1,d_2,d_3)in{pm1}^3.
]
Pass 11548 identifies the repository PSp linear subgroup as the even-sign
signed-permutation group.

Define
[
oxed{chi(d)=d_1d_2d_3in{pm1}.}
]

That one cubic sign is the old history orientation.

## Exact match to the frozen 108-edge cycle

The old certificate
`w33_20260924_history_invariant_cycle_orientation.json`
constructed its orientation by averaging a seed triangle over all (648)
PSp history symmetries and then taking a boundary.

In the repository's canonical edge ordering (i<j), with original
symmetric-matrix history labels (s_i), Pass 11549 finds

[
oxed{
c_{ij}
=
-chi!left(H_{m sym}(s_j-s_i)ight)
}
]

for **all 108 of 108 edges**.

The minus sign is only the old certificate's arbitrary choice of global
orientation.  The invariant line is one-dimensional, so the opposite sign is
the same orientation datum.

This eliminates the group average from the theorem: the cycle can be computed
edge-by-edge from the displacement alone.

## Why it is an orientation

Because the dimension is odd,
[
chi(-d)=-chi(d).
]
So (chi) is automatically antisymmetric under edge reversal.

The eight body diagonals split into two tetrahedra:
[
chi=+1:
quad
(+++) , (+--),(-+-),(--+),
]
[
chi=-1:
quad
(++-),(+-+),(-++),(---).
]

These are exactly the (4+4) oriented-null sectors that the older Fourier
analysis found abstractly.

## The 36 temporal triangles are a (4	imes9) foliation

The distance-three Hamming graph has a remarkable property here:

[
oxed{	ext{every null edge lies in exactly one triangle}.}
]

There are (108) edges, so there are
[
108/3=36
]
edge-disjoint temporal triangles.

Choose the four projective body-diagonal directions by the convention
[
chi(d)=+1.
]
For each direction there are
[
27/3=9
]
affine lines.  Hence
[
oxed{36=4	imes9.}
]

Orient each affine line
[
x	o x+d	o x+2d	o x
]
along its unique representative (d) with (chi(d)=+1).

The sum of those 36 directed triangle boundaries is exactly the (chi)
edge cycle.

The old producer used sorted vertices and coefficient (+1) on every one of
the 36 triangles; that convention gives the opposite global sign, exactly as
the 108-edge comparison predicts.

## Symmetry law in one equation

For a signed permutation (Min O(3,3)), let
[
sigma(M)=	ext{product of its three coordinate signs}.
]
Coordinate permutations do not change a coordinate product, so

[
oxed{
chi(Md)=sigma(M)chi(d).
}
]

Therefore:

- actual repo (PSp_B^{m lin}=kersigma) preserves the arrow;
- all (27) translations preserve it because it depends only on displacement;
- the affine PSp group of order (648) fixes it;
- the outer (sigma=-1) coset of size (648) negates it.

This is exactly the old orbit-average theorem, now proved locally.

## Connection to the corrected orientation square

Pass 11548 showed that the global history bit is (sigma), not the ordinary
orthogonal determinant.  Pass 11549 explains why:

[
sigma
]
is literally the representation character carried by the arrow cycle.

The outer element (-I) sends
[
dmapsto-d
]
and therefore
[
chimapsto-chi.
]

The inner (648	o324) bit, coordinate-permutation parity, does **not** flip
the cubic arrow character; it acts inside the orientation-preserving PSp
group.

## TOE relevance, with firewall intact

This is a meaningful simplification of the finite temporal architecture.

The history arrow is no longer an opaque invariant obtained by averaging
(648) transformations.  It is the parity of the three simultaneous ternary
coordinate changes on a null step.

Equivalently: the finite arrow distinguishes the two tetrahedral chiralities
inside the cube of eight null directions.

What this does **not** prove is equally important.  It does not derive:

- thermodynamic irreversibility;
- a low-entropy cosmological boundary condition;
- observed CPT dynamics;
- a continuum time orientation;
- Einstein causality or gravity.

It is the exact local formula for the repository's already-certified finite
history-orientation line.

## Reproducibility

Producer:
`analysis/w33_pass11549_history_arrow_cubic_character.py`

Certificate:
`data/PART_W33_PASS11549_HISTORY_ARROW_CUBIC_CHARACTER.json`

Upstream:
- old invariant-cycle certificate;
- Pass 11547 Hamming conjugacy;
- Pass 11548 corrected subgroup dictionary.

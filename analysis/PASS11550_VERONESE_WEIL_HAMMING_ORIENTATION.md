# Pass 11550 — The Veronese null sheet is the (-i) Weil sheet

**Status:** exact objectwise weld across the clock code, finite cone, spectral
clock and corrected Hamming arrow.

Three older constructions had the same (4+4) pattern but no closed formula
connecting all of them:

1. Pass 10946 mapped the four clock rays in (mathbb P^1(mathbf F_3)) to
   four rank-one null matrices by the symmetric-square Veronese map;
2. the September temporal Fourier analysis split the eight oriented rank-one
   null modes into Weil-phase orbits (4_{-i}+4_{+i});
3. Pass 11549 found the finite history arrow
   [
   chi(z)=operatorname{sgn}(z_1)operatorname{sgn}(z_2)operatorname{sgn}(z_3)
   ]
   on the Hamming null cube.

Pass 11550 proves these are the same two sheets.

## 1. The square/Veronese sheet

Take
[

u(x,y)=(x^2,xy,y^2).
]
For the four projective clock rays of Pass 10946,
[
[1:0],quad[0:1],quad[1:1],quad[1:2],
]
the four images are
[
(1,0,0),quad(0,0,1),quad(1,1,1),quad(1,2,1).
]

They are precisely the four rank-one null matrices already attached
coordinate-by-coordinate to the Hesse/M36/tetracode clock.

Map them into Pass 11547's Hamming coordinates.  Every one has
[
oxed{chi=-1.}
]

Their negatives form the other tetrahedron and all have
[
oxed{chi=+1.}
]

So the eight oriented rank-one null vectors split canonically as

[
oxed{
{
u(v)}_{[v]inmathbb P^1(mathbf F_3)}
;sqcup;
-{
u(v)}_{[v]inmathbb P^1(mathbf F_3)}.
}
]

Over (mathbf F_3), this is literally the square-scaled versus nonsquare-scaled
rank-one cone: the only nonzero square is (1), and the nonsquare is
(-1=2).

## 2. Exact weld to the old Weil phases

The old fixed-nine certificate stores all eight oriented rank-one symmetric
matrices and their normalized Weil phases.

Objectwise comparison gives

[
oxed{
	ext{Veronese sheet}
=
chi=-1
=
	ext{Weil phase }-i,
}
]
and

[
oxed{
-	ext{Veronese sheet}
=
chi=+1
=
	ext{Weil phase }+i.
}
]

For all eight modes,

[
oxed{	ext{Weil phase}=i,chi.}
]

This is an exact row-by-row identity, not a matching of orbit sizes.

## 3. Outer temporal reversal is nonsquare scaling

The outer history extension multiplies a symmetric matrix by
[
-1=2,
]
the unique nonsquare scalar of (mathbf F_3).

It therefore sends
[
Smapsto-S,
]
which swaps the two rank-one sheets.

In Hamming coordinates this sends
[
zmapsto-z,
qquad
chimapsto-chi.
]

On the old spectral certificate it sends
[
+ileftrightarrow-i,
]
i.e. complex conjugation of the normalized Weil phase.

Thus the same outer operation is simultaneously

- nonsquare scaling of the finite quadratic form;
- exchange of the two Veronese cone sheets;
- reversal of the Hamming cubic arrow;
- exchange of the two (4)-dimensional fixed Fourier sectors;
- complex conjugation of their (pm i) Weil phases.

This is the objectwise finite content behind the old statement that the outer
coset reverses temporal orientation.

## 4. Why the cubic character is natural for (W(D_3))

In the Hamming frame, actual repo
[
PSp_B^{m lin}=W(D_3)
]
is the even signed permutation group.

The polynomial
[
p(z)=z_1z_2z_3
]
is invariant under every even sign change and under coordinate permutations:
[
p(Mz)=p(z),
qquad Min W(D_3).
]

For the full signed-permutation group
[
O(3,3),
]
it transforms by the sign-product character
[
p(Mz)=sigma(M),p(z).
]

On the null cube, (p=1) is (chi=+1) and (p=2=-1) is (chi=-1).

This is exactly the degree-(3) type-(D_3) invariant.  Standard invariant
theory gives the type-(D_n) basic degrees
[
2,4,ldots,2n-2,n;
]
for (D_3) these are (2,3,4), and the degree-(3) coordinate product is the
odd-degree generator.  The repo-specific result is not that classical
invariant theory, but that **this generator is the already-certified finite
history orientation**.

## 5. Important distinction: cone orientation is not tetracode sign choice

Pass 10946 retained explicit oriented representatives in
(mathbb P^1(mathbf F_3)) in order to freeze the tetracode evaluation signs.

Pass 11550 does **not** remove that choice.

Indeed
[

u(v)=vv^{mathsf T}=
u(-v),
]
so the Veronese map forgets the sign of the projective representative.

What it *does* canonically select is the orientation sheet of the **rank-one
symmetric-matrix cone**: the Veronese/square sheet versus its nonsquare
negative.

Thus there are two related but distinct lifts:

- a code-coordinate lift on the (P^1) evaluation representatives;
- a temporal/cone lift on the oriented rank-one symmetric matrices.

They are compatible in Pass 10946, but one does not determine the other.

## TOE relevance

This collapses several formerly separate signs into one exact finite datum.

The finite history arrow, normalized Weil phase, square/nonsquare cone sheet
and PGSp outer reversal are now four languages for the same binary structure.

That is a substantial internal unification of the temporal architecture, but
it remains finite.  It does not derive a thermodynamic arrow, continuum CPT
theorem, weak-interaction chirality, or gravitational time orientation.

## Reproducibility

Producer:
`analysis/w33_pass11550_veronese_weil_hamming_orientation.py`

Certificate:
`data/PART_W33_PASS11550_VERONESE_WEIL_HAMMING_ORIENTATION.json`

Upstream exact owners:
- Pass 10946 clock-code-cone objectwise theorem;
- `w33_20260924_fixed9_chiral_null_fourplusfour.json`;
- Pass 11547 Hamming conjugacy;
- Pass 11548 corrected orientation square;
- Pass 11549 cubic history arrow.

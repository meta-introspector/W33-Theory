# Pass 11540 — The 27-event history chart is `VO(3,3)`, and its two orientation bits form a square

**Status:** exact finite theorem; standard classical-group identification; no continuum/gravity promotion.

## Executive result

The Bell-centered history chart already frozen in the repository is

[
mathcal H = operatorname{Sym}_2(mathbf F_3)cong mathbf F_3^3
]

with interval

[
q(t,x,y)=t^2-x^2-y^2
]

and null adjacency

[
usim vquadLongleftrightarrowquad u
e v,;q(u-v)=0.
]

That graph has a standard name:

[
oxed{mathcal H_{m null}cong VO(3,3)}
]

—the **parabolic affine orthogonal polar graph** in dimension three over (mathbf F_3).

This is not a new numerical spectrum calculation.  The repository already froze the
(27)-vertex, degree-(8), (108)-edge null graph, its Fourier spectrum, the
(1296/648) Bell-stabilizer actions, and the (648/324) line-orientation quotient.
The new content is the classical-group identification that puts those facts into one
canonical ladder and proves that the two sign reductions are **different characters**.

## Exact group ladder

The verifier exhaustively checks all (3^9=19683) ternary (3	imes3) matrices against (q).

It finds

[
|O(3,3)|=48,qquad |SO(3,3)|=24,qquad |Omega(3,3)|=12.
]

The four projective null directions are

[
Q(2,3)=mathbb P^1(mathbf F_3),
]

and (SO(3,3)) acts faithfully on them as all (24) permutations.  Thus, in the
standard rank-three exceptional isomorphisms,

[
oxed{SO(3,3)cong PGL(2,3)cong S_4}
]

and

[
oxed{Omega(3,3)cong PSL(2,3)cong A_4}.
]

Internally, the verifier identifies the (A_4) subgroup twice and finds the same (12)
matrices:

1. the even-permutation kernel of the action on the four null directions;
2. the derived subgroup ([SO(3,3),SO(3,3)]).

For odd (q), standard finite orthogonal-group theory identifies (Omega(n,q)) with the
kernel of the spinor norm on (SO(n,q)).  Hence the repository's old
(S_4	o A_4) orientation quotient has a standard name: **the spinor-norm quotient**.

References used for the nomenclature/classification:

- Sage/PassageMath, `AffineOrthogonalPolarGraph`: odd-dimensional `sign=None`
  is the parabolic graph (VO(d,q)).
- Magma Handbook, `Omega(n,q)`: for odd (q), (Omega(n,q)) is the kernel
  of the spinor norm on (SO(n,q)).
- Standard rank-three exceptional isomorphisms:
  (SO(3,q)cong PGL(2,q)) and (Omega(3,q)cong PSL(2,q)) for odd (q).

## The affine history symmetry

Translations by (mathbf F_3^3) give

[
oxed{
3^3{:}O(3,3)
supset
3^3{:}SO(3,3)
supset
3^3{:}Omega(3,3)
}
]

with orders

[
oxed{1296supset648supset324}.
]

All (27cdot48=1296) affine isometries are explicitly checked to preserve every
null-history neighborhood.

These are exactly the three orders already sitting in the temporal corpus:

- (1296): the full (PGSp) Bell-line action;
- (648): the (PSp) Bell-line action (3^3{:}S_4);
- (324): the existing oriented-line kernel (3^3{:}A_4).

So the older pieces were not three unrelated group-order coincidences.  They are the
orthogonal subgroup ladder of the same finite Minkowski chart.

## The real breakthrough: there are two independent orientation bits

The quotient from (1296) to (648) and the quotient from (648) to (324) are not the
same (C_2).

### Bit 1 — global history/CPT-like orientation

[
oxed{
(3^3{:}O)/(3^3{:}SO)cong C_2.
}
]

Its character is the determinant of the (3)-dimensional orthogonal linear part.

This matches the already-frozen outer Bell-stabilizer bit:

[
PGSp_B/PSp_Bcong C_2,
]

which reverses the unique global history-cycle orientation and is the same outer
unitary/antiunitary and half-spin/chirality bit used in the Forty Points paper.

### Bit 2 — null-frame/spinor orientation

[
oxed{
(3^3{:}SO)/(3^3{:}Omega)cong C_2.
}
]

Its character is the parity of the (S_4) permutation of the four projective null
directions; in standard orthogonal language this is the spinor-norm quotient.

This matches the older

[
3^3{:}S_4longrightarrow3^3{:}A_4
]

line-orientation quotient of order (648	o324).

### They are independent

The central inversion

[
-I_3in O(3,3)
]

has determinant (-1), but projectively fixes all four null directions.  It therefore flips
the first bit and not the second.

Exhaustive enumeration gives exactly (12) linear orthogonal transformations in each of
the four character sectors

[
(+,+),quad(+,-),quad(-,+),quad(-,-).
]

Hence

[
oxed{
(3^3{:}O(3,3))/(3^3{:}Omega(3,3))
cong C_2	imes C_2.
}
]

This is the **history orientation square**.

## Why this matters for the TOE programme

The four-set that the paper independently meets as

- the four projective null directions,
- the four Hesse/qutrit MUB frames,
- and the four tetracode coordinates

is now also the natural projective null conic on which

[
SO(3,3)cong PGL(2,3)cong S_4
]

acts.  Its even subgroup

[
Omega(3,3)cong PSL(2,3)cong A_4
]

is therefore not an arbitrary "even permutation" convention: it is the spinorially
oriented half of the finite orthogonal group.

That sharpens the clock–cone–code bridge.  The tetracode's four coordinates and the finite
light cone's four rays are not merely two (S_4)-sets of size four; they sit on the standard
rank-three orthogonal boundary (mathbb P^1(mathbf F_3)), and the (S_4/A_4) reduction is
the orthogonal/spinor reduction.

## Fourier compatibility

The same quadratic form is self-dual under the additive Fourier characters of
(mathbf F_3^3).  The already-frozen adjacency eigenvalues are recovered by quadratic
shell:

[
lambda(k)=
egin{cases}
8,&k=0,\
-1,&q(k)=0, k
e0,\
-4,&q(k)=1,\
2,&q(k)=2.
end{cases}
]

So the paper's (1+8) clock-fixed Fourier sector is exactly the trivial character plus the
nonzero dual null cone.  This is old repo content; Pass 11540 records it only to show that
the new (VO(3,3)) naming is objectwise compatible with the existing spectral clock.

## Firewalls

This pass **does not** claim:

- that (VO(3,3)) is continuum Minkowski space;
- that (O(3,3)) is the physical Lorentz group;
- that the finite null relation derives the measured speed of light;
- that either (C_2) chooses the observed weak chirality;
- that the affine graph supplies Einstein dynamics or gravity.

The theorem is a finite quadratic-geometry/classical-group identification.  Its value is
that it removes an ambiguity in the internal architecture: the temporal chart has two
different canonical orientation characters, and the existing corpus had already measured
both of them without naming the second one as the spinor norm.

## Reproducibility

Producer:

`analysis/w33_pass11540_history_vo33_orientation_square.py`

Frozen certificate:

`data/PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json`

Prior exact owners welded, not overwritten:

- `data/w33_20260924_history_bigcell_q43_compactification.json`
- `data/w33_20260924_history_invariant_cycle_orientation.json`
- `data/w33_20260924_null_history_spectral_clock.json`
- `data/PART_W33_PASS9741_9748_ORIENTATION_CHARACTER_WELD.json`

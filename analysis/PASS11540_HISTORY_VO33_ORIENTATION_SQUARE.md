# Pass 11540 — The 27-event history chart is `VO(3,3)`

> **Correction (Pass 11548).** The graph and group-order theorem below survives.
> The original version incorrectly identified the repository's (648)-element
> PSp Bell stabilizer with the affine (SO(3,3)) subgroup from equality of
> orders. Objectwise reconstruction in Pass 11548 shows instead
> [
> PSp_B^{m lin}=kersigma=W(D_3)cong S_4,
> ]
> where (sigma) is the product of the three coordinate signs in the
> Pass-11547 Hamming frame.  (SO(3,3)=kerdet) is a *different* order-(24)
> (S_4).  Their intersection is
> (A_4=Omega(3,3)).  The global (1296	o648) history-orientation bit is
> (sigma), not the orthogonal determinant.  See
> `PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.md`.

**Status:** exact finite theorem, corrected subgroup dictionary; no
continuum/gravity promotion.

## Exact graph theorem

The Bell-centered history chart is
[
mathcal H=operatorname{Sym}_2(mathbf F_3)congmathbf F_3^3
]
with quadratic interval
[
q(t,x,y)=t^2-x^2-y^2
]
and adjacency
[
usim viff u
e v, q(u-v)=0.
]

This graph is the parabolic affine orthogonal polar graph
[
oxed{mathcal H_{m null}cong VO(3,3)}.
]

It has
[
27	ext{ vertices},qquad 108	ext{ edges},qquad k=8
]
and spectrum
[
8^1,quad2^{12},quad(-1)^8,quad(-4)^6.
]

Pass 11547 later strengthens this by an explicit linear conjugacy to the
distance-three graph of the ternary Hamming scheme (H(3,3)).

## Classical orthogonal group

Exhaustive enumeration of all (3^9) ternary (3	imes3) matrices gives
[
|O(3,3)|=48,qquad |SO(3,3)|=24,qquad|Omega(3,3)|=12.
]

The four projective null directions carry the standard rank-three action
[
SO(3,3)cong PGL(2,3)cong S_4,
]
with
[
Omega(3,3)cong PSL(2,3)cong A_4.
]

The affine orthogonal groups therefore have orders
[
1296,qquad648,qquad324.
]

**Important:** the middle affine (SO) subgroup has the same order (648) as
the repository PSp Bell stabilizer, but they are not the same subgroup.

## Corrected repository stabilizer ladder

Pass 11548 reconstructs the actual repository action
[
Smapsto ASA^{mathsf T},
qquad
Ain GL(2,3)/{pm I},
]
and conjugates it to the Hamming coordinates of Pass 11547.

There
[
O(3,3)=C_2^3{:}S_3
]
is the signed-permutation group. Define
[
sigma=	ext{product of the three coordinate signs},
qquad
pi=	ext{coordinate-permutation parity}.
]
Then
[
det=sigmapi.
]

The exact subgroup dictionary is
[
oxed{
PGSp_B^{m lin}=O(3,3),
qquad
PSp_B^{m lin}=kersigma=W(D_3)cong S_4,
}
]
and
[
oxed{
PSp_B^{m lin}cap SO(3,3)=A_4=Omega(3,3).
}
]

Thus the repository's affine ladder is
[
1296
supset
648=3^3{:}W(D_3)
supset
324=3^3{:}A_4.
]

## Corrected orientation square

The two repo-relevant bits are

[
oxed{sigma}
]
for the global (1296	o648) history/outer bit, and

[
oxed{pi}
]
for the inner (648	o324) line-orientation bit.

The ordinary orthogonal determinant is the third nontrivial character
[
oxed{det=sigmapi}.
]

Hence
[
oxed{O(3,3)/A_4cong C_2	imes C_2.}
]

This is the corrected history orientation square.

In the Hamming model the eight oriented null vectors are the cube corners
((pm1,pm1,pm1)).  The actual PSp subgroup (kersigma) preserves their
coordinate product, splitting them into two tetrahedra (4+4); the full
PGSp/O group fuses all eight.  This exactly matches the independently frozen
old (4+4) chiral-null result.

## Spinor-norm firewall

The standard statement
[
Omega(3,3)=ker(	ext{spinor norm}:SO(3,3)	o C_2)
]
remains correct.

The withdrawn shortcut was to identify the repository
(648	o324) quotient itself with (SO	oOmega).  The repo quotient is
instead
[
W(D_3)	o A_4
]
through (S_4) permutation parity.  The same (A_4) also equals
(Omega(3,3)) because it is the intersection of the two distinct
order-(24) (S_4) subgroups.

## Firewalls

This pass does **not** claim:

- (VO(3,3)) is continuum Minkowski space;
- (O(3,3)) is the physical Lorentz group;
- the finite null relation derives the measured speed of light;
- either finite sign character selects observed weak chirality;
- the graph supplies Einstein dynamics or gravity.

## Reproducibility

Primary producer:
`analysis/w33_pass11540_history_vo33_orientation_square.py`

Corrected certificate:
`data/PART_W33_PASS11540_HISTORY_VO33_ORIENTATION_SQUARE.json`

Superseding character audit:
- `analysis/w33_pass11548_history_orientation_square_correction.py`
- `data/PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json`
- `analysis/PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.md`

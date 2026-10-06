# Pass 11548 — Correction: the PSp history subgroup is (W(D_3)), not (SO(3,3))

**Status:** exact objectwise correction of one character assignment in Pass 11540.

Pass 11540 correctly identified the 27-event graph as (VO(3,3)), correctly
computed
[
|O(3,3)|=48,qquad |SO(3,3)|=24,qquad |Omega(3,3)|=12,
]
and correctly recovered the affine orders
[
1296supset648supset324.
]
It then made one inference from equal orders that was too fast: it identified
the repository's (648)-element PSp Bell stabilizer with
(3^3{:}SO(3,3)).

Pass 11547's Hamming coordinates make it possible to decide this objectwise.
The identification is **wrong**.  This pass reconstructs the actual PSp action
from the old producer itself and replaces the order-only inference.

## 1. Reconstruct the PSp subgroup from the source action

The old history producer acts on
[
S=egin{pmatrix}a&b\b&cend{pmatrix}inoperatorname{Sym}_2(mathbf F_3)
]
by
[
Slongmapsto A S A^{mathsf T},
qquad
Ain GL(2,3)/{pm I},
]
giving (24) linear maps and (27) translations.

Conjugate those (24) maps by Pass 11547's Hamming map
[
H_{m sym}(a,b,c)
=
(2a+b+c,;2a+2b+c,;2a+2c).
]
Every one becomes a signed permutation matrix.

Write a signed permutation as a coordinate permutation (piin S_3)
together with three signs (epsilon_i=pm1). Define
[
sigma=epsilon_1epsilon_2epsilon_3,
qquad
arepsilon=operatorname{sgn}(pi).
]
Then
[
det=sigmaarepsilon.
]

Exact enumeration gives
[
oxed{
PSp_B^{m lin}
=
kersigma
=
W(D_3),
}
]
the even-signed permutation group of order (24), isomorphic to (S_4).

It is **not**
[
SO(3,3)=kerdet.
]

## 2. The corrected (C_2	imes C_2) square

In Hamming coordinates
[
O(3,3)cong C_2^3{:}S_3.
]
Its two most transparent independent characters are

[
sigma=	ext{product of coordinate signs},
qquad
arepsilon=	ext{coordinate-permutation parity}.
]

The ordinary determinant is their product:
[
oxed{det=sigmaarepsilon.}
]

Each of the four ((sigma,arepsilon)) sectors contains exactly (12)
elements.  Therefore
[
O(3,3)/A_4cong C_2	imes C_2.
]

The three index-two kernels are:

- (kersigma=W(D_3)cong S_4): the **actual repository PSp** linear subgroup;
- (ker(sigmaarepsilon)=SO(3,3)cong S_4): a different order-(24) subgroup;
- (kerarepsiloncong C_2	imes A_4): the third index-two subgroup.

Their relevant common kernel is
[
oxed{
kersigmacapkerarepsilon
=
PSp_B^{m lin}cap SO(3,3)
=
A_4
=
Omega(3,3).
}
]

The verifier independently finds this same (A_4) as the derived subgroup of
both (PSp_B^{m lin}cong S_4) and (SO(3,3)cong S_4).

## 3. Which bit is the arrow-of-time bit?

The existing repo theorem says the (1296)-element full Bell stabilizer has
an inner (648) that preserves the global history orientation cycle and an
outer (648) that reverses it.

The corrected character is
[
oxed{sigma,}
]
not the (3	imes3) orthogonal determinant.

The outer multiplier (epsilon=2=-1) in the old history producer becomes
[
-I_3
]
in Hamming coordinates.  It has

[
sigma=-1,qquad
arepsilon=+1,qquad
det=-1.
]

Thus
[
oxed{
PGSp_B/PSp_B
quadleftrightarrowquad
O(3,3)/kersigma.
}
]

This is the already-frozen history-cycle / unitary-vs-antiunitary /
half-spin-chirality outer bit.

## 4. The old (4+4) null-mode split becomes two cube tetrahedra

In Hamming coordinates the eight oriented null vectors are exactly
[
(pm1,pm1,pm1).
]

Their product
[
chi(d)=d_1d_2d_3in{pm1}
]
splits them into
[
4+4.
]

The actual PSp subgroup (kersigma) preserves (chi), so the two orbits
are

[
chi=+1:
quad
(+++) , (+--),(-+-),(--+),
]

and

[
chi=-1:
quad
(++-),(+-+),(-++),(---).
]

These are the two tetrahedra inside the cube.

The full (O(3,3)) fuses them, and (-I) exchanges them.  This is exactly the
older repo result that PSp splits the eight oriented null Fourier modes as
(4+4), while the PGSp extension fuses all eight.

So the correction does more than repair a label: it gives a closed-form
coordinate meaning to the old chiral-null split.

## 5. What is the (648	o324) bit?

Inside the actual PSp subgroup (kersigma), the second quotient is
[
oxed{
kersigma/A_4cong C_2,
}
]
and its character is
[
oxed{arepsilon,}
]
the parity of the coordinate permutation.

The same parity is the sign of the induced permutation on the four
projective null directions.  Hence this is precisely the old
[
3^3{:}S_4longrightarrow3^3{:}A_4
]
line-orientation quotient.

The orthogonal determinant restricted to PSp happens to equal
(arepsilon), because (sigma=+1) there.  But globally determinant and
history orientation are different characters.

## 6. Spinor-norm firewall

The standard finite-orthogonal statement remains true:
[
Omega(3,3)=ker(	ext{spinor norm}:SO(3,3)	o C_2).
]

What is withdrawn is the shortcut
[
	ext{repo }648	o324
=
SO	oOmega.
]

The correct statement is:

- repo (648	o324) is (W(D_3)	o A_4) via (S_4) sign;
- independently (SO(3,3)	oOmega(3,3)=A_4) is the spinor-norm quotient;
- the same (A_4) sits at the intersection of the two distinct (S_4)'s.

## 7. Corrected orientation dictionary

The exact finite square is therefore

[
egin{array}{c|c|c}
	ext{character} & 	ext{kernel} & 	ext{repo role}\ hline
sigma & W(D_3)cong S_4 & 	ext{global history / antiunitary / chirality bit}\
arepsilon & C_2	imes A_4 & 	ext{independent coordinate-permutation bit}\
det=sigmaarepsilon & SO(3,3)cong S_4 & 	ext{orthogonal determinant bit}
end{array}
]

and the inner (324) kernel is
[
3^3{:}A_4.
]

## Scope

This is an exact finite subgroup correction.  It does not alter the
(VO(3,3)) graph theorem, the Hamming conjugacy, the clock factorization, or
the finite group orders.  It adds no continuum Lorentz, physical CPT,
weak-chirality-selection, or gravity claim.

## Reproducibility

Producer:
`analysis/w33_pass11548_history_orientation_square_correction.py`

Certificate:
`data/PART_W33_PASS11548_HISTORY_ORIENTATION_SQUARE_CORRECTION.json`

Primary old sources rechecked objectwise:
- `analysis/w33_20260924_history_bigcell_q43_compactification.py`
- `analysis/w33_20260924_history_invariant_cycle_orientation.py`
- `analysis/w33_20260924_fixed9_chiral_null_fourplusfour.py`
- Pass 11547 Hamming conjugacy.

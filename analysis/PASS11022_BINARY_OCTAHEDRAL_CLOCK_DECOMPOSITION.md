# Pass 11022 — The exact signed clock carrier splits 12 + 12

Producer: `analysis/w33_pass11022_binary_octahedral_clock_decomposition.py`
Certificate: `data/w33_pass11022_binary_octahedral_clock_decomposition.json`
Regression: `tests/test_w33_pass11022_binary_octahedral_clock_decomposition.py`

> **Terminology correction (Pass 11025).** The numerical decomposition below is
> unchanged, but (GL_2(3)) is the plus Schur cover (2^+S_4), not the
> non-isomorphic binary-octahedral minus cover (2^-S_4). Accordingly the
> classical (SU(2))/ADE McKay identification does not apply to this group.

Pass 11021 identified the minimal exact signed clock carrier as the 24
noncentral coordinates of the current (H_{27}) gauge.

This pass determines the actual complex representation carried by those
24 coordinates under the exact split (GL_2(3)) monomial symmetry.
## Exact character

Using the eight conjugacy classes of the concrete matrix group (GL_2(3)),
with class sizes

[
1,8,1,8,6,6,6,12,
]

the 24-dimensional signed representation has character

[
oxed{chi_{24}=(24,0,0,6,0,0,0,2)}.
]

The executable table is independently matched against GAP's
`CharacterTable(GL(2,3))` and verifies its own row orthogonality.
With irreducibles named by their projective or central behavior, the exact
decomposition is

[
oxed{
V_{24}
=
2,mathbf1
oplus det
oplus (3_{
m std}!otimes!det)
oplus 2,3_{
m std}
oplus 3,4_{
m spin}.
}
]

The quotient two-dimensional irrep and both faithful two-dimensional spinor
irreps occur with multiplicity zero.

The pullback (3_{
m std}) is identified internally from the action on the
four projective clock labels (mathbf P^1(mathbf F_3)), not by a name match.
## The central involution

Let (z=-Iin GL_2(3)). Because it is central, its (pm1) eigenspaces are
group-invariant. The exact trace is zero on (V_{24}), hence

[
oxed{
V_{24}=V_+oplus V_-,
qquad
dim V_+=dim V_-=12.
}
]

Their characters are

[
chi_+=(12,3,12,3,0,0,0,2),
]
[
chi_-=(12,-3,-12,3,0,0,0,0).
]
The irreducible decomposition sharpens this completely:

[
oxed{
V_+
=
2,mathbf1
oplusdet
oplus(3_{
m std}!otimes!det)
oplus2,3_{
m std},
}
]

while

[
oxed{
V_-=3,4_{
m spin}.
}
]

Thus every irreducible in the plus sector descends through
(GL_2(3)/{pm I}=S_4), whereas the entire minus sector is carried by three
copies of one faithful four-dimensional representation on which the central
double-cover element acts as (-1).
## Every clock fibre already contains both halves

Each of the four projective clock directions is a six-coordinate fibre.
The central involution acts inside every fibre with

[
oxed{6=3_+oplus3_-}.
]

This is not merely a global dimension balance.

Take the four coarse fibre indicators from Pass 11021 and project them using

[
P_pm=
rac{1}{2}(1pm z).
]

Closing either projected set under the exact signed group gives rank 12.
Even more strongly, start only from the three coarse clock augmentation
generators (v_1-v_4,v_2-v_4,v_3-v_4): their plus projections span all of
(V_+), and their minus projections span all of (V_-).
So the three-dimensional clock order parameter is genuinely a coarse seed for
both exact twelve-dimensional sectors.

## Cover-type firewall

The exact group is (GL_2(3)=2^+S_4). It has 13 nonidentity involutions.
That alone rules out a faithful embedding in (SL_2(mathbb C)), because a
finite subgroup of (SL_2(mathbb C)) has only one nontrivial involution,
(-I). Binary octahedral is the other, non-isomorphic Schur cover
(2^-S_4).

The repository result here is therefore strictly a decomposition in the
complex representation ring of the exact (GL_2(3)) action.

## Boundary

"Spinorial" here is a representation-theoretic word: the central element of
the double cover acts as (-1). It does **not** identify these twelve modes
with physical fermions, Lorentz spinors, or particle species.

Pass 11023 originally reported an affine-(E_7) tensor graph because its
order-eight cyclotomic character values had been simplified incorrectly from
(pm isqrt2) to (pmsqrt2). Pass 11025 corrects that arithmetic: the
faithful 2D tensor quiver is directed and non-symmetric, and its underlying
undirected graph has 11 edges rather than the seven-edge affine-(E_7) tree.
No (E_7) gauge theory, spacetime symmetry, or dynamics is inferred.

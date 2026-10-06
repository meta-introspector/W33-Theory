# Pass 11541 — The finite history chart is a self-dual 3-class quadratic-shell scheme

**Status:** exact finite theorem; prior-art-aware; no continuum promotion.

Pass 11540 identifies the 27-event null graph as the parabolic affine polar graph
(VO(3,3)) and separates its two orientation characters.  Pass 11541 now closes
the **entire four-way separation algebra** of the same chart.

Let
[
V=mathbf F_3^3,qquad q(t,x,y)=t^2-x^2-y^2.
]
For a pair of histories (x,y), classify the difference (y-x) by
[
R_0:0,qquad
R_1:q=0
e y-x,qquad
R_2:q=1,qquad
R_3:q=2.
]

The relation valencies are exactly
[
oxed{(1,8,6,12)}.
]
These are the paper's coincident/null/two anisotropic history shells.

## Complete Bose--Mesner closure

The verifier checks every triple count
[
p_{ij}^{,k}
=
#{z:(x,z)in R_i,;(z,y)in R_j}
]
for every ordered pair ((x,y)) in each relation (R_k), and proves that it
depends only on (i,j,k).  Hence the four relations form a commutative
3-class translation association scheme.

The full intersection tensor is frozen in
`data/PART_W33_PASS11541_HISTORY_QUADRATIC_SHELL_SCHEME.json`.

## Exact Fourier eigenmatrix

The additive characters of (mathbf F_3^3) diagonalize all four adjacency
relations.  Ordering both primal and dual shells as (0, q=0
e0, q=1, q=2)
gives
[
oxed{
P=
egin{pmatrix}
1&8&6&12\
1&-1&-3&3\
1&-4&3&0\
1&2&0&-3
end{pmatrix}.
}
]

The null-graph column is therefore
[
(8,-1,-4,2)
]
with multiplicities
[
(1,8,6,12),
]
recovering the already-frozen null-history spectrum
(8^1,(-1)^8,(-4)^6,2^{12}).

## The stronger result: exact formal self-duality

Direct multiplication gives
[
oxed{P^2=27I_4.}
]

Thus, in this shell ordering,
[
oxed{Q=P}
]
and the primitive multiplicities are the same numbers as the primal valencies:
[
oxed{(1,8,6,12)}.
]

So the spatial shell census and the Fourier-shell census are not merely the
same multiset by accident.  The whole shell algebra is formally self-dual.

This is the finite harmonic object sitting behind the paper's spectral clock.
The order-three clock used only one column of this table; the association
scheme shows that all three nontrivial interval types belong to one closed
commutative algebra.

## Relation to the literature

Translation association schemes attached to quadratic forms and symmetric
bilinear forms are standard.  In odd characteristic the quadratic-form and
symmetric-bilinear-form schemes are known to be isomorphic/formally self-dual;
see, e.g., work of Wang--Ma--Ho and the modern treatment by Kai-Uwe Schmidt,
*Quadratic and symmetric bilinear forms over finite fields and their
association schemes*, Algebraic Combinatorics 3 (2020), 161--189.

The repo-specific content here is narrower and explicit: the small
(mathbf F_3) **value-shell fusion** relevant to the Forty Points history
chart, its exact (4	imes4) eigenmatrix and intersection tensor, and its weld
to the existing temporal-clock certificate.

## TOE relevance without overclaim

The paper's finite interval already separated histories into
[
1+8+6+12=27.
]
Pass 11541 shows that this split is simultaneously:

- a translation-distance partition;
- a Bose--Mesner algebra;
- a Fourier spectral partition;
- formally self-dual.

That is a stronger statement than "the counts match."  Position-like and
character/momentum-like descriptions are two bases of the **same exact finite
association scheme**.

But this is still finite harmonic analysis.  It does not establish continuum
position--momentum duality, a physical Lorentzian field theory, Einstein
dynamics, or any measured coupling.

## Reproducibility

Producer:
`analysis/w33_pass11541_history_quadratic_shell_scheme.py`

Certificate:
`data/PART_W33_PASS11541_HISTORY_QUADRATIC_SHELL_SCHEME.json`

Upstream owners:
- Pass 11540 (VO(3,3)) orientation square;
- `w33_20260924_null_history_spectral_clock.py`;
- `w33_20260924_history_bigcell_q43_compactification.py`.

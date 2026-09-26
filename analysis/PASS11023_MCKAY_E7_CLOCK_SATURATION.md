# Pass 11023 — The clock carrier has an exact affine-E7 McKay placement

Producer: `analysis/w33_pass11023_mckay_e7_clock_saturation.py`
Certificate: `data/w33_pass11023_mckay_e7_clock_saturation.json`
Regression: `tests/test_w33_pass11023_mckay_e7_clock_saturation.py`

Pass 11022 decomposed the exact 24-dimensional signed clock carrier under
(GL_2(3)), the binary-octahedral double cover of (S_4).

This pass reconstructs the McKay graph directly from that same exact character
table and places the clock multiplicity vector on it.
## Reconstruct the graph

Choose either faithful two-dimensional spinor (S). For every irreducible
representation (
ho_i), decompose

[
Sotimes
ho_i.
]

The resulting adjacency graph has eight nodes, seven simple edges, is connected,
and has degree sequence

[
1,1,1,2,2,2,2,3.
]

Removing the trivial node leaves the seven-node finite (E_7) tree, with the
four-dimensional spinorial node as the trivalent branch.
For one choice of the defining spinor, the edges are

[
mathbf1-2_a-3_t-4_s-3-2_b-det
]

together with the branch

[
2_q-4_s.
]

The other faithful two-dimensional spinor gives the determinant-twisted
version of the same graph.

The irreducible dimension vector is

[
d=(1,1,2,2,2,3,3,4)
]

and the executable calculation verifies the affine McKay relation

[
oxed{A d=2d}.
]
## Place the actual clock carrier

Pass 11022 gives the multiplicity vector

[
m=(2,1,0,0,0,1,2,3)
]

in the node order

[
(mathbf1,det,2_q,2_a,2_b,3_t,3,4_s).
]

Let (n_{
m nonlin}) be the indicator of the six irreducibles of dimension
greater than one:

[
n_{
m nonlin}=(0,0,1,1,1,1,1,1).
]

Then for **either** defining spinor,

[
oxed{A m=3n_{
m nonlin}}.
]
Equivalently,

[
oxed{
Sotimes V_{24}
=
3(2_qoplus2_aoplus2_boplus3_toplus3oplus4_s).
}
]

So tensoring the actual clock carrier by one defining spinor populates every
non-one-dimensional McKay node with exactly the same multiplicity three, while
both one-dimensional leaves vanish.

This is much stronger than merely observing that the group belongs to the
binary-octahedral/ADE family.
## The 12 + 12 split saturates opposite bipartitions

Write

[
V_{24}=V_+oplus V_-
]

for the central-(-I) eigenspaces of Pass 11022.

The plus sector contains only representations descending to (S_4). The minus
sector is (3,4_s). Tensoring by (S) flips central parity and gives

[
oxed{
Sotimes V_+
=
3(2_aoplus2_boplus4_s),
}
]

[
oxed{
Sotimes V_-
=
3(2_qoplus3_toplus3).
}
]
Thus each twelve-dimensional half lands on **all three** nontrivial nodes of
the opposite McKay bipartition, uniformly with multiplicity three.

A second exact graph identity follows:

[
A,n_{
m nonlin}=d-n_{
m nonlin},
]

hence

[
oxed{A(A+I)m=3d}.
]

This gives a compact graph-theoretic characterization of the clock
multiplicity vector relative to the affine-(E_7) dimension vector.
## Prior-art boundary

The McKay correspondence between the binary octahedral group and affine
(E_7) is classical. The repository increment is the placement of the
**already-derived signed W33/H27 clock carrier** in that representation ring,
including the exact saturation identities above.

## Boundary

No (E_7) Lie-algebra gauge symmetry has been derived here. The McKay graph is
a representation graph of the finite group (GL_2(3)); its ADE label does not
supply gauge bosons, a continuum action, masses, or couplings.

Likewise the tensor product by a defining spinor is an exact operation in the
finite representation ring, not yet a physical interaction vertex.

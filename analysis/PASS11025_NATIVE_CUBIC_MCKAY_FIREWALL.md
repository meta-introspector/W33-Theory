# Pass 11025 — Schur-cover correction and native-cubic tensor firewall

Producer: `analysis/w33_pass11025_native_cubic_mckay_firewall.py`
Certificate: `data/w33_pass11025_native_cubic_mckay_firewall.json`
Regression: `tests/test_w33_pass11025_native_cubic_mckay_firewall.py`

This pass corrects a real error in the original interpretation of Passes
11022–11023 while preserving the parts that actually survive exact checking.

The group is

[
oxed{GL_2(3)=2^+S_4},
]

not the binary-octahedral minus cover (2^-S_4). The concrete group has
13 nonidentity involutions. A faithful finite subgroup of
(SL_2(mathbb C)) can have only the single nontrivial involution (-I), so
the exact clock group cannot be the binary-octahedral (SU(2)) group.
The character audit also finds the arithmetic source of the mistaken graph.
For the two faithful 2D characters the order-eight values are

[
oxed{pm isqrt2},
]

not (pmsqrt2). Equivalently every irreducible two-dimensional character
has

[
Lambda^2(2)=det,
]

rather than the trivial determinant required for an (SL_2) embedding.

With the cyclotomic values restored, tensoring by one faithful 2D irrep gives
a directed, non-symmetric quiver with 14 arrows; the conjugate faithful irrep
gives the transpose quiver. Its underlying undirected graph has 11 edges, not
the seven-edge affine-(E_7) tree.
The most interesting representation identity from Pass 11023 nevertheless
survives exactly. For either faithful 2D irrep (S),

[
oxed{
Sotimes V_{24}
=
3(2_qoplus2_aoplus2_boplus3_toplus3oplus4).
}
]

So the saturation was real even though the ADE label was not.

The second half of this pass asks whether the corrected finite symmetry/quiver
determines the native (E_6) cubic. It does not:

[
oxed{dimoperatorname{Sym}^3(V_{24})^{GL_2(3)}=71}.
]

The central split resolves this as

[
71=23+48,
]
where

[
dimoperatorname{Sym}^3(V_+)^G=23,qquad
dim(V_+otimesoperatorname{Sym}^2V_-)^G=48,
]

while both odd-central-parity cubic sectors vanish.

The actual noncentral restriction of the signed (E_6) cubic is far more
specific. Its 32 noncentral triads expand in a central-eigenbasis into exactly

[
oxed{16;(+,+,+)+48;(+,-,-)}
]

nonzero monomials, with no odd-parity terms. After removing the common
(1/(2sqrt2)) basis factor, every surviving coefficient has magnitude two.

Therefore the corrected (GL_2(3)) tensor quiver constrains parity and
representation content but does not determine the interaction. The native
(E_6) incidence/sign tensor is indispensable additional structure.
This also resets the categorical boundary correctly. The classical ADE McKay
correspondence applies to finite subgroups of (SL_2(mathbb C)), with binary
octahedral giving affine (E_7). The exact clock group belongs instead to the
broader finite-(GL_2) setting, where generalized McKay/reconstruction
algebras exist but the ADE identification is not available.

A richer quiver-with-potential construction using the extra (E_6) tensor is
still open. What is now ruled out is the simpler claim that an affine-(E_7)
McKay graph or the finite symmetry alone determines the committed cubic.

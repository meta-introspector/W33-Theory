# Pass 11023 — corrected GL2(3) directed tensor saturation

Producer: `analysis/w33_pass11023_mckay_e7_clock_saturation.py`
Certificate: `data/w33_pass11023_mckay_e7_clock_saturation.json`
Regression: `tests/test_w33_pass11023_mckay_e7_clock_saturation.py`

The filename is retained for provenance. Pass 11025 corrects the original
McKay/E7 interpretation.

Pass 11022 decomposes the exact 24-dimensional signed clock carrier under
(GL_2(3)). The GAP character table has two faithful two-dimensional
characters whose order-eight values are

[
pm isqrt2,
]

not (pmsqrt2). With those cyclotomic values restored, the faithful tensor
matrix is not symmetric.
Let (A) denote tensoring irreducibles by one faithful 2D irrep and (B)
tensoring by its conjugate. The exact result is

[
oxed{B=A^T,qquad A
e A^T.}
]

Each quiver has 14 directed edges and is strongly connected. The underlying
undirected graph has 11 edges, so it is **not** the seven-edge affine-(E_7)
tree.

For (A), the arrows are

[
egin{aligned}
1&	o2_a, & det&	o2_b, & 2_q&	o4,\
2_a&	odet+3, & 2_b&	o1+3_t,\
3_t&	o2_a+4, & 3&	o2_b+4,\
4&	o2_q+3_t+3.
end{aligned}
]
The irreducible dimension vector remains

[
d=(1,1,2,2,2,3,3,4),
]

and the ordinary tensor-dimension identity holds:

[
oxed{Ad=Bd=2d}.
]

The exact clock multiplicity vector is still

[
m=(2,1,0,0,0,1,2,3).
]

Remarkably, the strongest saturation identity from the original pass survives
the correction unchanged. For **either** faithful 2D irrep (S),

[
oxed{
Sotimes V_{24}
=
3(2_qoplus2_aoplus2_boplus3_toplus3oplus4).
}
]
The central (12+12) split also survives:

[
Sotimes V_+
=3(2_aoplus2_boplus4),
]

[
Sotimes V_-
=3(2_qoplus3_toplus3).
]

Likewise

[
A(A+I)m=3d.
]

These are exact representation-ring identities of (GL_2(3)). They are not
classical ADE McKay statements.

The group itself is (GL_2(3)=2^+S_4), SmallGroup(48,29), with 13
nonidentity involutions. Binary octahedral is the non-isomorphic minus cover
(2^-S_4) and is the group occurring in the classical (SU(2)) affine-(E_7)
McKay correspondence. The exact clock saturation found here therefore stands
on its own, without an (E_7) group identification.

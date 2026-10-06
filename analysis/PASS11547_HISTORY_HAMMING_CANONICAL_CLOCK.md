# Pass 11547 — Finite Minkowski history is exactly the ternary Hamming cube

**Status:** exact finite theorem, with explicit conjugacy and clock factorization.  
**Boundary:** no continuum-space-time or physical-subsystem promotion.

Passes 11540--11541 identified the 27-event history chart with (VO(3,3)) and
closed its four quadratic shells as a formally self-dual association scheme.
There is a stronger description.

## 1. The interval is Hamming weight mod 3

Write the paper's finite interval as
[
q(t,x,y)=t^2-x^2-y^2.
]
Over (mathbf F_3), define
[
H(t,x,y)=(u,v,w)=(x+y,;x-y,;t).
]
Then, identically,
[
q(t,x,y)=u^2+v^2+w^2.
]
But in (mathbf F_3), every nonzero element has square (1). Therefore
[
oxed{q(t,x,y)=operatorname{wt}_H(u,v,w)pmod 3.}
]

Because a length-three ternary word has Hamming weight (0,1,2,) or (3),
the quadratic value plus the identity relation recovers the *entire* Hamming
distance:
[
oxed{
egin{array}{c|c}
	ext{history separation} & 	ext{Hamming distance}\ hline
x=y & 0\
q=+1 & 1\
q=-1=2 & 2\
q=0, x
e y & 3
end{array}}
]
Thus the paper's finite labels "timelike / spacelike / null" are, objectwise,
the distance-(1/2/3) relations of the ternary Hamming scheme after this
linear coordinate change.

In the repository's symmetric-matrix coordinates
[
q(a,b,c)=ac-b^2,
]
one direct version is
[
(a,b,c)mapsto
(2a+b+c,;2a+2b+c,;2a+2c).
]

## 2. Pass 11541 is literally (H(3,3))

The ordinary Hamming scheme (H(3,3)) has valencies
[
(1,6,12,8)
]
for distances (0,1,2,3).  Reordering relations as
[
(0,3,1,2)
]
gives exactly the Pass-11541 shell valencies
[
(1,8,6,12).
]

Its q-ary Krawtchouk eigenmatrix is
[
P_H=
egin{pmatrix}
1&6&12&8\
1&3&0&-4\
1&0&-3&2\
1&-3&3&-1
end{pmatrix}.
]
Applying the same ((0,3,1,2)) row/column permutation gives
[
egin{pmatrix}
1&8&6&12\
1&-1&-3&3\
1&-4&3&0\
1&2&0&-3
end{pmatrix},
]
exactly the eigenmatrix frozen in Pass 11541.

So "formally self-dual like a Hamming object" can be strengthened to

[
oxed{	ext{the history shell scheme is the ternary Hamming scheme }H(3,3).}
]

## 3. The null graph is (K_3	imes K_3	imes K_3)

Null separation is Hamming distance three: all three coordinates change.
Hence the null adjacency matrix factors as
[
oxed{
A_{m null}=(J_3-I_3)^{otimes3}.
}
]

This makes its spectrum immediate.  Since (J_3-I_3) has eigenvalues
(2,-1,-1), tensoring three copies gives
[
8^1,quad(-4)^6,quad2^{12},quad(-1)^8,
]
which is exactly the frozen history spectrum.

The null graph is therefore the direct/tensor product of three (K_3)'s,
the distance-three relation of (H(3,3)).

Even more strongly, the null graph alone reconstructs the other Hamming
relations:
[
24A_1=-A_3^3+9A_3^2+18A_3-64I,
]
[
12A_2=A_3^3-3A_3^2-24A_3+16I.
]
So (A_{m null}=A_3) generates the entire four-dimensional Bose--Mesner
algebra.  Any automorphism of the null graph automatically preserves the
ordinary Hamming graph.

## 4. The (1296) history symmetries are the full Hamming automorphism group

The standard automorphism group of (H(3,3)) is
[
S_3wr S_3,
]
with order
[
6^3cdot6=1296.
]

For the alphabet (mathbf F_3), every one-symbol permutation is affine:
[
zmapsto pm z+b.
]
Therefore
[
S_3wr S_3
=
mathbf F_3^3times(C_2^3times S_3).
]
In Hamming coordinates the linear stabilizer is exactly the signed
permutation group of order (48), i.e.
[
O(3,3)cong C_2^3times S_3.
]
Hence
[
oxed{
operatorname{Aut}(mathcal H_{m null})
=
operatorname{Aut}(H(3,3))
=
S_3wr S_3
=
3^3{:}O(3,3),
}
]
recovering the (1296) group of Pass 11540 as the **full** null-history
automorphism group, not merely a natural subgroup.

A recent independent preprint by S. M. Mirafzal, *The automorphism groups of
the Hamming-like graphs* (2026), studies the general direct product
(M(n,m)=K_m	imescdots	imes K_m) and proves
(operatorname{Aut}M(n,m)=operatorname{Aut}H(n,m)).  Pass 11547 does not
claim that general theorem; here the equality also follows internally because
(A_1) is a polynomial in (A_3).

## 5. The interval picks a canonical three-trit frame

Among the non-null projective lines of the quadratic space there are exactly
four mutually orthogonal three-line frames, reproducing Pass 11182.  Their
norm patterns are
[
1	imes(1,1,1),qquad3	imes(1,2,2).
]
There is therefore a **unique** orthogonal projective frame all of whose axes
have (q=1).

It is precisely the preimage of the three standard Hamming coordinate axes.
Thus (q) canonically selects an unordered decomposition
[
mathbf F_3^3=L_1oplus L_2oplus L_3
]
into three ternary address coordinates.  What remains unfixed is exactly what
the Hamming automorphism group should leave unfixed: permutation of the three
coordinates and independent affine relabeling of each three-symbol alphabet.

This is a finite-address/mereology theorem.  It does **not** say Nature has
three independently physical qutrit subsystems.

## 6. The spectral clock becomes three identical one-qutrit clocks

For Fourier Hamming weight (w=0,1,2,3), the distance-three adjacency
eigenvalue is
[
lambda_w=2^{3-w}(-1)^w,
]
so the null Laplacian eigenvalues are
[
0, 12, 6, 9.
]

At the paper's clock time (t_*=2pi/9),
[
U=e^{-it_*L}
]
has phases
[
1, omega^2, omega, 1
=
omega^{-w}.
]
Each nonzero Fourier coordinate contributes one factor (omega^2). Hence
[
oxed{
U=Cotimes Cotimes C,
qquad
C=F_3^daggeroperatorname{diag}(1,omega^2,omega^2)F_3
=e^{-2pi i p^2/3}.
}
]

This identifies the canonical Hamming frame with one of the four
position/momentum-aligned local clock splittings found by Pass 11182.

There is no conflict with Pass 11166, which found one-trit exchange in the
paper's light-cone coordinate factorization.  The Hamming transform is a
global linear change of tensor factorization; it diagonalizes the quadratic
form and therefore changes what "local" means.

## 7. Why (3) is maximally sharp over (mathbf F_3)

For any (n),
[
sum_{i=1}^n z_i^2
equiv
operatorname{wt}_H(z)pmod3.
]
For (nle3), the positive distances (1,ldots,n) have distinct residues
modulo (3).  At (n=4), distances (1) and (4) both have residue (1).

Therefore
[
oxed{
n=3	ext{ is the largest ternary dimension in which this quadratic value,
plus the identity relation, resolves every Hamming distance.}
}
]

This is a genuine (q=3), dimension-(3) rigidity statement.  It is not an
explanation of the observed dimension of physical space-time.

## Prior-art / interpretation firewall

Standard material:
- Hamming association schemes and q-ary Krawtchouk eigenmatrices;
- (operatorname{Aut}H(n,m)=S_mwr S_n);
- direct products of complete graphs as distance-(n) Hamming relations.

Repo-specific result in this pass:
- the explicit Forty-Points history/Hamming conjugacy;
- the identification of Pass 11541 with literal (H(3,3));
- the null-graph polynomial reconstruction;
- the weld to Pass 11540's (1296) orthogonal symmetry;
- the unique all-(q=1) Hamming frame and Pass 11182;
- the exact three-factor spectral-clock decomposition.

No continuum Lorentz group, physical light cone, Einstein equation, measured
speed of light, or fundamental tensor factorization is inferred.

## Reproducibility

Producer:
`analysis/w33_pass11547_history_hamming_canonical_clock.py`

Certificate:
`data/PART_W33_PASS11547_HISTORY_HAMMING_CANONICAL_CLOCK.json`

Upstream exact owners:
- Pass 11166 — spectral clock as a Clifford propagator;
- Pass 11182 — clock mereology;
- Pass 11540 — (VO(3,3)) and the orientation square;
- Pass 11541 — quadratic-shell association scheme.

# Passes 11054–11060 — finite cubic-weld breakthrough

The previous frontier said that the E6 cubic mechanism was only a tangent-level compiler. This packet sharpens that boundary substantially.

## Finite gate closure

For every matter background (v), the cubic Jacobian (D_v) is exactly real skew-symmetric. Hence

[
U_v(t)=exp(tD_v)in SO(81).
]

More importantly for an algebraic certificate, the Cayley gate

[
C_t=(I+tD_v)(I-tD_v)^{-1}
]

is defined for every real (t), because the spectrum of a real skew matrix is purely imaginary. For (t
eq0),

[
C_t-I=2tD_v(I-tD_v)^{-1},
]

so

[
operatorname{Im}(C_t-I)=operatorname{Im}(D_v),qquad
ker(C_t-I)=ker D_v.
]

For both diagonal FI welds the Jacobian rank is 66 and its projection onto the dressed (S2+L) quotient has rank 54. Therefore the same full 54-dimensional symmetry-changing capacity survives at finite amplitude. The old "tangent only" linear obstruction is gone. What remains is physical pulse synthesis, not existence of a finite unitary map.

## Local structure

The center-plus-external generator is supported on a 20-regular Cayley graph of

[
K=H_{27}	imes C_3.
]

Its twenty connection elements form ten inverse pairs of order-three elements. Each pair defines a rank-54 skew factor consisting of 27 disjoint triangles. Thus one factor exponential splits exactly into 27 independent three-mode rotations.

The ten weighted factors do not commute pairwise, so their product is not the exact exponential of the sum. Their kernel geometry is nevertheless rigid:

- symplectic-zero factor pairs share a 9D kernel;
- symplectic-nonzero factor pairs share a 3D kernel;
- all ten factors share exactly one ray, the diagonal background itself;
- the full sum has a 15D kernel.

Hence fourteen of the fifteen global dark directions are **interference dark modes**: they are created by cancellation between active local factors, not by being locally inactive.

## Fourier structure

The diagonal background is itself only seven-sparse in the frozen (K)-Fourier basis. For the (c+p) orientation,

[
v_+
=
2L(0,0,0)
+rac{omega-1}{3}sum_{r=0}^2S1(1,r,0)
+rac{omega^2-1}{3}sum_{r=0}^2S2(2,r,0).
]

The (c-p) orientation swaps the external Fourier labels (1leftrightarrow2).

Even more strongly, the complete 15D dark kernel is transverse to all three 27D Fourier sectors (S1,S2,L): it contains no pure-sector vector and projects with rank 15 onto every sector. The dark space is therefore a genuinely mixed triple-graph subspace.

## Orientation

External phase reversal (R:(h,p)mapsto(h,-p)) gives

[
D_-=-R D_+R^T,
]

and consequently

[
C_-(t)=R,C_+(-t),R^T=R,C_+(t)^T R^T.
]

The two FI orientations are exact inverse-conjugate finite gates. This is a kinematic equivalence; it does not select a physical chirality.

## Current boundary

The finite representation/compiler problem is now closed algebraically through an explicit orthogonal gate family. The sharp remaining frontier is to realize that gate with allowed microscopic interactions while controlling noncommuting triangle-factor synthesis, leakage, loss, timing, and fault tolerance.

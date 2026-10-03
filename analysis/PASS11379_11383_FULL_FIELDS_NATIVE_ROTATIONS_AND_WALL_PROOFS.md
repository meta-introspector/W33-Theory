# Passes 11379–11383: full fields, native rotations and wall proofs

Reservation `8265b1cd1`; producer `w33_pass11379_11383_full_frontier.py`;
certificate `../data/w33_pass11379_11383_full_frontier.json`;
regressions `../tests/test_w33_pass11379_11383_full_frontier.py`.
The calculations advance five previously named boundaries. Four use exact
symbolic/rational witnesses; the native cut-algebra calculation uses the stored
floating-point tensor certificate and is explicitly numerical. The supplied
actions do not yet predict observed masses, couplings or vacuum energy.

## 11379: a complete normal Hessian on the three-ray field space

Let the columns of $X$ be three arbitrary complex three-vectors: **18 real
fields**, with no fixed phase ansatz. Put $g=X^\dagger X$ and
$B=g_{12}g_{23}g_{31}$. Choose the supplied positive sum-of-squares potential

$$
V=\sum_i(g_{ii}-1)^2+
(|g_{12}|^2-1/2)^2+(|g_{13}|^2-1/3)^2+
(|g_{23}|^2-1/3)^2+(\operatorname{Re}B-1/6)^2.
$$

It is CP even, invariant under a common $U(3)$ and independent column phases,
and coercive because its norm terms control all columns. At the prior frame
$u=(1,0,0)$, $v=(1,1,0)/\sqrt2$, $w=(1,i,1)/\sqrt3$, all seven squares vanish.
The exact $7\times18$ constraint Jacobian has positive leading Gram minors
(stored in the certificate). Therefore $H_V=2J^TJ$ has **seven positive normal
directions**. The twelve common-unitary/column-phase tangent generators have
rank eleven and are annihilated by $J$. They exhaust the kernel: no additional
flat ray directions survive.

The zero constraints also give $|B|^2=1/18$ and $\operatorname{Re}B=1/6$,
hence $B=(1\pm i)/6$. The Gram determinant is $1/6>0$. Column phases remove
two overlap phases; the remaining loop sign distinguishes exactly **two
CP-conjugate minimum orbits**. Equal positive-definite Gram matrices determine
frames up to a common unitary. For the prior supplied weights $(1,2,4)$ and
$(2,3,5)$, $\operatorname{Tr}[A,B_Y]^3/i=\pm1$.

This closes the full-ray normal-stability boundary of
[11374](PASS11374_11378_FIVE_TOE_FRONTS.md), while citing earlier CP toy
selection in [11286](w33_pass11286_misaligned_family_vacua.py) and the
[11326 flavor operator](w33_pass11326_noncommuting_flavor_operator.py).
It does **not** derive the seven targets, weights or degree-12 operator from
the native E6 action. The invariant's physical interpretation follows
[Jarlskog](https://cds.cern.ch/record/161795).

## 11380: compress the native nonanalytic cuts, not arbitrary counterterms

Replay the actual 22-field quotient Hessian, cubic tensor, 11-field Weyl mass
and Yukawa tensor from [11287](w33_pass11287_quotient_self_energy_matrix.py).
There are eight Goldstone directions and fourteen normal directions. Crucially,
sum cut dyads **within equal internal mass pairs** before counting matrices;
individual dyads depend on the choice of degenerate mass eigenvectors.

The 21 mass-pair groups supply 25 nonzero real/imaginary normal coefficient
matrices. After Frobenius normalization their span has numerical dimension
**seven**, with its seventh singular value above 0.1 and the eighth below
$10^{-10}$. Their generated real matrix algebra has dimension **eight** and
center dimension **five**. A generic center element resolves orthogonal
projectors of ranks **6, 4, 2, 1, 1**. The rank-six, rank-two and two rank-one
blocks are scalar; the rank-four block carries **two copies** of the full real
$2\times2$ matrix algebra. The producer constructs a normal-coordinate
intertwiner taking all seven channel matrices on this block to
$a_{2\times2}\otimes I_2$, with a numerical residual. This part is
noncommuting, so replacing every channel by diagonal weights loses information.

An independent orthogonal rotation of all 22 scalar coordinates, transforming
the internal cubic legs and external Yukawa leg consistently, reproduces the
same normal span. The certificate stores the normal embedding, seven channel
matrices and five central projectors, with residuals. At $s=0.002,0.02,0.2$
the absorptive normal matrix is positive definite, of rank fourteen.

This suggests a practical matching architecture: four scalar dispersive
functions and one symmetric two-by-two matrix function can describe the native
cut sector. It does **not** shrink the unrestricted 105-entry local hard block
of [11337](w33_pass11337_quotient_hard_ward_columns.py). Symmetry of the cuts
alone does not impose that symmetry on unknown UV counterterms. Neither exact
representation identification nor a full renormalized pole spectrum is claimed.
Both the identity and tree normal Hessian lie in the same seven-dimensional
span to numerical tolerance. A concrete absorptive-truncation denominator
$D=(0.02+0.003i)I-10^{-4}H_{\rm normal}+i\rho(0.02)$ is inverted by four
scalar reciprocals and one $2\times2$ inverse repeated twice; multiplying by
the original $14\times14$ denominator checks the result. This is a usable
reduced propagator architecture, with unknown real hard matching still absent.

## 11381: the actual native pair action fails lapse linearity after rotations

Use **all 80 vertices and 160 edges** of the prior native Levi graph, with
one frame per vertex. The original pair interaction is

$$
I_{ij}=\tfrac12[\det E_i\operatorname{Tr}(E_i^{-1}E_j)
+\det E_j\operatorname{Tr}(E_j^{-1}E_i)].
$$

Set boosts and shifts to zero at zero momentum and take
$E_i=\operatorname{diag}(N_i,R(\theta_i)e_i,1)$: a two-dimensional active
spatial plane and one spectator spatial direction inside the four-dimensional
vielbein. Exact matrix contraction gives

$$
I_{ij}=\frac{N_i+N_j}{2}
[a_{ij}\cos(\theta_j-\theta_i)+b_{ij}\sin(\theta_j-\theta_i)
+\det e_i+\det e_j],
$$

where $a=\operatorname{Tr}(e_j\operatorname{adj}e_i)$ and
$b=\operatorname{Tr}(J_2e_j\operatorname{adj}e_i)$.
On the native fundamental eight-cycle $(0,40,1,44,4,53,13,41)$, choose

$$
e_i=\begin{pmatrix}1+x_i/10&y_i/10\\y_i/10&1-x_i/10\end{pmatrix},
\quad (x_i,y_i)=(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1).
$$

All other frames are the identity. These are positive definite, rational
coframes. Only the cycle's eight $b_e$ are nonzero, each $\pm1/50$;
the actual incidence matrix satisfies $Bb=0$ exactly. Thus $N=1,\theta=0$
is a stationary rotational branch, even with every other native edge retained.
Fix the common rotation, write $D=B$ with its vertex-zero row removed,
and $U=|B|$. At this point,

$$
K_\theta=-D\operatorname{diag}(a)D^T<0,\qquad
G=\tfrac12U\operatorname{diag}(b)D^T,\qquad
H_N=-G K_\theta^{-1}G^T.
$$

The exact rational rank of $G$ is **six**. Consequently the reduced lapse
Hessian has exact rank six and is nonzero. The planar shift/boost Jacobian,
with common variables fixed at vertex zero, is an integer matrix after scaling
by 2000; reduction modulo prime 1000003 proves its full rank 158. Its edge
weights are $(\det(e_i)e_j+\det(e_j)e_i)/2$, positive definite, giving an
independent graph-Laplacian invertibility proof. Zero boosts
and shifts remain a stationary branch as planar rotations and lapses vary.
An independent finite-difference re-solve of rotations agrees with the
predicted lapse curvature. Direct regressions recover the pair expression from
the original four-dimensional determinant/trace action.

The mechanism is **rotational frustration**, not independent connection
holonomy. A tree allows each edge angle to extremize separately, yielding
$(N_i+N_j)\sqrt{a_e^2+b_e^2}/2$, a lapse-linear expression. Cycles allow a
nonzero divergence-free residual edge torque. Altering lapses changes its
weights and forces the stationary vertex rotations to move.
The earlier [11338 collinear boost theorem](w33_pass11338_native_boost_constraint_slice.py)
and [11376 flat-transport count](PASS11374_11378_FIVE_TOE_FRONTS.md) remain true;
neither implies general native lapse linearity. This supplies a native witness
of the danger described by [de Rham–Tolley](https://arxiv.org/abs/1505.01450),
not a new general multigravity no-go theorem. A full Dirac count and propagation
analysis are still necessary; rank six is **not** a count of six ghosts.

## 11382: exact positivity of the coupled wound-wall radial vector block

Use the supplied quantized winding family of
[11341](w33_pass11341_winding_lifted_wall.py) and the exact mixed zero of
[11346](w33_pass11346_winding_wall_mixed_vector_modes.py):
$1<q\le4/3$, $h=4(q^2-q)/5$, $c=q/5$, $b_w=5/4$,
$\rho=q(4-3q)/10$, and $F=C/b-\rho b^2/3-q^2/(2b^2)-h$.
The cap is positive in its interior by the endpoint polynomial check in the
certificate. For one torus axis the full radial form is

$$
Q=\int_1^{b_w}\left[F\chi'^2+\frac{c^2b^4}{2}v'^2
+2qc\chi v'+\frac{hc^2b^2}{F}v^2\right]db.
$$

The prior exact zero $v_0=-F/(qcb^2)$ is node-free in the interior and solves
$-(b^4v_0')'+2(hb^2-q^2)v_0/F=0$. Integration by parts and its ground-state
identity give the **exact** factorization

$$
Q=\int F(\chi'-qcv/F)^2db
+\frac{c^2}{2}\int b^4v_0^2[(v/v_0)']^2db\ge0.
$$

The boundary terms vanish at a regular pole ($v=O(b-1)$), and for either wall
parity: odd $\chi(b_w)=0,v'(b_w)=0$; even $\chi'(b_w)=0,v(b_w)=0$.
The odd kernel is precisely $(1/b-4/5,v_0)$; the even kernel is constant
$\chi$, $v=0$. Two axes give four radial vector zeros. This upgrades the
previous winding-wall numerical Ritz check to an exact full radial-block
positivity statement. Other harmonics, scalar breathing, wall bending and
nonlinear integrability of these vector zeros remain open.

**Correction to 11377:** write the wound scalar $\theta=y+a$ and metric fiber
$dy+V$. Under $y'=y+\lambda$, $a'=a-\lambda$, $V'=V-d\lambda$;
$da-V$ is invariant. Constant axion shifts can therefore be undone by constant
torus translations. With the dynamical metric and no boundary frame/charge
fixing the translation, they are not two additional physical moduli. The
fixed-background positive axion principal block remains correct.

## 11383: an actual fixed-coupling scale variation

At $q=4/3,h=16/45,\rho=0$, take $g(L)=L^2g_0$ while holding the coordinate
Maxwell flux, axion winding and **input couplings** fixed. With $K=M_{\rm Pl}^2$
the Euclidean action uses $-K\int R/2$, a GHY term $-K\int K_{\rm extr}$
on **each cap**, positive Maxwell and axion terms, and one explicit wall term.
Do not also count distributional wall curvature.

Divide by the common angular factor $(2\pi)^3$. The exact background integrals
are bulk EH $-2/3$, axion $+2/3$, Maxwell $4/3$, two-cap GHY $-2$, and wall
$4/3$. EH/GHY/axion scale as $L^2$, Maxwell as $L^0$, wall as $L^3$:

$$
\frac{S_E(L)}{(2\pi)^3}=\frac43-2L^2+\frac43L^3,
\quad S_E'(1)=0,\quad \frac{S_E''(1)}{(2\pi)^3}=4>0.
$$

Dropping GHY destroys stationarity. This establishes positive curvature along
the **global uniform-scale direction** for the fixed supplied action. It is
not positivity of nonuniform conformal modes or a total Morse index.

**Units correction to 11378:** repository $\kappa^2$ means the EH coefficient
$K$, not Newton's coupling. The restored relation is
$T_{\rm phys}=(64/75)K_{\rm phys}/L$ and
$K_{\rm phys}^2R_{\rm phys}/T_{\rm phys}^2=5/8$. With compatible other background
inputs, the junction fixes $L$ from supplied $K_{\rm phys}/T_{\rm phys}$.
W33 does not yet derive those dimensional inputs or observed vacuum energy.

## External checks used and open boundaries

Primary [Jarlskog](https://cds.cern.ch/record/161795) supports the CP-invariant
interpretation; [de Rham–Tolley](https://arxiv.org/abs/1505.01450) explicitly
explains the lapse dependence induced by auxiliary Lorentz rotations. These
are prior literature, not novelty claims. The elementary Gram-orbit and
ground-state factorization proofs are written out above and independently
tested. The tensor-algebra compression is a numerical result of the native
certificate, not an assumed symmetry theorem. The five steps leave substantial
physical gaps: UV-supported flavor targets and masses; hard 1PI matching;
a consistent gravitational constraint architecture; the remaining coupled wall
sectors; and a derivation of dimensional input parameters.

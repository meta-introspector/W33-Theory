# Pass 11218 corrected — cubic rank 78+3 and an abelian centralizer

> **Correction (Pass 11255).** The exact rank and centralizer computations in
> this pass remain valid, but the displayed three-dimensional kernel is **not**
> a Vinberg Cartan subspace. Its first stored basis vector is nonzero and
> satisfies $(18\,\mathrm{ad}\,v)^4\ne0$ and
> $(18\,\mathrm{ad}\,v)^5=0$. It is therefore nilpotent rather than
> semisimple. The former identification of these three vectors with $G_{26}$
> quotient coordinates is withdrawn.

## What the exact witness still proves

In the canonical coordinates $n=3i+a$, the stored nine-coordinate background
has

\[
\operatorname{rank}D_v=78,\qquad \dim\ker D_v=3.
\]

Exact rational nullspace computation gives three integral kernel vectors.
Their pairwise brackets vanish and the background is the second basis vector,
so the kernel is exactly a three-dimensional abelian centralizer of this
background. Deleting flat coordinates $1,2,16$ leaves a $78\times78$
principal minor with

\[
\det=2^{76}3^{24}5^2 7^2
 =\left(2^{38}3^{12}\cdot5\cdot7\right)^2\ne0.
\]

This exact minor certifies rank 78. The Vinberg classification independently
says that the order-three $E_8$ grading with fixed algebra
$E_6+A_2$ has theta rank three. Together these facts retain the generic
rank-78 conclusion. They do not make every regular centralizer a Cartan
subspace; semisimplicity is the missing condition that Pass 11255 tests.

## The decisive regression audit

Pass 11255 builds the full repository $E_8$ adjoint matrix of the first kernel
vector. After the exact integral scaling $N=18\,\mathrm{ad}(v)$, its powers
have respectively 469, 171, 16, 1, and 0 nonzero entries. Thus its nilpotency
index is five.

The same pass factors the fixed $78\times78$ restricted principal Pfaffian on
the old kernel coordinates $(x,y,z)$:

\[
-\frac{5112848641047920640}{5}
 H_6(y,z)^3 Q_2(y,z)^{10}(96x+5y-120z),
\]

where

\[
\begin{aligned}
H_6={}&y^6+18y^5z-4320y^3z^3-77760y^2z^4\\
     &-653184yz^5-2239488z^6,\\
Q_2={}&y^2+24yz+72z^2.
\end{aligned}
\]

The identity is checked exactly on the 820-point triangular interpolation grid
for total degree at most 39. Its multiplicities $(6,3),(2,10),(1,1)$ differ
sharply from the degree-33 Jacobian and 21 reflection hyperplanes of the
standard $G_{26}$ invariants. The Pfaffian is therefore not a disguised
$G_{26}$ reflection Jacobian.

## What the Vinberg classification still supplies

Reeder, Levy, Yu, and Gross classify the grading in Table 21, $E_8$ row 3b, of
*Gradings of positive rank on simple Lie algebras*. The abstract Cartan
subspace has dimension three, little Weyl group Shephard–Todd $G_{26}$, and
basic invariant degrees $6,12,18$. Pass 11256 constructs those standard
invariants directly. What remains open is an explicit **semisimple** embedding
of their coordinates into the repository's signed $3\otimes27$ tensor, plus a
physical kinetic metric.

The distinction from `data/w33_extended_clifford_g26_no_go.json` remains: the
genuine $G_{26}$ action belongs to an abstract Cartan subspace and is not the
full retained-phase extended qutrit Clifford group.

## Boundary

This pass certifies a regular rank-78 background, a three-dimensional abelian
centralizer, and an exact principal-minor witness. It does not certify the
displayed kernel as semisimple, as a Cartan subspace, or as physical vacuum
coordinates. Passes 11255–11259 preserve that boundary in every potential,
Hessian, CP, and supertrace calculation.

Evidence:

- `analysis/w33_pass11218_e8_cubic_cartan_kernel.py`
- `data/w33_pass11218_e8_cubic_cartan_kernel.json`
- `analysis/w33_pass11255_restricted_pfaffian_cartan_audit.py`
- `data/w33_pass11255_restricted_pfaffian_cartan_audit.json`
- `tests/test_w33_pass11218_e8_cubic_cartan_kernel.py`
- [Reeder–Levy–Yu–Gross, Table 21](https://abel.math.harvard.edu/~gross/preprints/PosRank.pdf)

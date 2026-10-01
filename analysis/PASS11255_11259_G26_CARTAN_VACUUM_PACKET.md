# Passes 11255–11259 — the Cartan correction and a scoped G26 vacuum sandbox

This packet executes five proposed vacuum-dynamics calculations and changes the
frontier at the first step. The three-dimensional kernel exhibited in Pass
11218 is an exact abelian centralizer, but it is not a Cartan subspace. The
remaining four calculations are therefore built in standard $G_{26}$
coordinates and are explicitly conditional until a semisimple embedding into
the repository's signed $3\otimes27$ tensor is constructed.

## 11255 — the degree-39 Pfaffian factors, and the old slice is nilpotent

Let $(x,y,z)$ multiply the three stored Pass 11218 kernel vectors. The fixed
$78\times78$ principal Pfaffian is exactly

\[
-\frac{5112848641047920640}{5}
H_6(y,z)^3Q_2(y,z)^{10}(96x+5y-120z),
\]

with

\[
\begin{aligned}
H_6={}&y^6+18y^5z-4320y^3z^3-77760y^2z^4\\
     &-653184yz^5-2239488z^6,\\
Q_2={}&y^2+24yz+72z^2.
\end{aligned}
\]

Evaluation on the 820-point triangular grid $(x,y,1)$, $x,y\ge0$ and
$x+y\le39$, proves the total-degree-39 identity over $\mathbb Q$.

The decisive test uses the full repository $E_8$ bracket. For the first stored
kernel vector, $N=18\,\mathrm{ad}(v)$ is an exact integral $248\times248$
matrix. Its powers have 469, 171, 16, 1, and 0 nonzero entries. Hence
$N^4\ne0$ and $N^5=0$: the vector is nonzero nilpotent. A Vinberg Cartan
subspace consists of commuting semisimple elements, so the displayed kernel is
not one.

This retains the exact rank 78, three-dimensional commuting centralizer, and
principal-minor witness. It withdraws the claim that these particular vectors
are $G_{26}$ quotient coordinates. The factor degrees and multiplicities
$(6,3),(2,10),(1,1)$ also rule out identifying this Pfaffian with the genuine
$G_{26}$ reflection Jacobian.

## 11256 — a genuine standard-coordinate G26 potential

The literature realizes $G_{26}=G_{25}\rtimes\mathbb Z_2$. The standard basic
invariants have degrees $6,12,18$, with the degree-18 invariant equal to the
square of the $G_{25}$ degree-9 invariant. The producer implements those
polynomials and factors their Jacobian. Its degree is 33 and it resolves into
21 complex reflection hyperplanes with the expected reflection
multiplicities.

At the regular standard-coordinate point $p=(1,2,3)$,

\[
(u_6,u_{12},u_{18})(p)=(-1716,1389996,11957764),
\]

and the Jacobian determinant is
$-4033247385424158720$. The controlled orbit-distance potential

\[
V=\sum_i\left(\frac{u_i}{u_i(p)}-1\right)^2
\]

has exact Hessian $2J_{\rm norm}^{T}J_{\rm norm}$ and eigenvalues

\[
2.01443132984,\quad62.5328323663,\quad1269.92337680.
\]

This constructs an invariant mathematical potential with a positive local
minimum. It does not choose physical coefficients or supply a kinetic metric.

## 11257 — the 1+10+16 pullback is a useful counterfeit

To expose what the missing embedding controls, the producer forms the unique
minimum-Frobenius-norm quadratic form in the span of the old kernel matrix $K$:

\[
M=K(K^TK)^{-1}H(K^TK)^{-1}K^T,
\qquad K^TMK=H.
\]

The compression residual is $1.32\times10^{-11}$ and $M$ has rank three.
Every one of the 27 choices of a Schläfli reference line was then decomposed
into its $1+10+16$ self/intersecting/skew blocks. All 27 block-trace profiles
are distinct. The one-block trace ranges from 0 to 6.0881, the ten-block trace
from 0.1180 to 10.0503, and the sixteen-block trace from 2.1540 to 11.6841.

This reference dependence, together with the nilpotent kernel, blocks a mass
reading. A physical pullback needs a semisimple Cartan embedding, a kinetic
metric, and either a selected $SO(10)$ line or a proof of reference
independence.

## 11258 — conditional cubic phase produces a CP-odd target unitary

Normalize $(1,2,3)$ as a qutrit control state $|v\rangle$ and use

\[
U=(T\otimes T)\operatorname{SUM},\qquad
T=\operatorname{diag}(1,\zeta_9,\zeta_9^{-1}).
\]

Postselecting the same control gives $K_v=\langle v|U|v\rangle$. Its singular
values are

\[
0.815724995,\quad0.810401739,\quad0.421712882.
\]

The polar unitary $W$ has

\[
J=\operatorname{Im}(W_{00}W_{11}W_{01}^*W_{10}^*)
 =0.00952765713467.
\]

Thirty-two random row/column rephasings preserve $J$ to
$5.3\times10^{-18}$; complex conjugation reverses its sign. SUM alone and a
target-only $T$ give zero, as does the symmetric control $(1,1,1)$. This is an
explicit postselected CP interface. It is not a CKM matrix because its control
coordinates have not been embedded into the W33 matter carrier.

## 11259 — discriminant zero modes do not supply CP or vacuum cancellation

For an orbit-centred unnormalized potential
$V_p=\sum_i(u_i-u_i(p))^2$,

\[
\operatorname{Hess}_p(V_p)=2J(p)^TJ(p).
\]

Thus reflection ramification produces exact Hessian zero modes. The audit
checks five singular points as well as the regular target. It also separates
that statement from CP: $(1,1,2)$ lies on the reflection discriminant but has
conditional $J=0.001625567511\ne0$.

Finally, all six nontrivial boson/fermion signings of the three normalized
Hessian modes were exhausted. No signing cancels any supertrace moment from
one through six. The smallest proxy

\[
\frac12\left|\sum_i(-1)^{F_i}\sqrt{\lambda_i}\right|
\]

is 13.1544536404. This is a negative result for the controlled three-mode
model, not a cosmological-constant calculation: no complete spectrum or
statistics assignment has been derived.

## Frontier after the five computations

The abstract rank-three $G_{26}$ classification remains intact. The missing
object is now sharply defined: an explicit commuting **semisimple**
three-plane in the signed $3\otimes27$ bracket and its coordinate map to the
standard reflection representation. Until that map exists, invariant
potentials can be studied exactly, but their pullback to masses and mixing is
conditional. The CP circuit supplies a separate, viable mechanism whose
physical use depends on the same embedding. The three-mode supertrace does not
hide a vacuum-energy cancellation.

Evidence:

- `analysis/w33_pass11255_restricted_pfaffian_cartan_audit.py`
- `analysis/w33_pass11256_g26_invariant_potential.py`
- `analysis/w33_pass11257_trinification_hessian_firewall.py`
- `analysis/w33_pass11258_g26_cubic_cp_interface.py`
- `analysis/w33_pass11259_g26_discriminant_supertrace.py`
- `data/w33_pass11255_restricted_pfaffian_cartan_audit.json`
- `data/w33_pass11256_g26_invariant_potential.json`
- `data/w33_pass11257_trinification_hessian_firewall.json`
- `data/w33_pass11258_g26_cubic_cp_interface.json`
- `data/w33_pass11259_g26_discriminant_supertrace.json`
- [Reeder–Levy–Yu–Gross, Table 21](https://abel.math.harvard.edu/~gross/preprints/PosRank.pdf)
- [JHEP 08 (2022) 264, Appendix A](https://scoap3-prod-backend.s3.cern.ch/media/files/72291/10.1007/JHEP08(2022)264_a.pdf)

# Passes 11384–11388: execute the five native-dynamics directions

Reservation `82b94e570`. Producer `w33_pass11384_11388_native_dynamics.py`;
packet `../data/w33_pass11384_11388_native_dynamics.json`; compact native input
`../data/w33_pass11384_native_inputs.json`; nine independent regressions in
`../tests/test_w33_pass11384_11388_native_dynamics.py`.
This executes all five directions from 11379–11383. It does not claim that
their outstanding physical problems are solved. In particular, native invariant
couplings, EH scales and UV matching remain inputs.
The CP completion treats both compact E6 and family SU3 as gauged; its moments
are part of a new supplied potential. The hard-matching audit concerns the
previous E6 quotient action, and the hub, wall and transmutation calculations
name their own actions. They are separate tested components awaiting a common
UV construction. The scalar CP order parameter is not a compiled quark
Yukawa invariant or an observed CKM prediction.

## 11384 — native invariant CP selection on all 162 real fields

The three-ray potential in 11379 cannot simply be relabelled an E6 model.
Its individual column overlaps are not invariant under the native family
$SU(3)$ mixing the columns. Instead use the **existing global signed-cubic
circuits** $I_6,I_{12},I_{18}$ of
[the global covariants](w33_20261001_global_e6_cartan_covariants.py) and
[the degree-18 completion](w33_20261001_degree18_phase_completion.py).
Let $\mu_a=\Phi^\dagger H_a\Phi$ be the 86 compact $E_6\times SU(3)$ moments.
The new supplied real-coupling potential on the whole complex 81-space is

$$
V=\frac12\sum_a\mu_a^2+|I_6|^4-|I_6|^2+(\operatorname{Re}I_6)^2
+\frac12\operatorname{Re}I_6
+|I_{12}-I_6^2|^2+|I_{18}-I_6^3|^2.
$$

Completing squares gives

$$
V+\frac5{16}=\frac12\sum_a\mu_a^2+
(|I_6|^2-1/2)^2+(\operatorname{Re}I_6+1/4)^2
+|I_{12}-I_6^2|^2+|I_{18}-I_6^3|^2\ge0.
$$

Consequently its invariant minima are **derived from the displayed real
couplings**, rather than prescribed ray overlap targets:

$$
I_6=(-1\pm i\sqrt7)/4,\qquad I_{12}=I_6^2,\qquad I_{18}=I_6^3.
$$

The native degree-six invariant's imaginary part is CP odd and gauge invariant.
The common phase symmetry remaining after its real linear term is $Z_6$,
which fixes $I_6$; it cannot identify these two branches. This differs from
the earlier balanced T selector, whose real invariant triple identified its
conjugate by the finite gauge quotient. No exhaustive classification of
generalized CP outside the named gauge/phase group is asserted here.

The producer discovers and Newton-refines a regular point in the prior exact
D-flat Cartan plane, then evaluates the actual global 81-field circuits.
The native compact-generator orbit has rank **78**. The complete canonical
162-real-field Hessian has **84 positive normal eigenvalues**, with smallest
about **0.0295**, and 78 gauge zeros. This uses analytic moment/invariant
Jacobians, not a six-field Hessian promoted to 162 fields. The certificate
stores the actual field and Cartan coordinates, the full normal spectrum and
gauge-Hessian residual. Regressions check finite gauge covariance and analytic
gradient differences independently. The compact input archives the prior 45
signed triads and 78 integer generators with source hashes; all 156 defining
lower/dual tensor-generator checks are replayed exactly.

This is a native invariant **completion**, not a derivation of its real EFT
coefficients or of observed Yukawas. The high-degree operators need a UV
origin. Full-field stability is numerical, not an interval certificate. Earlier
flavor selectors and invariant circuits retain ownership; this adds their
native CP-breaking coupling and full-field witness.

## 11385 — distinguish hard matching from commuting with cuts

The prior numerical cut algebra is $\mathbb R^4\oplus M_2(\mathbb R)$,
acting on scalar blocks of ranks 6, 2, 1, 1 and on a rank-four block
$\mathbb R^2\otimes\mathbb R^2$. Its noncommuting factor acts on the first
two-space; the second is its multiplicity. Four different matching hypotheses
have four different dimensions:

| Assumption on a symmetric 14-normal hard matrix | Free entries |
| --- | ---: |
| Stationary Goldstone Ward columns only | 105 |
| Preserve the five central projectors | 36 |
| Commute with every cut-algebra matrix | 29 |
| Be invariant under the full orthogonal commutant | 7 |

The third condition is **not** the fourth. On the doubled block, commuting
with cuts allows $I_2\otimes b$, whereas the protecting multiplicity symmetry
permits $a\otimes I_2$. The full real commutant has dimension 46. Its
orthogonal group is $O(6)\times O(2)\times O(1)\times O(1)\times O(2)$;
invariance under it leaves four scalar coefficients and one symmetric
two-by-two matrix. The producer embeds an actual random element of this
group in native normal coordinates and checks it commutes with all seven cuts.

An explicit normal counterterm connecting two central blocks breaks this
protection while leaving stationary Goldstone columns zero. This proves why
[11337's Ward completion](w33_pass11337_quotient_hard_ward_columns.py) alone
cannot constrain hard matching to seven structures. The displayed compact
group protects cuts; whether the complete interacting UV theory has it must
be established before using it to predict poles. No unrestricted hard1PI
calculation is silently replaced by these symmetry assumptions.

## 11386 — a concrete acyclic spin-two architecture with native cyclic matter

The native pair-vielbein rotational counterexample and the older obstruction
to a symmetry-invariant spanning tree make a direct tree truncation unsuitable.
An explicit alternative adds **one invariant hub metric** $g_0$ to 80 site
metrics $g_i$. Couple only $(g_0,g_i)$ through Hassan–Rosen bimetric potentials,
with equal spoke coefficients and individual EH terms. The spin-two graph is
the 81-node star, with 80 edges and no cycle. Every native site automorphism
extends by fixing the hub; no native symmetry is lost by selecting a root site.

Retain the original 160 native edges as **internal matter couplings**, with
80 real scalars all minimally coupled to the hub metric:

$$
S_m=-\frac12\int\sqrt{-g_0}\left[
\sum_i(\nabla\chi_i)^2+g^2\sigma^2\chi^TL_{\rm native}\chi\right].
$$

This separates the two graphs explicitly. Coupling matter to multiple site
metrics would require another ghost audit; it is not part of this design.
The equal-spoke linear species Laplacian has spectrum $0^1,1^{79},81^1$.
The native internal matter Laplacian retains
$0^1,8^1,(4-\sqrt6)^{24},(4+\sqrt6)^{24},4^{30}$.
The known regular tree-bimetric branch has one massless plus 80 massive
spin-two fields, conditionally **402 polarizations**, not the old native
80-species spectrum. Parameters must give a positive Fierz–Pauli coefficient
for the displayed linear graph. The graph matrices are checked directly.

Constraint preservation uses published tree-multimetric results, reviewed in
[Flinckman–Hassan 2026](https://arxiv.org/html/2604.07625v1), and the earlier
[minimal matter constraint construction](w33_pass11296_covariant_matter_constraints.py).
This is a changed, fully specified architecture with supplied continuum,
hub, EH scales and couplings. It is neither an independent new Dirac theorem
nor an emergent-spacetime or observed-spectrum derivation.

## 11387 — fixed-flux bulk relaxation exposes a scalar wall negative direction

Fix **all action inputs** $K=1,h=16/45,T=64/75,\rho=0$, the three angular
periods, and magnetic Dirac flux. Vary the physical torus radius $u$ at each
regular pole and radius $v$ at the wall. The full regular bulk cap remains

$$
F(b)=C/b-q^2/(2b^2)-h,\quad
C=q^2/(2u)+hu,\quad c=q^2/(4u^3)-h/(2u).
$$

It satisfies $F(u)=0$, $F'(u)=2c$ and both bulk Ricci equations. The fixed
flux condition $q(1/u-1/v)/c=1$ determines the positive branch

$$
q=2u^2(1-u/v)+\sqrt{4u^4(1-u/v)^2+2hu^2}.
$$

Thus **q, C and c all relax**; this is not the old frozen-bulk wall shift.
With both cap GHY terms and one wall term, the exact on-shell-cap action is

$$
\frac{S(u,v)}{(2\pi)^3}
=q-\frac{v^2F'(v)+4vF(v)}c+\frac{T\sqrt{F(v)}v^2}c.
$$

Bulk EH and winding action cancel at $\rho=0$; fixed-flux Maxwell contributes
$q$. At $(u,v)=(1,5/4)$ the action is $2/3$ and **both derivatives vanish**.
The exact bulk-relaxed Hessian is

$$
H=\begin{pmatrix}400/49&-1536/245\\-1536/245&16064/3675\end{pmatrix},
\quad\det H=-13312/3675.
$$

Its eigenvalues are approximately $-0.282617,12.817039$. In particular, the
regular fixed-flux path $u=L,v=5L/4$ has exact second derivative
**$-100/147$**, verified by an independent finite difference. Nearby caps
remain smooth and positive. This means the Euclidean saddle **is not a local
minimum on this real regular-cap family**.
Independent quadrature of the original four-dimensional Ricci, Maxwell,
winding, two GHY and wall densities at the saddle and two nearby caps matches
the reduced action to better than $10^{-12}$; the negative direction does not
depend on assuming the reduced formula.

This path changes the cap profile at fixed flux and couplings; it is not the
pure metric homothety of 11383, whose curvature is positive. The exact radial
vector factorization of 11382 also remains true. Neither restricted positive
block proves full stability. The two-dimensional reduction underlying the
calculation includes the smooth-pole term and reduced GHY cancellation, and
is recorded in the certificate. A full gauge/conformal-contour treatment,
nonconstant wall bending, determinant and Lorentzian continuation are still
needed. The [primary negative-mode literature](https://arxiv.org/abs/1210.4740)
explains why Euclidean negative directions and their physical interpretations
must be kept distinct; no nucleation rate or Lorentzian instability is inferred.

## 11388 — native dimensional transmutation meets gravitational stationarity

For the native matter action above, $M_\chi^2=g^2\sigma^2L_{\rm native}$.
The graph's known traces are $\operatorname{Tr}L=320$ and
$\operatorname{Tr}L^2=1600$; these are graph arithmetic, not new physical
predictions. Applying the standard [Coleman–Weinberg calculation](https://doi.org/10.1103/PhysRevD.7.1888)
to 80 real scalars gives the logarithmic coefficient

$$
V=\sigma^4[A+B\log(\sigma^2/\mu^2)],\qquad B=25g^4/\pi^2>0.
$$

Its flat-space minimum is $v=\mu\exp[-1/4-A/(2B)]$, with
$m_\sigma^2=8Bv^2$ and $V(v)=-Bv^4/2$.
Thus $m_\sigma^2=-16V(v)/v^2$: in this truncation a nonzero stable scale does
not give zero vacuum energy. The renormalized quartic boundary condition A
still supplies the dimensional-transmutation integration constant.

A second, stronger test adds only the scale-invariant induced-EH coefficient
$F(\sigma)=\xi\sigma^2$, with constant $\xi>0$. Constant-field gravitational
vacua extremize the Einstein-frame potential $V_E=V/F^2$, not V alone:

$$
\frac{dV_E}{d\sigma}=\frac{2B}{\xi^2\sigma}>0.
$$

So the flat-space minimum fails gravitational stationarity in this named
one-loop constant-$\xi$ truncation. Running nonminimal couplings, extra
sectors, curvature operators or higher loops could change the result.
A supplied bare EH coefficient would change it while reintroducing a scale
input. This is an obstruction to a specific proposed scale mechanism, not
a general theorem forbidding radiative gravity.

## Interpretation

The most concrete advances are the native full-field CP-breaking witness and
the exact bulk-relaxed wall negative direction. Hard matching now has an
explicit protection criterion; an acyclic spin-two/native-matter design names
every metric coupling; dimensional transmutation has been tested against the
gravitational equation rather than assumed to solve it. Observed masses,
mixing, UV coefficients, dimensional inputs and vacuum selection remain open.

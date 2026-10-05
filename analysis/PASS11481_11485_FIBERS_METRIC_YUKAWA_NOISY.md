# Passes11481–11485: finite vacuum-shape fibers, coupled metric dynamics and noisy inference

Reservation `ef90a0dc7` was pushed before computing. Remote intake was
`6bc4ff928`, a refresh of the formula universe incorporating11476–11480;
no new physical model was introduced by that commit. The five requested
frontiers have constructed experiments, certificates and independent checks.
Their larger targets are not all solved. In particular, an exact stationary
moduli theorem and a spatially local chiral measure remain open.

Producer: `analysis/w33_pass11481_11485_fibers_metric_yukawa_noisy.py`.
Certificate: `data/w33_pass11481_11485_fibers_metric_yukawa_noisy.json`.
Regressions: `tests/test_w33_pass11481_11485_fibers_metric_yukawa_noisy.py`.
The JSON binds its three source certificates by canonical JSON SHA256.

## Ownership and novelty audit

11466 owns the numerical native vacuum and four soft modes.11476 owns its
numerical unitary Albert isotope and the reduction to two3-by3 complex
matrices.11389 owns the actual periodic spatial cover.11471 owns the native
160-edge phase masks.11480 owns full correlated two-block Steane recovery.
These objects are reused rather than reidentified from their dimensions.

A deeper search found that **11271 already owns the identical-Weyl Yukawa
zero and its(27,bar6) Higgs repair**.11281 owns the renormalizable one-sextet
selection obstruction.11483 cites both and transports their interface to the
actual11466 numerical frame. It does not claim a new representation rule.
Relevant files read in full include
`w33_pass11271_chiral_symmetric_yukawa_completion.py` and
`w33_pass11281_family_sextet_selection.py`; the report11270–11274 supplies
the canonical SM interface and the changed beta-function inventory.

Result searches covered the index, recent reports, source scripts,
`w33_paper.tex` and `docs/index.html`, including exchange symmetry, family
sextets, common-left deformations, hollowisation, noisy syndromes, Regge
constraints and local-current reconstruction. Classical implicit-function,
unitary level-set, ADM/DeWitt, local cohomology and Bayesian inference facts
retain their existing ownership; no absolute novelty claim is made.

##11481: integrate two physical potential-value fibers; keep stationarity separate

Let the native fixed-stratum fields be Phi=E A and Psi=E B, with E from11476.
A common-left transformation A->U A, B->U B, U inSU3, preserves both
determinants, norms, A-transpose B-conjugate, and family moments. The remaining
frame moment contributions depend on the diagonals of A A-dagger and
B B-dagger. Consequently preserving both diagonals preserves the entire
restricted potential value, not merely a few low-degree invariants.

Use eight traceless Hermitian generators L_a and four independent diagonal
constraints. Their derivative is

    C[:,a] = (diag(i[L_a,A A†])[:2], diag(i[L_a,B B†])[:2]).

At the stored numerical vacuum the four singular values are
1.0893563,0.6598792,0.00015737,0.000086377. Its four-dimensional kernel has
two directions inside the internal gauge orbit. Removing the ten internal
gauge directions leaves singular values0.8692187,0.6340446 and two values
below1.1e-12: **two physical common-left fiber tangents**. Their projection
into the prior four-dimensional soft space is independently retained.

For each physical tangent, integrate U=exp(i[t L_free+sum y_a L_normal])
by solving the four nonlinear diagonal constraints, for t=+-0.001,+-0.003.
Eight finite controls have diagonal residuals below5.8e-15 and full native
potential changes below1.6e-14. Their nongauge chords are nonzero.

There is also an exact Gaussian-integer pair independent of the vacuum.
Exact ranks of C, the gauge tangent matrix, and the joined gauge/fiber matrix
are4,10,12. An explicit rational kernel verifies C N=0. This proves that the
regular locus with two physical fiber directions is nonempty. At any exact
rank-four point the implicit-function theorem supplies a local level manifold;
it does not certify that our numerical stationary point is exact.

**Crucial distinction:** full native gradient norms on the finite curves grow
to1.25e-8, from the prior4.4e-10 stationary residual. Flat potential-value
curves are not automatically curves of stationary minima. This packet does
not label all four soft modes exact moduli, nor infer masses from their tiny
numerical Hessian eigenvalues. An analytic common-left exponential-path Hessian splits as

    H_ab = sum_fields (delta_a mu) dot (delta_b mu)
           + sum_fields mu dot delta_a delta_b mu.

The positive Gram term has four nonzero eigenvalues:7.4594e-9,8.2564e-9,
0.395442,0.435577. Its kernel includes two left-diagonal gauge directions and
two physical fiber directions. The moment-curvature remainder has norm6.13e-10;
left-sensitive moment norms are about3e-6. Thus the four previously reported
soft directions have a concrete candidate decomposition: two potential-value
fiber directions and two weak restoring directions. The tiny numerical
stationarity defect and remainder keep exact vacuum classification open.
These generator-coordinate eigenvalues are not canonically normalized
physical masses.

[Damm–Fassbender](https://arxiv.org/abs/1910.08813), section2.4, supplies the
classical simultaneous-hollowisation connection for two Hermitian matrices;
it is a useful geometric route, not a W33 discovery.

##11482: a supplied covariant metric action with actual coupled evolution

Use the160 oriented edges and harmonic displacements h_e of the actual
80-vertex four-regular Levi cover. Introduce a homogeneous SPD3-by3 metric
G(t), a lapse N(t), scalar values u_v(t), and a reference cell volume c.
This is an explicit model choice. Put v=c sqrt(detG), Q=G^{-1} dotG and

    K = v ||dotu||²/2 + v[Tr(Q²)-(TrQ)²]/8,
    V = (1/2) sum_e v (u_b-u_a)²/(h_e^T G h_e),
    L = K/N - N V.

The spring action is a graph model on the inherited cover, not a derived
Einstein spatial-curvature term. Metric and scalar kinetics are coupled
through v and the actual edge lengths. Constants and units are supplied.

Under a constant orientation-preserving coordinate transformation J,
h->J h, G->J^{-T}GJ^{-1}, c->det(J)c, the action is invariant. Under a time
reparameterization t'=alpha t, velocities and lapse divide by alpha and the
integrated action is unchanged. A second independently chosen spatial basis
checks the covariance in regression tests. This is affine covariance and
one time-reparameterization symmetry, not full local diffeomorphism symmetry.

The lapse equation is **K+N² V=0**. Positive scalar/shear kinetics alone cannot
satisfy it for a nonconstant scalar configuration. The deliberately naive
positive-kinetic trial has lapse force-43.9621954. Including the negative
conformal DeWitt direction gives a genuine finite constraint solution; this
is not a claim that its sign is a physical ghost or that constraints are solved
in an Einstein theory.

Canonical momenta p=v dotu and
P=v[G^{-1}dotG G^{-1}-Tr(Q)G^{-1}]/4 yield

    H = ||p||²/(2v) + (2/v){Tr[(G P)²]-(TrG P)²/2} + V.

The producer evolves the coupled Hamilton equations at N=1, with metric force
including both spring backreaction and volume dependence. Five stored states
through time0.002 preserve H=0 within3.3e-11 and keep the smallest metric
eigenvalue positive. Independent energy derivatives check all four canonical
blocks: G,P,u,p. This is an actual inhomogeneous scalar trajectory coupled to
six homogeneous metric components; no inhomogeneous metric, local constraint
algebra, W33-selected action, Newton constant or continuum gravity is derived.
See [Bahr–Dittrich](https://arxiv.org/abs/0905.1670) for why discrete constraint
symmetry must be checked rather than inferred from continuum terminology.

##11483: prior-owned Yukawa rule, applied to the current vacuum and its fibers

For one identical left-Weyl species in(27,3), the Lorentz-contracted fermion
bilinear is symmetric in the combined species indices. The native E6 cubic
d_ABC is symmetric; the family epsilon_ijk is antisymmetric. Thus

    Sym_{Ai,Bj}(d_ABC epsilon_ijk)=0.

The stored tensor verifies this directly. This is the **11271 result**, not a
new obstruction. Distinct fermion species, a different E6 tensor, or a larger
family-Higgs representation can evade it; the statement is scoped to this
renormalizable self-Yukawa and this field inventory.

Reuse11271's H_C^{ij} in(27,bar6), with H symmetric in i,j. On the actual
numerical frame construct H_C^{ij}=sum_k(Phi_Ck) S_k^{ij}, using three declared
symmetric family spurions. The additional S must transform if covariance is
required: there is no invariant triplet-to-bar6 map hiding in this formula.
Then Y_{Ai,Bj}=d_ABC H_C^{ij} is an explicitly stored symmetric81-by81 matrix.
It has rank81 at this declared background. This masses the entire native81;
it is **not** the earlier canonical SM-preserving rank30 exotic interface.

The all81 singular values feed a supplied MSbar Weyl fermion loop at mu=10:

    V1 = -sum_r m_r^4[log(m_r²/mu²)-3/2]/(32 pi²).

Zero singular values, if present, use the continuous zero-mass limit.
No bare mass is added to singular values. The same portal is evaluated on all
eight tree-level potential-value fibers, with loop shifts retained separately.
The loop shifts reach3.07e-6 on the first fiber and6.70e-6 on the second,
while the native tree shifts stay below1.6e-14. Central first derivatives
at step0.001 are0.00096653 and-0.00221703 in the declared fiber parameters;
step0.003 independently controls the derivative. Thus this particular fixed
portal produces a shape tadpole despite tree-level potential-value flatness.
A fixed spurion can lift otherwise flat shape directions, but its coefficients
and finite matching conditions remain free; this does not predict masses or
select a physical vacuum. This81-Weyl inventory is distinct from11478's876
omitted Dirac triplets; their loop and beta-function certificates cannot be
silently combined.11271 already counts the UV cost of its extra scalar fields.

For the established E6/family-Higgs approach see
[Stech–Tavartkiladze](https://doi.org/10.1103/PhysRevD.77.076009).

##11484: support-limited anomaly primitive attempt; locality still fails to close

Use the same actual4D overlap projector and native hypercharges1,-4,2,-3,6,
with multiplicities6,3,3,2,1. Compute the analytic64-link derivative of the
summed16-site anomaly density. Its columns sum to zero and its gauge-direction
variation J B vanishes at numerical precision, with B the oriented gauge
boundary matrix. A separate random direction replays the overlap derivative.

For each perturbed link, solve B^T flow_l=-J_l using only links within graph
radius r of that link's initial site. This constructs an explicit finite-support
candidate for the differentiated anomaly-divergence primitive. The controls are:

| Radius | Maximum divergence defect | Gauge-derivative defect |
|---|---:|---:|
|0|6.20e-9|0|
|1|4.32e-9|1.12e-8|
|2|2.93e-9|7.27e-9|
|3|1.46e-9|3.91e-9|
|4|1.16e-15|1.32e-14|

The full-radius flow recovers the differentiated Ward equation. Smaller balls
retain measurable defects. A vanishing radius-zero gauge defect is vacuous:
that candidate is zero and fails the divergence equation. L2 has no room to
establish infinite-volume exponential tails. A background derivative flow is
also not a globally integrated nonlinear current.

[Lüscher](https://arxiv.org/html/hep-lat/9811032v2), equation5.8, reconstructs
an integrable local current from a **local gauge-invariant anomaly primitive**
and a radial curvature term, with finite-volume corrections. Our truncated
flows do not yet supply that primitive. The previous Coulomb construction
passes finite-patch Ward/integrability but does not replace this locality step.
Axial lattice covariance, curvature matching, nonlinear primitive integration,
nonzero flux sectors and the non-Abelian theory remain open. This is a retained
construction attempt with failed controls, not a completed chiral measure.

##11485: noisy syndrome inference with unconditional erasures and native leakage

Apply the untwirled native routed relay marginal at damping0.01 to a full
seven-qubit Steane block. Retain all64 syndrome-conditioned logical images and
the true syndrome probabilities for a maximally mixed encoded state. Add an
explicit independent bit-flip probability e on the six reported syndrome bits.
Compare the raw reported syndrome with maximum-a-posteriori inference using
that actual channel-dependent syndrome distribution.

A declared ideal final code-space verification accepts only a correctly
identified syndrome; rejected branches are stored as an orthogonal erasure
state. This verification is an extra supplied resource, not a simulated noisy
circuit. The resulting normalized2-to3 Choi matrix retains failure probability
and is trace preserving and positive. No postselection renormalization is used.
MAP maximizes acceptance under the stated logical prior; it need not maximize
quantum fidelity for every input or every channel.

| Six-bit readout error | Raw acceptance | MAP acceptance | Raw unconditional entanglement fidelity | MAP fidelity |
|---|---:|---:|---:|---:|
|0|1|1|0.973697|0.973697|
|0.001|0.994015|0.994015|0.967869|0.967869|
|0.01|0.941480|0.945530|0.916716|0.921238|
|0.05|0.735092|0.875305|0.715757|0.865219|

This is an explicit adaptive marginal experiment;11480's full correlated
two-block channel is preserved as the prior owner and is not renamed an
independent-channel approximation here. Adapting a joint correlated noisy
recovery circuit remains a separate task. See
[Fletcher–Shor–Win](https://arxiv.org/abs/0710.1052) for the established
channel-adapted error-correction setting.

Separately, actual11471 grouped phase masks on160 edge modes preserve the
12-cell while mixing its logical and leakage sectors. A0.1-radian phase on a
stored group gives mean logical leakage0.00124876353. The actual12-by12 unitary,
logical subchannel, group and noiseless inverse-echo residual are stored and
rebuilt independently. The exact inverse cancels this coherent error in the
ideal control model; it does not handle unknown phases, dissipative leakage,
noisy verification or fault-tolerant hardware thresholds.

## Validation and remaining boundary

The five producer sections PASS means their stated controls replay, including
negative controls. Eleven independent regressions pass in107.15seconds, including exact generic
fiber ranks, independent portal-loop replay and full-native weak-restoring
derivatives. The prior combined suite also passed all16 owned/parallel tests. Syntax and the forced-arithmetic
scan/selftest pass; result-index refresh and full batch intake accompany
publication. Their completed outputs are recorded in session notes. None of these checks
upgrades a numerical point, supplied metric/spurion, finite flow or ideal
verification into a solved TOE.

# Passes11280–11284: smaller Higgs fields, family selection, decay and metric/flux architectures

Reservation ff3447be2 followed the reviewed inventory refresh d5bb013b4. This
packet executes all five independent targets left by Passes11275–11279.
Per Continuity decision-ef36a1ef-df21-4675-a3ff-0f7ddaf468ec, that earlier
packet was published as21afbc093 with its unresolved physical boundaries.
The five scoped PASS certificates below establish constructed objects and their
checks; they do not establish a complete theory of everything.

## 11280 — factor the Higgs selector instead of buying two End27 fields

The prior signed-cubic selector N=d(v0)†d(v1) has rank5 at the canonical SM
point. Use three complex27-by5 fields U,V,W and an added auxiliary U(5) gauge
factor. Replace the old degree8 selector by squared constraints

    UV†=N, U†U=V†V=I5, ΛW=AU, (A+I)W=6U/Λ.

Together with the old Higgs Gram, cubic and AΦ constraints, these form a
nonnegative degree-four potential. A acts in the E6 representation; all three
matrix fields transform on the left under E6 and on the right under U(5).
At the reference U is an orthonormal basis of im N, V=N†U and W=AU/Λ.
The scalar count falls from3102 to996 real components. Eliminating the new
linearized constraints recovers the old selector differential. Its kernel is
66 E6 directions plus25 auxiliary unitary directions:91 gauge directions and
905 positive normal directions. This is an analytic elimination argument,
not a numerical996-dimensional Hessian. The desired SM stabilizer survives
as a diagonal subgroup; gauge matching requires additional couplings.

The combined E6 one-loop coefficient is13, versus−80 for the previous matrix
lift. The added SU(5) coefficient is29/6; its U(1), with unit scalar charge,
has coefficient−135. Thus E6 and auxiliary SU(5) remain asymptotically free. The existing family
SU(3) coefficient remains−79/3, and the auxiliary Abelian UV completion
also remains open. Rank5 is minimal only
for this factorization architecture. Global zero-locus equivalence, global
minimality, full gauge spectrum and quantum vacuum stability are not proved.
Producer: `w33_pass11280_factor_higgs_mediator.py`.

## 11281 — a family vacuum, a tree-level obstruction and a missed bilinear

For a complex symmetric family sextet F, a nonnegative invariant potential

    Σ(k=1..3)|Tr[(FF†)^k]−c_k|² + |det F−6|²,
    (c1,c2,c3)=(14,98,794),

selects the SU(3) orbit of diag(1,2,3). Newton identities fix the squared
Takagi singular values and the determinant fixes the remaining phase.
Its real12-field differential has exact rank4 and kernel the eight-dimensional
SU(3) orbit: modular rank lower bounds plus integer JT=0 give the rank
sandwich. This degree12 EFT selects an orbit, but its three targets are inputs,
not observed quark/lepton mass predictions.

There is also an obstruction for a single-sextet tree renormalizable potential.
The independent invariants are Tr S, (Tr S)^2, Tr S² and det F+h.c., S=FF†.
For positive singular values xi, stationary pair subtraction gives

    (xi²−xj²)[2κ+μ_eff xk/(xi xj)]=0.

Three distinct nonzero singular values force κ=μ_eff=0, leaving shape flat
rather than an isolated hierarchy. Quantum logarithms or extra fields are
outside this tree-level result.

The earlier pure-Higgs cubic spectator proposal vanishes on SM singlets.
A different contraction does not: F_eff^{ij}=v^A H†_A^{ij}. With neutral
v=e0 and H=e0⊗diag(1,2,3), this is precisely the full-rank family sextet.
Replacing F in χχ Sym(F³) gives a dimension9 operator suppressed by M_UV^5,
versus dimension6 for elementary F. If elementary F is retained, the
matching square ||ΛF−vH†||² has degree4. This is a new named contraction,
not a claim that the old vanishing cubic was wrong. A single aligned sextet
still gives commuting up/down Gram matrices and no CKM mixing. Misalignment
and absolute scales remain open.
Producer: `w33_pass11281_family_sextet_selection.py`.

## 11282 — all internal modulus cuts expose a radial width

The external sector computed is radial only, not the full22-by22 self-energy
matrix. The declared canonical fixed-slice EFT includes all22 internal real
scalars and11 Weyl modes, using the previous full-modulus mass profiles.
Vertices are partial_r M_scalar²/√2 and −8M_fermion/√2. Eight global
Goldstones have the Ward vertex M_radial²/√2. Their open decay cut alone gives

    Γ(h→8GG)=M_radial³/(8π r0²), r0=1 here.

This absorptive channel is independent of local subtraction counterterms.
Twice-subtracted spectral dispersion at spacelike s0=−m² avoids the
massless zero-momentum logarithm and declares the mass/kinetic matching.
For m=.01 the strict one-loop squared pole is approximately
.009382185308374077−i .0000333836529982515, with Γ/M=.00360625881154.
The pole is complex: the earlier real single-scalar result omitted open
channels. Its real shift used a different subtraction prescription and must
not be compared as if the schemes were identical. The independent integration
split changes the computed real self-energy by8.13e−20; this is numerical
consistency, not an interval error bound. Nonlinear quotient Kähler vertices,
heavy vectors, full external mixing and physical UV matching remain open.
Producer: `w33_pass11282_all_modulus_radial_self_energy.py`.

## 11283 — what extra fields actually allow full local metric variation?

A positive compatible symplectic metric G satisfies GJG=J. Pick a G-unit
vector u, set v=−JGu, wu=Gu,wv=Gv and P=wu wu^T+wv wv^T. The explicit map

    g=e^(2ω){e^(−χ)(−wu wu^T+wv wv^T)+e^χ(G−P)}

has Lorentz signature and determinant−e^(8ω). At G=I,u=e0, the integer
Jacobian has rank8 without χ,ω, rank9 with χ, and rank10 with both.
Two modular primes independently check these ranks. Thus the naive reflected
compatible metric misses two local variations; a relative plane scale and an
overall scale repair the local chart. Eleven source fields contain one local
redundancy. A second positive realization and symplectic frame change check
signature, determinant and covariance.

An explicitly declared real continuum base and Einstein–Hilbert action can
then vary all ten metric components locally. This does not derive the base,
Lorentz signature, Newton constant, cosmological constant or a ghost-free
coupled completion from the finite geometry. The positive compatible G target
has positive sigma-model kinetic metric; that fact alone does not settle the
whole coupled theory.

Legacy `w33_BREAKTHROUGH_AdS4_Siegel.py` and
`w33_BREAKTHROUGH_366_spacetime_emergence_Minkowski.py` claim that finite
Sp(4,3) embeds in real Sp(4,R) and forces spacetime. Reading finite-field
entries as real numbers does not define such an embedding. Those claims are
not adopted here; the conflict was surfaced to the user and the files remain
untouched pending their preference. The current paper's discrete AdS passage
already scopes its correspondence as an analogy. A rank-ten Jordan real form
elsewhere in the paper is a different object from this four-dimensional
metric chart.
Producer: `w33_pass11283_symplectic_metric_field.py`.

## 11284 — a literal closed history with one top flux

The W33 Levi graph has80 vertices and160 flags; its actual integer incidence
matrix has rank79, hence cycle rank81. Add an oriented handlebody thickening
and identity double: M3=#81(S1×S2). Add a history circle to get closed
M4=M3×S1. Its Betti numbers are(1,82,162,82,1). This is an explicit extra
topological construction, not the old clique complex or an intrinsic finite
spacetime.

In the standard minimal product CW model with N history ticks, D4 has the
circle incidence block and zero spatial block. Rank D4=N−1 and its kernel is
the all-ones top cycle. Therefore F=D4^T A3+(Q/N)1 has fixed total flux Q;
exact gauge variations leave Q unchanged. Multiple history lengths check the
integer rank and flux identity. The old clique-history H4=0 obstruction stands.

With added linear sequestering functions and fixed fluxes, the declared action
imposes Vol=Q/μ4 and <R>=−2μ4 Qhat/(M²Q), giving
ΔΛ=−κ²μ4 Qhat/(2M²Q). The conditional metric equation removes a constant
stress shift at fixed flux: κ²G=T−g<Tr T>/4−ΔΛg. This supplies a closed-sector
constraint architecture; it does not predict the flux ratio or the observed
cosmological constant. Compact gauge quantization would be another input.
Handlebody gluing, continuum dynamics, couplings and convergence remain open.
Producer: `w33_pass11284_closed_history_flux.py`.

## Corpus and primary-source checks

Result-specific searches covered RESULTS_INDEX, the papers, the newest site
cards, recent analysis markdown, Python and certificates: rank-five cubic
support, factor inventories996/905, Takagi selection, Goldstone width,
compatible metrics and handlebody/history flux. Owners11270–11279 are reused
explicitly; searches by bare integers were noisy and were narrowed to formulas
and contractions. No claim of novelty rests on a topic-only search.

Canonical E6 condensate spectra are published prior art in
[Goh et al.](https://arxiv.org/html/2505.07931v1), whose near-SUSY/global-family
model must not be conflated with the added gauge/Higgs architectures here.
The pole calculation uses the distinction between curvature and momentum
self-energy emphasized by [Martin](https://arxiv.org/abs/hep-ph/0502168).
The compatible symplectic target is the standard Siegel geometry;
[Ohsawa](https://arxiv.org/abs/1504.03963) supplies the primary context.
The handlebody identity double is the classical construction described in
[Birman–Johnson–Putman](https://academicweb.nd.edu/~andyp/papers/SymplecticHeegaard.pdf).
The added four-form architecture is checked against
[local vacuum-energy sequestering](https://arxiv.org/html/1505.01492v2).
These citations support background and scope, not the claim that finite W33
alone supplies the new continuum inputs.

## Validation

Five focused pytest regressions passed in65.59s. They independently check
finite E6/U5 covariance, family congruence and the mixed bilinear, threshold
and Goldstone width identities, a second metric realization/prime, and
variable-length closed-history flux invariance. Five producer certificates
report scoped PASS. Paper/PDF and final source/certificate intake are recorded
in the publication receipt after completion.

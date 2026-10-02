# Passes11275–11279: local Higgs alignment, spectator economy, poles and geometry tests

Reserved and pushed in b8ffb3b99 after reviewing and pulling c3f8f4f4c, the
formula-universe inventory refresh. All five requested investigations have
runnable producers and scoped certificates. Successful witness execution does
not mean their larger physical targets are all solved.

## 11275 — a polynomial with an isolated SM orbit locally

The earlier Pass11271 configuration is V=(e0,e1), A=Y in the canonical
signed27 representation. Per Continuity decision-b65c46f8, its exact
SM stabilizer and exotic-mass interface were established, while the selecting
potential was an imposed nonpolynomial orbit distance. This pass supplies
an explicitly gauge-invariant polynomial alternative, with the adjoint
reference rescaled to6Y. Take two complex27 Higgs vectors and a real adjoint:

    V_H=||V†V-I₂||² + Σ_i≤j||d(v_i,v_j)||² + ||AV||²
        + ||(A²+A−6I) N||²_F,  N=d(v0)† d(v1).

The key connection uses the already-certified cubic tensors, rather than
imposing a spectrum on the entire27. At the reference N has rank5 and image
indices17,18,19,20,22. Those are three charge2 and two charge−3 eigenvectors
of6Y, so(A−2I)(A+3I)N=0. Because N transforms as U N U†, this constraint is
E6 covariant even as both Higgs vectors vary. Its square has field degree8;
the other three terms have degree at most4. Appropriate dimensionful
coefficients and cutoff powers are understood. Every positive choice of
coefficients gives the same local kernel. This is a polynomial Higgs EFT,
not a renormalizable result.

There are186 real fields. The exact integer443-by186 constraint differential
has rank120 and its kernel is precisely the66-dimensional compact E6 gauge
orbit. The certificate uses a rank sandwich, not floating eigenvalues: modulo101
gives rank J≥120 and rank T≥66; integer JT=0 gives the matching upper bounds.
A second prime103 independently checks these ranks. The differential includes
both adjoint variations and the full variation of N through the two Higgs
vectors. Treating N as a fixed spurion would answer a different question.

Consequently the Hessian is positive in all120 normal directions. A transverse
slice and the implicit-function theorem isolate this gauge orbit locally among
zero-energy minima. This is stronger than merely naming a desired stabilizer:
there is an actual polynomial with no additional local flat directions.
Disconnected distant minima, economical alignment, natural coefficients and
selected vev scales remain open. A finite E6 rotation verifies covariance.

There is also a concrete renormalizable scalar lift, rather than a promise to
find one. Add two complex End(27)=1+78+650 scalars X,Z and replace the final
term by three nonnegative squared constraints:

    ||ΛX−N||² + ||ΛZ−AX||² + ||(A+aI)Z+(b/Λ)X||²,
    a=1, b=−6 in reference units.

Each constraint is at most quadratic in fields, so the full scalar potential
has degree at most4. At a zero, X=N/Λ and Z=AN/Λ²; the third equation is exactly
(A²+aA+bI)N=0. The first two linearized equations uniquely determineδX,δZ,
so the same66-dimensional gauge kernel survives in3102 real fields, with3036
positive normal directions. This is an analytic constraint-elimination proof,
not a numerical3102-dimensional mass calculation. Finite coefficients preserve
the zero locus; integrating out the added fields does not reproduce the original
off-shell potential exactly.

The UV price is severe and explicitly counted: T(End27)=2·27·T27=162 for each
complex scalar. The two extra fields change the combined E6 coefficient from28
to−80. Thus a power-counting-renormalizable scalar alignment exists, but this
large lift is not an asymptotically free or established viable UV model.
Finding a smaller mediator inventory is a substantive remaining target.

The roots, couplings and extra fields are inputs. The standard E6 breaking
problem is prior art; see [the explicit renormalizable supersymmetric model](https://arxiv.org/abs/1504.00904),
whose much larger field inventory must not be conflated with this real-adjoint,
non-supersymmetric polynomial construction.

## 11276 — replace28 by a composite of6

The anomaly-free chiral interface uses Weyl(27,3)+(1,10bar), scalarH(27,6bar)
and previously an elementary family28 scalarΣ. ReplaceΣ by an E6-neutral
family sextet F(1,6):

    Σ_eff^{ijklmn}=Σ_15 pairings F^{ab}F^{cd}F^{ef},
    L_mass=(c/M_UV²) χ_ijk χ_lmn Σ_eff^{ijklmn} + h.c.

This dimension6 operator gives mass scale c v_F³/M_UV². At F=I the exact
10-by10 symmetric mass matrix has determinant592704 and is positive definite.
For every invertible complex symmetric F, congruence by Sym³ of an invertible
3-by3 matrix preserves rank10. A nondiagonal positive F and a singular control
are tested independently. This is a representation-theoretic construction;
it still needs a mediator UV completion and a potential selecting F.

A tempting even smaller construction fails exactly: d restricted to the two
canonical SM-neutral27 directions is zero. Thus the E6 singlet made from three
H† cannot give this spectator mass at the SM-preserving Higgs point. An
independent E6-neutral sextet avoids that zero.

The actual family one-loop coefficient changes from−93/2 to−79/3; it improves
but remains non-asymptotically free. With H plusF, bE6=32. Adding the two
family-neutral complex27s and the real78 of11275 gives bE6=28; using a complex78
instead gives26. These inventories must not be mixed. The family Landau scale
and the spectator mass still depend on boundary couplings andUV input.

## 11277 — compute the momentum bubble rather than rename a curvature

Declare a separate single-real-scalar radial EFT, not the full eleven-complex-
modulus model:

    L=½(∂h)²−V(1+h/√2), V(r)=m²[81/(49r^14)−27/(7r^6)].

Set the renormalized tadpole tozero and curvature mass to V_hh at r=1.
Subtract momentum-independent loop terms; retain the momentum-dependent bubble:

    B(z)=−∫₀¹ log[1−z x(1−x)] dx,
    Γ_M(s)=s−M_curv² + g₃² B(s/M_curv²)/(32π²).

The sign follows directly by twice varying½Tr log of the Euclidean fluctuation
operator and Wick rotating. Below threshold, an independent closed form is
B(z)=2−2√(4/z−1) arctan[1/√(4/z−1)]. It agrees with direct quadrature and
B′(0)=1/6. The local kinetic coefficient is fixed toone; the bubble induces
additional momentum dependence. This is a stated renormalization prescription,
not an arbitrary fitted wavefunction factor.

For m=.01, M_curv²=.009257142857, g₃=−.150553135240. The strict one-loop pole is
.009243779639; the Dyson root of this one-loop inverse propagator is
.009243801274, with residue.998381370470 and residual below1.2e−16.
The squared pole is about0.144% below the curvature. The difference between
strict and iterated values is selected higher-order content, not two-loop
accuracy. The pole is below the two-particle threshold and real in this EFT.

The other modulus directions, Goldstones, Weyl fermions and gauge thresholds
are excluded, so this is not the full physical-modulus pole spectrum or an
observed particle mass. It supplies a working momentum-pole calculation and
regression target for that larger task. General pole-mass methodology is
published prior art: [Martin's self-energy and pole-mass treatment](https://arxiv.org/abs/hep-ph/0502168).

## 11278 — a literal refinement tower, with a metric obstruction

Starting from the symplectic4-space behind W33, lift to(Z/3^k)^4. Primitive
projective points have count40·3^{3(k−1)}. The level2 witness constructs all1080
points modulo9 and checks every reduction fiber has27 points. Its orthogonality
graph is116-regular. The complete59,049-element first congruence kernel consists
of I+3X with JX symmetric; all are explicitly checked modulo9.

The group-order formula is51840·3^{10(k−1)}, from smooth symplectic reduction.
This classical fact is not claimed as new. The repo already owns the base
symplectic-basis torsor and the dual-number adjoint kernel in Pass4937. The
mixed-characteristic extension is not asserted to split or to equal the
previous dual-number group simply because its order agrees.

The inverse limit is3-adic, not automatically the Archimedean torus used by
the previous Wilson test. In fact, an integral symplectic shear fixes e0 and
sends e2 to e2+e0. Invariance TᵀGT=G forces G00=0; no positive definite real
quadratic metric can be invariant under all such lifts. Choosing a Euclidean
metric therefore requires a frame, a symmetry reduction or a dynamical
symmetry-breaking field. The four module coordinates do not establish a
Lorentzian spacetime or select the physical dimension.

## 11279 — test the actual history complex before adding a flux constraint

The full W33 clique complex has cells(40,240,160,40). Its topology is prior art:
Pass1448 and the September21 ledger show it collapses to a bouquet of81 circles;
Pass1944 rejects an earlier topological flux interpretation. We cite these
owners rather than rediscovering their homology.

Here the explicit new candidate is its product with a periodic three-tick
history circle. Its120 four-cells have a600-by120 boundary matrix D4. The
D3⊗I block is injective because every triangular face belongs to only one
W33 tetrahedron. Hence D4 is injective over the integers, not only at the
checked prime101. This persists for the independently tested two- and four-tick
circles. In particular H4=0 and real H^4=0: there is no intrinsic closed
fundamental four-cycle carrying a topological flux.

For the named discrete topological action

    S_top=⟨σ(Λ),F⟩, F=D4ᵀ A3+F_background,

free variation of A3 requires D4σ(Λ)=0, which forcesσ(Λ)=0 cellwise here.
It fails to produce the nonzero rigid multiplier of the continuum mechanism.
The all-ones four-volume chain has nonzero boundary. Fixing appropriate
boundary3-forms, passing to relative cohomology or choosing a closed oriented
four-complex changes the problem and introduces additional physical data.

The [local covariant sequestering action](https://arxiv.org/html/1505.01492v2)
uses additional four-form sectors and retains a flux-dependent integration
constant. It is an external conditional mechanism; neither that constant nor
the requisite closed geometry is selected by the construction tested here.
This obstruction applies to this explicit W33 clique history, not every
possible emergent spacetime.

## Corpus and validation boundary

Searches covered the result index, earlier tensor/Higgs producers, signed27
certificates, the base symplectic torsor, dual-number and Hjelmslev results,
clique topology and flux retractions, the site, and the relevant source passages
of w33_paper.tex and its included body. Historical numericalGeV identities and
the spectral action on a supplied manifold do not provide a UV matching input
for the present poles. They are not adopted as derivations of absolute masses.

Five producers and five JSON certificates are paired with independent regression
controls. The report, newest site card and scoped paper ledger distinguish exact
local algebra, numerical propagator results and negative geometry tests. The
full TOE questions of selected scales, observed Yukawa matrices, complete
physical poles, emergent spacetime and cosmological constant remain open.

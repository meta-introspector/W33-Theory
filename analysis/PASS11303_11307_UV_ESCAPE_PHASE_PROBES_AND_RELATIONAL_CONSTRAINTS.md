# Passes11303–11307: UV escape, phase probes, total-jet boundaries and relational constraints

Reservation `bbe2c12a6` follows reviewed parallel inventory `ea5f87cca`.
Parallel reservation `e128b140e` was reviewed and integrated during the work.
All five producers write scoped PASS certificates. Nine independent regressions
pass in247.00s. These are explicit mathematical and conditional EFT results,
not a completed TOE or predictions of observed parameters.

## 11303 — One scalar suffices for two channels; the unchanged54 has a UV barrier

Pass11298's rank-one obstruction concerns one **fermion pair**, not every
mediator. A single complex scalar H in(27,bar6) admits the source

```
J = sum_ab mu_ab v_a F_b†,
V = ||M H + J/M||²,
L_Y = y d_ijk psi_i psi_j H_k + h.c.
```

Here mu has mass dimension one. Eliminating H gives
`C_ab = -y mu_ab/M²`, of rank two for invertible mu. Direct signed-E6 matching
checks the actual tensor, not just representation dimensions. The positive
square includes mandatory cross quartics. The E6/SO10 gauge coefficients are
(8,8); gauged family remains non-AF. Sparse contractions on324 real H components
and81 Weyls give `16pi² beta_y/y = 25y²-52gE²-8gF²` for this active Yukawa.
Full frame/mediator/portal running is not supplied by this equation.

The54 obstruction can be strengthened beyond the isolated fixed-ray test in
Pass11298. Let x=lambda1/gSO², y=lambda2/gSO² and tau integrate gSO²/(16pi²).
At bSO=8, the two pure-S ratio polynomials f1,f2 obey the exact identity

```
z = x + 3y/10,
f1 + 3f2/10 = (16580/159) z² -104z +36 + (56x-15y)²/159.
min_z [(16580/159) z² -104z +36] = 41736/4145 >0.
```

The trace-fourth-power controls r=1/10 and73/90 occur on explicit traceless
unit54 backgrounds. Weights23/32 and9/32 give r=3/10. Every other real scalar
block and mixed block adds nonnegative trace-of-Hessian-square terms to these
backgrounds. The unchanged inventory has no SO10-charged Weyls and therefore
no S Yukawa or negative fermion box. Consequently arbitrary scalar portals
cannot remove this differential inequality. The positive quadratic Riccati
comparison blows up at finite tau, while tau diverges toward the UV for bSO=8.
This excludes complete weak-coupling one-loop asymptotic freedom of the
unchanged inventory, including scalar portals. It excludes neither changed
matter nor asymptotic safety, strong dynamics or gravitational completions.

**Constructed matter escape.** Add eight E6/family-singlet Weyls in the real
SO10 vector10; two have identical S Yukawas and six have zero S Yukawa.
Real representations introduce no perturbative gauge anomaly. Then

```
bSO = 8/3,
16pi² beta_y/y = (41/5)y² -27gSO²,
y²/gSO² = 365/123.
```

Direct54-tensor contractions verify the Yukawa coefficient. The quartics gain
`4 nactive y² lambda_i` and beta_lambda2 gains `-4 nactive y4`.
A bounded ray has

```
lambda1/gSO² = 0.31514549843936723,
lambda2/gSO² = -0.12087307962152648,
minimum normalized quartic = 0.2171040005241291,
quartic ratio stability eigenvalues = -18.73307907, -88.95647030.
```

Thus the isolated gauge/Yukawa/two-quartic S block admits an asymptotically
free ray with two attractive quartic ratios and a tuned Yukawa ratio.
Exact Sturm screens find no real rays for total1–7 Weyls within the scanned
identical-active-Yukawa family. This is not a minimum over arbitrary matter.
The full frame/family/H portals, protection of spectator zeros, masses and
phenomenological decoupling of the added charged fields remain open.

**Alternative representation.** With K²=-I and [K,S]=0, L=SK is a real
antisymmetric45 and S=-KL. The alignment constraints become
`[K,L]=0`, `L²+KL+6I=0`, `Tr(KL)=0`, `AU=iUL`, together with `UK=-iU`.
They are quadratic, retain the111-dimensional orbit kernel under the local
invertible map, and leave1245 normal directions among1356 real fields.
SO10 covariance is checked in a second basis. One H gives bSO=26/3.
Full two-adjoint/portal UV stability is still open.

Conventions: [Luo–Wang–Xiao, eqs19 and33](https://arxiv.org/abs/hep-ph/0211440).
Context, not ownership of this54 result:
[Einhorn–Jones on complete AF with adjoint/fundamental scalars](https://arxiv.org/abs/1705.00751).

## 11304 — SM-singlet loops can see the four missing phases

Pass11299 proved that the canonical E6 determinant sees only Su+Sd. Instead
add two vectorlike E6-singlet family3+bar3 Weyl pairs, with real input couplings:

```
M_k = [[y(Fu†+k Fd†), M I], [M I,0]], k=1,2.
```

These interactions preserve the E6/SM stabilizer and are CP even. The six
invariants `Tr[((Fu+k Fd)(Fu+k Fd)†)^p]`, k=1,2 and p=1,2,3, have Jacobian
rank four on the old four Takagi-phase flats at epsilon=.18,.22,.30,.40.
The smallest scaled singular value at.22 is.01317404176. Using only p=1,2
missed a direction. The actual twelve-Weyl Coleman–Weinberg determinant
changes under a deformation preserving both old Gram matrices. An independent
SU3 basis change preserves its spectrum. Two vectorlike pairs are anomaly
free, cost4/3 family beta units and no E6/SO10 units.

This supplies named renormalizable probes and phase identifiability. It does
not establish a stationary vacuum with positive Hessian, radiative spontaneous
CP breaking, an observed mixing phase or a generated hierarchy. Masses,
couplings and finite matching terms remain input.

## 11305 — Matrix vector cuts work; scalar IR prevents an unqualified total match

The construction uses the **earlier11287/11290 condensate EFT**, not the
new11303 Higgs inventory. It retains all70 massive E6 vectors and22 scalar
coordinates. Actual derivatives of the vector mass matrix give VV vertices;
antisymmetric scalar-representation generators give VS derivative vertices.
Physical cut matrices are sums of positive outer products with phase-space
weights. For vector masses squared a,b at external invariant t, the weights
are proportional to

```
VV: beta [2+(t-a-b)²/(4ab)]/(32pi),
VS: beta [(t-a-b)²-4ab]/(16pi a).
```

Unequal VV species carry the appropriate factor two. Independent radial VV
replay gives .02169596627046357 against .021695966270463427 from the earlier
width formula. Three zero-momentum subtractions define the vector dispersion
matching. Adding it to11295's scalar/Weyl spacelike-matched bubbles gives14
conditional leading massive pole entries; active massive-block vector cuts
are below their on-shell thresholds.

**Goldstone correction:** invariant vector eigenvalues do not force every
mass-matrix derivative or VV/VS vertex to vanish. Quotient tangents include
E6 compensators rotating vector eigenvectors; both couplings can be nonzero.
The vector contribution Sigma(0)=0 is a subtraction condition, not a vertex
selection rule. The certificate's historical field name
`Goldstone_VV_projection_error` stores the nonzero coupling norm, not a failed
assertion of exact vanishing.

The actual combined scalar/Weyl/vector first jet is computed. Regulated
Hermite spectral polynomials match the first two derivatives at all vacuum
nodes. However the scalar Goldstone first mass variations produce a positive
semidefinite **rank11 log(delta) Hessian coefficient**. The vector Gram-kernel
regularity argument does not transfer to the scalar sector. A regulator-free
C² polynomial counterfunctional cannot be claimed at this tree vacuum.
Ward-consistent Goldstone resummation and full UV matching remain necessary.
The formal total spectral counterfunctional assumes covariant full mass maps;
full nonlinear scalar-background derivatives are not implemented. The regulated
jet audit and conditional dispersion scheme are separate calculations, not a
single completed renormalization prescription or physical pole prediction.

Relevant primary work:
[Martin–Patel](https://arxiv.org/abs/1808.07615),
[Braathen–Goodsell](https://arxiv.org/abs/1609.06977),
[Goodsell–Paßehr](https://arxiv.org/abs/1910.02094).

## 11306 — An enlarged native-graph constraint algebra closes exactly

Use the quadratic part of11301's80-site/160-edge Levi oscillator Hamiltonians
H_i=.5 z^T K_i z. Add80 canonical clocks (t_i,P_i), and define

```
A_i=J K_i,
U(t)=ordered product exp(t_i A_i),
B_i=(partial_i U)U^-1,
C_i=P_i + .5 z^T(-J B_i)z.
```

The right Maurer–Cartan identity
`partial_i B_j-partial_j B_i-[B_i,B_j]=0` gives `{C_i,C_j}=0`.
The clock momentum Jacobian is I80, so the constraints are independent;
Omega=sum c_i C_i is a nilpotent classical BRST charge. The relational
observable is z0=U^-1z. Dressed two-clock paths agree to1.12e-16 while the
undressed paths differ by7.45e-4 in the witness. The symplectic error is1.12e-16.

Phase space320 reduces to160 dimensions: the original80 oscillators remain.
This is an explicit parametrized graph theory, **not gravity**, no graviton
polarizations or derived spacetime. Ordering and clock choices are inputs;
quartic site interactions and a Lorentzian3D refinement are not implemented.
The Maurer–Cartan construction is standard; its native-graph realization is
what is built here. Context:
[parametrized field theory](https://arxiv.org/abs/1011.2463).

**Parallel intake control.** Pass11312's exact finite-depth probabilities are
prior-owned and not rederived here. Its inference from Haar-measure-zero
reversibility to *every sufficiently long magic circuit* is false: T^9=I,
and identity Clifford insertions give reversible circuits of arbitrary depth9m.
Measure-zero alone also does not prove random-walk convergence or an exact
geometric rate. Pass11309's evolving certificate initially lacked its new
affine-codimension-two controls; do not accept a producer/prose result without
its matching certificate. The parallel files are preserved, and this packet
records its own explicit counterexample.

The prepublication parallel regression run gave **7 passed,1 failed in86.66s**:
11309 raises a KeyError for `nonviolating_frames_affine_subspace_dims`.
Its exact count is retained as an asserted certificate result, not newly
replayed here. For11308, the fast cyclotomic/stabilizer identity replays; its
stronger sample bound is not promoted to a theorem. For11310, fast exact group
counts replay18 versus36 order-nine elements; the GAP scripts are also run
independently. For11311, 6000 sampled F9-perfect three-magic-leg words all
violate; that is **sample evidence**, not an exhaustive universal law. Its
non-F9-perfect one-leg certificate reports430/6000=7.1667%, whereas the prose
says6.8%. None of these parallel files is changed or bundled with this packet.
The additional three-qutrit Weyl-flag producer was read: its remaining
24,879,312 layer-scrambling maps are explicitly still open, so its witness
must not be called a certified global optimum. The older qubit-cycle-index
handover was read and kept as prior work rather than rederived.

## 11307 — Integer membrane flux retains a reduced radial negative direction

Add a compact Maxwell membrane form distinct from both sequestering forms.
The reduced Euclidean O4 two-cap action includes GHY terms, wall tension,
and fixed global volume/curvature/Maxwell fluxes. Off shell, the integrated
junction curvature is `6(cplus+cminus) Area/R`; substituting its on-shell
Israel value before varying would give the wrong constrained Hessian.

With q=.4, T=.4, quantize `q integral Fmem = 2pi N`. Sector N=24 has
`Qmem=376.991118430775`, replacing the previous unquantized value377.1856663.
Solving all four stationarity equations gives

```
R=1.7468742236732588, Lambda=.10084255247763585,
kappa²=1.0002578270145999, fminus=.799459922342953.
```

The bare radial second derivative is -134.2591735076. Eliminating Lambda,
kappa² and fminus by their fixed-flux equations gives the Schur curvature
**-204.7110165626**. An independent re-solve at neighboring radii confirms it.
This proves one negative collective-radius direction in the declared reduced
model, not the full fluctuation index, conformal-factor treatment or rate.
A common bare-energy shift is absorbed by Lambda without changing the gradient.

A degree-one collapse map from the genus81 history4-manifold to S4 transfers
its integral H4 flux class. The native three-tick top cycle has unit-chain
(1,1,1) and admits the same period with an exact cochain fluctuation. This is
only a topological flux transfer: no nonsingular round-cap metric or instanton
has been pulled back. Compactness, sectors, charge, tension and geometry remain
inputs; the observed cosmological constant is not predicted.

Primary context:
[sequestered vacuum decay](https://arxiv.org/abs/1604.04000),
[negative-mode subtleties](https://arxiv.org/abs/1210.4740).

## Validation and ownership

Prior owners11287,11290,11293–11302 remain cited in the producers/certificates.
Result-specific corpus/index searches preceded new claims. The nine regressions
check alternative bases, actual tensor contractions, independent spectral
normalizations, exact polynomial identities and constrained numerical response.
Full intake, refreshed index and rebuilt paper are recorded in session notes.

Publication validation: five scoped certificates PASS; nine owned regressions
passed; producers/tests compile; refreshed12-file intake clean with no
rediscovery collisions or forced arithmetic; site card unique; rebuilt paper
PDF835.46KiB. Independent GAP runs confirm SmallGroup(81,9) for the magic group
and SmallGroup(81,7) for both chamber signs. The initial parallel regression failure is recorded historically; the
publication-time integration section distinguishes its later resolution.

## Publication-time parallel integration correction

Remote8f8fb86bd supplies the previously missing11309 certificate field;
the earlier7/8 regression result refers to the prepublication snapshot.
All local cloud copies were backed up under
`/tmp/w33_cloud_11308_11312_before_merge` before integration; normalized
text comparisons show no source differences, and the remote certificate
includes the completed affine-subspace controls. The additional0096da599
formula-universe inventory update was reviewed. The shared ledger merges
both packets, and the binary PDF is rebuilt from the combined sources.

A sentence-level audit finds a more important11309 coverage limitation:
the producer decides all81 frames only for the6480 classes flagged by
0,e1,...,e4. For the other45360 classes it tests those five frames only.
The remaining frames are certified by that shortcut **only if the global
affine-complement law is proved independently**. The300-class sample and
full flagged-class check do not prove it for unflagged classes. Thus the
reported223/2430 is a conditional census result, not an independently
exhaustive4,199,040-Clifford certificate. The merged paper now states this
condition explicitly. Its11311 wording is also corrected to sampled6000/6000
and the actual430/6000 non-F9-perfect one-leg count. Cloud producers and
reports remain preserved; this audit does not replace their numerical data.

After integration, all eight cloud regressions pass in89.21s. This resolves
the certificate-field mismatch, not the separate unflagged-class coverage gap.

The merged18-file cloud intake is clean (no rediscovery collisions or forced
arithmetic); the corrected combined PDF builds at839.72KiB. Intake checks and
passing regressions do not discharge the affine-coverage proof obligation.

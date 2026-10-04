# Passes11438–11442: finite-field obstruction, curved stationarity and deterministic erasure recovery

Reservation: `afa5402d3`; original `2fbe9aa28` released after the later parallel `2c25c12cd` namespace collision. Producer:
[w33_pass11438_11442_coercive_curved_recovery.py](w33_pass11438_11442_coercive_curved_recovery.py).
Certificate:
[w33_pass11438_11442_coercive_curved_recovery.json](../data/w33_pass11438_11442_coercive_curved_recovery.json).
Independent regressions:
[test_w33_pass11438_11442_coercive_curved_recovery.py](../tests/test_w33_pass11438_11442_coercive_curved_recovery.py).

These execute all five requested directions as named models and tests. They do
not complete the five ultimate physical constructions. The main gains are an
exact obstruction to the old finite-stiffness extension, a change of gravity
building blocks that satisfies every interior equation, and explicit quantum
recovery channels. Changed actions, supplied renormalization conditions and
uncompiled operations remain visible.

## Intake and ownership

Integrated parallel `3daf05782` and inventory `57af3f9d8` before reservation.
Read all five new11418–11422 reports and the affected paper material: J8 remains
a conjectured complete witness; the pseudo-reflection J6 blind spot has a
twelfth-order scaling mechanism; the magic-axis law remains partially proved;
four-qutrit sampling is not a census; the apparent constant0.72 depth rate was
withdrawn. None is used to predict a Standard Model parameter here.

Searches included `RESULTS_INDEX.md`, `docs/index.html`, `w33_paper.tex`, recent
reports, scripts and certificates, then result phrases for the normalization
runaway, curved subdivision and complete erasure map. These are bounded searches,
not a claim to have read the entire corpus. Prior owners are11384 native CP,
11428 stiff alignment,11429 local measure,11430 normal-stationary blocking,
11431 threshold forces,11423 ideal erasure ML and11432 flagged extraction.
Generic curved Regge, Abelian chiral measures, renormalization and leakage
reduction are prior literature. This packet makes their maps and limits explicit
for the existing supplied models.

##11438: why the stiff alignment cannot simply be made finite

The11428 map uses
`C=Phi^T Psi*/g0`, `g0=tr(Phi^T Phi*)/3`,
`U=2I+(C+Cdagger)/2`, `B=I+CCdagger`, and
`q=Im tr[U,B]^3`. On the native positive-flux witness, shrink only the first
field: `Phi=t Phi0`, `Psi=Psi0`, with positive real `t`.

Exactly, `C=C0/t`, `q=q0/t^9` and `Im I6=t^6 chi0`. Thus
`V_align=-kappa chi0 q0/t^3`. Both native scalar potentials stay finite as
`t->0`: moments remain zero and the invariant relations scale homogeneously.
Every finite positive stiffness therefore leaves the complete old action
**unbounded below**. A global finite-stiffness vacuum and its Hessian do not
exist for that extension. This does not retract the compact fixed-orbit result.

A separately declared replacement uses
`g=(||Phi||F²+||Psi||F²)/2+sigma²`, the same pin formulas, the bounded tilt
`-kappa chi q/(1+chi²)`, and
`rho (||Phi||F²+||Psi||F²)²`, with supplied positive `sigma,rho,kappa`.
AM–GM gives `||C||<=1`; hence `||[U,B]||<=2`, `|q|<=24`, and the tilt has
absolute value at most `12 kappa`. The native square-completed potentials give

`V >= -5/8-12 kappa+rho (||Phi||F²+||Psi||F²)²`.

This continuous gauge-invariant CP-even action is coercive on all324 real
field coordinates and has a finite global minimizer. The below-bound trial witness described next also forces every global minimizer
to have nonzero native CP order and pin CP flux.
The stored native seed gives nonzero bounded flux `3.635115310781e-6`.
The common compact gauge orbit at the pair has rank86. Two separate native
orbits have total dimension156, leaving70 relative directions. A single-field
84-normal Hessian cannot certify stability of this two-field system.
An extra all324-field optimization supplies a stationary candidate and its
complete324x324 Hessian. The reproducible analytic-gradient helper is
[w33_pass11438_finite_native_model.py](w33_pass11438_finite_native_model.py).
Coefficients are explicitly supplied: `kappa=10000,rho=1e-6,sigma=.1`.
Every state with a vanishing alignment term has energy at least `-5/8`,
because the native potentials obey their separate bounds and the added norm
term is nonnegative. The trial energy below this bound forces nonzero native
CP order and pin CP flux at every global minimizer of this supplied action.
This uses the trial energy; it does not identify the candidate as the minimizer
or settle the soft-mode stability. CP here means conjugation composed with the
named compact gauge group, not a classification of all generalizedCP symmetries.
Action value is `-2.935743239533`, gradient norm `1.32e-6`, and common-gauge
Hessian residual `6.04e-10`. The first field has
`I6=-.1162882825+.7047959636i`, with cross-pin CP flux `.0007214784773`.
The complete gauge-normal Hessian has234 positive modes (smallest clearly
positive `.40612573`) and four near-zero modes from about `-4.13e-7` to
`3.74e-7`. Smaller-step directional controls are about `1e-7` to `2.3e-7`.
Thus the four modes remain unresolved; the candidate is not called an isolated
stable vacuum or the global minimizer. Positive unrelaxed line responses in
those directions scale approximately quartically, but relaxation of the massive
fields and coupled quartic terms could change that conclusion. Those analyses,
stronger stationarity certification and physical coefficient selection remain
open.
The full common compact orbit also leaves no continuous stabilizer. This candidate is not an unbroken Standard Model gauge vacuum; an embedding or a modified symmetry-breaking pattern is still needed.
Its rational saturation and high-degree native potential are a supplied EFT,
not a fundamental renormalizable theory.

##11439: the admissible Abelian branch has a precise size obstruction

One-family integer charges are `[1,-4,2,-3,6,0]`, with multiplicities
`[6,3,3,2,1,1]`. Linear and cubic sums vanish. Absolute odd-charge multiplicities
are6 at1 and2 at3, both even. Three families preserve these conditions.
The imported all-sector Abelian construction also needs admissibility,
locality and a sufficiently large lattice; the charge test alone supplies no
implemented measure current. See [Luscher, hep-lat/9811032,
sections2.3 and5](https://arxiv.org/pdf/hep-lat/9811032).

At the stored background scaled by `1e-4`, the largest charged plaquette phase
is `0.004169882804`, below the sufficient bound `1/30`. Charged Wilson gaps and
projector ranks are independently replayed. Under this common bound an L=3
lattice cannot support nonzero magnetic flux. A unit flux requires
`6*2pi/L² < 1/30`, so the smallest possible integer L is34.

Explicit periodic links are
`U0(x,y)=exp(-2pi i y/L²)` and
`U1(x,y)=exp(2pi i x/L)` on the `y=L-1` seam, unity elsewhere.
Every plaquette has phase `2pi/L²`; the plane flux is exactly one unit.
All six charged representations pass, with maximum phase `0.032611688446`.
L=33 fails the bound. This constructs a nontrivial admissible background, not
its four-dimensional fermion measure. L=34 is a flux size bound, not proof of
the separate localization-length hypothesis. Non-Abelian reconstruction,
sector weights, an implemented current and a native continuum limit remain open.

##11440: all centroid equations close for a changed curved action

The flat-simplex15-radial model still leaves full gradient norm about
`9.29e-6` after normal elimination. Direct root attempts did not reduce it;
a subsequent center solve left the positive-simplex domain. Neither is a
stationary gravity result.

Use geodesic simplices of supplied sectional curvature `k=.001` instead:

`S_k = sum_internal A_k delta_k + sum_boundary A_k psi_k + 3k sum V4_k`.

The volume coefficient is conventionally `Lambda=3k` for comparison with11430.
The map is explicit: `Gij=cos(sqrt(k)*ell_ij)`, `X=Cholesky(G)`,
`c=normalize(sum wi Xi)`, and `ri=acos(Xi dot c)/sqrt(k)`.
Positive weights select an interior point. Each coarse simplex subdivides1->5,
giving the same15 fine simplices and15 radial variables as before.

Every interior hinge lies in the original constant-curvature simplex and has
zero deficit. The curved Schlaefli identity `sum A dtheta=3k dV` cancels the
volume derivative, leaving `dS_k=sum delta dA`. Hinges involving an interior
vertex have zero deficit; boundary-only areas have no radial derivatives.
**All15 radial equations therefore vanish**, including the twelve centroid
equations. Two distinct positive-weight placements pass independently; tests
use spherical cosine-law triangle areas in place of the producer's inverse
Gram formulation. Orthogonal changes of the spherical coordinates preserve
the constructed lengths.

The complete15-radial Hessian has three large positive normal eigenvalues near
`5415.365,5505.059,7702.802` at the equal-weight placement. Twelve displacement
null directions follow from the exact subdivision identity. Two difference
steps show the small spurious eigenvalues shrinking quadratically; the largest
at step `1e-5` is about `3.8e-4`, compared with normal modes above5000.
This is numerical normal stability with analytically identified null directions,
not interval-certified eigenvalues.

This changes the action and building blocks. It does not solve the old flat
equations. The stationary fine action equals the coarse curved action because
interior deficits vanish, boundary angles add and geometric volumes partition.
This local subdivision statement supplies no globally perfect4D gravity action,
Lorentzian evolution, selected curvature or cosmological constant prediction.
No numerical4-volume integration is claimed. The mechanism and identities
are prior [Bahr–Dittrich0907.4323, AppendixA](https://arxiv.org/pdf/0907.4323).

##11441: a running stationary slice, with its inputs exposed

Use the actual two486-state quark spectra and480 link radial/480 phase modes,
with `phi>.8` so scalar logarithms stay on the non-tachyonic slice. Choose three
renormalization conditions explicitly: `phi*=.81`, radial curvature100 and
vacuum energy0. Add `c0+c2 phi²+c4 phi⁴` to the full slice threshold plus tree
potential. The supplied conditions yield coefficients approximately
`[571430.9650,187077.0375,3679.2955]`.

Their one-loop running is fixed by the slice supertrace:
`dc/dln(mu)=coefficients of STr M4(phi)/(32pi²)`.
Fitting the even quartic and checking off-grid points gives beta coefficients
approximately `[-367816.1201,-87558.9557,-1020.7096]`.
At scales `.5,1,2`, the total potential is invariant within `1e-7` on the
tested nearby points. Central-difference force is `0.00047672` against the old
force of order300000, and curvature is `100.00658`.

This demonstrates a self-consistent locally stable radial slice under imposed
renormalization conditions and fixed spectral couplings. It derives neither
the chosen modulus nor its curvature or energy. A different constant counterterm
changes the vacuum energy without changing stationarity. Full joint pin/link/
Majorana equations, gauge-completed counterterms away from this slice, fundamental
coupling running, perturbative control and the995-field Hessian remain open.
The standard threshold construction is prior
[Martin, hep-ph/0111209](https://arxiv.org/abs/hep-ph/0111209).

##11442: deterministic recovery is a channel, not block rejection

For known erased physical sites, replace each by `|0>` with Kraus maps
`K_b=|0><b|` on that site and identity elsewhere. In the explicit Steane
encoding V, every single- or double-erasure support satisfies
`Vdag K_bdag K_c V=delta_bc I/r`, with `r=2^|E|`.

Set `R_b=sqrt(r) V(K_b V)dag`. Their adjoint-product sum is the projector onto
the noisy code subspace. Add `|0L><v_j|` for an orthonormal basis of its
complement to make the recovery trace-preserving on the full128-space.
Composition with the nonselective reset channel is exactly identity on the
logical qubit. All28 supports pass entanglement-fidelity and trace-preservation
checks. Seven three-site supports carry a minimum-weight logical operator,
showing that arbitrary three erasures cannot be included in that guarantee.
Prior11423 already counted ideal erasure ML; this supplies an explicit complete
recovery map for that model. Leakage reduction as a general technique is prior
[Aliferis–Terhal, quant-ph/0511065](https://arxiv.org/abs/quant-ph/0511065).

Routes use80 independent copies of the prior12-dimensional native cell, indexed
by A80, with graph-neighbor joint couplings supplied as architecture inputs.
Data copies0–6, syndrome7 and flag8 need no degree-eight hub: shortest paths have
remote CNOT costs `[19,7,19,19,7,19,19,19]`, including forward and reverse SWAPs
and three CNOTs per SWAP. Each compiled CNOT costs107635 old-model ticks.
Exhaustive path-basis tests restore every spectator and perform the remote CNOT.

This is a routing recipe and exact deterministic recovery channel. It does not
construct80 independent cells inside one cover, synthesize the recovery channel
with native gates, or implement native location detection/reset/preparation/
readout. Routing introduces many noisy gates and changes the fault model;
the old1714-case single-fault certificate cannot be reused as a routed or leakage
threshold proof. During-pulse noise and leakage propagation remain open.

## Validation and scope

Five producer sections PASS. Thirteen independent regression functions PASS.
Source certificates are hashed canonically. Full batch intake, rediscovery,
forced-arithmetic and syntax checks accompany publication. Numerical tolerances
above describe finite models, not observed physics or a complete TOE.

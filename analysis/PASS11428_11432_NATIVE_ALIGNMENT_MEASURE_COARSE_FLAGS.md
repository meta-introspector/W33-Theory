# Passes11428–11432: native CP reaches a pin, local chiral measures, curved coarse responses and flagged native correction

Reservation `8e8bfdcb3`. Producer:
`analysis/w33_pass11428_11432_native_alignment_measure_coarse_flags.py`.
Certificate: `data/w33_pass11428_11432_native_alignment_measure_coarse_flags.json`.
Regressions: `tests/test_w33_pass11428_11432_native_alignment_measure_coarse_flags.py`.
This executes the five directions following11423–11427 as concrete calculations.
It does not complete their outstanding physical problems or solve the TOE.

The main constructive connection is a map from **two actual native81-component
fields** to full-rank family pins, transferring11384's native invariant CP
order into an endpoint-matched cubic CP invariant. The extra field and the
stiff-orbit approximation are explicit. Independently, an actual native CNOT
is attached to a flagged circuit that satisfies the stated one-fault contract.
The measure, gravity and vacuum calculations expose the remaining obligations
rather than replacing them by named-but-unbuilt objects.

## Intake and prior ownership

GitKraken fetched and integrated `fdccffd97`, an inventory-only update, before
reservation. The remote had no later scientific packet at the final review.
Result-oriented searches checked RESULTS_INDEX, existing paper/index material,
recent reports, Python producers and certificates. The following owners matter:

* [11325–11329](PASS11325_11329_COMPLEX_PHASE_MIXING_QUOTIENT_AND_WALL.md),
  especially11326 and `w33_pass11326_noncommuting_flavor_operator.py`, already
  constructs an engineered CP-even noncommuting flavour vacuum and its complete
  normal Hessian. Generic noncommuting flavour selection is **not new here**.
* [11384–11388](PASS11384_11388_NATIVE_CP_HARD_MATCHING_GRAVITY_AND_RELAXED_WALL.md)
  owns native invariant CP selection on162 real fields, its actual Cartan seed,
  moments and signed E6 generators. It supplies our two stiff vacuum orbits.
* `w33_20261001_global_e6_cartan_covariants.py` and
  `w33_20261001_degree18_phase_completion.py` own the global native invariants;
  their actual values on both new fields are replayed in tests.
* [11403–11407](PASS11403_11407_NATIVE_BIMODULE_COLLECTIVE_GATES.md) owns the
  actual colour/family slots, anomaly family and conditional modular gates.
* [11408–11412](PASS11408_11412_CLOCK_SPURION_SEAM_TICKS.md) and
  [11413–11417](PASS11413_11417_PIN_REGGE_VACUUM_NOISE.md) own native endpoint
  matching, FCC periods, finite-tick calibration and the shortened Hadamard.
* [11423–11427](PASS11423_11427_FLUX_INDEX_CURVED_UV_NOISE.md) owns the
 19-coordinate pin potential, actual flux-index background, curved mesh with
  flat inserted stars, full declared heavy spectrum, during-Hadamard noise and
  the bare-ancilla hook failure. The present packet uses its stored witnesses.
* [BT4081–BT4088](BT4081_BT4088_deep_physics.md), especially4084, already uses
  external4D overlap fermions. Neither their construction nor a supplied
  spacetime is relabelled a native physical-chirality derivation.
* `exploration/PART_CCCCIV_W33_CSS_STEANE_LIFT.py` owns the prior Steane lift.

D-flat scalar Grams, Schur's lemma, CP invariants, overlap measure geometry,
stationary Schur complements, fermion determinants and flag correction are
established methods. Bounded corpus searches found the owners above; they
are not worldwide novelty proofs. Full source hashes are canonical JSON.

##11428 — the isotropic route fails; a second native field transmits CP

Write the actual native field as a27-by3 matrix Phi. Family D-flatness forces
`Phi^T Phi* = g0 I`. This consequence of the existing moment equations is
prior mathematics. The exact integer Cartan embedding gives more structure:
all coefficients of this Gram are scalar, and the native symmetric cubic

`T_ijk = d_abc Phi_ai Phi_bj Phi_ck`

has cyclic X invariance and only Z-neutral coefficients. These identities are
checked on its full integer coefficient tensor, not only at one numerical
vacuum. Therefore Hermitian family covariants in the named algebra generated
by T and invariant contractions commute with the irreducible family clocks
X,Z; their commutant is scalar. In particular the first two cubic contractions
K and L are scalar. This is a restriction on the named algebra, **not a
classification of every possible native covariant**.

Introduce another native field Psi on the same11384 vacuum orbit. The fields
transform under the **same** compact E6 and family SU3. Their relative E6
orientation is physical once common transformations are quotiented. Define

```
C = Phi^T Psi* / g0,       g0 = tr(Phi^T Phi*)/3
U = 2I + (C+Cdagger)/2
B = I + C Cdagger
q = Im tr[U,B]^3
chi = Im I6(Phi) = sqrt(7)/4 on the stored branch
V_align = -kappa chi q,    kappa > 0.
```

C,U,B transform by common family conjugation and are E6 singlets. The action
is CP even: chi and q both reverse under conjugating the native fields.
It transfers **existing native CP** into a family mass invariant rather than
claiming a new generic mechanism for spontaneous CP breaking.

Both fields have scalar Gram g0 I, so C is a contraction, U has eigenvalues
at least1, and B has eigenvalues at least1. Both pins are full rank without
prescribed individual entries or eigenvalue targets. A Lie-orbit ascent using
all78 compact E6 directions gives q=0.07509685822590392. The actual signed
invariants I6,I12,I18 and all86 moments are replayed on both fields. A common
E6/family basis transformation checks covariance independently.

The relative orbit is compact. Since a positive q witness exists, minimizing
this named alignment action globally on that orbit has nonzero q. The stored
ascent point is a numerical witness, not a certified unique global optimizer.
The native seed and its full-field stability retain11384's numerical boundary.
Finite stiffness and a324-real-field Hessian for the enlarged action are not
computed. The added field, kappa, singlet offsets and a unit pin scale are
supplied. This is a stiff-orbit EFT construction, not a renormalizable UV model.

The pins feed the existing actual endpoint propagator:
`Yu=D_R U D_R`, `Yd=D_R B D_R` at eta=.08. Its squared-Gram CP cubic is about
8.96974e-23; smallness follows the supplied path hierarchy. The actual
rank-nine colour carrier commutes with all eight native colour generators.
That makes the slot map concrete; it does not identify observed particles or
predict CKM parameters. The endpoint map includes a supplied family spurion
D_R; a common family basis change must conjugate D_R along with the pins.
Deriving its alignment from covariant microscopic attachments remains open.
Extending the old fixed Majorana source to this fully
family-covariant model is a separate obligation, explained in11431.

##11429 — an actual local Weyl measure frame and its nonzero curvature

On a supplied L=3 periodic4D gauge lattice with antiperiodic temporal spin
boundary conditions, the Wilson sign defines a rank162 Weyl projector in324
spinor dimensions. Use integer hypercharges `1,-4,2,-3,6,0` with multiplicities
`6,3,3,2,1,1`, one Standard Model-plus-neutrino family template. Their first and
third charge moments vanish exactly; this anomaly family is prior11403.
Repeating three families multiplies the measure curvature by three; the table
below records one family, not a derivation of physical generations.

The producer computes

`F12 = i Tr P[dP1,dP2]`

using the divided-difference derivative of the actual Wilson sign. Charged
projectors pass finite-difference and site-dependent gauge-realization tests.
For the stored weak-link background and tangent directions:

| Background epsilon | One-family summed measure curvature |
| --- | ---: |
| .01 | -1.06431695e-3 |
| .005 | -1.33839046e-4 |
| .0025 | -1.67544630e-5 |

Linear response cancels, while the finite residual scales approximately as
 epsilon cubed. **Anomaly cancellation does not require a flat determinant
bundle connection.** This nonzero curvature is not itself a residual gauge
anomaly. A measure connection must have the appropriate curvature and gauge
properties rather than simply declaring it zero.

A smooth local section is actually constructed on a contractible patch:

`F(A)=P(A)E [Edagger P(A)E]^-1/2`.

The square-loop determinant holonomy is approximately-2.27298e-10; it agrees
with the independently computed local curvature at the loop's midpoint, with
orientation convention accounted for. The loop area is1e-6. This frame also
exists for anomalous charges on a sufficiently small patch, so its existence
is not a proof of anomaly-free global gauge invariance.

The required global local gauge-covariant current, topological-sector gluing,
non-Abelian measure and physical mirror selection remain unbuilt. The
[primary Abelian construction by Luescher](https://arxiv.org/abs/hep-lat/9811032)
shows what an exact measure theorem must establish; it is not supplied by our
local polar frame or by the prior background index.

##11430 — curved fine stars and a constructed stationary-normal coarse action

Use11425's actual native-period boundary lengths and its15-simplex refinement.
Add a supplied cosmological term:

`S = S_Regge - Lambda sum4-simplex volumes`.

Analytic Schlaefli/area derivatives and Gram-determinant volume derivatives
produce the full radial gradient. At fixed centroid coordinates, solve the
three radial-normal equations. For Lambda=.001 this gives local fine-star
curvature about0.0019747, whereas the previous inserted stars were flat.

The constructed blocked action depends on a boundary length and on the twelve
centroid coordinates; three radial-normal variables have been eliminated by
stationarity. The remaining full radial gradient norm is about9.29e-6, although
all three normal equations pass. The centroid equations are **not removed by
calling them gauge**. Center probes also give nonzero responses.

For the tested boundary length, the stationary Schur Hessian is-25.2255479;
an independent re-solved boundary finite difference gives-25.2256893. The
negative value is a boundary response of this restricted Euclidean model,
not a Lorentzian particle mass, ghost count or physical instability theorem.
The elimination is a saddle calculation on the named slice, not minimization
over every fine degree of freedom.

This implements a specific coarse response and exposes exactly where full
curved constraints are still missing. It is not a perfect action or a full
stationary gravitational solution. The established coarse-graining programme
is described by [Bahr–Dittrich](https://arxiv.org/abs/0907.4323); its perfect-action
results are prior art, not consequences silently assigned to our4D mesh.

##11431 — small relative phases are lifted; the old link modulus is driven away

Keep11423's selected pin and its declared real Majorana source. The quark
singular spectra are unchanged under diagonal family rephasing of the pin.
The six-state neutrino determinant generally changes because the source is
held fixed. Computing tiny differences requires high precision throughout,
including subtraction and Hessian operations, not just the eigensolver.

The phase minimum found on this fixed-source torus is approximately
`(alpha,beta)=(3.06978030,3.08966561)` modulo pi, with potential shift
-5.92275e-13 relative to the original phase representative. Its two positive
canonical phase mass-squared values are approximately1.00989e-9 and1.37143e-9
in supplied model units. The normalization uses the actual induced kinetic
Gram of the pin's phase tangents. This is a numerical local minimum, not a
global interval certificate or full scalar-vacuum solution.

**Gauge control:** simultaneously transform the Majorana source as
`S -> G* S Gdagger`. The full symmetric neutrino mass matrix then transforms
by unitary congruence, keeping its Takagi spectrum and threshold unchanged.
Thus the fixed-source relative-phase lift is not a lift of gauged family
Goldstones in11428. Combining the two models requires a covariant dynamical
source and its own vacuum, not importing this phase potential unchanged.

For link phases, independently anchor each of the three node-phase profiles
at the pin. The incidence rank is79 per network, giving237 node-gradient
phase directions with unchanged full quark singular spectra. There remain
243 cycle-phase directions. These are spectral flat directions of the stated
one-loop determinant, **not** protection against nonperturbative gauge
anomalies or a complete quantum symmetry classification.

The common link-modulus slice includes both complete486-state quark spectra,
all480 link radial and480 link phase scalar eigenvalues, and the tree link
potential. Other declared thresholds are constant on this slice. At phi=.8,
its one-sided quantum force is approximately-306859.4272. The tree force is
zero there. Hence the old tree background is **not radiatively stationary**
on even this restricted joint slice. The calculation remains in the stated
one-loop scheme; it is not a controlled large-field extrapolation.

Full995-field stability, cycle-phase lifting, running couplings, finite
counterterm conditions and a protected cosmological constant remain open.
This is a connection between the actual declared heavy spectrum and vacuum
selection, not a physical vacuum-energy prediction.

##11432 — native CNOTs and a one-fault flagged extraction contract

The actual native Hadamard, a calibrated pi phase and an intercell density
phase compile a CNOT:

`CNOT = -(I tensor H)(Z tensor Z) CZ00 (I tensor H)`

up to global phases. The exact joint carrier has dimension144. The compiled
gate takes107635 ticks, phase-aligned operator error0.0131377, and average
coherent leakage8.21018e-5. A specified native edge-phase channel is applied
**at the compiled gate boundary**; it is not called noise at every internal
tick. Logical Pauli averaging is explicit, with errors placed after the
ideal CNOT to match the syndrome fault convention.

Replace the bare extraction ancilla by a syndrome ancilla and a flag ancilla.
Three flagged rounds use108 CNOTs. If any flag occurs, a complete additional
unflagged syndrome round gives132 total. The earlier fixed-round candidate
had six ambiguous histories; the conditional branch resolves them.
Under the single-fault contract, a flag already consumes the one allowed
fault, so the conditional follow-up is clean. Multi-fault simulation includes
its actual noisy gates and preparation/readout.

The certificate stores a decoder for476 measurement histories and checks1714
input/single-fault cases, including:

* all21 incoming single data Paulis with no circuit fault: output restored
  modulo stabilizers;
* all1620 single CNOT Pauli faults and72 representative preparation/readout
  flips on clean input: residual equivalent to at most one data Pauli.

This second condition is the relevant correctability condition; requiring
zero residual syndrome would incorrectly reject harmless late single errors.
Actual Steane code states independently check recovery, and a nine-qubit
statevector replays a flagged hook including ancilla measurements.
Flag correction is established prior art:
[Chao–Reichardt](https://arxiv.org/abs/1705.02329). This is a concrete compiled
model and certificate, not a new generic flag-correction theorem.

For16384 accepted-Pauli simulations per boundary-noise rate, with independent
preparation and readout flips1e-5:

| Boundary edge flip | Uncorrectable residuals | Estimated rate | Wilson95 interval | Conservative no-leak acceptance |
| --- | ---: | ---: | --- | ---: |
| 0 | 0 | 0 | [0,2.34418e-4] | .989221 |
| 1e-5 | 0 | 0 | [0,2.34418e-4] | .986613 |
| 1e-4 | 2 | 1.22070e-4 | [3.34760e-5,4.45025e-4] | .963450 |

Zero observed failures is not a zero-rate theorem. Unknown multi-fault
histories use the declared majority/minimum-weight fallback; it is not an
optimized multi-fault decoder. The leakage acceptance bound uses the maximum
132 compiled gates and is distinct from conditional Pauli correctness.
Any detected native leakage is **rejected**, not deterministically corrected.
During-gate CNOT noise, correlated drift, leakage-aware recovery, physical
reset/readout and a scalable threshold remain open.

## Validation and common-action boundary

All five producer fronts PASS; nineteen focused independent regressions cover
exact Cartan coefficients, global native invariants, second gauge/basis
realizations, physical colour slots, finite projector derivatives, local
holonomy, gravity gradients and Schur elimination, relative/gauge phase
controls, link-gradient rank, actual code states, a nine-qubit circuit and
the native after-gate channel. Full four-file corpus intake is clean: no
rediscovery collisions, no forced arithmetic and no certified-value
contradictions. Arithmetic selftest and compilation pass; RESULTS_INDEX was
regenerated.

These are explicitly distinct models. In particular, the fixed-source loop
phase calculation is not transplanted into the fully family-covariant native
alignment theory, and the supplied Regge Lambda is not identified with a
measured vacuum energy. The new maps make those compatibility questions
computable; a common physically selected microscopic action remains open.

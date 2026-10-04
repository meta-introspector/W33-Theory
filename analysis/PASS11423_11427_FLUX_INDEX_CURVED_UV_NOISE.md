# Passes11423–11427: spontaneous pinning, a flux index, curved gauge islands, heavy thresholds and noisy extraction

Reservation `69e0bfcd1`. Producer:
`analysis/w33_pass11423_11427_flux_index_curved_uv_noise.py`.
Certificate: `data/w33_pass11423_11427_flux_index_curved_uv_noise.json`.
Tests: `tests/test_w33_pass11423_11427_flux_index_curved_uv_noise.py`.
All five producer sections and seventeen independent regressions pass.
These execute the five follow-ups to11413–11417 as explicitly supplied models;
none establishes a physical TOE or measured masses, angles or vacuum energy.

## Intake and ownership

The initial Codex reservation `6c6bb42fe` collided with Claude's earlier
`dae280475` reservation for11418–11422. The remote rejection was respected:
`1ec47c9ab` integrated the formula inventory `ea3b54dc5`, and `69e0bfcd1`
released the conflicting claim and reserved11423–11427 before computation.
Claude's11418–11422 packet is not modified. GitKraken fetch and remote review
were repeated before publication; source hashes use sorted compact JSON.

Result-oriented searches checked the result index, paper/index regions and
existing analysis, scripts and certificates. Relevant earlier owners are:

* [11374–11378](PASS11374_11378_FIVE_TOE_FRONTS.md): CP-even phase selection
  on three rays, not a full-field scalar completion.
* [11384–11388](PASS11384_11388_NATIVE_CP_HARD_MATCHING_GRAVITY_AND_RELAXED_WALL.md):
  actual native E6 invariant CP selection on162 real coordinates,84 positive
  normal directions and78 gauge directions; also hard matching, hub gravity
  and negative wall modes. Its certificate was reread. Full-field spontaneous
  CP selection already exists. The present19-coordinate pin construction
  connects another supplied potential directly to the11413 family mass map;
  it does not supersede or rediscover that162-field construction.
* [11408–11412](PASS11408_11412_CLOCK_SPURION_SEAM_TICKS.md): native Yukawa
  arrow tensors, path protection assumptions and gate calibration.
* [11413–11417](PASS11413_11417_PIN_REGGE_VACUUM_NOISE.md): the full-family
  common-pin map, native FCC periods, declared low-energy determinant,
  shortened29091-tick Hadamard and actual Steane code states.
* [BT4081–BT4088](BT4081_BT4088_deep_physics.md), specifically4084:
  an external4D overlap/Ginsparg–Wilson construction with a vectorlike native
  fiber. External spacetime and the overlap mechanism already exist here.
* [11390–11394](PASS11390_11394_CONTEXT_MATTER_CLOCK.md): native80-context
  carriers, family/colour slots and the finite chirality boundary.
* `exploration/PART_CCCCIV_W33_CSS_STEANE_LIFT.py`: prior Steane lift.

Generic spontaneous CP, overlap fermions, Regge subdivisions, vectorlike
messenger matching, one-loop determinants and the Steane code are established
methods. This packet supplies concrete compatible maps and failure controls,
not external novelty claims for those methods. Searches are bounded corpus
checks, not a proof of worldwide novelty.

##11423 — real couplings select a full-rank CP-breaking family pin

Let U,B be Hermitian3-by3 matrices and chi a real pseudoscalar. In canonical
Hermitian coordinates there are19 real fields. The actual potential is

```
Vup = ||offdiag U||² + ||U³-6U²+11U-6I||²
      + (tr U-6)² + (tr U²-14)²
Vdown = ||diag B-a||² + sum_cycle(|Bij|²-r²)²
        + (chi²-vchi²)² - g chi Im(B01 B12 B20)
a=3, r=1/4, vchi=1/2, g=1/10.
```

It is invariant under CP `(U,B,chi)->(U*,B*,-chi)` and under conjugation by
family cyclic shift X and diagonal qutrit clock Z. These clocks act on the
family multiplicity space of the stored colour carrier; they are not a claim
that the original H27 colour action has become a family action. The clock
axis and real spectral/radial coefficients are supplied. In particular the
polynomial chooses eigenvalues1,2,3; these are not predicted observed masses.
The up potential has degree six and is a cutoff scalar EFT.

Up minima are diagonal permutations of1,2,3. The trace-square term excludes
repeated root configurations. For the down sector, triangle phases align
with chi. Young's inequality bounds the mixed quartic by g/4 times the sum
of four fourth powers, proving coercivity for this g<4. Any chi=0 or any
zero-radius boundary has nonnegative energy; the nonzero branch has negative
energy. At an interior stationary point all three radii exceed r, and
`4 ui²(ui²-r²)=g |chi| u1u2u3` makes them equal, since the left side is
strictly increasing there. Writing w=u² gives the certificate's cubic:

```
(256-g⁴)w³ -(768r²+16vchi²g²)w²
 +(768r⁴+16vchi²g²r²)w -256r⁶ = 0.
```

An exact rational Sturm count gives one root above r². The positive branch
has u=0.2563388764360664 and chi=0.5008400781229808. A representative is
`U=diag(1,2,3)`, `B=3I+zX+z*X^T`, `z=u exp(-i pi/6)`.
The full19-coordinate numerical Hessian has17 positive directions, smallest
positive eigenvalue about0.0192577, and two diagonal-rephasing flat directions.
The pair invariant `Im tr[U,B]^3=0.20212716062708724` flips sign under CP.
It also excludes a common unitary conjugation identifying the two pairs.

Embedding B in the actual family/colour frame commutes with native colour.
The previously constructed endpoint map `Y=D_R B D_R` is full rank and has a
nonzero cubic CP invariant; its tiny magnitude follows the supplied path
parameters, not a measured CKM fit. Both scalar minima and the full mass map
are now named. Loop-stable pinning, the origin of the clock anisotropy, the
E6-to-pin map and physical parameter selection remain open.

##11424 — a nonzero background index and an actual rectangular Weyl map

On a supplied periodic4D Euclidean torus with L=4, uniform U1 fluxes are
inserted in the12 and34 planes. The1024-dimensional Wilson kernel at mass1
is explicitly diagonalized. The standard overlap construction is
`D=I+gamma5 sign(H_W)` and `hat_gamma5=gamma5(I-D)`.

| Fluxes | Overlap index | Wilson gap |
| --- | ---: | ---: |
| (1,1) | -1 | 0.6192073257 |
| (-1,1) | +1 | 0.6192073257 |
| (0,0) | 0 | 1 |

The Ginsparg–Wilson residual is below3e-13. In the first background,
`hat_P_minus` has rank513 and ordinary `P_plus` rank512; their actual
intertwining map `D hat_P_minus=P_plus D` has residual below2e-13.
A stored normalized zero mode has negative chirality and residual below3e-15.
Plaquettes, flux reversal, gauge rephasing and the zero mode are replayed.
Tensoring the actual native rank-nine family/colour frame gives index-9.

The overlap relation is prior [Neuberger](https://arxiv.org/abs/hep-lat/9801031),
and external4D overlap is prior4084 in this corpus. This adds an explicit
flux background, zero mode and rectangular map to the native rank-nine fiber.
A background Dirac index is not an anomaly-safe Standard Model Weyl measure,
a mirror-decoupling construction, emergent spacetime or Lorentzian dynamics.
No numerical index is relabelled a count of physical generations.

##11425 — global curvature coexists with locally flat translation directions

Three coarse4-simplices share an internal triangle. Boundary coordinates
use the actual common-index-five FCC spatial periods of11415 and a supplied
Euclidean time scale. Perturbing one boundary length by0.2 percent produces
internal hinge deficit0.005153998065082099. Each coarse simplex is split1-to5
by adding its centroid, giving15 simplices and15 radial lengths.

The three newly inserted vertex stars are individually flat. Moving each
inserted vertex through its coarse-simplex interior preserves all boundary
lengths, coarse curvature and the Regge action including its boundary term.
The twelve tested finite translations have action residual below3e-13.
Three5-by5 radial Hessians have four translation directions each; using the
Schlaefli identity and analytic area derivatives gives tangent residuals
8.9e-7,1.05e-5,3.03e-6. Independent simultaneous finite moves test the nonlinear
statement rather than relying only on those numerical Hessians.

An earlier second difference of the total action suffered cancellation;
replacing it by derivatives of the analytic radial gradient fixed the
numerical calculation without relaxing its criterion. Perturbing one radial
length normally creates local curvature, with maximum deficit about0.392.
This is a restricted locally-flat subdivision mechanism. It does not give
translation symmetry at generic curved moved stars or a full constraint
algebra. It is consistent with11415, which put curvature inside the moved
star, and with the generic broken-symmetry problem described by
[Bahr–Dittrich](https://arxiv.org/abs/0905.1670).

This suggests testing stationary coarse graining with native boundary data:
flat internal subdivisions exactly reproduce their coarse action. A perfect
4D action for arbitrary curvature is not constructed here.

##11426 — the complete declared messenger extension exposes the vacuum cost

The quark bridge is explicit: on each native vertex there are three families
of vectorlike Q,U,D fields with Standard Model representations
`(3,2,1/6)`, `(3,1,2/3)`, `(3,1,-1/3)`. At the pin, three F_u and three F_d
Dirac mediators have the matching right-singlet representations. The chain
is `Q -- Higgs --> F -- neutral family pin --> right network`. Its fermion
vertices have dimension four. All heavy fermions are vectorlike, so their
own gauge anomalies cancel; that says nothing about a future chiral measure.

For M=10, eta=.08 and Higgs v=1 the two483-by483 heavy blocks use
`K=M(I-eta A_Levi) tensor I3`, mediator mass M, Higgs coupling v/sqrt2 and
pin B. Attaching the three external left/right families gives two486-by486
full matrices. The exact Schur map is
`Yeff=-v/(sqrt2 M³) D_R B D_R`. All full singular values are stored;
light singular values are not equated to an unnormalized Schur matrix.

Three link networks on160 edges give480 complex neutral link scalars.
Their radial mass squared is0.256;480 phases are flat at tree level. Together
with19 pin coordinates,12 Majorana-source and four Higgs coordinates there
are995 real scalars, of which485 have zero tree mass. There are2898 heavy
Dirac fields,2919 total Dirac fields including quarks and charged leptons,
and six Majorana fields, for fermion real weight11688.

The neutrino EFT now uses the same selected pin, `Ynu=.01 B_down`, with the
prior real Majorana source. This dimension-five lepton-pinning vertex and
the degree-six up potential make this a declared cutoff extension, **not a
fundamental ultraviolet completion**. The full six-state seesaw is stored.
All fixed coefficients and sources are real; conjugating the chosen pin
conjugates the full fermion matrices. Opposite CP branches therefore have
identical singular spectra and one-loop thresholds. CP degeneracy supplies
no cancellation of vacuum energy.

In the stated MSbar/Landau fixed-mass calculation at mu=1,

```
STr M⁴ = -133981012.15279265
V1     = -692586.7185293995
```

The certificate enumerates every declared massive quadratic group and checks
`dV1/dln mu=-STr M⁴/(32pi²)`. These are supplied model units and background
thresholds, not an observed cosmological constant or a loop-minimized vacuum.
Link phases, scalar loop stability, physical running and the independent
vacuum counterterm remain unresolved. The standard determinant formula is
reviewed by [Martin](https://arxiv.org/abs/hep-ph/0111209).

##11427 — native gate noise, erasure flags and a failing extraction circuit

The actual shortened29091-tick native Hadamard is propagated as a density
channel. After each tick, both addressed160-edge coordinates independently
undergo phase flip with probability q. The joint invariant carrier has
complex dimension12, and exactly contains the logical seeds and noise axes.
It retains coherent control error and leakage. Explicit logical Pauli
averaging and ideal leakage flags define the Pauli/erasure recovery model;
those operations are stated assumptions, not properties silently assigned
to the unprocessed native channel.

Exact maximum-likelihood Steane recovery sums16384 Pauli patterns over all128
known erasure supports, grouping by syndrome and logical class. Actual code
states independently verify dephasing and erasure controls. The code and
correction mechanism are prior [Steane](https://arxiv.org/abs/quant-ph/9601029).

| Native edge flip q/tick | Uncoded infidelity | Ideal encoded failure |
| --- | ---: | ---: |
| 0 | 4.6191e-5 | 2.4328e-11 |
| 1e-8 | 4.8579e-4 | 1.3817e-6 |
| 1e-7 | 4.4323e-3 | 1.2342e-4 |
| 1e-6 | 4.2938e-2 | 1.0647e-2 |

Six independent noisy classical syndrome bits give
`Fencoded=(1-eps)^6 Fideal`. At q=0, readout eps=1e-5 makes correction worse
than leaving the native gate uncoded. Its one-shot break-even eps is about
7.70e-6. At q=1e-8 the same readout rate still improves fidelity. This is a
specified one-shot comparison, not a fault-tolerance threshold.

A separate **actual extraction circuit** uses seven data qubits, one reused
bare ancilla, six checks and24 CNOTs. It corrects all21 single input Paulis
when ideal. Enumerating all360 single-CNOT nonidentity two-qubit Pauli faults
finds291 malignant faults for the stated minimum-weight decoder. Under a
uniform15-Pauli gate fault model the first-order failure coefficient is19.4.
The stored correlated witness is replayed independently on an eight-qubit
statevector, including ancilla preparation and projective measurements.
Every fault is also classified against the actual code subspace.

Thus the bare-ancilla design is **not fault tolerant**. Noisy native CNOT
compilation, flagged/repeated extraction, leakage-aware readout, preparation
and multi-fault correlations remain open. Ideal Steane recovery is useful
but cannot be promoted to a physical universal correction machine.

## Validation and physical boundary

All five producer sections PASS; seventeen focused independent regressions
PASS. Canonical source hashes, exact rational root counting, gauge-transformed
kernels, full-field finite differences, nonlinear mesh moves, full mass
matching, code projectors and an eight-qubit circuit supply distinct checks.
Four-file corpus intake is clean: no rediscovery collisions, no forced
arithmetic findings and no certified-value contradictions. The arithmetic
guard selftest and compilation pass; RESULTS_INDEX was regenerated.

The useful connection is between a real-coupling scalar selector, an actual
full-rank mass map and its complete declared heavy threshold. Its limitations
are equally constructive: background index does not solve chiral matter,
local flatness does not solve curved gravity, equal CP-branch thresholds do
not solve vacuum energy, and ideal decoding does not fix ancilla hooks.
These separate maps still lack a common physically selected action.

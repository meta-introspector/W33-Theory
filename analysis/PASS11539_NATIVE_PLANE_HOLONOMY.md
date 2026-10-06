# Pass11539 — Computing through neutral composite channels

Reserved by `b6a492645`. This packet branches from the previous checklist into
one architecture: use the certified native vacuum to identify **gauge-neutral
composite control channels**, rather than regard motion of a gauge frame as a
physical gate. The exact result is a three-dimensional singlet pair space,
an explicit two-dimensional dark band inside it, and a conditional universal
control construction. It is not a physical processor or a solution of the TOE.

## What is actually established

The inputs are the signed native E6 cubic and 78 Hermitian generators from
`w33_pass11384_11388_native_dynamics.py`, and the two-condensate phase from
`PASS11526_11530_VALIDATED_DYNAMICS.md`. The latter uses two labelled fields,
Phi=e0⊗f0/√2 and Psi=e1⊗f0/√2, with a 28-dimensional unbroken Lie algebra.
This packet does not replace that phase or choose a new physical scale.

**Elementary obstruction.** Under all 28 surviving generators in native
27⊗3, the one-particle singlet space has dimension two. A rank-two band
contained entirely in that fixed two-dimensional space is constant. Therefore
projector holonomy alone cannot implement nontrivial gates while maintaining
that strict elementary-neutrality condition. Checking only SU5 misses the
residual family SU2 and diagonal U1 and is insufficient.

**Composite repair.** The internal antisymmetric pair representation
Λ²(27⊗3), dimension 3240, has exactly three connected-gauge singlets.
The exact calculation imposes six Cartan constraints (21 zero-weight pairs)
and 22 root-generator constraints (rank18). The resulting three vectors have
squared norms **1,10,10**. All 28 induced generators annihilate each vector,
not just its charge sum. The certificate stores every signed coefficient.
This is a finite representation-theoretic result about internal pair states;
physical spin, orbitals, statistics and binding must still be specified.

Let V:C³→Λ²81 have these normalized pair vectors as its columns. For an
independently controlled bright direction b with b†b=1, define the named map

    H_ref(b)=Δ[I−VV†+Vbb†V†],   Δ>0.
    P_dark(b)=V(I3−bb†)V†.

This Hamiltonian has spectrum **0 twice, Δ 3238 times** and commutes with all
28 unbroken generators. Its ground-band rank stays two as b traverses CP².
Four normal-direction curvature values span u(2), so the standard holonomic
control theorem gives full connected U(2) control under supplied adiabatic
bright-direction controls. Gauge generators themselves act trivially in the
three singlet coordinates; they cannot supply these controls.

The cubic and family alternating tensor make this frame explicit. They are
the existing graded E8 matter bracket from
`w33_e6_cubic_hybrid81_transport.py`, used here in adjoint form. The native
and archived canonical45 signed triads agree exactly; this is not a new tensor. For ordered
pair coordinates i f < j g, define

    B_(if wedge jg,mh)=d_ijm epsilon_fgh,   B†B=10I81.
    V(u,v)=[u wedge v, −B conjugate(v)/sqrt10, B conjugate(u)/sqrt10].

Here u,v are the unit-normalized directions of the two labelled condensates
on the certified orbit. At u=e0f0,v=e1f0 these are precisely the stored three
normalized pair vectors. B intertwines conjugate81 into Λ²81; invariance of
the native cubic and the SU3 alternating form supplies the full-group proof.
Consequently V transforms covariantly and fixes the orbit-section ambiguity,
including any disconnected stabilizer that fixes both condensates. This is a
concrete effective coupling map rather than just an existence prescription.
The columns' polynomial numerators have different condensate degrees; their
normalization and physical actuator coefficients still require matching.

The pair Casimirs expose an additional resource boundary:

| Channel | Native E6 pair C | Family pair CF | Full-group sector |
|---|---:|---:|---|
|0, ancilla|−1/9|2/3|Λ²27⊗Sym²3, dimension351×6|
|1, logical|−13/9|−4/3|conjugate27⊗conjugate3, dimension27×3|
|2, logical|−13/9|−4/3|conjugate27⊗conjugate3, dimension27×3|

The standard multiplicity-free decomposition also has a 351prime⊗conjugate3
sector, with no singlets in this phase. Full-gauge invariant Hamiltonians are
scalar on each irreducible sector and **cannot mix the ancilla with the two
logical channels**. Condensate-referenced controls are essential. Attractive
Casimir exchange alone selects a representation channel, not a bound neutral
qubit or its physical mass. The same cubic supplies an explicit way to write
those relative controls; it does not determine their accessible strengths.

## Declared computing process

At b=e0, use the last two normalized pair channels as |0_L> and |1_L>.
Prepare one internal pair per register, traverse closed bright-direction paths
for local gates, and return each register to b=e0 before applying the composite
exchange pulse. Readout uses the two logical-channel projectors. Initialization,
readout and the addressed pair interaction are supplied resources. Repeating
local gates and entangling pulses on connected registers gives the usual
universal circuit architecture; it does not derive those resources from W33.
Gauge neutrality removes charged excursions in this model, but the group acts
trivially on the three-channel multiplicity space and therefore **does not
protect its chosen degeneracy** against arbitrary neutral perturbations.

## Exact native exchange and the cubic connection

For trace Gram G_ab=tr(Ta Tb), define C=Σ_ab(G⁻¹)_ab Ta⊗Tb.
The full 729-dimensional native two-carrier calculation gives

    spectrum(C)={−13/9:27, −1/9:351, 2/9:351},
    SWAP27=(9C²+11C)/2−4I/9.

A separate replay using the native cubic gives D_(27i+j,k)=d_ijk and

    D†D=10I27,
    18C=I729+3 SWAP27−3DD†,
    SWAP27=6C+DD†−I729/3.

Thus the actual cubic used by the native Yukawa construction also identifies
an exchange channel in this declared control model. Casimir completeness and
E6 tensor calculus are classical; the contribution here is the exact signed
native map and its independent replay, not a new general Fierz identity.

For the family fundamental, trace-dual SU3 exchange obeys
SWAP3=CF+I/3. After regrouping native and family indices,
SWAP81=SWAP27⊗SWAP3. For pair registers in slots (1,2) and (3,4), their
whole-register swap is S13(81)S24(81), restricted to antisymmetric pairs.
It preserves the three-dimensional singlet register space and any aligned
logical two-dimensional subspace. The Hamiltonian proportional to this
**product** gives exp(−iπ S13S24/4), an encoded sqrtSWAP.
It is a supplied four-body interaction, not the sequence of two swaps and not
an interaction already derived from the native scalar action. Local U(2) plus
exchange has exact control-Lie rank15 (su4); the certificate stores a replayable
commutator basis. This is conditional universal computation, with its required
resources named rather than concealed.

## Why the elementary geometric route needed this correction

For unit-norm control vectors e0,e1, their joint native moment K has spectrum

    {5/18:2, 1/9:10, −1/18:10, −2/9:5}.
    P=36(K−I/9)(K+I/18)(K+2I/9).
    H0=5I/18−K; gap=1/6.

The physical scalar phase normalizes both vectors by 1/√2. Its moment and the
corresponding constructed gap are half these values. Neither normalization
chooses a physical energy scale or a mass.

The E6 ordered-frame orbit has dimension54. Its equal-weight Berry form has
rank50: the four missing directions are fibers of the moment/Slater-ray
quotient. With unequal positive weights 1,2 the rank is52, but the logical
band splits by (r0−r1)/6. These are first-order coherent-control statements,
not additional scalar modes of the original gauged model. That model already
quotients all scalar gauge-orbit directions.

For Q=I−P the elementary curvature is

    F_ab=P(Ta Q Tb−Tb Q Ta)P.

The native generator pairs (6,7), (6,74), (6,75), (10,11) span u2. Each witness
satisfies T³=T and PTP=0, allowing exact rectangular-loop formulas. This
provides a useful mathematical control witness, but moving the elementary
plane generally leaves the fixed unbroken singlet space. The composite
construction above supplies the missing neutral ancilla instead.

Two regression guards prevent a physics over-read:

- With D=∂−iA and A=−i(dot U)U†, a simultaneous pure gauge change satisfies
  DP=0 and U†DU=0. The apparent Berry connection cancels.
- Filling both levels gives Λ²(C²), dimension1; U acts only through det U.
  A filled two-level band is not a qubit. The composite proposal instead uses
  two particles in one logical pair state drawn from a three-channel pair space.

## Finite-time numerical witness and limits

An elementary 27-state rectangle, epsilon0.6 and supplied flattened gap1, was
integrated with DOP853 (rtol1e−10, atol1e−12):

| Total duration | Worst leakage probability | Geometric operator error |
|---|---:|---:|
|64|0.00417356|0.0468241|
|256|0.000843864|0.0114021|
|512|0.000227026|0.00583742|

Norm drift is below1.6e−11. These are floating convergence witnesses for the
elementary control model, **not** certified finite-time bounds, a composite
simulation, protection against noise, or a hardware threshold.

## Ownership, search and external checks

Existing owners include `2026-09-23_execute_all5_hybrid_lie_golay_calibration.md`
and its graded E8 bracket producer `w33_e6_cubic_hybrid81_transport.py`,
`2026-09-23_cubic_jacobian_rank_stratification.md`, native dynamics11384–11388, two-condensate
11526–11530, SU5/Higgs dictionaries11271, and wedge351/branching7097–7104.
The result-index and corpus searches included the results themselves:
`3240` with singlet/neutral qualifiers, `1,10,10`, rank2/three-singlet bands,
the specific swap/Casimir polynomial, and tripod/holonomic terms.
`docs/index.html`, the paper entrypoints and native owner reports were checked.
BT882 already treats finite subgroup holonomy on the matter graph; that is a
different object from this continuous degenerate-eigenbundle. No claim of
being first to use holonomy, coherent orbits, invariant d×epsilon, wedge351 or universal exchange
is made. The contribution is the fully specified native neutral-pair map and
its obstruction/repair within the certified phase.

Primary external checks used:
- Zanardi–Rasetti, [Holonomic Quantum Computation](https://arxiv.org/abs/quant-ph/9904011): the general degenerate-eigenbundle method.
- Duan–Cirac–Zoller, [Geometric Manipulation of Trapped Ions](https://arxiv.org/abs/quant-ph/0111086): established geometric-gate mechanisms; no ion implementation is inferred here.
- Deppisch, [E6Tensors](https://arxiv.org/abs/1605.05920): classical E6 tensor and invariant-interaction framework.
- Brylinski–Brylinski, [Universal quantum gates](https://arxiv.org/abs/quant-ph/0108062): universality context; this packet also supplies its own su4 replay.

The producer, certificate and focused regression suite accompany this report.
Open physical questions are pair binding, occupation, accessible relative
controls, physical control matching, four-body exchange synthesis and
error budgets. This does not determine masses, mixing, couplings, gravity,
dynamical spacetime or the cosmological constant.

## Intake advisory review

The intake harness reports no forced-arithmetic claims or certified-value
contradictions. Its rediscovery guard flags E6/E8+Singer compounds because
this producer cites the Ambrose–Singer holonomy theorem. The named candidate
release reports563–567,573–577,578–582,583–587, BT1251, BT808 and the two
September23 representation-physics reports were read in entirety. They concern
finite Singer normalizers/cycles or McKean–Singer cochain indices, not this
continuous eigenbundle. These candidates do not establish ownership of the
new neutral-pair application. The actual d×epsilon owner found by deeper
result search is credited above. Candidate warnings are retained as reviewed,
not silently described as absent.

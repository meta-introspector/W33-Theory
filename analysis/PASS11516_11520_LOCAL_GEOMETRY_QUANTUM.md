# Passes11516–11520: full vacuum normals, adjoint flavor, complete metric stress, local projector current and a quantum decoder

Reservation948b8e23d was pushed before computation, after reviewing740dcbf17
and the parallel11511–11515 reservation. All five preceding targets were
investigated. The interval CP-vacuum certificate, dynamical flavor selection,
physical spin identification, reference-free chiral measure and full faulty
hardware remain open. The models below are named separately; no common
physical vacuum is inferred between them.

Producer: `analysis/w33_pass11516_11520_local_geometry_quantum.py`.
Certificate: `data/w33_pass11516_11520_local_geometry_quantum.json`.
Regressions: `tests/test_w33_pass11516_11520_local_geometry_quantum.py`.
Five input certificates and the frozen producer are SHA256-bound. Each
producer function was executed; changed functions were rerun before the
five sections were assembled. Early exploratory outputs are not certificates.

## Ownership and external checks

11271 owns the canonical SM charges and sextet repair;11293 owns mixed-source
mediation;11326/11336/11343 own engineered and constrained flavor selection.
11438 owns the finite native potential;11481 owns the nongauge value fibers;
11493 owns the lower-energy CP stationary witness.11494/11508 own the native
graph metric and certified five-shape determinant.11324/11341 own the supplied
flux/winding framework.11429/11509 own the weak overlap and finite Green
primitive.11475/11480/11510 own Steane recovery, the correlated channel,
the classical MAP decoder and the complete flagged cell instrument.

Result searches included the actual insertion fractions, projector
commutator identity, metric-rank/trace blindness, two-hop paths and fidelity
formula, plus RESULTS_INDEX.md, current papers and docs/index.html. Existing
two-adjoint mediator work is a different construction, and is cited above.
The following general frameworks are established primary literature:

* [Antusch et al., GUT Yukawa Clebsches and adjoint operators](https://arxiv.org/abs/1405.6962).
  No novelty in adjoint flavor splitting generally is claimed.
* [Benzi–Rinelli, sparse spectral-projector decay](https://arxiv.org/abs/2110.11833).
  Gapped projector locality is prior art; here it replaces the inappropriate
  shrinking graph-Laplacian gap in the named relative-current construction.
* [Luscher, local anomaly cohomology](https://arxiv.org/abs/hep-lat/9808021)
  and [chiral reconstruction](https://arxiv.org/abs/hep-lat/9811032).
  Local transgression alone does not meet reconstruction/integrability.
* [Schumacher, entanglement fidelity](https://arxiv.org/abs/quant-ph/9606012)
  and [correlated stabilizer-code fidelity](https://arxiv.org/abs/1005.3374).
  Choi overlap and finite Bayesian decision optimization are standard tools.

##11516 — inspect the entire324-field normal problem

The earlier low-energy point was relaxed on a36-dimensional tripotent
stratum. Compute the full324-by324 numerical Hessian from the analytic
gradient at steps3e-5 and1e-5. Construct all86 common compact-gauge tangent
columns and project onto their orthogonal complement.

At the supplied witness, the numerical gauge rank is58, leaving266 normals.
Both scales find262 eigenvalues above1e-3. The normal spectrum starts with
four unresolved soft values, then approximately0.4252,0.4931,0.7583.
The smallest full eigenvalue moves from-2.53e-7 to-2.61e-8 as the difference
step shrinks. The minimum-direction energy replays increase at the tested
finite displacements. Thus there is no resolved transverse descending
direction in these controls, but neither positivity nor four exact moduli
follow from them.

This substantially enlarges the tested carrier. It does **not** supply an
interval stationary point or remove the weak/flat conditioning problem.
11481's two nongauge value fibers and weak pair retain ownership. A validated
normal/Morse–Bott treatment must distinguish genuine moduli from weak modes;
an isolated-root claim cannot be made by simply discarding small eigenvalues.

##11517 — the adjoint Ward identity fixes the first allowed split in a subclass

Use the actual11271 hypercharge and11507 neutral Hd component in the signed
E6 cubic. For identical fermions with the symmetric family sextet, the
single-adjoint insertion is symmetrized on the two fermion legs. Invariance
of d gives the exact identity

    d(A psi,psi,H)+d(psi,A psi,H)=-d(psi,psi,A H).

For the actual down and charged-lepton blocks, Yi+Yj=1/2 in both cases.
Therefore this single-adjoint correction cannot split their family relation.
Inserting the same adjoint once on each fermion instead gives

    down: Yi Yj=1/18; charged lepton: Yi Yj=-1/2.

The relative coefficient is exactly **-9**. The explicit covariant operator
is d_abc(A psi)^a_i(A psi)^b_j H^c_ij/Lambda², with a transforming E6 adjoint
and H in(27,bar6). Thus the dimension-six operator can split the restrictive
11507 relation while the tested dimension-five operator cannot.

An exact engineered control takes Yd=diag(1,2,3), Ye=diag(1/3,6,3):

    F0=(9Yd+Ye)/10, F2=9(Yd-Ye)/5,
    Yd=F0+F2/18, Ye=F0-F2/2.

The minimality statement is only within this same-adjoint, same-cubic,
symmetrized-fermion-leg subclass. Other E6 representations/operators have
not been excluded. The adjoint vev, sextets and Wilson coefficients are
supplied; these textures are a construction, not dynamical selection,
observed masses, CKM prediction or a new general GUT Clebsch theorem.

##11518 — exact blindness, a native-path repair and a declared full-metric minimum

Every inherited one-hop displacement has equal component magnitudes.
Consequently h_e^T D h_e=0 for every diagonal traceless D. As long as G and
G+D are positive definite, they have exactly the same ungrounded L(G).
The linear metric-observation map has rank4 and a two-dimensional kernel.
**No action depending only on that L can have an isolated full-metric minimum.**

This is compatible with11508's fixed-volume shape minimum: the determinant-one
constraint bends the otherwise blind directions, and the unfixed volume
force was explicitly nonzero there. It is not a retraction of that certificate.

Build324 actual oriented two-hop paths on the same80 vertices, with signed
displacements h_ab+h_bc. Adding those paths raises metric-observation rank
to6; a six-row integer witness has determinant2048. A nonorthogonal second
basis preserves ranks4 and6. The path coupling0.1 and the new determinant
inventory are inputs. Three six-variable searches yield lower anisotropic
controls than the isotropic saddle, but their tiny soft curvature is not
certified. They are not used as a full-volume theorem.

A separate explicit completion reuses the certified11508 determinant W,
adds a declared3-torus with integer3-form flux and three unit axion windings,
and takes

    V(G)=Lambda v+3/(2v)+(v/2)Tr(G^-1)+gamma W(G),
    v=sqrt(detG), gamma=1/1000.

Define the finite counterterm exactly by the rational matrix expression

    Lambda=1-(2gamma/3)Tr[(L+I)^-1 L-(L+4I)^-1 L], L=L11508(I).

This makes the full volume derivative vanish exactly at I. The base stress
has Hessian diag(6,1,3,1,1,1) in the stored logarithmic metric coordinates.
The native shape gradient and volume/shape cross terms vanish by its exact
sign-flip automorphisms and equal-component displacements.

For each of79 nonzero modes, r(lambda)=3lambda/[(4+lambda)(1+lambda)]<=1/3,
and |Wzz_mode|<=1/4. Hence

    Lambda >= 1-158gamma/9 > 0,
    Hessian V >= diag(6-237gamma/4,1,3,1,1,1) > 0.

The shape part also uses11508's rational positive-Hessian certificate.
Thus **I is an exact strict local minimum in all six metric directions**
of this supplied completion. Flux controls v->0, Lambda controls v->infinity,
winding controls the shape boundary, and W is bounded; a global minimum
exists, without a claim that I is the unique global one.

The added stress really depends on geometry beyond graph L. Its flux coupling,
windings and matched counterterm are not W33-derived physical constants.
The graph Dirac inventory remains unmatched to SM spinors. This gives a
constructed stabilization mechanism, not Einstein dynamics or the measured
cosmological constant.

##11519 — a gapped projector current avoids the Laplacian-gap obstruction

For every differentiable orthogonal projector,

    K=[Pdot,P], Pdot=[K,P],
    j_xy=2 Re tr(K_xy P_yx)=-j_yx,
    sum_y j_xy=tr_x Pdot.

Integrate this actual pair current along the weak-overlap homotopy tA,
with the existing anomaly-free charge inventory. Twenty-point quadrature
gives divergence defect about3.62e-14 and simultaneous gauge-conjugation
error below2e-18. These are numerical controls for the exact identity.

The locality argument uses the **fermion gap**, not inverse graph Laplacian.
For a_mu=1-cos p_mu>=0 the free Wilson operator satisfies

    H0²=1+2 sum_(mu<nu) a_mu a_nu >=1

at every torus size. The inherited componentwise weak-field bounds imply
a uniform homotopy gap g>0.565, and ||H||<=7<8. With M=8,

    sign H=(H/M) sum_k binom(2k,k)/4^k (I-H²/M²)^k,
    r=1-g²/M²<1.

Polynomial radii2k+1 and geometric derivative-tail bounds give
volume-independent exponential locality for P,Pdot and their pair current.
Routing pair currents along shortest lattice paths preserves exponential
decay up to polynomial factors. The certificate stores conservative tails.

**Covariance boundary:** the current is invariant under a time-independent
gauge transformation of the entire path, including its flat reference.
Fixing the reference to free and transforming only the endpoint is not that
theorem. A canonical reference-free gauge-invariant current and its measure
curvature/integrability are still unbuilt. Flux sectors are not covered.
This is an explicit local relative transgression, not a completed chiral
measure or a replacement for Luscher's construction.

##11520 — quantum fidelity changes the decoder

Compute all256 logical Pauli-transfer entries by sparse finite stabilizer
contraction, retaining the full native cross-block transfer. The zero-readout
standard decoder exactly reproduces11480's stored16-by16 decoded transfer;
that baseline belongs to11480. The evaluated normalized-Choi expression
for the seven-channel roundoff bound is below4e-15. This norm evaluation,
like the transfer and fidelity calculations, is floating, not interval-certified.

For the11510 complete flagged instrument, let K be its native logical-cell
subchannel and E=(Kdag K) tensor(Kdag K). For normalized logical Choi J,

    F_E=vec(E^T)^dag J vec(E^T)/4,
    F_flagged=(1-v)[(1-v)F_E+v F_(I-E)], v=.001.

The target is the ideal isometry U tensor U applied to two logical qubits
in the actual two12-cells. Accepted wrong physical syndrome sectors and
erasure flags are orthogonal to that target. Accepted140-dimensional cell
leakage is retained and contributes through I-E. Nothing is postselected.

The classical MAP decoder maximizes syndrome probability, whereas quantum
fidelity depends on the conditional logical error. Build all4096 conditional
logical channels and choose, for each12-bit report, both a syndrome pair
and one of16 logical Pauli corrections to maximize the complete-instrument
fidelity contribution. Actual decision maps and all branch scores are stored.
This is optimal within that finite Pauli decision class, not arbitrary recovery.

| readout bit error | classical MAP complete fidelity | quantum decoder complete fidelity |
|---|---|---|
|0|0.81484723|0.88136586|
|.001|0.80726323|0.87283058|
|.01|0.74531637|0.80440806|
|.05|0.63293425|0.67086392|

An explicit additional fault case inserts a two-qubit depolarizing layer
of strength.001 before ideal extraction at each of seven paired positions,
and independent postdecode logical Z faults of probability.002. Quantum
decoder complete fidelity is0.87626482 at zero readout error and0.79905167
at.01 readout error. This is a named fault layer, not the complete gate-level
syndrome circuit or a threshold theorem. Source, ancilla, extraction/reset
gates and general temporal correlations remain unmodeled.

## Validation boundary

Exact algebraic identities, source-bound numerical witnesses, explicit maps,
decision-class optimality and physical open questions are separated above.
The fifteen independent regressions, index/intake and final publication
receipt are recorded in Continuity and SESSION_NOTES.md. No paper claims
or parallel scientific files are changed by this packet.

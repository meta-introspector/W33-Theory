# Passes11521–11525: exact reductions, coupled constraints and an explicit flag circuit

Reservation6388c0ab3 was published before computation, after integrating the
catalog-only0936daaa4. Producer: `w33_pass11521_11525_exact_frontiers.py`.
Certificate: `../data/w33_pass11521_11525_exact_frontiers.json`.
Regressions: `../tests/test_w33_pass11521_11525_exact_frontiers.py`.

All five preceding targets were investigated, changing direction where an
assumption failed or an optimization saturated. An exact reduction is not
an exact physical vacuum. No measured mass, Einstein dynamics, executable
local chiral measure or fault-tolerance threshold is claimed here.

##11521 — the vacuum problem has an exact form smaller than324

11471 owns the signed integer frame and its Spin8 dictionary;11476 owns the
polydisc restrictions;11481 owns the generic value-fiber/rank witness.
Use the canonical frame `[0,17,26]`, with the third sign negative. Embedding
11493's two3x3 matrices there reproduces its low energy and tiny gradient;
the numerical vacuum-frame conjugator is unnecessary for this control.

Instead of introducing square roots to orthonormalize E6 generators,
invert their rational trace Gram matrix. On this frame their restrictions
are diagonal, and the exact row-moment matrix is

    C=(1/3)(I-J/3),
    V=V_U+(1/6)sum_(D=A,B)||diag(DD†)-tr(DD†)/3||².

Here V_U includes the family moment term, determinants, norm terms and
the cross-field CP term. It is invariant under common-left SU3:
det(UA)=detA, A†A is unchanged, and (UA)^T(UB)*=A^T B*.
The penalty and its entire first derivative vanish when both row Grams
are hollow (equal diagonal entries).

**Conditional stationary-fiber theorem:** if an exact reduced stationary
pair has both Grams hollow, every common-left SU3 transform retaining
hollow Grams is stationary as well. The canonical Spin8 fixes this
stratum and has no transverse fixed vector, so reduced stationarity
implies full stationarity. At the regular rank-four locus the fiber is
four-dimensional before the two diagonal gauge directions are removed.
11481's exact generic ranks establish the available two nongauge tangents;
they do not locate an exact stationary pair on that regular locus.

There is also an exact transverse reduction:

    H324=H36 ⊕ (Bv⊗I8) ⊕ (Bs⊗I8) ⊕ (Bc⊗I8),

with three symmetric12x12 matrices. Integer commutant equations have
rank63 on each8 and intertwiner equations rank64 between every pair.
Thus each real8 has scalar commutant and the three modules are inequivalent.
Schur's lemma proves the displayed form at every point of the fixed stratum.
Floating QR selects six minors; exact ranks modulo prime65521 certify the lower bounds,
with scalar-commutant and ambient-dimension upper bounds closing equality.
Every minor is stored. Numerical full-gradient Hessian blocks provide a
separate replay of the representation reduction, not certified masses.

Classical simultaneous hollowisation also gives the exact minimization
identity min_common-left V(UA,UB)=V_U(A,B): every pair of traceless Hermitian
Grams can be hollowed together. Thus the global infimum on the polydisc is
that of V_U. On its regular SVD chart, write A=e^(i phase)diag(a1,a2,a3)
and fix B12/B23 real positive with residual diagonal SU3. The quotient has
20 real coordinates. A70-digit stationary solve is retained as a numerical
control: gradient maximum6.02e-67, energy-2.93574323966686637 and
smallest numerical quotient-Hessian eigenvalue1.31551. This removes the
previous weak-coordinate conditioning from the stationary solve; it is
still not interval root isolation or a transverse mass certificate.
[Damm–Fassbender](https://arxiv.org/abs/1910.08813) and11481 retain
simultaneous hollowisation ownership.

This explains the correct next mathematical object: a stationary hollow
pair, its genuine fiber, and the reduced normal blocks. It does not turn
11493's approximate pair into an interval/Morse–Bott certificate. In
particular, the four unresolved normal eigenvalues are not all declared moduli.

##11522 — a smaller Higgs inventory has an exact symmetry boundary

Build the actual quadratic covariant

    H^c_ij=dbar^{cab} Phi*_ai Phi*_bj/M,

in(27,bar6), allowing the symmetric identical-Weyl Yukawa without an
elementary sextet. The native signed cubic's infinitesimal covariance is
checked coefficientwise. This is a different composite from11276's
Sym(F⊗F⊗F) spectator;11271/11281 retain the original sextet repair and
selection obstruction. No general composite-Higgs novelty is claimed.

The tempting shortcut nevertheless fails in the SM-preserving sector.
Every polynomial covariant of SM-invariant vevs is SM-invariant. More
directly, native hypercharge conservation requires Yc=-Ya-Yb. Two Y=0
vevs cannot populate the actual Hu/Hd components with Y=±1/2. Their
quadratic source coefficients vanish exactly. SM-preserving adjoint
insertions do not alter this conclusion.

For repeated identical fermion legs the true insertion algebra is the
symmetric polynomial ring Q[Yi+Yj,YiYj]. The first generator is fixed by
the Higgs charge; the second distinguishes1/18 from-1/2.11517 owns the
first allowed split, whose ratio is-9. More insertions alone do not supply
a vacuum-selection principle.11438's potential contains neither the
independent adjoint nor the sextets, so their vacuum equations and Wilson
coefficients cannot be inferred from its stationary point. A composite
can reduce inventory, but electroweak-breaking dynamics still has to be named.

##11523 — stabilizing a metric does not satisfy its static lapse equation

Keep the actual11518 determinant and vary the complete declared stress:

    V(G)=Lambda v+alpha/v+(beta v/2)tr G^-1+gamma W(G),
    v=sqrt(detG), G=e^z I.

AtI, stationarity sets

    Lambda=alpha-beta/2-(2gamma/3)Wz,
    V0=2alpha+beta+gamma(W0-2Wz/3).

For this inventory W0=sum log[(lambda+4)/(lambda+1)] and
Wz=sum3lambda/[(lambda+4)(lambda+1)] is **positive**. The exact inequality
log(1+x)>=x/(1+x), x=3/(lambda+1), proves W0>=Wz>=0. Therefore

    V0>=2alpha+beta+gamma W0/3>0

for positive flux/winding coefficients and gamma>=0. Simultaneously
requiring stationarity and V0=0 gives

    alpha=-beta/2-gamma W0/2+gamma Wz/3<0,

contradicting positive flux energy. The11518 point remains a strict local
metric minimum; this does not retract that theorem. It is incompatible
with a **static flat solution of11482's particular lapse action**, where
K+N²V=0 and zero velocities imply V=0. Time evolution, spatial curvature
or a different stress completion is required. This is not a theorem
excluding flux stabilization in other geometries or theories.

The native graph incidence Dirac also has a precise carrier boundary:
80 vertex states plus160 edge states,79 nonzero Laplacian modes and82
zero modes. Exactly,

    det(mI+iDgraph)=m^82 product_(j=1)^79(m²+lambda_j)
                 =m^80 det(m²I80+L).

An explicit supplied spin lift now names the map. Let spatial gamma
matrices obey the Clifford relations and set Ce=80 Gamma(h)/||h||² times
the endpoint difference, using actual integer h_times80. Gamma(h)²=||h||²I
proves Cdag C=L⊗I4 exactly. The960-dimensional lifted incidence Dirac has
328 zero modes and determinant m^320 det(m²I80+L)^4. A Lorentzian time
gamma anticommutes with its spatial gammas and squares to-I. This is a
constructed cochain spin lift, with additional fibers and surplus modes;
it is not the original fermion inventory or an identified SM spinor.

Its determinant does not specify a4D Clifford principal symbol, spin
structure, tetrad or spin connection. The paper's genuine Spin(1,9)
Albert matter construction is a different internal map. An explicit
site-spinor operator and its full determinant must be built before
matching these inventories. The established geometry method is described
by [Brower et al.](https://arxiv.org/abs/1610.08587); no invention of
lattice spin transport is claimed.

##11524 — the missing chiral-measure term is an exact one-form

The relative density current of11519 is not itself a measure connection.
For a projector P and its variation dP, P dP P=0 gives exactly

    tr P[[omega,P],dP]=-tr omega dP.

This is the gauge variation that the radial Berry term must cancel.
If k is a gauge-invariant, spatially local nonlinear primitive with
div k=anomaly, the required correction is

    delta_eta S(A), S(A)=integral_0^1 <A,k(tA)>dt,
    delta_eta S=integral_0^1 [<eta,k(tA)>+t<A,Dk(tA)[eta]>]dt.

It cancels the radial gauge variation while its exterior derivative is
zero, preserving the required curvature. This is the standard construction
of [Luscher, equation5.8 and theorems5.1–5.4](https://arxiv.org/pdf/hep-lat/9811032).
11438 already owns the charge-condition audit, and11484 already identifies
the local primitive as the missing executable object. These are cited,
not rediscovered. Locality, admissibility, sufficiently large finite volume
and sector corrections remain conditions of the imported construction.

The new exact projector/chain-rule checks isolate why the relative
current cannot simply be relabelled a completed measure. The nonlinear
local k and finite-volume correction still require implementation;
the Coulomb primitive's massless inverse does not discharge locality.

##11525 — certify recovery limits, then repair an actual extraction hook

11520 owns all conditional logical channels and the Pauli decoder. Reuse
its sparse contraction, preserving coherence. Allow any CPTP recovery
on a decoded four-dimensional logical fiber through the SDP

    max tr(QX), X>=0, tr_output X=I4;
    dual min trY, I4⊗Y>=Q.

Q includes the actual native-cell effects and verification flags.
An independent explicit-Kraus test fixes the Choi factor ordering.
This is the standard [Fletcher–Shor–Win recovery optimization](https://arxiv.org/abs/quant-ph/0606035),
not a new SDP theorem. Sixteen dominant-branch SDPs find no resolved
improvement over the prior Pauli choices. Their solver labels are not
certificates.

Store an exact Gaussian-rational Pauli primal and a feasible rational
dual for the leading branch. Positive exact LDL pivots of the dual slack
bound every CPTP recovery for that stored dyadic objective. Additionally,
all4096 Bell-basis objectives are stored as integer numerators. Exact
Gershgorin row bounds give scalar duals and a complete gain cap; tests
replay every integer bound. These exact optimization certificates apply
to the rounded objectives. The native channel contraction and its bridge
to those objectives remain floating, not interval-certified hardware.

The exact code symmetry tightens the all-branch result further. Four
three-bit syndrome vectors are classified under simultaneous GL3(2) by
their relation kernels:66 orbits cover all4096 branches. The group permutes
seven Hamming columns, fixes both Steane logical states and commutes with
the seven identical paired channels. Every orbit receives a full CPTP SDP
and an exact Gaussian-rational feasible dual, certified by positive LDL
pivots. Their weighted rounded-objective bounds are

    0.88136585817 <= optimal complete fidelity <= 0.88136929909,

so coherent recovery can improve these stored objectives by at most
3.45e-6. The exact result is for the66 representative dyadic objectives
and their orbit multiplicities; the physical coefficient bridge remains
floating. Native orbit-channel replays agree within1.12e-16. Independent
finite-field group enumeration checks every orbit, and every dual is replayed.

Branch-wise coherent recovery is not the only physical bottleneck.
In an actual weight-four Steane Z-check, data-control CNOTs target a
syndrome ancilla. A syndrome-ancilla Z fault after the second coupling
propagates to the last two data controls without flipping the ancilla's
Z readout. Later ideal X checks and the ordinary single-error decoder
complete a three-qubit Fano line: **logical Z**. A sole Bernoulli fault p
there gives entanglement fidelity1-p even with perfect remaining readout.

Construct the repair, rather than stopping at that failure: initialize
a flag ancilla in|+>, bracket all four couplings with two flag-control
CNOTs targeting the syndrome ancilla, and measure the flag inX.
The flag CNOTs commute with the data CNOTs and cancel in the fault-free
circuit. Enumerate all seven positions of a single syndrome-ancilla Z
fault. The flag outcome and subsequent X syndrome select stored
hook-aware corrections; every residual is exactly a Steane stabilizer.
Thus **this entire named fault family is corrected** by an actual circuit
and lookup table. Faulty flag gates, arbitrary CNOT faults, multiple
faults and faulty readout are not covered by that statement.

## Validation and retained boundaries

The producer separates exact rational identities/minors/dual bounds from
full-gradient numerical controls and imported existence theorems. The
certificate binds six canonical JSON inputs and the producer. Focused
regressions, result searches, corpus guard/intake and publication receipt
are recorded in Continuity. No paper or parallel scientific file is changed.

Final focused validation:19 tests passed in230.50s. Regenerate with NumPy,
SciPy, SymPy, mpmath and CVXPY/CLARABEL; certificate regressions do not need
CVXPY because they check stored exact duals independently. The producer's
five functions were executed separately, changed sections rerun, and the
certificate assembled atomically with its frozen producer/input hashes.

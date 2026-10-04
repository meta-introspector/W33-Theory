# Passes11413–11417: a full-family pin, native Regge patch, declared vacuum and counted noise

Reservation:77dc3d5fb; remote-inventory merge:e129a91b3.
Producer:`analysis/w33_pass11413_11417_pin_regge_vacuum_noise.py`.
Certificate:`data/w33_pass11413_11417_pin_regge_vacuum_noise.json`.
Independent regressions:`tests/test_w33_pass11413_11417_pin_regge_vacuum_noise.py`.
All five producer fronts pass. The packet is a set of conditional models,
not a solution of the TOE or a prediction of measured parameters.

## Intake, prior ownership and what is being tested

The latest remote change df91a4354 only refreshed the complete formula-search
inventory. GitKraken integrated it before the reservation was published.
This work follows11408–11412 and cites the earlier native carriers in
11389,11390–11394 and11403–11407. Existing paper/index sections on Yukawa
textures, seesaw, mirror selection, finite gravity and vacuum counterterms
were checked alongside result searches, not treated as physical predictions.

Prior owners include:

* `analysis/w33_pass235_yukawa_texture.py`: an earlier Froggatt–Nielsen
  texture ansatz; its all-one coefficient example is rankone, so it supplies
  hints rather than a verified three-mass model.
* `analysis/w33_wolfenstein_substrate.py`: numerical parameter identifications,
  without a constructed common microscopic mass map.
* `analysis/PASS11408_11412_CLOCK_SPURION_SEAM_TICKS.md`: common scalar-pin
  rankone obstruction, three-conserved-copy repair, stable supplied orbit
  mixing, native Majorana rank2 mask and the counted gate calibration.
* `analysis/PASS11270_11274_PHYSICAL_FRONTIERS.md`: E6/Higgs/chirality
  boundaries, hard loops and a nonvanishing vacuum supertrace.
* `analysis/BT986_edgewise_regge_sphere_convergence.md`: a2D curvature proxy,
  distinct from the4D length-Regge action constructed here.
* `analysis/w33_discrete_einstein_hilbert.py`: prior equilateral triangle
  counting/Euler-characteristic action, not this variable4D simplex patch.
* `analysis/w33_pass11225_flavour_from_symmetry.py` and
  `analysis/PASS11225_FLAVOUR_FROM_SYMMETRY.md`: the complete residual-symmetry
  scan found no Cabibbo or full PMNS pattern; its certificate was reread.
  `analysis/PASS11236_CUBIC_FLAVOUR.md` extends the restricted no-go to
  one-cubic-gate residual eigenbases. Our pin coefficients supply dynamics
  beyond those finite residual-symmetry scans; they do not refute them.
* `exploration/w33_levi_selector_amplitude_bridge.py`: prior rational
  visible/null selector-amplitude identities, including9/40 and3/80.
  This is different from the endpoint propagator and full pin matching.
  Its generated summary JSON and the old triangle-count JSON are not present
  in this worktree, so their source claims are not upgraded to certificates.
* Pass11301: finite pointwise Leibniz obstruction; Pass11323: native lapse
  Hessian obstruction. Regge geometry does not repeal either result.
* `exploration/PART_CCCCIV_W33_CSS_STEANE_LIFT.py`: the prior Steane lift.
  This packet tests an explicit ideal channel and actual code states;
  it does not claim a new code or a physical correction threshold.

The generic FN mechanism, Regge vertex-translation symmetry, one-loop
determinant formula and Steane code are established external results.
The concrete contribution is their explicitly scoped native-map witnesses
and the compatibility/criticality controls below. Corpus searches found
the prior scalar-pin repair, not this full-family pin construction;
this is a bounded search finding, not a universal novelty claim.

During validation, parallel commit caf59a13e added11369–11373 and its paper
corrections; GitKraken fast-forwarded those changes plus inventory0c234a443.
All five new reports were read, and their certificate keys and11371 census
were inspected.11369 refutes global J6 completeness on PU(3), while keeping
the finite word census.11370 gives a twisted-indicator degree law;11372
separates exact design rates from the empirical reversible-fraction decay;
11373 gives exact three-qutrit fractions. None changes the source certificates
used here. The parallel contribution is integrated, not overwritten.

`analysis/PASS11371_FLAVOUR_VS_ARROW.md` now supplies a direct tick-as-mixing
map and separates its transpose-even flavour J from the substrate arrow.
Our full-pin mass map is a different, supplied-source construction. Its
CP phase must likewise not be called the microscopic arrow of time.
AZ/Floquet spectral chiral pairing in11371 also differs from a nonzero
Lorentzian Weyl index: it cannot remove our balanced-grading obstruction.
The actual pin eigenframes obey the same transpose-even J identity.

##11413: mixing at the pin without shortening the protected paths

Let A be the actual80-vertex Levi adjacency from signed incidence D:
A=4I-DD^T. Put a3-dimensional family fibre on every node of two messenger
networks. For K=I-eta A and a full3-by3 family bridge B at vertex0, use

    M = [[K tensor I3, P0 tensor B], [0, K tensor I3]].

Attach light family j to vertex i_j and fibre j on each network. The
offdiagonal inverse block is exactly

    -(K^-1 P0 K^-1) tensor B,

so its three selected endpoint matrix elements give

    Y = D_R B D_R,   D_R=diag((K^-1)_(i_j,0)).

The full neutral-family messenger dimension is480. The independent test
solves that block directly. This dimension excludes colour multiplicity,
weak partners and any unspecified ultraviolet scalar sector.

Unlike a scalar common pin, this matrix has rank3 when B is invertible and
all three endpoint propagators are nonzero. It also allows up/down mixing
without separately conserved family-copy tags. Independent node phases act
uniformly on each family fibre. Mixing only at the pin therefore does not
create a shorter cross-network route: the analytic link-spurion degree of
entry ij remains at least d_i+d_j. This is conditional perturbative
protection by the declared architecture, not an anomaly-certified gauge
symmetry or a nonperturbative theorem about a completed UV theory.

Choose native endpoints [44,1,0], distances [3,2,0]. Each has one shortest
path to the pin. The supplied pin matrices are

    B_u = diag(1,2,3),
    B_d = diag(2,3,4)+z X+z* X^T,  z=.2+.15i,

where X is cyclic shift. B_d is Hermitian positive definite. Both matrices
embed into the native family/colour frame as C(Y tensor I3)C^dagger,
and commute with the actual native colour generators. The prior rectangular
Higgs/Yukawa traces remain gauge invariant for these family coefficients.
Pin entries, endpoints, link couplings and Higgs data are supplied inputs.

The propagators have leading powers [eta^3,eta^2,1]. For a positive pin
matrix with nonvanishing Schur pivots this implies mass orders
[eta^6,eta^4,1]. The down-sector leading coefficients are approximately
[1.9624293,2.984375,4]. Let

    S = B_d[0:2,0:2]-B_d[0:2,2] B_d[2,0:2]/B_d[2,2].

Then the leading three mixing-angle coefficients are

    |S01|/S11, |B12|/B22, |B02|/B22
      = [.0857527055,.0625,.0625],

with respective powers [eta,eta^2,eta^3]. This resolves the restricted
equal-offdiagonal competing-clock correlation from11408 within a different
explicit action. It does not determine a Cabibbo angle.

For these Hermitian Yukawas, their eigenframes also diagonalize their left
mass Grams. The triangle identity gives

    |J| = |Im(Y01 Y12 Y20)| / product_(i<j)(m_j-m_i).

Consequently |J| has order eta^6, with leading coefficient
|Im(B01 B12 B20)|/(B22^2 S11). Complex pin data explicitly break CP;
conjugating the source reverses J. No spontaneous or observed CP phase is
deduced. The producer and independent test check both operations.

### Criticality removes the hierarchy

The actual distance shells have sizes [1,4,12,36,27]. Their exact radial
resolvent, ordered by distance0 through4, is

    [1-18eta^2+36eta^4, eta(1-15eta^2),
     eta^2(1-12eta^2), eta^3, 4eta^4]
        / [(1-16eta^2)(1-6eta^2)].

The symbolic regression checks all five radial recurrences. All80 entries
are compared with the direct inverse at four couplings. As eta tends to1/4
from below, every radial endpoint/pin ratio tends to1. The Perron pole
screens the path hierarchy instead of enhancing it. Thus the weak-coupling
power counting cannot simply be extrapolated to eta≈.225 as a CKM fit.
Generic FN power counting is established; this native propagator supplies
a quantitative limit on that analogy.

##11414: a full-rank Majorana/Higgs minimum with its source exposed

Use a complex symmetric3-by3 SM-singlet field S and a Higgs doublet h.
Declare an explicit H27-breaking source S0=diag(1,2,3) and the potential

    V = lambda_H(h^dag h-v^2/2)^2 + t||S-S0||_F^2
      + lambda_S||SS^dag-S0 S0^dag||_F^2
      + kappa h^dag h ||S-S0||_F^2.

All coefficients are positive: [.13,.3,.2,.1], with v=1. Every term is
nonnegative, so S=S0 and h^dag h=v^2/2 is a global minimum. The source term
makes all12 real symmetric-field directions strictly stable. The Higgs has
one radial mass squared2 lambda_H v^2 and3 electroweak Goldstones.
The coefficient t is dimensionful in a dimensionful implementation; here
all numerical masses are expressed in one supplied reference scale.
Kinetic coordinates have canonical1/2 Euclidean normalization.

The certificate lists compensating character charges for all six symmetric
Majorana entries. The unbroken11409 mask only allowed the23 pairing; the
diagonal source explicitly supplies family-breaking channels. The source
is assumed neutral under the SM, using the prior bare-Majorana condition
h_charge=3q. A charged source would change that argument.

Lepton fibres are explicitly attached at the pin: Ynu=.01 B_d. They do not
inherit the quark distance hierarchy. This declared choice avoids claiming
an unresolvable light seesaw singular value: an earlier trial copied the
quark hierarchy and put the smallest mass below double-precision accuracy.
That failed independent regression is recorded in Continuity.

The full symmetric6-by6 matrix has blocks [0,vYnu/sqrt2; transpose,S0].
Its six positive Takagi masses are checked by SVD, its exact determinant
identity and its small-Dirac seesaw limit. This is an explicit conditional
vacuum, not a source-free origin of its three Majorana scales.

Chirality remains an obstruction: the12-by12 selfadjoint completion has
balanced grading and index zero. Stable internal mass data cannot select
Lorentzian Weyl roles or remove mirrors. The calculation retains this
negative result rather than calling a full-rank mass matrix a chiral theory.

##11415: nonlinear flat vertex motions on a native4D Regge patch

Use the index5 common FCC periods from11411 as three spatial boundary
edges, and add a supplied fourth, Euclidean time edge of length2. The
five boundary vertices define one4-simplex; an interior centroid subdivides
it into five4-simplices. Its5 radial edges are variable; boundary lengths
are fixed. The action includes both internal deficits and the Regge boundary
term:

    S_R = sum_internal A_t(2pi-sum theta_t)
        + sum_boundary A_t(pi-sum theta_t).

Actual finite interior displacements leave every internal deficit zero and
preserve the boundary action. The full Hessian at the flat point has four
near-zero numerical eigenvalues and one nonzero curvature eigenvalue,
approximately1232.10. Vertex-translation tangent vectors are stored.
Independent Cartesian face normals check every Gram-derived angle.
Finite-difference refinement decreases the tangent residual; the final
stored step is5e-5, with a separate3e-5 regression and8e-5 comparison.
The near-zero eigenvalues are diagnostics, not an exact numerical rank
certificate; the flat vertex family is the geometrical reason for them.

Off-flat radial-normal perturbations produce deficits and lift all four
flat tangent directions. These are off-shell probes, not curved solutions
or a ghost diagnosis. This reproduces the standard Regge symmetry boundary
on explicit native periods. Newton's constant, Lorentzian signature,
microscopic metric dynamics, refinement and the nonlinear constraint algebra
remain unbuilt. Flat symmetry alone does not establish nonlinear gravity.

##11416: the complete determinant of the declared low-energy model

The inventory is48 Weyl fields,12 real singlet scalar coordinates,4 real
Higgs coordinates and12 gauge generators. At the chosen stationary
background there are13 massive physical real scalars,3 Goldstones,
3 massive vector bosons and9 massless gauge bosons. The three Dirac up and
down masses have colour multiplicity3; charged leptons have3 Dirac masses;
the full seesaw has6 Majorana masses. Fermionic weight sums to96 real
degrees of freedom. Charged-lepton Yukawas and gauge couplings are supplied.

With Landau-gauge MSbar conventions and scale mu=1:

    V1 = sum n_a m_a^4 [ln(m_a^2/mu^2)-c_a]/(64pi^2),

where c=3/2 for scalars/fermions and5/6 for massive vectors. Zero-mass
Goldstones, Landau ghosts, gluons and photon give zero at this background.
The stored result is V1≈-0.30872745 and STr M^4≈-784.15631546 in supplied
reference units. The independent scale test verifies
dV1/dln(mu)=-STr M^4/(32pi^2), holding masses fixed.

This is complete for the declared low-energy quadratic spectrum. The
480-dimensional heavy messenger family blocks have been integrated out;
their matching thresholds, colour/weak completion and dynamical link-scalar
potential are not specified or counted. It is therefore not a complete
UV/TOE vacuum calculation. The nonzero supertrace and freely renormalized
vacuum counterterm leave the cosmological constant undetermined. Evaluating
one loop at a tree minimum does not prove a quantum-stable minimum.

##11417: remove phase windings, then test a specified correction channel

An ideal rankone projector kick is exactly2pi periodic in its phase.
Reduce the five prior Hadamard pulse angles to their principal values and
recalibrate the actual finite-penalty Floquet shifts. The exact160-edge
unitary replay now uses29091 ticks instead of140082, a4.8153-fold reduction.
The full encoded coherent error is.00796074, versus.00636689 previously;
the conservative dressing/rounding bound is.0499972. Tick savings come
with a slightly larger coherent error at this particular calibration.
No assumption of coherent-error cancellation is used.

Separately, declare independent logical Z errors during an idle interval
of the same tick count, followed by perfect Steane syndrome and recovery.
The aggregate idle probability is p=[1-(1-2q)^N]/2. Enumerating all128
Z-error patterns gives uncorrectable counts by weight
[0,0,21,7,28,0,7,1]. The independent test constructs the actual128-by2
code basis and checks Knill–Laflamme for all21 single-qubit Paulis plus
identity, rather than using the code parameters as a substitute for recovery.

At q=1e-8 and N=29091, p≈2.90825e-4 becomes logical Z failure≈1.77376e-6
under this ideal recovery. This does not include the compiled coherent
error, native during-gate edge noise, massive leakage, noisy syndrome
extraction or physical correction hardware. It supplies neither a physical
threshold nor a fault-tolerant universal architecture.

## External checks used

* Froggatt and Nielsen, hierarchy/angles/CP,1979,
  [CERN bibliographic record](https://cds.cern.ch/record/133050),
  [original DOI](https://doi.org/10.1016/0550-3213(79)90316-X).
  The supplied-link power mechanism is FN-type, not a new generic principle.
* Bahr and Dittrich,
  [(Broken) Gauge Symmetries and Constraints in Regge Calculus](https://arxiv.org/abs/0905.1670).
  Flat vertex symmetry and broken curved constraints are established.
* Martin, [general effective-potential conventions](https://arxiv.org/html/hep-ph/0111209v2),
  sections1 and3. The one-loop determinant and scheme constants are standard.
* Steane,
  [Multiple Particle Interference and Quantum Error Correction](https://arxiv.org/abs/quant-ph/9601029).
  Code/recovery theory is prior art; the new work is the explicit declared
  noise calculation beside the counted native gate.

## Remaining independent directions

1. Replace the supplied full-family pin matrices by a native invariant
   potential and determine whether its noncommuting minimum retains path protection.
2. Construct a chiral spacetime/domain-wall extension with an actual index
   and an explicit native-matter coupling, rather than balanced finite grading.
3. Test a coarse-grained/perfect Regge action on glued curved FCC patches
   for nonlinear gauge directions beyond the flat branch.
4. Complete the heavy messenger/link-scalar ultraviolet spectrum and match
   its vacuum thresholds before proposing any radiative cancellation mechanism.
5. Simulate native during-gate edge noise and leakage with noisy syndrome
   circuits, then test whether the tick reduction survives correction overhead.

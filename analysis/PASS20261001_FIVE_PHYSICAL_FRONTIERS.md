# Five physical frontiers: exact fields, null-mode completion, full loops, curved refinement, and anomalous scales

This packet executes the five research targets following 74b52332e. It provides
stronger mathematical interfaces and explicit models, not measured Standard
Model parameters or a solved theory of everything.

## 1. Exactly balanced fields and the universal mass law

The canonical signed27 x3 coordinate order is retained. Three integer columns
have entries one at flat indices Q0=(0,52,80), Q1=(5,37,75), Q2=(26,40,54).
The78 stored E6 matrices have exact integer entries. The producer checks
Q†Q=3I, every E6 cross moment Qj† B Qk=0, SU3 partial cross Grams
Vj†Vk=delta_jk I3, and all pairwise signed cubic brackets. Every field
Phi=Qq/sqrt(3) is therefore exactly D-flat with the standard norm ||q||².
Moment-zero/closed-orbit theory and the prior rank-three theta classification
identify this as a Cartan plane. This alternative plane does not make the
old large floating gauge transformation algebraic.

The mass matrices have support components: three3-blocks, nine6-blocks,
one18-block. Cycle3-blocks supply one zero and two masses |qi|.
Each6-block is bipartite with relative signed permutation U satisfying
U³=-I and trU=trU²=0. The18-block has commuting order-three U,V with
tr(U^a V^b)=9delta_(a,b),(0,0), the regular Z3 x Z3 character.

These integer identities prove the squared singular spectrum for every complex q:

- |qi|² twice each;
- |q0+omega^a q1+omega^b q2|²/3 twice for every a,b;
- |qi-omega^k qj|²/3 six times for each i<j,k;
- three zeros, with additional mirror degeneracies.

This is the prior12-MUB/twofold and9-SIC/sixfold mirror law, now certified
on the exact plane rather than inferred from five directions. Moments are
TrM†M=20||q||² and Tr(M†M)²=8||q||⁴.
Producer/certificate: w33_20261001_exact_cartan_mirror_completion.py/.json.

## 2. Three light masses without changing the heavy spectrum

The global invariants are now explicit polynomial contraction circuits. Form
X_a(alpha,beta)=d_aij Phi_i(alpha)Phi_j(beta), using the signed E6 cubic.
Contract three X tensors with the dual cubic and two SL3 volume forms to get
raw6. Form the ternary cubic F(alpha,beta,gamma)=d_ijk Phi_i(alpha)Phi_j(beta)Phi_k(gamma)
and its four-F Aronhold contraction S. On the unnormalized Q plane, exact
symbolic expansion gives raw6=12u6 and S=(152u6²-32u12)/5. For E=Q/sqrt(3),
the global normalized polynomials are

    I6=(9/4)raw6,
    I12=(152I6²-3645S)/32.

All156 lower/dual E6 generator checks and8 SL3 volume-form checks vanish
exactly. Their analytic chain-rule differentials work on every complex81
field, including off-plane/nonregular inputs. Independent finite differences
and complexified gauge transport test the evaluator. The contractions are
classical invariant theory; the signed-field normalization and executable
covector interface are the result here. No enormous expanded monomial table
or gauge-fixing oracle is required.
Producer/certificate: w33_20261001_global_e6_cartan_covariants.py/.json.

For two distinct fermion species, the named compact-gauge-invariant operators are:

    (c5/Lambda) (Phi†chi)(Phi†psi)
    (c6/Lambda^9) dI6[chi] dI6[psi]
    (c12/Lambda^21) dI12[chi] dI12[psi]

Spinor epsilon contractions and Hermitian conjugates are understood. Operator
dimensions5,13,25 determine cutoff powers1,9,21; these are nonsupersymmetric
EFT interactions, with no string or holomorphic permission assumed.

At unit T the three covectors conj(T),grad u6,grad u12 have exact cyclotomic
Hermitian Gram diag(1,16,100). The latter two norms belong to the earlier
global-selector packet. The covectors span the normal space, lie in kerM†,
and their adjoint inputs lie in kerM. All cross terms with the heavy matrix
vanish. The completed map has rank81, preserving all78 heavy singular masses.
The three new masses are |c5|r²/Lambda,16|c6|r^10/Lambda^9,
100|c12|r^22/Lambda^21. A finite control supplies .001,.002,.003.

This is a conditional hierarchy, not an observed-family identification or
prediction of coefficients/scales. The earlier dimension-five packet owns
the first lift; the two gradient lifts complete the other normal directions.
The global circuit and differential compiler is built and independently
checked. Full higher-operator loops and physical assignments remain open.

A positive full81-field classical potential is also now explicit:

    lambda(||Phi||²-r0²)²+alpha|I6(Phi)|²+beta|I12(Phi)|²+kappa sum mu_a(Phi)².

Its global minimum is zero. Every minimizer is moment-zero/semisimple, so
Cartan conjugacy and uniqueness of the compact moment-zero representative
reduce it to the prior exact72-ray intersection at fixed radius. These72
Cartan rays are little-Weyl-related **gauge representatives**, not72 distinct
physical vacua or prepared magic states. The common phase, chosen radius,
coefficients and quantum stability remain separate inputs/problems.

## 3. Full declared gauge/scalar/fermion loop: a stronger obstruction

Field content: one canonical complex81 scalar, V_D=(1/2)sum g_a²mu_a²,
E6 x SU3 gauge fields, two original81 Weyl multiplets and their two conjugates.
Conjugate pairing cancels perturbative gauge anomalies. A specified Z2
separates sectors and excludes direct cross-sector bare masses. No Standard
Model assignment is made.

Rational Hermitian bases and inverse trace Grams certify all nine coefficients in

    2sum_E6(TPhi)(TPhi)†+(1/3)sum_SU3(TPhi)(TPhi)†=M†M/3.

At g_SU3=g_E6/sqrt(6), gauge squared masses are g_E6²/3 times the78 nonzero
cubic squared masses; eight gauge generators remain massless. The moment
potential gives the same78 massive real scalars,78 Goldstone and6 Cartan zeros.

MSbar Landau-gauge weights are3,1,-8 for vector/scalar/paired-fermion modes;
subtraction constants are5/6,3/2,3/2. Exact moments make subtraction and
scale terms angularly constant:

    64pi² V1_angular=(4g^4/9-8|y|^4)F(q),
    F=sum multiplicity sigma^4 log(sigma²).

The previous outward interval Rayleigh witnesses establish both Hessian signs
at T. T is a saddle for every nonzero prefactor; g^4=18|y|^4 gives angular
flatness at one loop, not a selected minimum. This extends the fermion-only
diagnostic to the complete declared renormalizable field determinant.
The higher-dimensional completion coefficients are zero in this test;
their loops, other relative couplings and global/two-loop stability remain open.
Producer/certificate: w33_20261001_complete_gauge_scalar_fermion_loop.py/.json.

## 4. Actual W33 mass fiber on a discrete curved refinement

The supplied four-torus is discretized with a Hermitian nearest-neighbor
Wilson operator and conformal weights D(g)=bD0b. A discrete exponential
average fixes volume exactly at every N. The finite perturbation formula
is independently replayed against a full small Fourier block at nonzero
epsilon. Exact signed-permutation multiplicities cover the entire N^4 lattice.

Fixed-time refinements and Richardson controls are stored in the certificate.
Wilson heavy modes cause substantial coarse artifacts. Removing the Wilson
term gives a wrong-sign N96 curvature coefficient: the failed coarse doubling
check is retained and replaced by a finer sequence. Neither extrapolation
nor a finite sequence is an interval-certified limit theorem.

The fiber is the actual signed81 mass map with its earlier one-mode lift,
self-adjointly doubled to162. Wilson breaks external chiral anticommutation.
Use the internal grading instead:

    D_total=D_Wilson(g) tensor Gamma_F+I tensor D_F,
    Gamma_F=diag(I81,-I81), {Gamma_F,D_F}=0.

Squared operators separate exactly; the heat trace and curvature response
factor by the positive actual finite heat factor. No40/480-carrier conflation.
The topology, four dimensions and metric remain supplied, not emergent.
Producer/certificate: w33_20261001_wilson_gravity_refinement.py/.json.

## 5. Scale dynamics: a mechanism and its remaining inputs

With r=xd, Lambda_heat=alpha d, Lambda_EFT=beta d, V=d^4U(x),
a nonzero stationary point requires Uprime=U=0 and leaves d flat.
This familiar homogeneous scale obstruction is not claimed as new physics.

A reduced anomalous potential fixes both radial directions:

    V=lambda(r²-v²d²)²+B d^4[log(d²/mu²)-1/2],
    rho=sqrt(2)r is the canonical radial coordinate.

At r=vmu,d=mu the exact Hessian has leading minor4lambda v²mu² and
determinant32lambda Bv²mu^4, positive for positive parameters.
Real singlet spectators of mass hd can supply a positive logarithmic term;
a counterterm/renormalization condition supplies the displayed constant.
The full matter contribution must determine B_total. Positive invariant
selector terms divided by powers of d provide tree angular stiffness;
a stated norm bound gives local robustness, not global quantum selection.

The named scales now share a field, but alpha,beta,v, couplings and the
transmutation datum remain free. Vacuum energy=-Bmu^4/2, not a small
cosmological constant. The leading volume and two-derivative terms of a single positive heat action give
V_volume/M_Pl^4=144pi²/(kappa Theta_F); its common threshold factor alone
does not explain a tiny ratio while retaining the same Newton coefficient.
Producer/certificate: w33_20261001_anomalous_scale_completion.py/.json.

## External checks used

- [Vinberg theta groups](https://m.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=2123&wshow=paper): restriction/classification input.
- [Luty–Taylor](https://arxiv.org/abs/hep-th/9506098): moment-zero/closed-orbit interpretation.
- [Martin section3](https://arxiv.org/html/hep-ph/0111209v2): Landau/MSbar field weights and subtraction constants.
- [Spectral action](https://arxiv.org/abs/hep-th/9606001) and [scale invariance](https://arxiv.org/abs/hep-th/0512169): curvature/finite-factor and dilaton context.
- [Ottaviani, classical Aronhold invariant](https://arxiv.org/abs/0712.2527): classical invariant context, not a new-invariant claim.
- [Coleman–Weinberg](https://doi.org/10.1103/PhysRevD.7.1888): anomalous reduced mechanism, with transmutation scale supplied.

Prior owners: Pass11260/61 Cartan/classification, Pass11269 mirrors/designs,
the dated global-selector and D-flat mass packets, Pass4057/4169 Wilson
controls, BT1033/1130 geometry and product bookkeeping.

## Validation and intake

Five focused regression functions pass, including deliberate sign-corruption
rejection, rational gauge identity replay, independent finite-epsilon matrix
control and brute-force/compressed lattice comparison. The updated forty-points
paper builds with Tectonic; typesetting warnings remain. The original five-file intake
audit is clean; the added global-circuit packet is audited at release. Its alpha@64 candidate was checked by reading all three prior
files in full: E8 real-form root signs, PGSp outer characters and quadrangle
coclique searches. These do not assert the scale model here; alpha denotes a
free cutoff ratio, not a claimed code parameter.

At t=.2 the Wilson N192/384 Richardson response differs from the supplied
continuum benchmark by0.0406276 (absolute); finite-matrix replay error is
3.98e-10. The zero-Wilson finer extrapolation gives16.028049 continuum species,
with the wrong-sign coarse response explicitly retained. These are numerical
controls, not directed error intervals.

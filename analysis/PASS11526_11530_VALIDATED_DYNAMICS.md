# Passes 11526–11530: certified dynamics and a change of physical vacuum route

Reservation: `057834b58`. Producer: `analysis/w33_pass11526_11530_validated_dynamics.py`.
Certificate: `data/w33_pass11526_11530_validated_dynamics.json`.
Independent replays: `tests/test_w33_pass11526_11530_validated_dynamics.py`.

The main result is a physical screening result, followed by a separately supplied compatible scalar action. The old CP root is rigorously real, but its Spin8 stabilizer cannot host commuting color SU3 and weak SU2. The replacement action has an exact global minimum, a Spin10-compatible first stage, and a complete rational 324-field normal certificate. Neither phase derives the Standard Model, its measured parameters, or gravity. Five requested targets were pursued, plus symmetry screening, the coherent parent, and an explicit constrained spin history.

## 11526: an exact CP stationary point, and why it is not an SM vacuum

The unchanged 11438 native action, reduced in 11521 to its regular 20-real-coordinate quotient, is evaluated using second-order jets over outward-rounded FLINT complex balls at 256-bit precision. Coefficients and the selected preconditioner are exact rational inputs; floating inversion only proposes the preconditioner. A cube of radius 10^-25 around the stored rational center satisfies strict Krawczyk inclusion. The infinity norm of I-RH is about 1.244e-14 and the inclusion displacement is less than one hundredth of the radius. Interval Schur elimination proves all twenty Hessian pivots positive throughout the box. Thus this cube contains a unique exact stationary point, a strict quotient local minimum.

Its energy is enclosed near -2.93574323966686636886, CP-odd invariant chi near 0.70479596409305084138101, and commutator invariant q near 0.000721478471303784195649. These are dimensionless action quantities, not physical masses. The existing simultaneous hollowisation theorem (11481) and 11521 stationary-fiber transport give an exact stationary representative in the full 324-field inventory. Full transverse stability of this unchanged action is not asserted.

Before attempting to turn its normal modes into particle masses, we screened its surviving gauge symmetry. Exact native generators give a 30-dimensional frame normalizer, consisting of Spin8 and a two-dimensional diagonal torus. Distinct singular values of A force any combined E6/family compensation to be diagonal; the two certified nonzero off-diagonal entries of B force that diagonal to vanish. Hence the combined compact stabilizer has Lie algebra Spin8. The only faithful SU3 representations on a real eight-dimensional space are the realified fundamental plus two singlets, and the adjoint. Their orthogonal centralizers are respectively u1+so2 and zero. Neither contains a commuting su2. This rules out an SM embedding in this stabilizer, including embeddings different from the canonical one, under the stated compact E6 x SU3-family gauge premise.

An independent exact obstruction uses the canonical SM singlets already identified in 11271. Their native cubic is zero, as are the relevant composite invariants and first jets. The unchanged action there is the strictly positive homogeneous quartic

    V = ||mu(Phi)||^2/2 + ||mu(Psi)||^2/2
        + rho (||Phi||^2+||Psi||^2)^2.

Euler's identity x.grad(V)=4V excludes every nonzero stationary point on that SM-preserving slice. This is an obstruction to this action and compact model, not to W33, all E6 models, or a TOE. The mathematical CP certificate remains valid.

## Additional result: a coherent parent with all 324 modes certified

We explicitly change the supplied action. Set R=||Phi||^2 and use

    V_new = (R-1/2)^2 + (8R^2/9-||mu(Phi)||^2) + ||Psi||^2.

The native Cartan trace-Gram certificate gives squared weight norm 2/9 for every E6 fundamental weight. Convexity of the moment polytope bounds the E6 contribution by 2R^2/9; the SU3-family contribution is at most 2R^2/3. Therefore ||mu||^2 <= 8R^2/9, with equality on the highest-weight coherent orbit. Every term above is nonnegative. Phi=e0(native) tensor e0(family)/sqrt(2), Psi=0 is an exact global minimum, with a Spin10 E6 stabilizer compatible with a first stage toward the SM.

The full real-coordinate Hessian has this exact spectrum:

| Eigenvalue | Multiplicity |
|---|---:|
| 0 | 37 |
| 2/3 | 20 |
| 2 | 162 |
| 7/3 | 64 |
| 8/3 | 40 |
| 4 | 1 |

The exact tangent matrix has rank 37 and its image is precisely the Hessian kernel. All 287 normal modes are positive. The certificate stores rational matrices, not just an eigenvalue list; tests reconstruct the Hessian and its gauge kernel independently.

Coefficients, units and this selector are supplied. Further breaking to the SM, chiral matter dynamics, CP/family selection and observed masses remain open. Setting its minimum energy to zero by the displayed additive constant is not a cosmological-constant solution. The coherent-state extremality mechanism is classical; the result here is its native application and full normal certificate, not a new general coherent-state theorem.

## Additional exact bridge: a condensate-produced spinor selector

The same native weight certificate gives the normalized E6 moment operator M=K_Phi/R_Phi eigenvalues 2/9 (one), 1/18 (sixteen), and -1/9 (ten) on the coherent orbit. Here K_Phi=sum_ab (G^-1)_ab tr(Phi^dag T_a Phi) T_b, with G_ab=tr(T_a T_b) for the native E6 Hermitian generators; the family moment is excluded from this operator. Consequently the actual spectral projectors are

    Pi1  = 18(M-I/18)(M+I/9),
    Pi16 = -36(M-2I/9)(M+I/9),
    Pi10 = 18(M-2I/9)(M-I/18).

They are orthogonal idempotents of ranks 1,16,10 and transform by conjugation. This is standard Spin10 branching, instantiated as a field-produced operator rather than a supplied fixed projector. Projector identities hold on the coherent orbit; they are not asserted for arbitrary Phi.

It also names a further, separately supplied positive polynomial action. Write D_X=8R_X^2/9-||mu(X)||^2, and replace the first-stage Psi mass by

    V2 = (R_Phi-1/2)^2 + D_Phi + (R_Psi-1/2)^2 + D_Psi
         + ||(K_Phi-R_Phi I/18)Psi||^2
         + R_Phi R_Psi-tr(rhoF_Phi rhoF_Psi),
    rhoF_X=X^dag X.

Every term is nonnegative; the last term follows from the PSD density inequality tr(A B)<=tr(A)tr(B) and aligns the rank-one family factors. Phi=e0 tensor e0/sqrt(2), Psi=e1 tensor e0/sqrt(2) makes every term zero: the native e1 belongs to the moment's sixteen-dimensional sector and is itself coherent. Thus this action has an exact global minimum on an SM-preserving two-vector witness. Direct native-generator rank calculation gives a 24-dimensional E6 stabilizer for the two fixed vectors, consistent with the SU5 interface already owned by 11271. Further adjoint breaking to the SM is not new or completed here. This connects the compatible first stage to a concrete second stage without reusing the incompatible Spin8 CP pair.

This strengthened action has its **own** rational full 324-field Hessian certificate. It splits into 292 singleton blocks and sixteen two-dimensional blocks. The exact spectrum is:

| Eigenvalue | Multiplicity |
|---|---:|
| 0 | 58 |
| 1/36 | 12 |
| 2/3 | 30 |
| 49/72 | 10 |
| 2 | 4 |
| 10/3 | 104 |
| 241/72 | 24 |
| 11/3 | 60 |
| 265/72 | 20 |
| 4 | 2 |

The actual combined compact gauge tangent has rank 58 and spans the entire kernel; all 266 normal modes are positive. The certificate stores the moment-selector Jacobian, separate coherent-field Hessian diagonals, every connected rational block and gauge tangent. An independent regression differentiates the actual nonlinear potential along a mixed complex tangent, including off-family entangled components; these matter because an exploratory alignment Hessian initially omitted them. That incomplete calculation was corrected before publication.

The selector is degree six, with a dimensionful EFT coefficient. Its coefficients and global uniqueness remain unselected/unproved; the alignment fixes common family orientation but supplies no CP-breaking source or observed flavor. The first-stage 37-zero/287-positive certificate is not transferred: the 58-zero/266-positive result is calculated anew. The eighteenth regression checks the exact projectors, zero energy of the witness, and actual native stabilizer dimension; the nineteenth checks the second-action normals and a direct potential derivative. Existing polynomial SM selectors (11275 and 11293) are prior architectures; this construction does not claim the first selected SM Higgs orbit.

## 11527: a positive EW action and a full-rank mixed composite

For two doublets define S=||Hu||^2+||Hd||^2, D=||Hu||^2-||Hd||^2 and p=Hu^T epsilon Hd. A separately supplied positive portal is

    V_EW = lambda(S-2k sqrt(1+chi^2))^2 + kappa D^2
           + eta |Hu^dag Hd|^2 + nu |p-k(1+i chi)|^2.

All five coefficients are positive. The identity |p|^2=(S^2-D^2)/4-|Hu^dag Hd|^2 gives an exact neutral minimum with D=0, Hu^dag Hd=0 and p=k(1+i chi). Its energy is identically zero for every chi, so this engineered portal has no vacuum-energy backreaction on that coordinate. With t=k sqrt(1+chi^2), its eight real-coordinate Hessian eigenvalues are zero (three), 4 eta t (two), 16 kappa t, 4t(4lambda+nu), and 4nu t. A generic uncompensated portal would shift chi.

A single-copy minimal singlet/doublet ansatz produces rank-at-most-two family sources. With both existing scalar copies the typed covariant

    H^c_ij = dbar^{cab}(Phi*_ai Psi*_bj+Phi*_aj Psi*_bi)/(2M)

escapes that restriction. Explicit native-tensor witnesses yield symmetric full-rank Hu and Hd, determinants 1/4 and 2chi^2, nondegenerate Gram spectra, and

    Im tr[Hu Hu^dag, Hd Hd^dag]^3 = -69 chi(chi^2-19)/2.

Thus this operator can carry three-family rank and CP violation. Source conjugation reverses the CP branch and is checked explicitly. The family vectors and CP-dependent source are supplied, not dynamically selected. In particular, the incompatible Spin8 CP root cannot be transplanted into an SM vacuum as a physical CP source. The new coherent parent itself has chi=0. Earlier 11326 already owns an engineered noncommuting family vacuum.

A healthy auxiliary realization uses dimension-two Jmix=M Hmix and dimension-one X in (27,bar6): M^2||X-Jmix/M||^2. Integrating X out gives F F Jmix/M. This introduces a heavy sextet and matching coefficients; it is not a UV derivation of them. The healthy auxiliary completing-square architecture is already owned by 11293 (`w33_pass11293_yukawa_uv_budget.py`); its earlier source is v Fdagger, whereas the source here is the mixed native bilinear. Neither the general mediator method nor supplied family CP capacity is claimed as new.

## 11528: complete local first jets and a constrained spin history

On the supplied 80-site harmonic periodic cover, incident lifted displacement ranks are 52 rank-one and 28 rank-three sites. Nonbacktracking paths through radius two leave twelve rank-one sites; radius three gives rank three at all eighty sites. For each site the exact positive Gram matrix G_i=sum_p d_p d_p^T supplies

    nabla_i psi = sum_p G_i^-1 d_p (U_p psi_endpoint-psi_i).

This reproduces affine lifted first jets exactly. Rational changes of spatial basis preserve the reconstruction; the selected harmonic embedding itself remains an input. Path transports transform by endpoint Spin rotations. This is a finite-radius kinematic construction, not an established continuum chiral Dirac operator.

A separate 320-state Hermitian site-spin Hamiltonian has explicit Spin3 links, nontrivial cycle holonomy, symmetrized hopping and a supplied Wilson term. Rebuilding all couplings after local Spin rotations verifies covariance. Its one-hop principal symbol still has the deficiency above; the radius-three reconstruction is a repair target, not silently substituted as a completed operator.

The supplied homogeneous metric action has lapse constraint

    H = -p_a^2/(12a) + V(a) + psi^dag (K/a+m beta) psi = 0.

An expanding positive-energy packet evolves for 10^-4 time units with a increasing from 1 to approximately 1.00045227480427, constraint residual below 1.74e-12 and spinor norm error near 10^-14. This is a commuting-spinor mean-field minisuperspace history with an independent Spin3 connection. It does not supply fermion statistics, local Einstein constraints, a Levi-Civita identification or quantum gravity. The general DeWitt mechanism is already owned by 11482.

## 11529: a different nonlinear primitive, with locality kept qualified

11509 already constructs an actual nonlinear gauge-invariant finite Green primitive; see `exponential_ward()` in `w33_pass11506_11510_five_physics_targets.py`. We credit that owner rather than calling nonlinearity new. Here a rooted shortest-path incidence right inverse T gives B^T T=I-e_root 1^T and k(A)=T q(A), using the actual charge-weighted weak-overlap projector density for charges (1,-4,2,-3,6) with multiplicities (6,3,3,2,1).

The L=3 all-four-direction localized patch has a genuine density peak 2.30e-8 and nonlinear half-amplitude defect 2.70e-8. Ward and gauge-covariance residuals are below 1.86e-13; the sampled Wilson gap exceeds 0.9955. For a bounded supported defect and uniformly gapped interpolation, resolvent decay of q gives an exponentially decaying routed current tail up to a polynomial prefactor. This is a source-tail statement. The globally rooted routing does not establish two-variable functional locality for arbitrary dense backgrounds. Flux, torons, measure curvature and integrability remain open. Earlier owners also include 11479, 11484 and 11524.

## 11530: coherent faults without twirling

11432 owns the complete three-round flagged extraction instrument, conditional clean followup, and all 1714 single Pauli/preparation/readout cases. We replay every stored history and correction, then append an ideal terminal recovery as a mathematical decoding map. Every residue becomes logical identity. At a fixed faulty port and fully recorded branch, all Pauli basis faults therefore act as scalar logical identity. Since the sixteen two-qubit Paulis span M4, arbitrary single-port CPTP faults, including coherent faults, are corrected under the ideal-other-operations assumptions. This is an application of standard error-span correction, not a new general QEC theorem.

For coherent unitary errors ||F_j-I||<=delta at at most 132 CNOT slots, the measurement-dilated circuit expansion bounds the multi-fault remainder by b=(1+delta)^132-1-132delta. Zero- and single-fault terms disappear under projection orthogonal to the encoded Bell state. Hence

    1-F_e <= min(1,b^2).

At delta=10^-5 the exact rational envelope is approximately 7.4818137152e-13. The infidelity bound is O(delta^4); the diamond distance is not claimed to have that order. The all-port bound covers the specified unitary CNOT errors; single preparation/readout error-span results are separate. Native leakage, correlated noise, noisy final decoding and a physical threshold remain open.

## Checks and literature

All nine producer sections PASS. All nineteen independent regressions passed in 89.55 seconds, including producer/input binding and a direct nonlinear second-action potential derivative. After the final source-comment clarification, a fresh producer/input binding replay passed in 50.88 seconds. Final batch intake reports no rediscovery collisions, no forced-arithmetic findings and intake clean for all four science files. Every rooted route is independently checked to have the shortest periodic Manhattan length. A preliminary decrement-only route was rejected before publication because it could spoil source-tail locality. Reproduction uses Python, NumPy, SciPy, SymPy, mpmath and python-flint 0.9.0:

    OPENBLAS_NUM_THREADS=1 python analysis/w33_pass11526_11530_validated_dynamics.py
    PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 OPENBLAS_NUM_THREADS=1 python -m pytest --noconftest -q tests/test_w33_pass11526_11530_validated_dynamics.py

The certificate binds producer bytes and canonical source JSON hashes. No changes to the physics scorecard of Forty Points are justified by these results. Remote 11511–11515 was integrated during this packet: its cubed-phase test counterexample reinforces keeping word-level computational evidence separate from general unitary claims.

Primary external checks used:

- [Krawczyk inclusion notes](https://ww2.ii.uj.edu.pl/~zgliczyn/cap07/krawczyk.pdf) and [python-flint documentation](https://python-flint.readthedocs.io/en/latest/) for classical inclusion/ball arithmetic.
- [Barnum et al., generalized entanglement](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.032308) for extremal coherent-state purity; [Slansky's representation tables](https://doi.org/10.1016/0370-1573(81)90092-2) for representation context.
- [Luscher's local cohomology lectures](https://luscher.web.cern.ch/luscher/talks/Andrejewski2.pdf): gauge-invariant divergence alone is insufficient for a local measure.
- [Knill–Laflamme error correction](https://arxiv.org/abs/quant-ph/9604034) and [Chao–Reichardt flag fault tolerance](https://arxiv.org/abs/1912.09549) for the established error-span and flag framework.

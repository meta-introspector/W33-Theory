# Five physical follow-ups: phase, charges, EFT loops, metric dynamics and ultraviolet scales

This packet executes the five follow-ups to c944a5060 and adds a separate
nonperturbative branch. Exact statements, numerical experiments, and conditional
physical models are distinguished. It does not solve the theory of everything.

## 1. Global degree18 and a phase selector

The signed E6 tensor maps every81-field Phi to the ternary cubic
F(a)=d(Phi a,Phi a,Phi a). Let H_F be its Hessian determinant, and let
T=dS_F[H_F], the directional derivative of the existing Aronhold scalar.
This classical circuit has103 coefficient terms. Exact restriction to Q gives

    T|Q = (48384u6^3-13824u6u12)/5 +1492992u18.
    I18 = (19683T-48384I6^3/5+13824I6I12/5)/1492992.

Hence I18(Qq/sqrt3)=u18(q). The global analytic differential passes complex
finite differences and complexified gauge transport. Its invariance follows
from the defining E6/SL3 tensors and the classical Hessian covariance; no
expanded81-variable monomial table or gauge fixing is needed.

Adding gamma|I18-tau|², tau=-r0^18/729 real, lifts the common phase of the
previous positive selector. The canonical phase mass squared is
324gamma r0^34/729². All exp(18i theta)=1 representatives share the invariant
triple and are in one finite little-Weyl/compact gauge orbit. These are not
18 physical vacua. Conjugation stays in that gauge orbit, so the real selector
does not generate physical CP violation. The coefficient tau is supplied.

Producer: `w33_20261001_degree18_phase_completion.py`.
Prior owners: global I6/I12 circuit, balanced T-ray CP audit, and Pass11269.
The Hessian/Aronhold construction is classical invariant theory, not a new
classical invariant discovery.

## 2. The three light modes are not charged families

The actual integer86-generator action at Q(1,2,3) has rank78 and an
8-dimensional exact stabilizer. Every stabilizer generator annihilates all
three columns of Q. The previous exact gauge/mirror identities extend this
stabilizer to regular T. Therefore every completed normal mass mode is a
singlet of the connected unbroken gauge group. A group of dimension8 cannot
contain unbroken SU3xSU2xU1, whose dimension is12.

Separately, the declared two(27,3) Weyl multiplets plus their two conjugates
have zero chiral representation index under every subgroup restriction.
Gauge-preserving mass pairings cannot leave a net three-generation chiral
index from this content. A CKM/PMNS matrix cannot be extracted from singular
values without charged-current generators and particle embeddings.

This audits the current model; it does not change the standard27 branching
in `w33_e6_27_standard_model.py`. That branching is useful prior art, but it
is not a charge assignment for the three normal modes here.
Producer: `w33_20261001_chiral_mass_assignment_audit.py`.

A different chiral candidate avoids conjugating every E6 family:
(27,3)+(1,10bar). Explicit symmetric-cube weight traces give A10=27 and
T10=15/2, so the family anomaly27 cancels with just ten E6-neutral spectator
components. Its total Weyl dimension is91; bE6=35, bSU3=-15/2. This is only
an anomaly repair, with spectator masses and symmetry breaking open. The
old antisymmetric cubic mass map requires two distinct species and is not
inherited by a single81. Two such species plus two spectators instead have
six27 generations and different beta coefficients.

## 3. Full declared completed-EFT determinant: failed stability hypothesis

The numerical calculation now includes all162 real scalar directions,86
vector directions, and the conjugate-paired fermions with all three EFT lifts.
It uses the global scalar potential and analytic invariant gradients,
central-difference holomorphic Hessians, exact inverse-trace-Gram generator
normalization, and the Landau/MSbar determinant. At T there are78 massive
vectors,83 positive real scalar modes,79 scalar zeros, and81 massive fermion
singular pairs. The light masses remain.001,.002,.003 at the supplied inputs.
Independent scalar-Hessian step replay differs by2.23e-16 in eigenvalues.

The positive one-loop curvature hypothesis **failed**. At angular steps.001
and.0005 the extrapolated tree-plus-loop eigenvalues are approximately

    -0.332978, -0.332272, -0.000129542, -0.0000971522.

The step-Hessian difference norm is.0129196. Nearby off-shell scalar
Goldstone directions acquire negative squared masses; the determinant uses
its real part log|m²| and records every diagnostic. This is not evidence of
physical tachyons by itself, and the unresummed infrared-sensitive curvature
is not a smooth vacuum-stability certificate. Neither the failed positivity
hypothesis nor its numerical signs are promoted to a global quantum theorem.
Goldstone resummation/physical observables and finite EFT matching
counterterms remain necessary.

Producer: `w33_20261001_completed_eft_one_loop.py`.
The common phase and new I18 phase term are excluded from this determinant.
No coefficient tuning was used to erase the failed stability experiment.

## 4. Varied metric and a directed infinite-sum heat certificate

The continuum torus heat response at t=1/5 is now enclosed by directed
40-digit interval arithmetic and an analytic Gaussian tail bound:

    -49.334231662265182681340176257403881358
    <= response <=
    -49.334231662265182677950680076383923603.

These decimals are displays; the certificate stores **exact rational binary
endpoints**, avoiding inward display-rounding claims. Enclosure width is
3.3895e-18; analytic tail upper bound is1.69475e-18. Exact integer shell
multiplicities cover every4-momentum in[-24,24]^4. Absolute summands are
bounded by shifted Gaussian moments and a union of four coordinate tails.
This certifies the continuum benchmark, not Wilson lattice convergence.

The supplied metric can now be varied in the named action
integral sqrtg[rho_vac-M_Pl²R/2+cR²] plus matter. For the exact-isovolume
conformal mode cos(kx), the epsilon² coefficient is

    V[-3M_Pl² k²/2 +18c k^4].

The Einstein term has the standard Euclidean conformal instability; multiplying
by the positive W33 mass-fiber heat factor cannot cure its sign. A positive
R² term can stabilize these modes on a finite fixed torus with a stated
size-dependent bound, but this is not the complete metric Hessian, quantum
measure or Lorentzian continuation. The flat Einstein equation also requires
total rho_vac+V_min=0. Scalar vacuum energy cannot be ignored once the metric
is varied. Dimension, topology and action coefficients remain supplied.
Producer: `w33_20261001_metric_dynamics_and_heat_bound.py`.

## 5. Actual ultraviolet inventory and scale limits

For four Weyl81 multiplets and one complex81 scalar, standard long-root
normalization gives

    bE6 =17, bSU3 =-59/2,
    beta(g)=-b*g^3/(16pi²).

The previous trace-normalized gauge match corresponds to equal standard
couplings at one scale. Their inverse squared couplings separate by
-93log(mu/mu0)/(16pi²). E6 is asymptotically free; family SU3 has a Landau
pole. The mirror-matched spectrum is not an RG-invariant coupling relation.
Even one Weyl81 with the same scalar makes family SU3 infrared free.
Singlet spectators cannot change this gauge coefficient.

The E6 strong scale is a candidate dynamical threshold, not a derived
condensate. The anomalous reduced potential needs running quartic coefficient
dL/dlogmu=2B, but that does not fix its integration constant. An independent
vacuum-energy counterterm changes gravity without changing nongravitational
stationarity. Absolute scale and the cosmological constant remain open.
Producer: `w33_20261001_uv_running_scale_audit.py`.

## Additional nonperturbative connection

An external check found a closely related **different theory**: Goh et al.
study E6 with three27 chiral multiplets, global family SU3 and small anomaly
mediation. Their condensate denominator has degree18; they discuss a vacuum
class with nonzero degree18 invariant and vanishing degree6/12 invariants.
These conditions match our Cartan T class algebraically after normalization.
Their assumptions do not establish our nonsupersymmetric EFT dynamics.

The repository tensor action independently gives E6 stabilizer dimension14
on the entire generic SIC wall q=(1,1,z) over Q(z); generic q=(1,2,3) and the
MUB-only control q=(0,1,2) retain dimension8. Two intersection controls give28.
If the condensate denominator is an invariant degree18 polynomial vanishing
on these enhanced loci, the9 order2 SIC mirrors force even multiplicity and
therefore force the squared SIC product, u18, up to overall normalization.
This supplies an explicit interface to the compiled global I18.

The resulting conditional canonical-Cartan action is

    W=Lambda^9*(-729I18)^(-1/3), on the branch W(T)>0,
    V=sum|dW/dq|²-18m ReW.

Its radial restriction is36Lambda^18/r^14-18mLambda^9/r^6, with stationary
r^8=14Lambda^9/(3m) and canonical radial mass squared648m²/7.
The exact six-direction canonical Hessian has eigenvalues
(162/49,162/49,486/49,486/49,648/7,486/7)m². Its real angular block is
(162/49)[[2,1],[1,2]], with the opposite offdiagonal sign in the imaginary
block: equal-weight scalar eigenmodes and a mass-squared ratio3. These are
scalar-slice modes, not CKM/PMNS angles. This is not a full11-modulus result, a Standard Model vacuum,
a proved nonsupersymmetric continuation, or a cosmological-constant solution.
The strong scale and AMSB mass remain inputs.
Producers: `w33_20261001_chiral_decuplet_and_condensate_bridge.py` and
`w33_20261001_condensate_cartan_potential.py`.

## Prior-art and external checks

Result searches checked the invariant degrees, exact contraction coefficients,
chiral index, decuplet/anomaly combination, beta values and conformal-factor
language against the result index, repo scripts, JSON certificates, paper and
visible site. The classical branching/anomaly work and earlier curvature
benchmark are cited above; the degree18 Cartan invariant itself is prior art.

Primary references:

- [Classical ternary cubic invariants, chapter9](https://people.dimai.unifi.it/ottaviani/tesi/tesidr_maurizio.pdf).
- [Martin, one-loop field determinants](https://arxiv.org/abs/hep-ph/0111209).
- [Martin, Goldstone resummation](https://arxiv.org/abs/1406.2355).
- [General gauge-theory RGEs](https://arxiv.org/abs/1809.06797).
- [E6 chiral dynamics and its stated near-SUSY regime](https://arxiv.org/html/2505.07931v1).
- [Euclidean conformal-factor problem](https://arxiv.org/abs/hep-th/0103186).

The focused regressions independently replay homogeneity, phase curvature,
light-group charges, scalar/fermion spectra, rational heat endpoints, gauge
inventory, decuplet weights and a condensate potential displacement.

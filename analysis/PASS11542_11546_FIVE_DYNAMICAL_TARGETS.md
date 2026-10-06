# Passes 11542–11546: five dynamical targets, with an exact two-body gate

Reservation: `fbaf36b43`. All five requested investigations were executed.
The principal result is an exact finite-time encoded square-root-of-SWAP
using two-body terms on the actual Pass11539 neutral pair code. This is a
conditional Hamiltonian construction, not a physical TOE or a hardware demonstration.

| Target | Result | Physical input still required |
| --- | --- | --- |
| Binding and local controls | Statistics-aware pair channel, exact scattering boundary solution, covariant logical Pauli operators | Binding interaction, scattering lengths, spin sector and coefficients |
| Two-body entanglement | Exact closed cyclic blocks and a finite revival implementing sqrtSWAP | Addressed pair gaps, constituent swaps and calibrated strengths |
| Invariant flavor | Positive-pin valuation theorem, exact weak-basis CP polynomial and an indefinite-pin counterexample | Texture selector, epsilon, Higgs matching and measured scales |
| Causal continuum | Explicit oriented integer cover with stable wave dynamics, finite support and quadratic dispersion convergence | Selection of the cover and dynamics; four-dimensional gravity |
| Vacuum energy | Native scalar loop coefficient and exact fixed-flux shift criterion for prior sequestering | Full physical supertrace, gravitational sector, flux matching and residual CC |

## 11542: binding and direct logical control

For identical left Weyl fields, epsilon_ab psi_Ia psi_Jb is symmetric in
internal indices. Therefore the native antisymmetric internal pair belongs to
the spinor-symmetric (1,0) channel, or needs additional species/orbital data.
A nonrelativistic S-wave spin-triplet with antisymmetric internal state obeys
Fermi statistics. Its three polarizations can be a spectator subsystem only
when controls act identically on them. No healthy relativistic spin-one UV
completion follows from writing an auxiliary pair field.

The zero-range radial boundary u'(0)/u(0)=-1/a gives
u=exp(-r/a), E_B=-1/(2 mu a^2), for a>0. The reduced mass and scattering
length are matching inputs. A Casimir label does not determine binding.

Use the existing cubic map B: conjugate(81) -> Lambda^2(81), B†B=10I,
from Pass11539 and its cited September23 graded E8 bracket owner. For
Phi=f e0f0, Psi=f e1f0, put X=[-B conjugate(Psi), B conjugate(Phi)].
Then X†X=10 f^2 I. The covariant mass-squared perturbation X C X†,
with C any Hermitian 2x2 matrix, restricts to 10 f^2 C on the normalized
logical code. After a supplied nonrelativistic reduction this supplies all
local Pauli directions without routing through the ancilla.

For a canonical dimension-one auxiliary pair field the logical interaction
has dimension4, while ancilla mixing has dimension5. For a microscopic
dimension-three fermion-pair operator these become dimensions8 and9.
These are operator counts, not proof of renormalizability. The family
Cartan diag(-2,1,1) assigns ancilla charge -4 and logical charge +2;
each available conjugate condensate supplies +2. Nonzero mixing on this
specific phase requires at least three condensate factors. Extra charged
backgrounds can change this lower bound.

## 11543: exact exchange synthesis

Take the two actual signed logical pair vectors stored in Pass11539,
each of squared norm10. Their antisymmetric wavefunction matrices have
squared Frobenius norm20 and support30 one-particle labels.
On registers (1,2) and (3,4), let

    W = S13 + S14 + S23 + S24,
    H = Delta(Q12 + Q34) + J W.

Each S is an ordinary constituent swap; each Q is the complement of its
two-dimensional pair code. Thus every summand acts on two constituents.
The gap projectors are supplied interactions, not interactions derived from
the native scalar action.

Exact contraction on the native tensors gives

    P W P = (I + SWAP_L)/10,
    P W^2 P = (19/5)(I + SWAP_L),
    P W Q W P = (189/50)(I + SWAP_L).

The last expression warns against using a bare exchange pulse as though it
preserved the code. More strongly, on every native code tensor V_j,

    (W^2 + 2W)V_j = 4(V_j + V_swap(j)),
    P12 W V_j = P34 W V_j = P W V_j.

These identities close the entire cyclic subspace, including its leakage
directions. Permutations cannot leave the30-label support. The singlet is
dark. Each logical triplet has the exact block

    [ J/5             3 sqrt(21) J/5 ]
    [ 3 sqrt(21) J/5  2 Delta-11 J/5 ]

Its center is Delta-J and Omega^2=Delta^2-(12/5)Delta J+9J^2.
Choose

    J/Delta = (649606 - 801 sqrt(337411))/25672045
            ≈ 0.007180121832944005,
    (Delta-J)/Omega = 801/800,
    t = 400 pi/Omega.

Then sin(Omega t)=0: all bright amplitude returns exactly to zero.
The triplet phase is -i and singlet phase1, hence the encoded gate is
exp(-i pi/4) exp(-i pi SWAP_L/4). This is sqrtSWAP up to global phase.
Together with the preceding declared local controls, it is a conditional
universal encoded gate set. It replaces the explicit four-body pulse of
the earlier architecture with an exact two-body model. It does not select
the physical gap operator, binding, exchange actuator or error tolerance;
generic calibration perturbations spoil the exact revival.

Exact closure also permits a shorter stronger-exchange protocol:
J/Delta=(-14+sqrt(571))/25 ≈0.39582425162788165,
Omega/Delta=2(1-J/Delta), t=3pi/Omega ≈7.799698999970913/Delta.
It yields the same triplet phase -i and singlet phase1 with exact return.
At pi/Omega it implements inverse sqrtSWAP. This does not use a weak-coupling
approximation, and no time-optimality claim is made. The ratio801/800 of
the earlier pulse is a convenient calibration choice, not a derived physical
constant; its unrelated occurrence in the prior W33 Green kernel does not
establish a connection between the two constructions.

Exchange-based encoded universality is established literature, not our
general invention: [exchange-only constructions](https://arxiv.org/abs/quant-ph/0309002).
The contribution here is the actual native code contraction, closed block
and algebraic revival certificate.

## 11544: coefficient-independent flavor scope

Pass11413 owns the supplied texture Y=D B D, D=diag(epsilon^3,epsilon^2,1),
and its mass/mixing orders. For fixed positive definite B,
lambda_min(B) D^2 <= DBD <= lambda_max(B) D^2. Min-max proves mass
valuations (6,4,0). Generic nonzero Schur off-diagonals give mixing
valuations (1,2,3), so val(theta_ij)=val(m_i/m_j)/2. This is a relation
between powers, not an equality between measured masses and angles.

With Bu=diag(1,2,3) and

    Bd = [ 3    1+i  1-i ]
         [ 1-i  4    1   ]
         [ 1+i  1    5   ],

both pins are positive definite. The certificate stores the exact
Im tr[Hu,Hd]^3 polynomial, whose leading term is369360 epsilon^22.
The classical [Jarlskog invariant](https://cds.cern.ch/record/161795)
has valuation6 from J and8 from each squared-mass Vandermonde.
CP-conserving pins can make it vanish. An invertible indefinite pin with
B01=B10=B22=1 and all other entries zero instead has singular-value
orders(5,5,0): invertibility alone is insufficient. No physical texture,
epsilon, endpoint or absolute scale has been selected.

## 11545: explicit causal refinement

Passes11540/11541 own the27-event finite history geometry. The integer
representatives of its8 null steps are (±1,±1,0), (±1,0,±1).
Choose an added cover Z^3 -> F3^3 and orient the four t=+1 steps.
The finite undirected graph does not select this cover or orientation.
Microscopic reach is a Manhattan diamond rather than a circular cone.

The supplied stencil delta_t^2 phi=c^2 Laplacian_xy phi, c^2=1/3,
has symbol

    -omega^2 + c^2(kx^2+ky^2)
    + a^2[omega^4-c^2(kx^4+ky^4)]/12 + O(a^4).

All Fourier frequencies are real for2c^2<1. The exactly conserved energy
E=||delta phi||^2/2+c^2<grad(phi_next),grad(phi)>/2 satisfies
E >= (1/2-c^2)||delta phi||^2+c^2||grad(average phi)||^2/2.
Four meshes verify quadratic dispersion convergence; a65x65 grid retains
finite support for20ticks with energy drift below1e-12. This is a conditional
fixed-background2+1 continuum example, not4D gravity, a selected metric,
spin2 dynamics or Einstein constraints.

## 11546: scalar loop and exact fixed-data shift protection

The prior11526–11530 two-stage Hessian has58 zero and266 positive modes.
Its sum of multiplicity times squared dimensionless eigenvalue is
2218769/864. With separately declared canonical scalar matching M_i^2=m^2 h_i,
the scalar-only Coleman-Weinberg log coefficient is
2218769 m^4/(55296 pi^2), with scale derivative
-2218769 m^4/(27648 pi^2). Gauge, ghost and fermion contributions and
kinetic matching are missing: this is not the complete physical supertrace.

Pass11274 owns global sequestering;11284 and11289 own local flux/history
extensions. Apply their affine-sigma sector to this specific scalar constant
shift: T -> T-Cg, Lambda -> Lambda-C leaves
T-g<Tr T>/4-DeltaLambda g invariant while retaining nonconstant sources.
This imports a known mechanism, not a new W33 gravitational solution.
See the primary [global proposal](https://arxiv.org/abs/1309.6562) and
[local formulation](https://arxiv.org/abs/1505.01492).

At fixed geometry and flux Q, the constraint
Volume=sigma'(Lambda/mu^4)Q/mu^4 is invariant under arbitrary constant
shift only if sigma''=0. Linear sigma passes exactly; quadratic sigma
changes the constraint by C Q/(mu^4)^2. This is an exact fixed-data Ward
criterion, not a general no-go for nonlinear sequestering when geometry
adjusts or approximate radiative stability suffices. Subtracting the
average in the entire matter potential instead makes that potential
identically zero and erases its forces; subtraction belongs in the
constrained gravitational source. Graviton loops and observed residual
vacuum energy remain open.

## Reproduction and ownership

Run `analysis/w33_pass11542_11546_five_dynamical_targets.py`, then
`pytest --noconftest tests/test_w33_pass11542_11546_five_dynamical_targets.py`.
All five producer sections PASS; all12 independent regressions PASS.
The JSON binds upstream certificates and code by SHA256 (canonical parsed
JSON, raw LF source). Tests independently rebuild the native pair map,
permutation closure, full block exponential, CP polynomial, wave symbol
and fixed-flux shift. Classical statistics, scattering, Jarlskog, lattice
waves and sequestering are credited rather than claimed as new theories.
Exact-result searches include189/25,801/800,369360 and2218769 alongside
the prior control, flavor and vacuum owners. Publication still requires
the refreshed corpus intake and index guard.

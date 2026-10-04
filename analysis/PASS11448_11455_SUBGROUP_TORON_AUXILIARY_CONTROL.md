# Passes11448–11455: identify the algebra, retain the determinant, expose relay faults

Reservation `ec4cee553`. [Producer](w33_pass11448_11455_subgroup_toron_auxiliary_control.py),
[certificate](../data/w33_pass11448_11455_subgroup_toron_auxiliary_control.json),
[tests](../tests/test_w33_pass11448_11455_subgroup_toron_auxiliary_control.py).

Five requested frontiers and three additional mechanisms are tested. Native
representation, action, overlap operator, neutrino determinant and routed CNOT
model retain their earlier owners in11438–11447. Generic branching, Schur
factorization, Berry transport and quantum instruments are established methods.
The physical theory is not closed by these finite constructions.

##11448: the preserved algebra is su3 inside E6

The eight compact generators close to `2.54e-15`, have zero center, and have
family-generator component norm `1.79e-15`. A compact center-free Lie algebra
of dimension8 is su3 by compact semisimple classification. Its adjoint Casimir
is the identity in the supplied normalization. On the actual81 representation,
the Casimir has27 zero eigenvalues and54 eigenvalues4/9. Factoring out the
three-family identity gives nine E6 singlets and eighteen non-singlet states
with fundamental/antifundamental Casimir. Casimir alone does not distinguish
3 from conjugate3; no weight assignment is inferred. The commutator-kernel
map on the78 E6 generators additionally measures the centralizer.

This sharpens11443's unnamed subgroup into a concrete colour-type su3 embedding.
It still requires an explicit intertwiner to the physical colour carrier and
charges; abstract branching is not a Standard Model identification. Trinification
branching is established prior literature and prior repo material, not claimed
as a newly discovered Lie-theoretic result.

##11449: mixed soft response

The preserved-gauge candidate's full normal Hessian is split into four soft
and242 massive directions. Twelve soft directions, including eight mixed
combinations, are tested at amplitudes.02 and.04. The quadratic massive
response is eliminated using the positive normal inverse. The resulting
energy/h^4 ranges from approximately `-.000391` to `.005870`. Signs occur at
energies near the finite-gradient and finite-difference error scale. These
samples neither prove stability nor establish exact moduli; no positivity
claim is made from unrelaxed or isolated axis scans.

##11450: actual closed toron transport

On a supplied flat L2 torus, vary the temporal toron through2pi and construct
the full Wilson-overlap projector at each step. Antiperiodic spin boundary
conditions remain explicit. Frames are parallel-transported by the polar part
of successive overlaps. The endpoint is closed with the actual large gauge
transformation, not by identifying distinct frames numerically. Charges1,2,3,4,6
are checked at32 and64 steps; the Wilson gap remains about1 and endpoint
projector covariance is at machine precision. Flat toron loop phases are
numerically trivial in this fixture. This is a global-loop control in one
flat sector, not the missing local current or nonzero-flux sector gluing.

##11451: a matter contribution to geometric equations

The Lorentzian scalar Schur action now varies with boundary vertex coordinates.
For a declared nonconstant boundary scalar, four geometric forces are computed
at two derivative steps. Massless refined and coarse forces agree, and a
translation test checks a second realization. Massive matter has an explicit
nonzero geometric response after its interior scalar is eliminated. This
provides a named term for a coupled equation; the Lorentzian gravitational
boost-angle action and causal branch prescription still have to be constructed.

##11452: auxiliary is not synonymous with loop-free

Use the actual80-node kernel `K=10I-.81A` with eigenvalues6.76–13.24, family
identity, and three native endpoint attachments. The block Gaussian identity
is checked directly:

`det(full)=det(K_heavy)*det(light-Schur)`.

Error is `9.24e-14`. Integrating out ordinary Grassmann fields retains the
heavy determinant. Its actual phi derivative is `-8.153992258`, checked against
finite differences. Dividing out detK changes the integration measure; it
requires a justified compensator or a different microscopic theory. Thus
renaming the972 elementary triplets as computational/auxiliary states does
not repair11446's beta-function obstruction. This is a finite native-kernel
control, not a computation of the full spacetime gauge determinant.

##11453: a deliberately supplied joint slice

A six-state neutrino determinant is evaluated at high precision while varying
a Majorana scale and relative pin phases. Add the declared potential
`100*(s^2-1)^2` and use the previous locally renormalized link curvature100
with portal `.1*(phi-.81)*(s-1)`. The slice gives
`s=1.001626636`, linear link response `phi=.8099983734`, radial curvature
`797.51`, radial residual about `-2.06e-6`, and pin angles near
`(3.06978039,3.08966600)`. Tiny phase differences are evaluated separately,
not added to a large double-precision heavy vacuum constant.

The pin matrices, source orientation and coefficients are fixed inputs. The
link response is a quadratic approximation, and portal feedback into the
Majorana scale is omitted at order portal²/100. This is an explicitly scoped
joint slice, not the complete nonlinear vacuum or a parameter prediction.

##11454: the complete instrument and the reset environment

Construct every Fourier-outcome Kraus coefficient from the actual native clock
spectrum at history dimensions32,128,1024. Completeness checks the entire
measurement instrument, rather than its zero-phase acceptance alone. The
logical dark space is accepted without disturbance in the ideal instrument.
Controlled clock and readout remain primitives requiring synthesis.

Preserve the two-dimensional logical sector and reset all ten orthogonal
leakage states to one logical state. The ten leakage environment states must
be mutually orthogonal and orthogonal to the retained logical environment
state. Therefore the minimal environment dimension is **eleven**. An actual
Pkeep plus ten reset Kraus operators has independent rank11 and is TP.
Closed unitary control of the cell cannot perform this dissipative reset.

One additional twelve-state native cell has sufficient environment dimension.
The producer completes the reset isometry to a joint-cell unitary permutation,
then derives every environment-outcome Kraus operator. The logical sector is
preserved and every leakage basis input is reset; discarding the environment
produces the stated channel. This is an actual virtual-hardware target in the
adapted logical/leakage basis. The existing logical CNOT does not supply its
leakage-sector interactions or its pulse sequence.

##11455: relay reuse is a physical part of the decoder

For one SWAP-out/CNOT/SWAP-back route, enumerate all105 nonidentity Pauli faults
on its seven two-qubit gates and propagate them through the remaining circuit.
Seventy-two leave a relay error; eighteen reach all three sites. The circuit
still has the correct ideal remote CNOT, but a clean initial relay is not
a clean final relay after a fault. Resetting or tracking relays must enter the
complete schedule. The288/372 resource optimum remains valid for the supplied
ideal routing model; these counterexamples prohibit inferring a noisy threshold
from that count.

## External checks and boundaries

- [Kikukawa, domain-wall chiral gauge construction](https://arxiv.org/abs/hep-lat/0105032) gives an existing5D measure/current route; our toron test does not implement it.
- [Bajc and Susič, realistic renormalizable E6 model](https://arxiv.org/abs/1310.4422) provides established E6 model-building context; our subgroup certificate does not derive its spectrum or vacuum.
- [Aliferis and Terhal, leakage faults](https://arxiv.org/abs/quant-ph/0511065) motivates explicit reset/environment operations beyond abstract erasure correction.

Searches covered the result index, paper, index.html, prior reports, scripts
and certificates, then result phrases for toron holonomy, determinant removal,
Casimir branching and environment dimension. No exhaustive corpus reading or
physical completion is asserted. Source certificate is canonically hash-bound.

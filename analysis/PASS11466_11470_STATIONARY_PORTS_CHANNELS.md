# Passes11466–11470: a limiting vacuum stratum, coupled Lorentz dynamics and minimal ports

Reservation90957fb96. [Producer](w33_pass11466_11470_stationary_ports_channels.py),
[certificate](../data/w33_pass11466_11470_stationary_ports_channels.json),
[regressions](../tests/test_w33_pass11466_11470_stationary_ports_channels.py).

All five requested frontiers are pursued with actual maps or channels. They
remain physical research problems. These controls do not establish a TOE.
Prior owners11461–11465,11428,11443 and10956 are consumed rather than renamed.
The native paper and site already contain D4/Spin8, graph spectra, trinification
and Lorentzian matter; those topics and the imported methods are not new.
Searches for the exact port dimensions, discarded determinant, and relay-state
independence found no earlier implementation of the constructions below.

## 11466: the numerical limiting stratum is D4, rather than the cutoff-dependent count

Start from11461's photon-preserving324-real-coordinate candidate. Take the
28-dimensional kernel separated from its resolved orbit singular values and
construct its actual native generators. Their brackets close at roundoff.
Their quadratic Casimir has nine fixed complex directions in81, corresponding
to three in27. Project the fields to this space and refine its36-real-coordinate
stationary equations, discarding unresolved Hessian directions from the Newton
inverse. This avoids assigning enormous motion to noise-scale eigenvalues.

The adjoint representation has a four-dimensional regular Cartan centralizer,
24 nonzero roots, and equal root lengths. Compact rank-four simply-laced
classification identifies the constructed numerical algebra as D4. Generator
matrices and root coordinates are stored, with a second-basis regression.
This is more information than the number28. Pass10956 already owns an exact
Spin8 on another Albert carrier; an intertwiner to that exact realization is
still required. Global group covering and exact stationary symmetry are not
inferred from floating matrices.

The refined orbit has rank58 at1e−8,1e−9 and1e−10. Its normal Hessian therefore
has266 directions:262 resolved positive directions and four unresolved modes.
This replaces the earlier *interpretation* of rank71 at1e−8 with a constructed
limiting numerical stratum; all older stored cutoff counts remain unchanged.
Full stationarity, finite-difference gradient and lowest-mode Hessian replay
are checked. No interval certificate, exact vacuum selection, relaxed quartic
stability or observed Higgs spectrum is asserted. Additional unbroken gauge
fields also prevent treating this stratum as an established SM vacuum.

A further soft-mode experiment removes the ten gauge directions inside the
36-real-coordinate fixed space, leaving26 normal directions: four soft and22
massive. Relax the latter at signed amplitudes.06 and.12 in four axis and two
mixed soft directions. An unprotected frozen-Hessian iteration can diverge;
the published construction instead uses scaled L-BFGS with a line search and
records all24 energies, coordinates and actual gradient residuals. The sampled
energy differences approach numerical precision after relaxation. This is
evidence against inferring restoring quartic coefficients from unrelaxed scans,
not a proof of exactly flat moduli or uniform stability. Each soft constraint,
energy and projected residual is independently replayed.


## 11467: all native charges fit an admissible flux sector, but sectors cannot be traversed

The actual integer6Y charges from11461 are applied to the uniformL34 links.
Each charged plaquette meets the1/30 bound, and its magnetic flux equals its
integer charge. A random gauge transformation independently preserves these
results. Prior11443 owns admissible charged links; this packet connects them
to11461's actual branching and adds an explicit sector-transition control.

The previous proposed continuous flux-sector transition is not a valid target
inside the admissible field space. [Lüscher, Lemma7.3](https://arxiv.org/html/hep-lat/9811032)
already proves those sectors are disconnected. An explicit continuous path
from trivial to unit-flux links crosses inadmissible plaquettes. The task is to
construct a local current within each sector and choose relative sector phases,
not to glue sectors with an admissible path. A toron Berry connection is not
that spatially local current. The full chiral measure remains open.

## 11468: simultaneous interior Lorentzian gravity and scalar equations

Branch beyond the single embedded simplex. Implement the regular homogeneous
two-frustum gravity–scalar action of
[Jercher and Steinhaus](https://arxiv.org/abs/2312.11639), equations19,31 and44.
The real action contains hyperbolic gravitational angles and timelike trapezoid
rotation deficits. Solve both lapse equations, the interior scale equation and
the interior scalar equation simultaneously. Enforce the causal-domain bounds
throughout the search. Independent Schläfli lapse equations and finite
variations check the solution. The independent height check differentiates
the area in equation5 directly: its derivative is
`H(a+b)/(2 sqrt(H²−(a−b)²/2))`. The displayed first derivative in equation51
omits this factor2; the action, direct derivative and finite differences agree.

For suppliedG=1, matter coupling1 andLambda0, boundary spatial scales1 and1.2,
and first height.3, the compatible witness has interior scale about1.094871,
second height.395907, interior scalar.045400 and final scalar.090974.
All four forces are below1e−10; the scalar current is conserved and both
trapezoid deficits are nonzero. The final scalar is solved as boundary
compatibility data and then held fixed for every stationary variation. The
first height parametrizes the search but its equation is also enforced.

This is an executable solution of an imported restricted model. It is not a
new gravity theory, local inhomogeneous Regge solution, W33-derived spacetime,
or prediction ofG or the cosmological constant. Its value here is to replace
boundary-only forces with simultaneous interior dynamics on a nonflat model.

## 11469: exact compression preserves flavour transfer, but discarded fermions still count

On11428's actual80-node Levi adjacencyA, verify the integer polynomial
A(A²−6I)(A²−16I)=0. For source0 and flavour endpoints44,1,0, build the integer
Krylov columns[A^k e0,A^k ej],k=0…4. Exact rational restrictions have dimensions
8,8,5. An orthonormal realization is used only for the singular-value replay;
the invariant integer bases and rational matrices are the certificates.

Each heavy leg occupies21 coordinates rather than240. Including the central
three-state bridge and three light states gives48 per quark sector instead of
486. At three supplied link moduli, the full Schur flavour matrix matches the
reduced one. All full singular values are recovered by adjoining the discarded
spectra, so this is an exact invariant decomposition, beyond truncated moment
matching. [Krylov transfer reduction](https://arxiv.org/abs/2102.11915) is prior art.
Controllability and symmetric observability establish the minimal all-port
realization separately for each family leg; they do not establish a minimal
microscopic physical theory.

Across the two sectors,876 discarded Dirac states contribute

`det_dark = 10^340 (100 − 6 phi²)^268`.

Its logarithmic derivative is nonzero and explicitly checked. These directions
are invisible to the chosen flavour ports but remain visible to the scalar
background and gauge field. An exact basis change retains all972 species and
b0=−637. A deliberately different theory containing only the96 retained fermions
would haveb0=−53: substantially smaller but still not asymptotically free. Removing
the dark determinant changes the model. No UV completion or new observed mass
prediction follows. The numerical Coleman–Weinberg force is not identified with
the zero-dimensional logdet derivative.

## 11470: complete relay channels separate Pauli independence from dissipative memory

Each ideal seven-CNOT route equals an endpoint CNOT tensor relay identity.
Every complete Pauli fault trajectory, with any number of faults, can therefore
be written as an endpoint Pauli tensor a relay Pauli after this ideal unitary.
Partial trace removes the latter even for initially entangled data and relay.
Any *state-independent* classical distribution over trajectories preserves this
property, including classically correlated fault probabilities. Thus endpoint
output depends only on the initial data marginal, not relay residues.

Compute normalized endpoint Choi matrices with noise after **every** CNOT,
for relay0,1 and+. Stochastic two-qubit Pauli noise produces identical endpoint
channels; amplitude damping on the participating wires produces different ones.
All are completely positive and trace preserving. This gives a whole-channel
counterexample to extending the Pauli conclusion to general noise, rather than
counting isolated fault residues. The proof applies after trajectory probabilities
are specified; it does not assert independence when physical probabilities
depend on a hidden quantum bath or relay state.
[Randomized compiling](https://arxiv.org/abs/1512.01098) provides a relevant
external route toward Pauli noise under its assumptions, not a guarantee that
this hardware realizes the required channel. Native leakage reset controls,
full decoded recovery and a threshold remain unbuilt.

The intake guard's broad background matches were read in full. In particular,
[PASS409_RESERVATION.md](PASS409_RESERVATION.md) owns earlier compact-control
and contact-parabolic algebra, while
[w33_pass348_every_multiplicity_is_a_torsor.py](w33_pass348_every_multiplicity_is_a_torsor.py)
owns the named-object guard history and symmetry-selection boundary.
Those earlier results retain ownership; they do not supply the present
stationary-field, port-restriction or whole-channel witnesses.

The remaining background owners are
[2026-09-23_diagonal_weld_e8_lie_generation.md](2026-09-23_diagonal_weld_e8_lie_generation.md),
[BT889_fermion_content_of_27.md](BT889_fermion_content_of_27.md), and
[docs/index.html](../docs/index.html).
The weld and finite-shell decompositions are earlier results with different
constructed objects; the present five investigations do not supersede them.

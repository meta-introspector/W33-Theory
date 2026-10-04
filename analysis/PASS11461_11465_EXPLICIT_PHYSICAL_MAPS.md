# Passes11461–11465: native charges, a photon test and five explicit maps

Reservation `3e9cc0867`. [Producer](w33_pass11461_11465_explicit_physical_maps.py),
[certificate](../data/w33_pass11461_11465_explicit_physical_maps.json),
[independent tests](../tests/test_w33_pass11461_11465_explicit_physical_maps.py).

The five requested continuations are implemented as finite constructions and
negative controls. They do not close the five physical problems or the TOE.
Remote42aef6d5d changes only the formula inventory. Scientific prior owners
are11438–11455; standard E6 branching, determinant-line geometry, Regge calculus,
Wilsonian running and Stinespring control are established methods.

##11461: an explicit native electroweak map, followed by a vacuum obstruction

The actual colour-type su3 from11448 has a16-dimensional centralizer inside
E6. Construct its matrices, form the adjoint representation, and solve its
commutant. The two-dimensional centroid separates two commuting8-dimensional
ideals. Their Cartan subalgebras support an explicit chosen hypercharge
assignment. The certificate stores the colour, left, right, weak, hypercharge
and electric-charge matrices in the same native27 basis as the current action.

The full branching, using colour/weak activity and6Y, is

| colour | weak |6Y| states |
|---|---|---:|---:|
| active | doublet |1|6|
| active | singlet |−4|3|
| active | singlet |−2|3|
| active | singlet |2|6|
| singlet | doublet |−3|4|
| singlet | doublet |3|2|
| singlet | singlet |6|1|
| singlet | singlet |0|2|

This is the standard one-family plus vectorlike coloured/weak pairs and neutral
singlets pattern. The quadratic Casimir identifies colour activity; an additional cubic colour
Casimir distinguishes3 and conjugate3 in a chosen orientation. It has nine
eigenvalues each−1,0,+1. TheY=1/6 doublet andY=−1/3 singlet are3;
Y=−2/3 andY=+1/3 singlets are conjugate3. Weyl chirality remains an input. Gravity/Y,
Y-cubed, colour-squared/Y and weak-squared/Y traces vanish numerically. Choosing
which centralizer ideal is weak and its weight convention remains an input.
Standard trinification is **not** a new discovery. Prior11271 already supplies
a separate canonical SM-preserving signed27 Higgs configuration and an imposed
selecting potential; it is a different field configuration and action.

The current preserved candidate fails the photon test. Its four weak/Y gauge
Gram eigenvalues are approximately `.315122,.621147,.621147,1.078617`; the chosen
Qem has Gram norm-squared `1.364127029`. Thus the current candidate preserves
colour but breaks electric charge. This makes the older eight-generator
insufficiency explicit; it does not supersede11271's different Higgs model.

A second search constrains **this same coercive action** to the colour/Qem
fixed space:15 complex coordinates per field. A full normal Hessian and
smaller-step lowest-direction replay distinguish stationary existence from
transverse stability. The search gives energy−2.935743240 and gradient norm6.24e−7, with photon
Gram at roundoff. At the recorded1e−8 cutoff the numerical little algebra has dimension15,
center zero and negligible family projection: **su4 inside E6 at that cutoff**. This follows from compact classification
and its contained su3, not the dimension alone. The stored constrained fields preserve a colour/charge subgroup, but the
limiting stationary little group needs the rank-sensitivity check below.
The full253-dimensional normal Hessian has249 resolved positive eigenvalues
and four unresolved near-zero values, with no resolved negative eigenvalue.
A smaller-step replay and finite-amplitude lowest-direction scan are included.
These do not certify relaxed quartic stability or a global minimum. The branch
is higher in energy than the current colour-only candidate; it is a concrete
constrained stationary candidate of the existing supplied action, not physical
vacuum selection. Its energy matches11438 within numerical tolerance; a newly
discovered vacuum energy or a distinct connected component is not claimed.

**Rank-sensitivity audit:** at cutoffs1e−8,1e−7,1e−6,1e−5 the photon candidate
has ranks71,58,58,58; the colour candidate has78,65,58,58; the older candidate
has86,73,58,58. The complete singular values are stored. Thus the previously
recorded numerical counts reproduce at their cutoff, but do not certify the
exact limiting stationary stabilizer. Several apparent broken directions are
at solver-noise scale. Orbit rank58 suggests an enhanced28-dimensional little
algebra; its exact identity and stationary existence are **not** asserted.
The normal Hessian counts also depend on this orbit partition. A constructed
high-precision or symbolic stratum is required before interpreting those
numbers as a physical gauge spectrum.

##11462: actual transport in an admissible nonzero-flux sector

Use the4D lattice16x16x2x2 with unit magnetic flux in plane01, periodic spatial
and antiperiodic temporal spin boundary conditions. Its maximum plaquette
norm is below1/30. Translation invariance in spectator directions gives four
exact1024-dimensional Wilson blocks. No4D mode is discarded. An independent
direct4D assembly at smaller size verifies the projector block decomposition.

On a two-toron slice, calculate the positive Wilson projectors, determinant-line
curvature, polar-frame rectangle transport, and an explicit local Stiefel-chart
connection. The connection convention is `A=-i Tr E†dE`, so its curvature is
`-i Tr P[dP,dP]`. The opposite curvature convention is not mixed with the loop
phase. Two derivative steps and the actual rectangle phase check orientation.

This is **field-space chart locality**, not a proved spatially local
anomaly-cancelling fermion-measure current. Charge1 is tested; all six SM
integer charges require a larger admissible flux plane. One magnetic plane
has zero4D index. Global transition cocycles, double-flux sectors and the
Luscher current construction are not claimed implemented.

##11463: a Lorentzian gravity action coupled to the scalar force

Every triangle of the prior supplied simplex is timelike. Their normal
planes therefore have ordinary rotation angles. Use the upper Wick branch
`z=-1+i epsilon` and the explicitly named complex boundary action

`Sg=-i (sum_h sqrt(det G_h)/2*(pi-theta_h)/G - Lambda*sqrt(z)*Vcoordinate)`.

HereG=1 andLambda=.001 are supplied. The Lorentzian limit gives
`Sg=163.391856381`, while the prior scalar Schur action gives`.349962954`.
All20 coordinate derivatives of both terms are stored, with their sum.
The Schlaefli residual is below5e−10; a boost changes the action by about3e−9
at the finite regulator, and the translation-force identity is checked.

This supplies a gravity-plus-matter action and equations on a concrete causal
fixture. It does not solve stationary coupled geometry. The simplex has no
interior curvature; spacelike-hinge boost branches, causal-sector gluing,
physical Newton coupling, spacetime emergence and CC selection remain open.

##11464: nonlinear determinant backreaction and a finite-EFT cutoff boundary

Replace11453's quadratic link response with the actual full quark spectral
potential:486 singular values per sector, the prior local radial counterterms,
the actual six-state Majorana determinant, and the declared portal
`.1*(phi-.81)*(s-1)`. The joint two-variable stationary search gives roughly
`phi=.809994742`, `s=1.001626961`. Two derivative steps test the small forces;
the radial Hessian eigenvalues are about100.002 and797.511.

All972 triplet species remain in the elementary reading, so the high-energy
coefficient staysb0=−637. For the suppliedg(5)=.2, threshold running to cutoffs
20,100,1000 gives inverse bare coupling-squared approximately
19.3610,6.37655,−12.2000. A finite cutoff can have a positive bare kinetic term;
continuing past that boundary cannot. This is an explicit Wilsonian EFT
interpretation with retained determinant, **not a UV completion** or a
continuum limit. Coefficients, pin orientation and neutrino phase remain inputs;
orientation equations and the full field Hessian are still absent.

##11465: reset pulse target and a sharper relay-reuse result

The stored144-state unitary reset target decomposes into127 transpositions.
Each transposition is the product of a two-level X pulse and a diagonal phase
pulse, giving254 algebraic pulses. Matrix-exponential regressions check the
actual signs, not just a permutation list. These joint leakage-sector controls
are additional assumptions. Existing controls commuting with the system
logical projector cannot synthesize the target: its projector commutator norm
is `sqrt(20)`. Environment reinitialization remains an open-system resource.

The relay finding improves on11455's residue census without retracting its
valid105-fault calculation. A complete ideal SWAP-out/CNOT/SWAP-back route is
exactly endpoint CNOT tensor relay identity. Propagate all105 first-route
Pauli faults through a second ideal route and compare perfect relay reset:
72 relay residues persist, but **zero endpoint Paulis change due to reset**.
Therefore a persistent relay error alone does not establish harmful ideal reuse
or a mandatory cleanup cost. Correlated/state-dependent noisy gates and a full
fault-tolerant schedule still need their own channel audit.

## External checks and prior ownership

- [Babu, Bajc and Susic, E6 symmetry breaking](https://arxiv.org/abs/2305.16398):
  trinification embeddings and symmetry-breaking choices are established;
  abstract subgroup embedding does not select the physical hypercharge.
- [Luscher, anomaly-free Abelian lattice theories](https://arxiv.org/abs/hep-lat/9811032):
  the external all-sector construction requires more than a chart connection.
- [Borissova and Dittrich, Lorentzian Regge calculus](https://arxiv.org/abs/2303.07367):
  complex area/angle branches and boundary action conventions.
- [11270–11274 canonical Higgs bridge](PASS11270_11274_PHYSICAL_FRONTIERS.md),
  [11438–11442 supplied coercive action](PASS11438_11442_COERCIVE_CURVED_RECOVERY.md),
  [11443–11447 preserved fields](PASS11443_11447_PRESERVED_GAUGE_LORENTZ_RECOVERY.md),
  [11448–11455 prior maps](PASS11448_11455_SUBGROUP_TORON_AUXILIARY_CONTROL.md).

Result searches included the actual numerical witnesses, native branching,
centralizer centroid, reset pulse count and photon requirement. The effective
paper entryw33_paper.tex forwards to its body and shared frontier ledger; its
existing physical-parameter and CC boundaries remain applicable. No new
physical constant, observed mass, mixing angle or TOE solution is asserted.

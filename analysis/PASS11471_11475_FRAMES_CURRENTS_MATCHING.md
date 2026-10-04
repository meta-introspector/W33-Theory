# Passes11471–11475: native frames, measure integrability and matched recovery

Reservation `c630384c2`. [Producer](w33_pass11471_11475_frames_currents_matching.py),
[certificate](../data/w33_pass11471_11475_frames_currents_matching.json),
[independent regressions](../tests/test_w33_pass11471_11475_frames_currents_matching.py).

All five requested branches were investigated. The strongest additions are an
integer native frame dictionary, explicit loop operators for the discarded
fermions, and a concrete leakage-control map with decoded dissipative channels.
The results do not complete gravity, vacuum selection or a TOE. The remaining
maps are stated below, rather than inferred from dimensions or solver flags.

## 11471: construct the native frame before identifying the soft modes

Pass11384 owns the actual signed cubic and integer E6 generators. Fix the three
coordinate vectors of its first nonzero cubic triad, indices `[0,17,26]`.
The exact integer fixing equations have a28-dimensional complex E6 kernel.
Store every integer coefficient and every resulting27-by27 generator. Bracket
closure follows algebraically: the parent algebra closes and a bracket of two
operators fixing the frame also fixes the frame. Regressions independently
replay all brackets. Transpose closure is additionally checked numerically.

The generator support decomposes into three fixed singleton coordinates and
three invariant eight-coordinate blocks. All45 actual signed cubic triads
become one frame term, three sets of four frame/eight pairings, and32
cross-eight contractions. This supplies an objectwise native dictionary with
the shape of the usual triality decomposition. It is not a new discovery of
E6 branching. Pass10956 already constructs exact Spin8 in the clock Albert
algebra; a further signed Gaussian-rational map now relates that carrier to the canonical native compact frame algebra. The map from Pass11466's numerical stationary frame remains unbuilt.

The map is built from10950's actual real/imaginary paired-root basis labels,
with native coordinate sign flips `[1,4,6,7,26]`. Its inverse has denominator2.
Conjugate all28 native generators, clear both inverse and Jordan-product
denominators, and verify the full derivation equation
`R(x circle y)=(Rx) circle y+x circle (Ry)` in integer arithmetic for all
basis pairs. Every transformed generator fixes the three actual clock
idempotents. Compact native combinations become real and span28 dimensions.
This identifies the canonical compact frame algebra with the exact clock
frame stabilizer already owned by10956; it upgrades a dimension comparison
to a verified carrier map. A further coefficientwise integer check proves full cubic equivalence: `N_native(T y)=N_Albert(y)`. It does not select a physical E6 real form or supply the numerical-vacuum conjugator. A naive
all-positive signed-frame basis failed the derivation test and was rejected.

The determinant comparison uses the actual Jordan trace and its powers:
`6N(y)=Tr(y)^3-3Tr(y)Tr(y circle y)+2Tr(y circle(y circle y))`.
Multiplying by12 clears all denominators in the paired-root basis. The whole
27-by27-by27 coefficient tensor equals the transformed native signed tensor
exactly, including zero coefficients. This gives a complex cubic-carrier
isomorphism and a compact frame-algebra dictionary, while keeping real-form
selection separate from an emergent physical Lorentzian spacetime.

Probe the four numerical soft vectors with gauge-invariant field-Gram spectra,
cross-Gram singular values and its complex trace. The largest differential
singular value is approximately3.42e-7. Those observables do not resolve the
four modes at the1e-6 cutoff. This **does not establish** that the modes are
pure gauge or exact moduli: the chosen observables need not separate orbits,
and second-order changes can survive a vanishing differential. The report
retains the prior unresolved-mode classification and the supplied potential.

## 11472: an actual integrable patch current, with locality kept precise

Pass11428 owns the analytic divided-difference derivative of the overlap
projector and its Weyl-bundle curvature. On a small L2 four-dimensional lattice,
use two individual link variations around a stored weak nonzero background.
All charged plaquette angles are bounded below1/30 at radial and derivative endpoints; their linear dependence bounds the entire homotopy in the trivial admissible sector. The integer hypercharges `[1,-4,2,-3,6]` have multiplicities `[6,3,3,2,1]`;
both their linear and cubic anomaly sums vanish. Neutral species do not affect
the gauge connection.

Let F(x,y) be the actual multiplet curvature. Define

```
j_x(x,y) = -y integral_0^1 t F(tx,ty) dt
j_y(x,y) =  x integral_0^1 t F(tx,ty) dt.
```

Differentiation gives `partial_x j_y - partial_y j_x = F(x,y)`: the integrand
is the derivative of `t^2 F(tx,ty)`. An eight-point quadrature and an independent
finite-difference curl agree to 1.08e-14 at the stored point.
The measured curvature is small, -1.31e-09. Vanishing anomaly sums do not
force finite-background curvature to vanish.

This is the standard radial homotopy applied to a previously constructed
native-charge overlap bundle, not a new reconstruction theorem. It constructs
a measure connection **on a contractible parameter patch**. No exponential
spatial locality estimate or gauge Ward identity is proved. Lüscher's
spatially local, gauge-covariant reconstruction conditions, relative sector
phases and a non-Abelian completion remain open. In particular, this patch
must not be called a completed local chiral measure.

## 11473: test inhomogeneous dynamics against the prior spatial constructions

The repo already owns spatial wave actions in
[w33_pass4041_4048_eight_physics_expansion.py](w33_pass4041_4048_eight_physics_expansion.py)
and a stronger periodic-cover construction in
[w33_pass11389_parabolic_spatial_cover.py](w33_pass11389_parabolic_spatial_cover.py).
The latter derives FCC cover data after selecting a parabolic, with supplied
wave dynamics. This packet does not replace those maps with a novelty claim.

Use the actual80-site W33 Levi adjacency A and `L=4I-A`. Map each Levi vertex
to a spatial field coordinate and each incidence edge to an action term;
supply an independent discrete time coordinate. The declared action is

```
S = sum_n [ ||q_(n+1)-q_n||^2/(2 dt) - dt q_n^T L q_n/2 ].
```

Every interior site/time equation is
`(q_(n+1)-2q_n+q_(n-1))/dt^2 + L q_n = 0`.
A generic inhomogeneous initial field is evolved for100 steps at dt=.5.
The exact conserved discrete energy is
`||q_(n+1)-q_n||^2/(2dt^2)+q_n^T L q_(n+1)/2`.
Its drift is approximately1.33e-15. A localized initial perturbation travels
at most one Levi edge per tick. This is a support statement, not a measured
relativistic light speed.

The dispersion is `4 sin^2(omega dt/2)=dt^2 lambda`. Since the maximum
Laplacian eigenvalue is8, strict nonzero-mode stability requires
`dt < 1/sqrt(2)`. At equality the extremal repeated root can cause secular
growth; a zero-mode uniform velocity is separately allowed. The dt=.8 control
has a characteristic root outside the unit circle. Stable evolution is not
evidence of emergent continuum Lorentz symmetry, three physical dimensions,
a W33-to-Regge metric or gravitational backreaction. The homogeneous
frustum equations remain owned by11466 and Jercher–Steinhaus.

## 11474: restore what port compression removes from loop calculations

Pass11466 owns the exact port restriction from486 to48 states per quark
sector, and its438 dark states per sector. Across both sectors the actual
masses and multiplicities are

```
m_-=10-sqrt(6) phi : 268
m_0=10            : 340
m_+=10+sqrt(6) phi : 268.
```

These876 species are Dirac colour fundamentals under the declared model.
Their Grassmann determinant and gauge contributions are not removed by an
orthogonal change of basis. The96 retained states reproduce port amplitudes,
but an equivalent effective action also needs the omitted loops.

Compute their one-loop MSbar potential, with supplied Nc=3 and scale mu=10:

```
V_dark(phi) = -Nc/(16 pi^2) sum_a n_a m_a(phi)^4
                        [ log(m_a(phi)^2/mu^2) - 3/2 ].
```

At phi=.81 its derivative is approximately17720.3736263 in the supplied
model units. An independent derivative replay verifies it. Its Taylor
coefficients at phi=0 include approximately12219.3347473 phi^2,
-977.5467798 phi^4 and1.46632017 phi^6. The constant and lower-dimensional
terms depend on counterterms; they are not a vacuum prediction. The initial
replay detected a precision downgrade inside automatic differentiation;
that evaluator was corrected before a certificate was written.

For Euclidean gauge polarization, explicitly define the convention

```
Delta Pi(Q^2) = T(R)/(2 pi^2) sum_a n_a integral_0^1 dx x(1-x)
                       log[1+Q^2 x(1-x)/m_a^2],   T(R)=1/2.
```

This is the subtracted scalar polarization with the gauge coupling factored
out. Its expansion is `T(R)/(2pi^2) [Q^2 sum(n/m^2)/30 -
Q^4 sum(n/m^4)/280 + ...]`. Independent adaptive integration checks the
quadrature at Q=.01,.1,1. These are explicit omitted scalar/gauge operators,
not a UV completion of the96-state model. Original high-energy species and
the lack of asymptotic freedom remain; observed masses, coupling boundary
conditions, physical units and scalar counterterms are not supplied by W33.

## 11475: native leakage controls and full decoded dissipative channels

The stored invariant frame embeds an actual12-state cell in160 edge modes;
its logical dark space has dimension2. The support components of its projector
partition the160 modes into groups of sizes
`[1,1,2,3,9,3,9,6,18,27,27,54]`.
A diagonal edge phase constant on each group preserves the cell. The stored
masks name **every edge address**, and their compressed matrices mix logical
and leakage sectors. Generic independent edge phases need not preserve the
cell, so projecting them without checking the outside-cell action would have
been an invalid control construction.

In the joint mask basis, the supplied native drift `H_point+.123 H_line`
has connected transitions. Assign mask energies `2^j`. Every positive
transition gap is distinct. Polynomial filtering in the adjoint diagonal
control isolates individual drift transitions; their real/imaginary
quadratures and commutators on a connected graph generate u12. Thus grouped
phase addressing supplies a concrete one-cell leakage-control resource.
It is an additional control assumption, not a demonstrated device capability.
Two-cell entangling hardware, environment resetting, calibration bounds and
fault-tolerant pulses remain open.

Independently, take the actual seven-CNOT dissipative relay channel from11466.
Undo the ideal endpoint CNOT, feed the other endpoint a maximally mixed state,
and trace that endpoint. This defines a declared one-qubit marginal noise
channel. Use seven independent copies on a Steane block, and perform full
encoding, six-syndrome recovery and decoding without Pauli twirling.
All64 recovery branches and arbitrary multiple damping events are included.
Tests verify every single-qubit Pauli correction and complete CPTP output.

| Damping per active wire/gate | Relay | Bare infidelity | Decoded infidelity |
| --- | --- | --- | --- |
| .01 | 0 | .04139 | .02630 |
| .01 | 1 | .06013 | .05091 |
| .03 | 0 | .11789 | .16045 |
| .03 | 1 | .16716 | .25601 |
| .12 | 0 | .37724 | .63763 |
| .12 | 1 | .48366 | .67728 |

Infidelity means entanglement infidelity against the identity after removal
of the ideal CNOT. The useful weak-noise range and high-noise degradation are
actual channel results, not a fault-tolerance threshold. Readout/recovery are
ideal, and taking independent marginal copies discards endpoint/block
correlations. Full correlated two-block decoded CNOT remains open.

## External checks used and ownership

- [Lüscher, abelian chiral gauge theories](https://arxiv.org/abs/hep-lat/9811032):
  reconstruction requires more than a smooth determinant-line patch.
- [Henning–Lu–Murayama, effective-field-theory matching](https://arxiv.org/html/1412.1837v2):
  fermion determinants and gauge-covariant loop operators; method ownership
  remains external. The potential and polarization conventions above are explicit.
- [Steane, multiple-particle interference and quantum error correction](https://arxiv.org/abs/quant-ph/9601029):
  code construction and standard single-error correction are prior results.
- [Jercher–Steinhaus, Lorentzian discrete cosmology](https://arxiv.org/abs/2312.11639):
  the prior11466 frustum action remains an imported method.
- [11466–11470](PASS11466_11470_STATIONARY_PORTS_CHANNELS.md),
  [11428–11432](PASS11428_11432_NATIVE_ALIGNMENT_MEASURE_COARSE_FLAGS.md),
  [11443–11447](PASS11443_11447_PRESERVED_GAUGE_LORENTZ_RECOVERY.md),
  native11384 and exact clock10956 retain their earlier constructions.

Validation includes source binding, integer fixing/brackets, an independent
homotopy identity, discrete energy/stability replay, independent loop
integration/differentiation, all single-qubit corrections, and logical Choi
positivity/trace preservation. No paper's physical conclusion is changed by
these finite controls.

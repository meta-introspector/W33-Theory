# Passes11298–11302: five executed frontier investigations

Reservation `4edf7dca1` follows the reviewed parallel formula-universe update
`660a5e06c`. These investigations use the actual signed E6 tensor and native
Levi incidence, rather than replacing them with dimension counts. Their
certificates label the scope of PASS. None asserts a solved TOE.

## 11298 — Running exposes a channel-rank bottleneck

For the economical two-frame/54 inventory of Pass11293, one heavy vectorlike
pair has active Yukawas on243 Weyl components and66 real scalar components.
Direct sparse contractions in the real-scalar/two-component-Weyl convention give

```
16π² βy1 = y1(40y1² + y2² − 52gE² − 8gF²)
16π² βy2 = y2(5y1² + 29y2² − 52gE² − 8gF²).
```

The gauge inventory includes both sextets and the anomaly spectator even
though only one active Yukawa channel is contracted. Its gauge coefficients
are `(bE6,bSO10,bFamily)=(2,8,−68/3)`, with `βg=−bg³/(16π²)`.
The five coupled gauge/Yukawa flows are integrated; these are dimensionless
toy initial conditions. For global family, the two nonzero fixed Yukawa ratios
are `y1²/gE²=40/33`, `y2²/gE²=50/33` in this truncation.

Crucially, a single pair gives the channel map
`C_ab=y1_a y2_b/M`, of rank at most one. It cannot independently implement
the rank-two diagonal map required for `d(vu)Fu + d(vd)Fd`. A second pair
costs twelve additional E6 beta units and gives `bE6=−10`. Two scalar
`(27,bar6)` mediators instead retain `bE6=2`, but their full flow is open.
This limits the use of the one-channel AF budget; it does not invalidate the
earlier single-channel construction.

The real SO10 symmetric-traceless54 necessarily admits two quartics,
`V=λ1(TrS²)²/4+λ2TrS4/4`. Their pure-S contributions are

```
16π² βλ1 = 124λ1² + (224/5)λ1λ2 + (159/50)λ2² −120gSO²λ1 +18gSO4
16π² βλ2 = 24λ1λ2 + (127/5)λ2² −120gSO²λ2 +60gSO4.
```

Independent full54-mode Hessian traces check the scalar terms. The exact
identity `Σi<j(si−sj)^4=10TrS4+3(TrS²)²` checks the gauge terms. The
isolated pure-S fixed-ray equations have no real solution, certified by a
Sturm count on their resultant. This is a statement about that restricted
sector, not a proof against the full model. K,S,U,V,A and family portals add
counterterms: the old positive-square coefficient ansatz is not RG closed.
**A complete two-channel quartic/Yukawa UV analysis remains open.**

Source conventions checked against equations19,33 of
[Luo–Wang–Xiao](https://arxiv.org/abs/hep-ph/0211440).

## 11299 — An E6 determinant cannot supply the missing flavor phases

The actual `d(e0),d(e1)` graph consists of five signed three-node paths, each
with one Fu edge and one Fd edge, and twelve unused vertices. Removing edge
signs by endpoint phases gives the exact block map `B=(Fu,Fd)`. Thus

```
spec(M†M) = 51 zeros ⊔ ten copies of spec(Su+Sd),
Su=Fu Fu†, Sd=Fd Fd†,
Tr((M†M)²) = 10Tr(Su²)+10Tr(Sd²)+20Tr(Su Sd).
```

The coefficient20 is fixed by the E6 support. Every spectral function of
this mass matrix, including its one-loop determinant, depends only on
`Su+Sd`. It is blind to relative Takagi phases that keep both Grams fixed.
Consequently this sector cannot generate the four phase locks or the
hierarchy coefficients imposed in Pass11294. Random complex symmetric
matrices and an independent phase deformation check the entire spectrum.
This result is restricted to the named two canonical Higgs directions and
their mass map; it is not a statement about all E6 interactions or all loops.
The old51 chiral zero modes are cited from Passes11271/11280, not claimed new.

## 11300 — An invariant vector counterterm with a declared EFT scheme

The actual vector Gram matrix at `gE=.5,gF=0` has70 positive eigenvalues in
six groups. A degree13 Hermite polynomial P matches the first and second
derivatives of the declared vector Coleman–Weinberg function at all six
positive eigenvalues, with `P′(0)=0`. The counterfunctional
`−3Tr P(X(q)/scale)/(64π²)` is invariant because X transforms by conjugation.
First trace jets use f′; second jets use f″ and divided differences of f′.
Those match even for noncommuting variations. At a Gram kernel the first
variation restricted to the kernel vanishes; emerging masses are O(q²),
so `x² log x` contributes no first or second field jet there.

This uses the earlier condensate EFT of Pass11290, not the enlarged
Higgs/vector spectrum of Pass11293. It cancels the **vector-sector** local
tadpole/Hessian jet. It requires
higher-degree EFT operators and a declared Landau/MS convention. It is not
a renormalizable UV counterterm basis, nor the full scalar/Weyl tadpole.
The residual constant is stored and not identified with the observed CC.

The radial two-vector cut is also integrated with three subtractions at
zero. All its thresholds are closed at the chosen scalar tree mass, so it
adds zero width and a leading squared-pole contribution approximately
`−1.2868583e−6`. The subtraction coefficients through s² are chosen inputs;
this conditional contribution is not a mass prediction or full pole matrix.
Prior owners: Passes11290/11295. Matching context:
[Goodsell–Paßehr](https://arxiv.org/abs/1910.02094),
[Goldstone resummation](https://arxiv.org/abs/1609.06977).

## 11301 — Native graph brackets generate missing edge constraints

On the native80-vertex/160-edge Levi graph, the named canonical scalar
Hamiltonian has the exact bracket

```
{H[N],H[M]} = Σedges (Ni Mj−Mi Nj) Dij,
Dij=(pi+pj)(φj−φi)/2.
```

The potential cancels. These mixed momentum/field edge currents are absent
from the site-only constraint ansatz; its pp/φφ quadratic terms cannot span
them with field-independent lapse coefficients. A symbolic single-edge
control and native incidence check this equation. A separate periodic-cycle
regulator has a nonzero continuum control with errors
`2.50e−4,5.35e−5,1.28e−5,3.17e−6` at16,32,64,128 sites.

The classical finite-algebra fact `e_i²=e_i ⇒ D(e_i)=2e_iD(e_i) ⇒ D=0`
also excludes a nonzero exact Leibniz derivative on all pointwise vertex
fields. This motivates deformed products, restricted fields or refinement;
it is not a no-go for every discrete gravity architecture. Related lattice
obstruction: [Kato–Sakamoto–So](https://arxiv.org/abs/0810.2360).
No gravitational edge variables, closed enlarged constraint algebra, real3D
refinement or W33-derived spacetime is supplied here. Pass11296's continuum
ADM calculation remains a separate conditional construction.

## 11302 — Separate membrane charge from sequestering flux

In linear local sequestering, variation of Λ fixes
`Fseq=μ4 vol` pointwise. Its extensive Qseq fixes volume. The same form
cannot also be the jumping Maxwell field in Pass11297. The named completion
adds a **distinct** Maxwell four-form Fmem and its charged membrane, while
retaining topological Fseq,Fhat. This is an explicit additional field.

A closed Euclidean two-cap S4 laboratory solves three equations for radius,
Λ and κ²: the Israel junction, total volume `Qseq/μ4`, and average-curvature
flux condition. Distributional wall curvature `3T area/κ²` is included.
Input flux sectors are generated from a reference solution and recovered
from perturbed starts. Residual is below1.1e−14. Common bare shifts are
canceled by `Λ→Λ−shift` with unchanged geometry and flux sectors.

This uses the established sequestering/vacuum-decay mechanism, not a new
mechanism: [Kaloper–Padilla–Stefanyszyn](https://arxiv.org/abs/1604.04000),
[local sequestering action](https://arxiv.org/abs/1505.01492).
The added cap geometry is not the prior genus81 W33 history. Charge, tension,
offset, flux sectors, negative modes, nucleation rates and the observed CC
remain inputs or open questions. Prior owners: Passes11284/11297.

## Corpus and validation

Result searches covered the coefficients, summed-Gram spectral law,
Hermite-vector construction and separate Maxwell/sequestering forms across
the results index, analysis, JSON, papers and site. Standard algebra,
general RGE formulas and sequestering are cited rather than claimed novel.
Five producers/certificates, independent focused regressions, the site,
scoped paper/ledger and rebuilt PDF accompany this packet. Continuity logs
record the reservation, execution and corrections; unrelated work is preserved.

Final validation: seven independent regressions passed in76.48s; the updated
vector-cut scaling control passed separately in82.81s. All current sources
compile. Twelve-file intake is clean with no guard collisions or forced
arithmetic; the results index is refreshed. The newest site card is unique
and the revised paper PDF built successfully (830.80KiB; typesetting warnings).

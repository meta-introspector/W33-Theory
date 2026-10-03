# Pass11389: a spatial cover from native harmonic periods

Reservation `4093ef96e`. This deliberately branches away from the previous
CP / hard matching / hub / wall / transmutation program.

**Result.** Select one W33 line and take its native unipotent radical
`N=C3^3`, rather than inserting a Cartesian spatial lattice. The `N`-fixed
Levi harmonic cycles define a saturated rank-three period quotient. Its
connected infinite cover has primitive FCC periods and a single isotropic
acoustic band. The construction also works for the point-side `H27` radical.
The selected context, equal edge weights and wave evolution are assumptions.
This is a concrete spatial-propagation candidate, not a derivation of gravity,
fermions, physical units or the selection of our space-time.

Producer: [w33_pass11389_parabolic_spatial_cover.py](w33_pass11389_parabolic_spatial_cover.py).
Witness: [stored integer voltages and harmonic representatives](../data/w33_pass11389_parabolic_spatial_cover.json).
Regression: [nine independent checks](../tests/test_w33_pass11389_parabolic_spatial_cover.py).

## Prior ownership and external checks

The Levi graph's 80 vertices, 160 edges and first homology rank81 are existing
results. [BT861](bt861_code_register_is_steinberg.py) owns its irreducible
Steinberg representation. The paper's `forty_points/sec04_time.tex` and
`sec09_physics.tex` already use that representation; we do not claim it anew.
[Pass11070](w33_pass11070_dual_parabolic_flag_diamond.py) owns the distinction
between the line radical `C3^3` and point radical `H27`, and their parabolics.
[Pass5677's prior-art audit](PASS5675_5682_external_prior_art.md) owns the
existing Levi voltage-cover tower and explicitly keeps physical speed open.
The earlier Cartesian-product continuum in
[BT4049–4056](BT4049_BT4056_five_front_outside_box.md) supplies its number of
directions separately. Here the period rank is computed from a selected
native subgroup and an explicit homology quotient.

Harmonic realizations of abelian graph covers are classical. The primary
reference is Kotani–Sunada, *Standard realizations of crystal lattices via
harmonic maps*, Transactions AMS353 (2001),1–20,
[DOI](https://doi.org/10.1090/S0002-9947-00-02632-5).
The framework is used here, not claimed as new. The rank-three harmonic
quotient, saturated W33 voltages, native FCC metric and its full-cell acoustic
corrections are this packet's constructed connection.
The classic diamond lattice is an abelian cover of the four-edge dipole;
this is explicitly described by Sunada in *Crystals That Nature Might Miss
Creating*, Notices AMS55 (2008),208–215
([primary author article mirror](https://www.fichier-pdf.fr/2013/11/08/tx080200208p/tx080200208p.pdf)).
our full cover has80 sites per primitive cell and is not asserted to be that
same abstract graph. The dipole appears after taking the radical quotient
and suppressing its degree-two vertices.

Before computation, searches checked `RESULTS_INDEX.md`, the paper, index and
analysis for the resulting period / harmonic / parabolic / FCC connection;
after computation they checked the actual metric entries, raw period index
and acoustic coefficient. Existing FCC integer numerology and unrelated
`256000` counts are not this period map. No universal novelty claim is made.

## 1. A period map, rather than a dimension assigned by hand

Use symplectic coordinates `(e1,e2,f1,f2)` overF3. For the line
`L=span(e1,e2)`, the radical consists of

    n(B) = [[I2,B],[0,I2]],       B=B^T overF3.

The producer constructs all27 matrices, verifies symplecticity and closure,
and computes their literal action on the points, lines and160 oriented
incidence edges. Orient every edge point to line, and call its boundary
matrix `D` (80×160). Let `O` have one indicator column for every edge orbit.
Then invariant harmonic flows are precisely

    H_N = O ker(D O).

There are16 edge orbits and14 vertex orbits. The quotient consists of four
internally disjoint length-four paths between two hubs. Each path has orbit
conductances27,9,3,1 and resistance

    1/27 + 1/9 + 1/3 + 1 = 40/27.

An invariant flow is specified by four branch currents with their sum zero.
This proves its dimension is three without reading dimension off an intended
space-time. It is also the elementary cell version of the rank-one Levi
Steinberg restriction. The full81-cycle space remains intact before this
quotient is selected.

An exact integer-scaled harmonic basis `C` has Gram

    C^T C = [[2160,1080,1080],
             [1080,2160,1080],
             [1080,1080,2160]].

Raw integer flows alone do not guarantee a connected Z3 cover. Tree-gauge
them: solve the tree potential `phi`, making
`P=C-D^T phi` zero on every tree edge. Chord rows of `P` are the actual cycle
periods. Their integer image has column-Hermite basis

    Lambda = [[160,120,120],[0,40,0],[0,0,40]].

Define `V=P Lambda^(-T)` and `H=C Lambda^(-T)`. The stored `V` is integral,
zero on the tree, and its chord voltages generate all ofZ3. Its infinite
derived graph has vertices `(a,z)`, `a=0,...,79`, `z in Z3`, and edges

    (a,z) -- (b,z+V_e).

Consequently the cover is connected. Independent regressions build the
actual reductions modulo2 and3 (640 and2160 vertices), checking connectivity
and degreefour at every vertex. `V-H` is an exact potential difference, so
both yield the same cycle periods and Bloch spectra.

The actual harmonic realization is also stored: with rational cell positions
`r_a` in primitive period coordinates, `H-V=D^T r`. There are only14 distinct
positions modulo the deck lattice among the80 graph vertices, with
multiplicities `1^5,3^4,9^4,27^1`. The realization is therefore **noninjective**.
Coincident vertices remain distinct internal states with the original graph
couplings. This is a periodic network with internal multiplicity, not an
asserted80-site embedded crystal with every vertex at a separate point.

The point radical is independently constructed by

    [[1,a,c,b],[0,1,b,0],[0,0,1,0],[0,0,-a,1]].

Its center has orderthree, while the line radical's center has order27.
Both have a rank-three fixed harmonic space and the same reported period
metric. Their subgroup actions and voltage maps are separate objects.

## 2. Native FCC metric and residual symmetry

At the equal-weight cell, define the acoustic tensor

    K = H^T H / 80
      = [[297/12800,-27/1600,-27/1600],
         [-27/1600,27/1600,27/3200],
         [-27/1600,27/3200,27/1600]].

These entries depend on the primitive period basis. The basis-independent
lattice statement uses the unimodular change

    U = [[-2,0,1],[1,1,-1],[2,-1,0]],       |det U|=1.

Exactly,

    U^T K^-1 U = (3200/27) [[2,-1,0],[-1,2,-1],[0,-1,2]].

The last matrix is the A3 root Gram, so the period metric is FCC up to one
overall scale. This identifies a lattice with a named map and an integral
basis change, not by matching a coordination number.

The actual GL2(F3) Levi action on these periods has24 distinct matrices after
its central scalar is removed. Their character has squared normone and
their invariant metric isK. Thus the residual S4 three-space is irreducible.
Its metric is unique up to scale; in a Cartesian FCC frame the acoustic
tensor is exactly `(27/3200) I3`. It has reflections as well as rotations;
it is not an asserted Lorentz group.

A symplectic exchange moves the selected line to a second line and carries
the full harmonic witness to another divergence-free realization with the
same Gram. Basis / drawing dependence is checked explicitly.

There is a useful full-symmetry boundary. A connected free-abelian regular
Levi cover is a quotient of `H1(Levi,Z)`. If its kernel is invariant under
the fullPSp(4,3), its rational quotient is a quotient of the earlier
irreducible81-dimensional Steinberg module. A nonzero such quotient must
therefore have dimension81. A rank-three cover cannot preserve that full
kernel symmetry. Selecting a parabolic is a real symmetry choice, not a
coordinate change. Parabolic automorphisms lift to cover automorphisms;
no split semidirect-product space group is asserted.

## 3. Actual propagation and its first anisotropic correction

The80×80 Bloch Laplacian has diagonal4 and, for an edge `a -> b`, entries

    L_ab(k) = -exp(i k.V_e),       L_ba(k)=conjugate(L_ab(k)).

It is Hermitian and positive semidefinite by its edge-energy factorization.
At zero momentum the constant band is simple and the next band is separated
by the prior native gap `4-sqrt(6)`. Harmonic gauge makes the first derivative
annihilate the constant vector. Consequently its lowest band obeys

    lambda0(k) = k^T K k + O(|k|^4).

In the Cartesian FCC momentum `p`, the producer computes the next term from
the full80-site graph, including relaxation into all optical modes:

    lambda0(p) = (27/3200)|p|^2
      - (223119/409600000) sum_i p_i^4
      - (367119/204800000) sum_{i<j} p_i^2 p_j^2
      + O(|p|^6).

For a direction `p`, put `t_e=E_e.p` and
`v2_a=1/2 sum_{e incident a} t_e^2`. The exact fourth coefficient is

    -sum_e t_e^4/(12*80) - v2^T L0^+ v2/80.

The second term is optical-mode relaxation; using only a frozen constant
mode would miss it. Directions `(1,0,0)`, `(1,1,1)`, `(1,1,0)` and `(1,2,3)`
verify the exact polynomial, while independent generic-momentum spectra
check the expansion numerically. The quartic term is not rotationally
invariant: Lorentz-like isotropy is a long-wavelength property here.

Supply the positive wave action

    S = 1/2 integral dt [sum_v dot(q_v)^2 - sum_edges(q_v-q_w)^2].

The single acoustic mode then has `omega^2=lambda0`, giving a three-spatial-
dimensional massless scalar wave equation at long wavelengths, with a time
coordinate supplied by this action. Its dimensionless squared speed in the
specified FCC coordinates is27/3200. Neither its SI speed nor a physical
lattice spacing is determined. This supplies a3+1 wave kinematics candidate,
not a proof that the physical universe chooses it.

Quantizing this particular free system is straightforward. With a positive
mass regulator, `Omega=sqrt(L+m^2)` is positive. Its Euclidean time covariance
`exp(-Omega |t|)/(2 Omega)` is reflection positive because, at positive times,
the reflected quadratic form is

    || (2 Omega)^(-1/2) sum_j exp(-Omega t_j) f_j ||^2 >= 0.

Thus the named model has a unitary free quantum realization. This standard
argument does not supply chiral fermions, interactions or Einstein dynamics.

## 4. An early test of geometric flexibility

Even this route needs a gravity boundary before promising curvature. Let all
160 conductances vary around their common value. Since the harmonic cell
corrector minimizes the Dirichlet energy, the envelope derivative is

    delta K = (1/80) sum_e delta w_e E_e E_e^T.

The exact160×6 response matrix has rankfour. Every displacement outer
product has equal Cartesian diagonal entries, so the two trace-free
diagonal metric directions are missing at first order. The certificate
stores their null tensors. Varying these conductances therefore cannot
produce an arbitrary local spatial metric near the symmetric cell at
linear order. Extra native channels, a different cell geometry, nonlinear
effects or a different geometric mechanism would be needed.

## 5. The field parameter also counts harmonic directions

This has a transparent generalization. For a selected line in `W(3,q)`, write
outside points as `(x,y)` with `y != 0`. The symmetric matrices of the line
radical act by `x -> x+B y`. For every direction in `P1(Fq)` this sweeps one
q²-point orbit. Points on the selected line are fixed. Lines transverse to
the selected one are the q³ symmetric graphs, a single regular radical
orbit. The q other lines through each selected point form one q-line orbit.

Thus there are `3q+5` vertex orbits and `4(q+1)` edge orbits: exactly q+1
length-four paths with edge conductances `q³,q²,q,1`. Their harmonic currents
sum to zero, giving **q independent directions** over characteristiczero.
This orbit argument works over every finite field. The underlying
Steinberg/Levi restriction is classical; its use here is an explicit period
quotient for spatial propagation.

An [independent prime-field builder](w33_pass11389_prime_period_family.py)
checks the stronger *integral* period statement at q=2,3,5. Put
`S=q³+q²+q+1`. In these exact witnesses,

    C_q^T C_q = q³ S (I_q+J_q),
    period image = S (I_q+J_q) Z^q.

The unimodular relation between this image and its Hermite basis proves the
primitive metric is the A_q root metric up to scale. The Cartesian acoustic
coefficient in root-lattice units is

    q³ / (2 S²).

The witnesses give A2, A3 and A5 respectively. The integral law has not been
proved here for every prime power; extension fields are not implemented.
The native ternary case has three directions because q=3 **after** the
selected-line harmonic quotient is chosen. This does not establish a
dynamical preference for ternary geometry or for that quotient.

## 6. Native H27 qutrits are present, but their original bands are flat

The coincident internal states can be analysed without assigning particles.
On the point-side H27 cover let `P_N` average the27 actual cell permutations.
It has rank14 and averages only within vertex orbits; set `Q=I-P_N`, rank66.
The character on the center is `(80,17,17)`. Exact cyclotomic character sums
give the complex80-state decomposition

    14 trivial copies
      + 3 copies of each of the8 nontrivial linear characters
      + 7 copies of each of the2 conjugate3-dimensional H27 irreps.

The two nontrivial central eigenspaces therefore have rank21 each; the
central-trivial eigenspace has rank38. These are module dimensions, not
fermion counts or Standard Model assignments.

The producer builds an actual intertwiner. For native H27 matrices `X,Y,Z`,
choose the central-omega projector and the `Y=1` projector. Their joint image
has dimensionseven. An orthonormal basis `S` of that image defines

    T = [S, X S, X² S],
    T^dagger T = I21,
    T^dagger L(k) T = I3 tensor h7(k).

The certificate stores the complex seed basis, group permutations and
operator residuals, so the map can be replayed. This is a genuine qutrit
factor in the named model, not a dimension-only identification.

But its spatial kinetics need repair. Every nontrivial radical
representation vanishes on singleton vertex orbits. Removing those vertices
from the radical quotient leaves a **tree**. The harmonic edge phases on
that tree have an exact vertex potential `psi`; all omitted edge blocks
annihilateQ. Hence for every momentum

    L(k) Q = U(k) L(0) U(k)^dagger Q,
    U(k)=diag(exp(-i k.psi_a)).

All66 internal bands are therefore flat, not just at sampled momenta. Their
exact spectrum is `(4-sqrt6)^20, 4^26, (4+sqrt6)^20`. Within either qutrit
sector the seven-mode spectrum is `(4-sqrt6)^2, 4^3, (4+sqrt6)^2`.
The14-dimensional invariant sector carries the dispersive acoustic band.
Native internal qutrit modes do not carry propagating information across
this space under the unmodified operator.

## 7. An explicit common-speed completion

There is a concrete extra coupling that repairs this. For the six positive
FCC root directions `(1,+/-1,0)`, `(1,0,+/-1)`, `(0,1,+/-1)`, define

    ell_FCC(p) = sum_roots [2-2 cos(p.root)] = 4|p|²+O(|p|⁴),
    L_completed(p) = L_native(p) + kappa ell_FCC(p) Q.

Q is a projector on internal states at coincident positions. Their stored
harmonic positions differ by integer periods within each orbit, so this
coupling has a local realization: hopping between neighboring FCC cells,
acting on the internal fiber. It is an additional channel, not a polynomial
in the native Laplacian; spectral functions of the latter would keep the
original flat bands flat.

The additional operator is positive semidefinite, preserves the native
radical and its normalizer, and leaves `P_N` untouched. Internal bands acquire
quadratic speed `4 kappa`. Requiring the same leading light cone as the native
acoustic mode fixes

    kappa = 27/12800.

That **common-speed condition is supplied**; symmetry alone does not fix
kappa. It needs no observed mass or speed input, but is still a model choice.
This completion is a positive free wave operator whose internal modes have
the same leading kinetic cone as its acoustic field. Neither their gaps nor
that cone have been calibrated to physical particles or SI units.

The qutrit factor remains a spectator to propagation:

    T^dagger L_completed(p) T = I3 tensor [h7(p)+kappa ell_FCC(p) I7].

Consequently every supplied logical gate `G3 tensor I7` commutes with this
kinetic operator. The regression checks a nontrivial arbitrary phase gate
as well as the factorization and hopping. This gives a mathematical
architecture for carrying an internal qutrit through the constructed space.
It does not derive access to arbitrary gates, entangling gates, fault tolerance
or a universal quantum computer from the native graph.

## 8. A gapless alternative and a concrete protected zero

The optical gaps above are orderone in lattice units. Under a continuum
refinement with fixed physical propagation speed, they generically grow as
inverse lattice spacing and decouple. They must not be read as finite
continuum particle masses without another mechanism.

The native bipartite adjacency gives an alternative. Define
`A(p)=4I-L_native(p)` and grade points by+1 and lines by-1, a matrixGamma.
Exactly `Gamma A+A Gamma=0`. The following **different supplied** positive
wave operator retains the acoustic band but changes the internal action:

    L_gapless = P_N L_native P_N + Q A² Q + (27/12800) ell_FCC Q.

It uses at most two native steps before the already named additional FCC
hopping. At zero momentum its exact internal zero projector is

    P_zero = Q - (L0-4I)² Q / 6,
    P_zero²=P_zero,       rank P_zero=26.

Inside either nontrivial central H27 sector, nine zero states remain: three
on the point side and six on the line side. Dividing by the actual qutrit
factor gives a3-by4 bipartite multiplicity operator, with index-1. The
dimension imbalance forces at least one zero qutrit copy under any
H27-equivariant bipartite first-order deformation. The native block has rank
two, so it happens to have three zero copies; two are accidental.

The producer constructs the actual rectangular block from the stored qutrit
intertwiner. Its singular values are `sqrt6,sqrt6,0`. Adding a rank-one map
between its left and right null vectors with supplied coefficientepsilon
raises the rank tothree. The squared seven-mode gaps become

    0, epsilon², epsilon², 6,6,6,6.

The epsilon=1/10 witness is replayed independently. One zero qutrit survives,
and a paired small gap can be produced without lifting the heavy pairs.
No value ofepsilon is predicted. The common-speed hopping turns the
remaining zero modes into dispersive massless wave modes in this action.

The generic rectangular-index mechanism is prior mathematics and appears
already in [the August flat-band audit](2026-08-29_PHYSICS_OUTSIDE_BOX_FLATBAND_SYNTHETIC_KOOPMAN.md).
The new object here is the native H27 intertwiner and its actual3-by4 block,
zero projector and spatial kinetic completion.
This point/line grading is **not Lorentz fermion chirality**. Index protection
holds within the stated first-order bipartite family; generic bosonic mass
counterterms are not automatically forbidden. Chiral fermion dynamics,
anomalies, radiative protection and observed flavor require separate maps
and checks.

## What this changes

There is now a specified route from finite W33 topology to connected spatial
translations, a primitive metric, a propagation operator and a free quantum
model. It does not require installing a3D Cartesian product first. The
selected parabolic still chooses the three-period quotient; a dynamical
reason to select it is open. The construction filters future spacetime
ideas through explicit period maps, dispersion and metric-response tests.
It neither replaces the internal E6 matter map nor derives its fermion
spectrum, masses, couplings, gravity or vacuum energy.

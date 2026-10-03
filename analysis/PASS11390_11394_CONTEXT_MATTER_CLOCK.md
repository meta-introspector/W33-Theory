# Passes11390–11394 — selection, matter, an interaction portal, chart holonomy and a local clock

Reservation86ebd9e68; follows Pass11389. All five requested directions are
executed in `w33_pass11390_11394_context_matter_clock.py` and its JSON witness.
Generic Potts ordering, CAR quantization, Wilson/overlap fermions, Schur
elimination and Grover/Szegedy walks are established tools, not new claims.

The substantive results here concern the **actual W33 period frames and H27
incidence carrier**: a shared native-edge null projector, a calculable virtual
path portal, singular versus integral chart transitions, their full FCC
holonomy group, and the unitary clock on the connected native voltage cover.
None of these establishes a complete theory of observed physics.

## Prior work and intake

Pass11389 owns the connected rank-three Levi voltage cover, primitive FCC
metric, actual H27 qutrit intertwiner, flat internal bands, and rectangular
3×4 multiplicity block. BT861 owns the rank81 irreducible Steinberg homology
and the obstruction to preserving a rank-three cover under the full group.
Pass4084 (`BT4081_BT4088_deep_physics.md`) already constructs overlap fermions
and audits their doublers on an external four-dimensional lattice. Pass4057–4064
records the vectorlike/Wilson boundary. Pass3325 and3348 already compile
Szegedy walks; their generic mapping is not rediscovered here.

The new remote commit8f379c08b updates the formula-search inventory and was
integrated before reservation. Result searches included 1/27 harmonic
matches, rank-one chart overlaps, common-left kernels, 2160/1080 triangle
counts, and the virtual-path mechanism. The generic square-root-three
coefficient appears elsewhere, including SU(3) structure constants; only its
specific projector identity below is asserted for this construction.
Pass4669 already builds a forty-to-eighty orientation cover of F4/triality moduli, and Pass4814 builds a twenty-seven-to-fifty-four cover. The orientation-double-cover construction itself is prior art here; no identification with those different carriers is assumed.
The native chamber transport and isoclinic 24+24 spaces in Pass4324–4334 are
other carriers, not these forty rank-three spaces.

##11390 — a selector that distinguishes quantum mixing from ordering

The forty isotropic lines are the possible selected parabolics. Their native
intersection adjacency A has exact spectrum

    12^1, 2^24, (-4)^15,
    A^2=8I-2A+4J.

A single selector with Hamiltonian -hA, h>0, has a **unique uniform ground
state** and gap10h. Native tunnelling alone therefore does not select a
spatial context.

An explicit many-cell proposal on a supplied connected spatial host is

    H_selector=-h sum_x A_x + J sum_<xy> [1-delta(L_x,L_y)].

It is equivariant under simultaneous permutation of all context labels by
PSp(4,3). At h=0 and J>0, its connected-host ground configurations are exactly
the forty constant labels. Any domain wall costs J per mismatched host bond.
An independent five-bond chain with distinct fixed endpoint contexts has
minimum energy J. Unordered label pairs comprise240 intersecting and540
disjoint wall types.

This supplies a genuine operator and defects; it does not derive J/h,
prove a finite-h phase transition, or assert symmetry breaking in a finite
symmetric system. A fixed spatial host in this ansatz is an input. A coupled
host-and-selector theory still needs an action that actually varies the cover.
The discrete vacuum set has domain walls; it alone does not produce continuous
vortices, textures or an Einstein metric.

##11391 — fermions need CAR and a spacetime spinor factor

A point/line grading index is not Lorentz chirality. We explicitly add a
four-component spatial Dirac spinor and impose canonical anticommutation
relations on its Fock space. Four Jordan–Wigner modes act on a16-dimensional
Fock space; their CAR are checked directly.

Let F be Pass11389's primitive FCC basis, theta=F^T p the primitive torus
momenta. The naive Weyl symbol is

    h_W(theta)=sigma dot [F^(-T) sin(theta)].

It has eight zeros, theta_j in {0,pi}. Their Jacobian charges are
sign(det(F^(-T)))*(-1)^(number of pi coordinates), which sum to zero.
Tensoring with the native H27 qutrit does not remove these doublers.

A supplied Wilson Dirac symbol

    h_D=alpha dot [F^(-T)sin(theta)]
        + beta*r*sum_j[1-cos(theta_j)]

has one light **Dirac** corner for a massless internal carrier; other corners
have Wilson mass2rn, n=1,2,3. The resulting representation is vectorlike.
For a supplied common U(1) charge q its cubic and mixed gravitational anomaly
coefficients cancel as3(q^3-q^3)=0 and3(q-q)=0. That cancellation is conditional
on the supplied representation, not a derivation of Standard Model charges.
H27's dimension-three fiber does not establish a continuous SU(3) gauge group.

The actual operator I3 tensor [alpha dot d tensor I7 + beta tensor (S^dagger A_native S + Wilson_mass I7)] has dimension84 in one central sector. Its Clifford square and Hermiticity are checked, using the stored native intertwiner rather than proposed masses. Native dimensionless gaps remain cutoff-unit data.

The 7-dimensional native multiplicity matrix has three zero copies before a
rank lift and one afterwards. These can be internal labels of Dirac fields;
the rectangular index protects a graph kernel, not a chiral Lorentz species.
The operator names every factor rather than calling an eigenvector a fermion.
Exact overlap chirality can use the prior Pass4084 construction; it is not
claimed anew or confused with a chiral gauge-theory solution.

##11392 — the native edges cannot supply epsilon, but an explicit portal can

Let I be the40×40 point-to-line incidence matrix and Z_p,Z_l the central H27
permutations on the two sides. Define

    C_p=(2I40-Z_p-Z_p^2)/3,
    C_l=(2I40-Z_l-Z_l^2)/3,
    P=(I40-I I^T/6) C_p,
    Q=(I40-I^T I/6) C_l.

Here C_p,C_l project onto both conjugate nontrivial central sectors. Exact
integer arithmetic proves P^2=P, rank(P)=6; Q^2=Q, rank(Q)=12.
In one central qutrit irrep these correspond to one point kernel copy and two
line kernel copies. The JSON stores18P and18Q as integer matrices.

For **every one of the16 H27 edge orbits**, let E_o be its40×40 incidence
matrix. The actual matrices obey

    P E_o=0.

Consequently every H27-invariant weighting of the original incidence edges
still annihilates the same point carrier. Its3×4 qutrit multiplicity block has
rank at most2. This is an exact common-kernel statement for all weights, not a
random-weight survey. Functions of native adjacency do not bypass its null
space. It does not forbid symmetry-breaking edge interactions, extra fields,
or arbitrary many-body mechanisms.

A nonincident point/line pair, vertices1 and47 in the reproducible ordering,
has a27-element H27 orbit. Let C be that orbit's incidence matrix. Exact
arithmetic gives

    P C Q C^T P = (3/4) P.

Thus the first-order small singular-value coefficient per qutrit fiber is
**sqrt(3)/2**. This coefficient has an exact carrier identity behind it.

We realize C through an actual three-edge Levi path, then take its27 H27
translates. Introduce two auxiliary modes along each path: one at its internal
line vertex and one at its internal point vertex. Their54×54 heavy operator
is a direct sum of [[0,M],[M,0]]. The light endpoints couple to them with
y_left,y_right. All three couplings follow native incidence edges, while the
auxiliary modes have distinct orbit labels; this is an **extended carrier**,
not a reweighting of the original80 states.

The exact energy-zero Schur complement is

    delta I_eff = -(y_left*y_right/M) C.

Its diagonal self-energy vanishes because the heavy block is bipartite.
The resulting leading small singular value is

    epsilon_eff = (sqrt(3)/2) |y_left*y_right/M| + O(g^2),
    g=y_left*y_right/M.

For equal y and the supplied M=sqrt6, independent scans converge to0.8660254
for epsilon_eff/(y^2/M). At y=.1 the effective small singular value is
0.00353552654. The producer computes the effective operator rather than
inserting the previous artificial SVD rank-one perturbation.

This is a calculable quadratic portal mechanism; **y and M remain inputs**.
It does not predict an observed fermion mass or independently generate the
small endpoint coupling. The zero-energy Schur complement is not the exact
finite-energy light dispersion. Wave-field positivity and quantum
renormalization require their own action and analysis. Neither an arbitrary
auxiliary heavy mass nor a cutoff gap is declared a measured particle mass.

##11393 — chart changes produce a lattice bundle with full FCC holonomy

Generate all forty selected-line frames H_L by actual symplectic
transvections of the W33 Levi graph. All lie in the same rank81 harmonic
carrier and have a common Gram G. Every pair is checked exactly.
In orthonormal frames the overlap singular values are

| Context relation | Pairs | Singular values | Rank |
|---|---:|---|---:|
| Intersecting lines |240|1/3,0,0|1|
| Disjoint lines |540|1/27,1/27,1/27|3|

Intersecting charts lose two directions under projection. Completing a polar
map there would introduce an arbitrary choice. Disjoint pairs instead have
a canonical transport

    R_{j<-i}=27 G^(-1) H_j^T H_i,
    R^T G R=G.

All540 undirected transport matrices are integral unimodular. Their
contragredient R^(-T) acts on the primitive deck lattice. This therefore builds
an actual rank-three integer lattice bundle over the **disjoint-context
graph**, not just a suggestive overlap plot.

For all3240 disjoint context triangles, the loop matrix W satisfies W^2=I:

| Trace(W) | det(W) | Triangles | Cartesian type |
|---:|---:|---:|---|
|-1|+1|2160|half turn|
|+1|-1|1080|reflection|

The complete fundamental-loop closure at one context has order48. Conjugating
the deck action into Pass11389's Cartesian FCC coordinates gives exactly all
signed permutation matrices: the full octahedral point group O_h. This is
larger than the24-element native line-Levi image. Individual edge transports
comprise24 matrices but are not that native subgroup in these coordinates;
loop closure, rather than a dimension count, establishes the48-element result.

The negative-determinant triangles obstruct a globally oriented choice of
frames: their sign cannot be removed by independently changing vertex frames.
The determinant orientation double cover has80 context states, is connected,
and restricts holonomy to the24 proper rotations. This builds an explicit
orientation repair, not a spin lift or a spacetime field equation.

All40 native symplectic transvections act by exact integral changes of frame at all40 contexts (1,600 maps checked). Their determinant lifts act transitively on these80 oriented contexts; transport covariance is checked. Thus full PSp acts on this family of rank-three period lattices, although no one rank-three cover kernel is invariant. This is a constructed way to retain covariance while allowing the context to be a degree of freedom, consistent with BT861. It does not construct the coupled host-selector dynamics.

For comparison, transitions built by choosing one symplectic representative
per context telescope to identity around every loop. The nontrivial connection
above comes from **projected harmonic transport**, not a choice of unrelated
coordinate representatives. Its holonomies are invariant up to conjugation
under independent changes of local basis.

The bundle is over a finite context graph. It does not yet identify context
loops with spatial plaquettes, glue intersecting covers, produce a continuous
metric or derive Einstein dynamics. A finite nonzero holonomy is not by itself
a gravitational solution. The missing map is now explicit: an evolving
selector configuration would have to pull this bundle back to physical host
bonds and specify what happens on singular intersecting transitions.

##11394 — a local unitary update on the native connected cover

The160 incidence edges have320 oriented arcs. Let d map each outgoing arc
at a vertex to1/2 of that vertex state, so d d^dagger=I80. Define the coin
C=2d^dagger d-I320. The voltage shift S reverses each arc while translating by
its signed integral voltage. Then

    S^2=I, C^2=I, U=S C is unitary,
    d S d^dagger = I-L_voltage/4.

This is a concrete discrete update on the **actual rank-three connected
cover**. Each coin is local at its vertex and each shift crosses one native
edge. After n ticks support lies within n graph edges. With Pass11389's
Cartesian harmonic edge displacements, the maximum squared step length is
2187/6400, giving a strict Euclidean support bound n*sqrt(2187/6400).
That bound concerns geometric vertex positions, which have internal
coincidences; it does not discard coincident states.

The standard Grover spectral mapping gives

    cos(omega)=1-lambda_Levi/4,
    omega^2=(27/6400)|p|^2+O(|p|^4).

The producer checks unitarity, discriminant identity, spectral eigenvectors
at nonzero momentum, and four-tick support. A bipartite partner occurs at
quasienergy pi. The native internal flat bands remain flat: a change from
continuous to discrete evolution does not create their spatial propagation.

A tick duration measured in seconds remains an input. The acoustic cone and
strict support cone need not coincide. Exact Lorentz invariance, a derived
physical speed of light, quantum universal controls, and the extension of
this walk to the added FCC fiber interactions are not asserted.

## Compatibility boundary between these branches

The chart bundle and primary clock above use the line-selected cover; the H27 fermion/portal carrier comes from the point-selected cover at zero Bloch momentum. Pass11389 separately constructs both with the same quadratic metric. Equality of that metric is not a point-to-line intertwiner. These five constructions are not yet one interacting spacetime/matter action, and no unbuilt identification of the two selected parabolics is used.

## External checks used

- Neuberger, *Exactly massless quarks on the lattice*,
  https://arxiv.org/abs/hep-lat/9707022 — established vectorlike overlap construction.
- Neuberger, *More about exactly massless quarks on the lattice*,
  https://arxiv.org/abs/hep-lat/9801031 — exact Ginsparg–Wilson identity.
- Szegedy, *Spectra of Quantized Walks and a sqrt(delta epsilon) rule*,
  https://arxiv.org/abs/quant-ph/0401053 — reflection-based walk spectral framework.
- Higuchi, Segawa, Suzuki, *Spectral mapping theorem of an abstract quantum walk*,
  https://arxiv.org/abs/1506.06457 — discriminant/unitary spectral mapping.

- Sjöqvist, Kult, Åberg, *Manifestations of quantum holonomy in interferometry*,
  https://arxiv.org/abs/quant-ph/0607198 — discrete subspace holonomies are established mathematics; the W33 integral-bundle and full-loop group computation above are the present specialization.

## Validation contract

The producer must replay all five sections to PASS. Ten independent regressions
reconstruct the native integer projectors and all edge orbit blocks, the
virtual-path Schur operator, sample chart overlaps and stored group closure,
the orientation cover, CAR/doublers, and voltage-walk eigenvectors/locality.
No generic low-dimensional numerical fit is treated as proof of a continuum
physical prediction. Papers are not amended with this supplied dynamics.

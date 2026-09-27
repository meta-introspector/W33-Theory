# Passes 11070–11075 — oriented chambers, dual parabolics, and the local-E8 overlap nerve

This packet cross-checks the building/oriented-chamber proposal against the older exact repository certificates and separates two orientation notions that initially looked identical.

## 1. The dual 648s form a parabolic diamond, not an equality

For an incident W33 flag p<L,

    G_p ∩ G_L = G_{p<L} = B,
    |G_p| = |G_L| = 648,
    |B| = 162.

The two order-27 normal radicals are different:

    O3(G_p) = H27,
    O3(G_L) = F3^3.

Inside the flag Sylow-3 group they satisfy

    |H27 ∩ F3^3| = 9,
    <H27,F3^3> = U81.

Thus the chamber update group is literally the weld of the noncommutative point-state torsor and the commuting line-program torsor.

## 2. Oriented chambers are a canonical double cover

The flag Borel is the Sylow-3 normalizer B=N_G(U81), with |B:U81|=2. Hence

    G/U81 -> G/B

is a two-sheeted homogeneous cover with 320 oriented chambers over the 160 W33 flags.

On one four-point line, choosing the omitted point leaves exactly two nonzero covectors +ell and -ell. Evaluating them on the four projective points gives the opposite pair of weight-three tetracode words. Therefore the eight oriented chambers over one line are the eight nonzero tetracode words after a local coordinate/sign choice.

## 3. Critical correction: chamber orientation is not cubic-tick chirality

The internal Borel deck involution may be represented by T=diag(1,-1,1,-1) in Sp(4,3). It swaps +ell <-> -ell, but on C2 root coordinates acts

    (a,b,c,d) -> (-a,b,-c,d).

So the highest-root coordinate d is fixed and the explicit H27 cocycle kappa is invariant.

By contrast, the outer symplectic similitude S=diag(1,-1,-1,1), with S^T J S=-J, acts

    (a,b,c,d) -> (-a,-b,c,-d)

and sends kappa to -kappa. It reverses the highest-root tick.

Both involutions swap the local tetracode sign, but only the outer PGSp involution reverses the cubic cocycle/tick. Therefore there are two distinct Z2 data: the internal chamber-sheet sign and the outer tick/chirality sign.

## 4. Forty local tetracode/E8 charts have an exact 81-cycle overlap topology

Use the 40 W33 lines as local-chart labels. Four charts meet at each W33 point. Declaring those four labels to span a simplex gives a 3-dimensional overlap nerve with

    f = (40,240,160,40).

Its boundary ranks are (39,120,40), so

    H0 = Z,
    H1 = Z^81,
    H2 = H3 = 0.

Distinct point tetrahedra intersect only in chart vertices. After barycentric subdivision each tetrahedron collapses onto the four incidence spokes from its point barycenter, leaving exactly the W33 Levi graph. Thus the chart nerve has the Levi graph homotopy type.

This supplies a clean new meaning to the recurring 81: it is the first homology rank of the local exceptional-chart overlap atlas.

## 5. Exact local E8 versus the open global amalgam

The existing W33 tetracode producer gives an exact rank-8, 240-root, reflection-closed E8 root system by the standard four-A2 Eisenstein lift. W33 is transitive on its 40 lines, so every line supports a coordinate-equivalent local construction.

But forty local charts do not yet imply one global E8. The next exact object is the amalgam diagram:

    40 local E8 charts
    160 A2 point-incidences
    240 pairwise chart overlaps
    81 independent overlap cycles.

We now need explicit A2 -> E8 embeddings on every incidence and must compute their universal completion. It may close to one E8, collapse to a quotient, require a cocycle/higher gauge field, or become infinite. The topology is frozen; the algebraic transition maps are the missing datum.

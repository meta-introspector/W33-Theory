# 2026-10-01 — class-three 54D transducer, declared Cartan grade-pair metric, and orientation bit

Certificates:

- `data/w33_20261001_u81_54d_cocycle_transducer.json`
- `data/w33_20261001_cartan_kinetic_wall_vacuum.json`
- `data/w33_20261001_u81_antiunitary_orientation.json`

## 1. The 54D chamber transducer survives the nonsplit weld

Pass 11062 found four single-foliation rank-54 transducers in the split K81
carrier: plus factors 7/8 and minus factors 6/9. Passes 11067--11069 later
constructed the explicit H27 two-cocycle and the nonsplit U81 laws s=+1,-1.

For all four sparse factors, their generators remain order three for s=0,+1,-1.
Using the common H27 quotient, each split 3-cycle is transported by keeping the
same H27 right-coset representative, the same central seed, and the same
position in the generator cycle. All eight nonsplit replays (four factors
times two orientations) remain 27 disjoint triangles.

At both split primes 103 and 109 every transported weighted factor has:

- raw rank 54;
- S2+L quotient rank 54;
- S2 rank 27;
- L rank 27;
- the same 54 pivot columns.

Because transport is permutation conjugation, raw rank and triangle structure
are exact over characteristic zero; the two nonzero split-prime minors certify
the Q(omega) quotient transversality.

There is a stronger simplification for the four perfect factors. Their H27
generators all have first coordinate A=0, while the certified cocycle is

    kappa((a,b,c),(A,B,C)) = A^2 b - 2 A c.

Therefore kappa(h,g)=0 for every carrier label h. Their right actions are
literally identical for s=0,+1,-1: no nontrivial transport permutation is
needed. In fact factors 6,7,8,9 are exactly the cocycle-blind four, and the
perfect subsets are plus 7/8 and minus 6/9.

For any nonzero real Cayley parameter t, Im(C_t(F)-I)=Im(F), so the same
finite one-layer 27-block gate retains the full 54-dimensional retyping image.
This closes the algebraic transducer-existence question. It does not yet supply
an optical Hamiltonian or hardware implementation.

## 2. The E8 bracket supplies a reconstructed Cartan grade-pair metric

On the exact semisimple Cartan basis q1,q2,q3, the frozen E8 bracket gives

`G_ij = Tr(ad(q_i in g1) ad(q_j in g2))`

equal to

`[[590,20,0],[20,980,0],[0,0,300]]`.

Its principal minors are 590, 577800, and 173340000, so it is positive definite
after the declared g1<->g2 grade-exchange identification. The adjoints are
integer reconstructions from numerical bracket coordinates; the trace identity
is not independently certified by exact bracket arithmetic. The local Hessian
calculation below is exact conditional on the displayed Gram matrix.
## 3. Equal mirror repulsion selects a CP-symmetric local vacuum

Let F=C3*N9*S9, the product of all 21 mirror forms up to a nonzero scalar.
The zero-fit projective functional

`Phi = |F|^2 / (q^dag G q)^21`

has an exact strict local maximum at the rational Cartan ray [0:0:1]. The real
and imaginary tangent Hessians of log Phi are both negative definite.

Under the exact Cartan-to-qutrit map this ray becomes
`[1, zeta9, zeta9^-1]` up to phase. Feeding it into the certified cubic
postselection interface gives Jarlskog magnitude below 1e-14: the simplest
maximally-regular vacuum is CP-conserving.

Thus this kinetic candidate plus equal mirror repulsion does not by itself
choose an orientation.
## 4. The missing orientation datum is already present in U81

In the qutrit Schrodinger representation, complex conjugation acts on H27 as

`tau(a,b,c)=(-a,b,-c)`.

This is an H27 automorphism and leaves the cocycle kappa invariant. Extending
it to the chamber coordinate by `Psi(h,d)=(tau(h),-d)` gives exactly

`Psi(x star_s y) = Psi(x) star_{-s} Psi(y)`.

All 13,122 element-pair identities were checked. The class-three central tick
also changes sign with s. Therefore the two nonsplit chamber orientations are
an antiunitary-conjugate pair.

This explicit mechanism supplies the odd
datum the symmetric Cartan potential lacks. It still does not energetically
choose one sign or identify the finite antiunitary with observed Standard-Model
CP.

# Passes 11570–11579 — executable exceptional gauge / curvature / overlap frontier

This packet executes the five post-11569 frontier attacks plus five auxiliary probes.

## 11570 — the old F4 Standard-Model intersection is now executable

Using the committed Albert multiplication table, the nine basis directions
`[0,9,26,10,11,1,2,24,25]` close under the Jordan product and form an explicit real `H3(C)` subalgebra.

The full compact F4 derivation algebra has dimension 52. Solving the exact preservation equations for this 9D subalgebra gives a 16D compact semisimple stabilizer. Its centroid has two 8D eigenspaces, so the stabilizer is
[
su(3)\oplus su(3).
]

Adding the exact primitive-idempotent condition defining the committed Spin(9) leaves dimension 12. That algebra has a one-dimensional center and 11D derived algebra; the derived centroid splits 8+3. Therefore
[
\boxed{(su(3)\oplus su(3))\cap spin(9)=su(3)\oplus su(2)\oplus u(1).}
]

This is no longer a dimension mnemonic or literature citation; it is solved directly on the repository's executable Albert matrices.

## 11571 — the unique Cartan connection produces actual curvature

The exact local Cartan solver from Pass11563 was applied to a nonuniform rational periodic frame on `Z3^3`. All 81 site/plaquette curvature matrices are nonzero. The natural scalar contraction has 11 distinct rational values and total
[
\sum_x R(x)=-\frac{360725}{209088}.
]

So the connection is not a pure-gauge artifact on the test frame.

## 11572 — Einstein-Hilbert remains conditional

The kinematic ingredients now exist: frame, unique torsion-free metric-compatible connection, curvature, and a Wilson-Dirac operator. What is still missing is a proof that the variable-frame lattice spectral action converges to the smooth Laplace-type spin operator strongly enough to inherit the standard scalar-curvature heat coefficient.

Therefore the packet does **not** promote older graph-curvature coincidences to Einstein gravity. The scalar-curvature term is conditional on proving the smooth variable-frame Dirac limit.

## 11573 / 11578 — Hamming overlap with gauge flux

A 3D overlap operator was built from the Hamming axis Wilson kernel and a doubled Clifford grading. With free links and with one unit of U(1) flux through the xy torus:

- the Wilson kernel is Hermitian after grading;
- the Ginsparg-Wilson residual stays below (5.5\times10^{-14});
- the index remains zero to numerical precision.

This is the correct odd-dimensional firewall: the event 3-torus supports a GW structure, but a genuine chiral index requires the prior even-dimensional extension.

## 11574 — the gauge selector canonically picks the event plane

On the exact Spin(9) vector 9, the 8D su(3) ideal fixes a unique 3D subspace. The 3D su(2) ideal fixes the complementary 6D subspace, while the u(1) center has eigenvalues
[
0^{\times3},\quad (-4i)^{\times3},\quad (+4i)^{\times3}.
]

Thus the same H3(C) selector that produces the gauge algebra also canonically selects a 3+6 vector splitting.

The caveat is crucial: the derived su(2) is exactly the rotation algebra of that selected 3-plane. If the plane is interpreted as physical spacetime, weak isospin has been conflated with spacetime rotations.

## 11575 — Hesse coupling does not yet determine masses

The commutant of the exact 12D subgroup on the Peirce-16 has real dimension 4. The computed basis exhibits invariant rank-12 and rank-4 sectors with commuting complex structures, consistent with a `C ⊕ C` real commutant.

So symmetry still allows more than one independent mass channel. A Hesse/vacuum scalar by itself does not uniquely determine a fermion mass matrix.

## 11576 — anomaly audit

For the unprojected real Peirce-16 carrier, the u(1) center has weights
[
(-6i)^2,(+6i)^2,(-2i)^6,(+2i)^6.
]
Hence `Tr Q = Tr Q^3 = 0`. Direct evaluation of all 364 symmetric cubic generator triples gives zero.

The carrier is therefore perturbatively anomaly-free because it is real/vectorlike. That is a firewall, not yet a derivation of the chiral Standard Model.

## 11577 — exact U(3) Pati-Salam bridge

The canonical event Spin(3) is the 3D su(2) ideal. Its centralizer inside the exact 12D algebra is
[
su(3)\oplus u(1)\cong u(3),
]
dimension 9.

This gives an exact color-plus-lepton bridge, but again shows that retaining a physical internal weak SU(2) requires a separate copy or selector.

## 11579 — common stabilizer chain

[
E_6 \to F_4 \to su(3)\oplus su(3)
\to su(3)\oplus su(2)\oplus u(1)
\to u(3)
]

where the arrows respectively impose the Hesse trace line, the explicit H3(C), the primitive Albert idempotent, and commutation with the canonical event Spin(3).

The major new positive theorem is the executable Standard-Model-shaped Lie algebra. The major new negative theorem is equally important: **the same selector cannot simultaneously make that su(2) both spacetime rotations and an independent internal weak gauge factor.**

## Physical boundary

Nothing here proves the observed Standard Model as a chiral quantum field theory or derives Einstein gravity, measured couplings, Yukawas, masses, or vacuum selection. The exact finite algebraic structure is now substantially tighter, and the remaining missing selectors are sharply identified.

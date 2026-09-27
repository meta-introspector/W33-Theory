# Pass 11023 — corrected GL₂(3) directed tensor saturation

Producer: `analysis/w33_pass11023_mckay_e7_clock_saturation.py`
Certificate: `data/w33_pass11023_mckay_e7_clock_saturation.json`
Regression: `tests/test_w33_pass11023_mckay_e7_clock_saturation.py`

The filename is retained for provenance. Pass 11025 corrects the original
McKay/E₇ interpretation.

The GAP character table has two faithful two-dimensional characters whose
order-eight values are ±i√2, not ±√2. With those cyclotomic values restored,
tensoring by a faithful 2D irrep is represented by a directed, non-symmetric
matrix A. Tensoring by the conjugate faithful irrep gives B = Aᵀ.

Exact graph facts:

- A ≠ Aᵀ.
- Each faithful tensor quiver has 14 directed edges.
- B = Aᵀ.
- Each quiver is strongly connected.
- The underlying undirected graph has 11 edges, not the seven-edge affine-E₇ tree.

The irreducible dimension vector is

d = (1, 1, 2, 2, 2, 3, 3, 4),

and ordinary tensor dimension gives A d = B d = 2 d.
## Clock saturation that survives the correction

The signed clock multiplicity vector is

m = (2, 1, 0, 0, 0, 1, 2, 3).

For either faithful 2D irrep S, the exact representation-ring identity is

S ⊗ V₂₄ = 3(2q ⊕ 2a ⊕ 2b ⊕ 3t ⊕ 3 ⊕ 4).

The central 12 + 12 split sharpens this to

S ⊗ V₊ = 3(2a ⊕ 2b ⊕ 4),

S ⊗ V₋ = 3(2q ⊕ 3t ⊕ 3).

The second-order identity A(A + I)m = 3d also survives exactly.

These are genuine GL₂(3) representation-ring identities. They do not define
a classical ADE McKay graph and do not imply E₇ field content, an E₇ gauge
symmetry, masses, couplings, or a continuum interaction.

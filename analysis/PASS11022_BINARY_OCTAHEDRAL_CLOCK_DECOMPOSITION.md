# Pass 11022 — exact signed clock carrier splits 12 + 12

Producer: `analysis/w33_pass11022_binary_octahedral_clock_decomposition.py`
Certificate: `data/w33_pass11022_binary_octahedral_clock_decomposition.json`
Regression: `tests/test_w33_pass11022_binary_octahedral_clock_decomposition.py`

> Correction recorded by Pass 11025: the exact group is GL₂(3) = 2⁺S₄,
> the plus Schur cover. It is not the non-isomorphic binary-octahedral
> minus cover 2⁻S₄, so the classical SU(2)/ADE McKay identification does
> not apply here.

Pass 11021 identified the minimal exact signed clock carrier with the
24 noncentral coordinates in the current H₂₇ gauge. Its exact character on
the eight GL₂(3) classes is

χ₂₄ = (24, 0, 0, 6, 0, 0, 0, 2).

The complex representation decomposes as

V₂₄ = 2·1 ⊕ det ⊕ (3std ⊗ det) ⊕ 2·3std ⊕ 3·4spin.

The quotient 2D irrep and both faithful 2D irreps occur with multiplicity zero.
## Central involution

For z = −I, the exact trace on V₂₄ is zero, giving the canonical split

V₂₄ = V₊ ⊕ V₋,   dim V₊ = dim V₋ = 12.

The even sector is

V₊ = 2·1 ⊕ det ⊕ (3std ⊗ det) ⊕ 2·3std,

while the odd sector is exactly

V₋ = 3·4spin.

Every six-coordinate projective clock fibre resolves as 3₊ + 3₋.
Projecting the four fibre indicators with P₊ = (1 + z)/2 and
P₋ = (1 − z)/2 and closing under the exact group gives rank 12 in each
sector. The same is true starting only from the three coarse clock
augmentation generators.

## Cover firewall

The concrete GL₂(3) has 13 nonidentity involutions. A faithful finite
subgroup of SL₂(C) has only one nontrivial involution, −I, so this group
cannot be the binary-octahedral SU(2) group. The exact character audit also
finds Λ²ρ = det for every irreducible 2D representation ρ, not the trivial
determinant required for an SL₂ embedding.

“Spinorial” therefore means only that central −I acts as −1 on the 4D irrep.
No physical fermion, Lorentz-spinor, ADE, or E₇ gauge interpretation is
asserted.

# Pass 11025 — Schur-cover correction and native-cubic tensor firewall

Producer: `analysis/w33_pass11025_native_cubic_mckay_firewall.py`
Certificate: `data/w33_pass11025_native_cubic_mckay_firewall.json`
Regression: `tests/test_w33_pass11025_native_cubic_mckay_firewall.py`

This pass corrects a real group-identification error while preserving the
exact finite-representation results that survive audit.

The exact group is GL₂(3) = 2⁺S₄, not binary octahedral 2⁻S₄. It has
13 nonidentity involutions, already excluding a faithful SL₂(C) embedding.
The faithful 2D characters take order-eight values ±i√2, and every
irreducible 2D representation has exterior square Λ²ρ = det.

After restoring those cyclotomic values, the faithful tensor quiver is
directed and non-symmetric: 14 arrows, with an 11-edge underlying undirected
graph. The conjugate faithful irrep gives the transpose quiver.

Crucially, the strongest clock identity survives:

S ⊗ V₂₄ = 3(2q ⊕ 2a ⊕ 2b ⊕ 3t ⊕ 3 ⊕ 4)

for either faithful 2D irrep S.
## Native cubic audit

The exact invariant-cubic space is large:

dim Sym³(V₂₄)ᴳ = 71.

Under V₂₄ = V₊ ⊕ V₋, the invariant sectors are

- dim Sym³(V₊)ᴳ = 23,
- dim (V₊ ⊗ Sym²(V₋))ᴳ = 48,
- dim Sym³(V₋)ᴳ = 0,
- dim (V₋ ⊗ Sym²(V₊))ᴳ = 0.

So central parity permits only an even number of minus fields, but finite
symmetry leaves 71 independent invariant symmetric cubics.

The committed E₆ tensor is far more specific. Of its 45 signed triads,
32 are wholly noncentral. In a central-eigenbasis those 32 triads expand to
exactly 64 nonzero monomials:

16 of type (+,+,+) and 48 of type (+,−,−).

No odd-parity monomial survives. After removing the common basis factor,
every surviving integer coefficient has magnitude two.

Therefore the corrected GL₂(3) tensor quiver constrains representation content
and parity but does not determine the interaction. The native E₆ incidence
and sign tensor is indispensable additional data. A richer finite-GL₂
reconstruction algebra or quiver-with-potential model remains open.

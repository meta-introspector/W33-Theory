# Pass 11140 — the zero-radius end of the neutral exit is the U(16) non-supersymmetric heterotic string

Producer: `analysis/w33_pass11140_zero_radius_u16.py`
Scan: `analysis/w33_pass11140_scan_zero_radius_roots.py`
Frozen: `data/w33_pass11140_zero_radius_roots_21.json`
Regression: `tests/test_w33_pass11140_zero_radius_u16.py`

As the neutral Wilson-line torus shrinks, a state stays massless only if its momentum on that torus cancels exactly. The
conditions are:
* l = π + N·W with l² = 2;
* the Wilson-line offset π·W + N·W²/2 is integral;
* the NS-sector Witten projection holds.

Validation:
* N = 0 gives the 224 roots of so(16) ⊕ so(16).
* The Witten-twisted sector gives no vectors at N = 0 (its 256 norm-2 states are the (128,1) + (1,128) fermions).
* Three candidate phase rules for winding states agree.
* Without the offset condition the set is **not** reflection-closed, so the condition is required.

**Result (21/21, other Wilson-line torus decompactified): 240 roots, rank 15, one irreducible component = A₁₅ = su(16)**,
plus u(1). This is the charge lattice Υ⁽¹⁾₁₆ ~ su₁₆ ⊕ u₁ of Fraiman–Graña–Parra De Freitas–Sethi (arXiv:2307.13745,
eq. 3.15): the **U(16) non-supersymmetric heterotic string**. It has 2¹ = **2 tachyons**, exactly the two norm-1 states
(one complex species) found independently at zero radius in Pass 11134.

With the other torus kept finite, the algebra breaks to products of su(n) that depend on the model, for example su(8)⊕su(4)⊕su(4), su(7)⊕su(5)⊕su(4) or su(10)⊕su(4)⊕su(2).

**Correction.** Pass 11128 suggested the vector class of the SO(32) or SO(16)×E8 strings. The T-dual string is the U(16)
one; this resolves the identification Pass 11134 left open.

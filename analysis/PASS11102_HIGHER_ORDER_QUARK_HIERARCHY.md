# Pass 11102 — the scalar vacuum does not split charm from up: an AG(3,3) point reflection protects the degeneracy at all orders

Producer: `analysis/w33_pass11102_higher_order_quark_hierarchy.py` (needs the full field dumps; results frozen)
Certificate: `data/w33_pass11102_higher_order_quark_hierarchy.json`
Regression: `tests/test_w33_pass11102_higher_order_quark_hierarchy.py`

Pass 11101 left 12 SO(16)×SO(16) A8 models that pass every test so far:
* tachyon-free, three generations, anomaly-free;
* fractional charges removable (Pass 11097, T* = 3);
* a single heavy top and bottom from a localized Higgs.

At tree level charm and up are degenerate. The three generations sit at the three θ-fixed points of one SU(3) torus.
Those points form an affine line of AG(3,3), the space the 27 fixed points of T⁶/Z3 make up, and the space-group
rule is collinearity in it. The point reflection of that line about the top's fixed point swaps the other two.

**The question.** The Standard-Model-neutral scalar VEVs that give the fractional fermions their masses need not
respect that reflection. Can they split m_c from m_u at higher order?

**Method.** Every Yukawa entry q_i ū_j H is extended to all orders, q_i ū_j H·φⁿ, using the same hidden-singlet VEV
scalars as Pass 11097:
* the order is minimal; n = 0 is cubic with exact rules, and n ≥ 1 uses H-momentum mod 3 (integer programs of
  `w33_so16_selection_rules.py`);
* cubic entries with different fixed points carry the geometric factor ε_geo^(tori differing);
* exponents are read off with ε_VEV = ε and ε_geo = ε^g, for g = ½, 1, 2.

## Result

| survivors | Higgs doublets with a tree-level coupling (up + down) | singular-value exponents | block reflection-symmetric |
|---|---|---|---|
| 12 | 72 | **(0, g, g)** in all 72, for every g | **72 / 72** |

For example, in model 2 (up, one Higgs) the order matrix is

    top          : cubic, O(1)
    charm-up block: off-diagonal cubic, one torus apart (e^{-area}); diagonal order-2 VEV terms

and it is symmetric under exchanging the two non-top generations. The same holds in every case. **The reflection is a
symmetry of all the selection rules together with the whole VEV set, at every order.**

## Reading

In the 12 survivors **m_c/m_u and m_s/m_d are O(1) parametrically**. The only way to a hierarchy is a fine-tuned
cancellation in the determinant of the 2×2 block:

  a d ε_VEV⁴ − b c ε_geo² ≈ 0.

The observed m_c/m_u ≈ 600 needs such a tuning. The observed m_s/m_d ≈ 20 is at the edge of what O(1)
coefficients can give. The up-sector hierarchy is therefore the obstruction that remains for the W(3,3)
SO(16)×SO(16) models, alongside the universal positive vacuum energy of Pass 11099.

Where the protection comes from. The reflection x ↦ −x of the Z3 torus about a fixed point normalises the point
group. Space-group selection rules respect it. Because charm and up differ only by the fixed point along that torus,
it also maps their quantum numbers into each other in every case found. A way out needs one of:
* a vacuum that breaks the reflection through fields whose charges are not reflection-symmetric; no such hidden-singlet
  VEV scalar exists in these models, since all 72 cases are symmetric;
* non-Abelian flavour structure from a hidden group;
* a different geometry (Pass 11100 found no tachyon-free Standard Model outside T⁶/Z3).

Scope. The pass computes orders and parametric exponents with random O(1) coefficients. The actual coupling constants
are not computed; they could at best produce an O(1), non-parametric splitting.

# Pass 11334 — the F₉-linear three-leg law, exhaustively: magic on three legs always breaks the arrow

Producer: `analysis/w33_pass11334_f9_three_leg_law.py`
Certificate: `data/w33_pass11334_f9_three_leg_law.json`
Regression: `tests/test_w33_pass11330_11334.py`

**Status before.** Pass 11311 found 6000/6000 sampled violations. Codex's audit correctly classed this as sample
evidence, not a universal law.

**Reduction (exact).**
* T^e = Z^{(e − e mod 3)/3}·T^{e mod 3}, because T³ = Z.
* The Z-type Paulis on the legs only change the Pauli frame of V = W(a)·V_M: on outputs directly, on inputs via
  V Z^k = W(M·z^k)·V. Every frame is enumerated anyway.
* So the 64 × 81 × (4 × 6³ × 3) = **13 436 928** three-leg dressings of the 64 F₉-linear perfect classes reduce to
  64 × 81 × 4 × 2³ = **165 888** cases (T^{1 or 2} on three legs, nothing on the fourth). Each was decided exactly
  (Pass 11252).

| magic legs | violating / cases |
|---|---|
| 3 | **165 888 / 165 888** |
| 1 | **3456 / 41 472 = 1/12** |

**Results.**
* **Theorem (by exhaustion):** every F₉-linear perfect two-qutrit tick, in every Pauli frame, with cubic phases on any
  three of its four legs, breaks substrate time reversal.
* With one magic leg the violating fraction is exactly **1/12**, the same as for a single cubic gate after a uniform
  one-qutrit Clifford (Pass 11266).
* For other perfect ticks the three-leg fraction is about 91% (Pass 11311, sampled), so the law is specific to
  F₉-linearity.

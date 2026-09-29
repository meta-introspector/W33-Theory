# Pass 11161 — perfect gates are the maximal scramblers; information is never concentrated; forgetting the partner is the arrow of time

Producer: `analysis/w33_pass11161_scrambling_and_arrow.py`
Regression: `tests/test_w33_pass11161_scrambling_and_arrow.py`

For S ∈ Sp(4,3), the Choi-state mutual informations (in units of log 3) are:
* I(A_in : A_out) = 2 − rank S_BA;
* I(A_in : B_out) = 2 − rank S_AA;
* I(A_in : A_out B_out) = 2.

So **I₃(A_in : A_out : B_out) = 2 − rank S_AA − rank S_BA.** These formulas were checked against brute-force entropies of
real Choi states.

* **Scrambling census.**

| I₃ | gates |
|---|---|
| 0 | 1152 (local and swap-type) |
| −1 | 36 864 |
| **−2** (the minimum) | **13 824 = the perfect gates** |

* **I₃ ≤ 0 for every two-qutrit Clifford gate.** This is a corollary of the determinant law det S_AA + det S_BA = 1: the
  two ranks can never both be deficient. Information about A's past is never concentrated in one output, and every
  entangling Clifford gate scrambles.
* **The arrow of time.** When S_BA is invertible, as for every perfect gate, and the partner B starts maximally mixed,
  the channel ρ ↦ Tr_B[U(ρ ⊗ I/3)U†] on A is **completely depolarising** (checked to 10⁻¹⁵).
  * So one tick with a partner, followed by forgetting the partner, reproduces exactly the paper's Theorem 4.7 arrow: a
    Pauli twirl equals complete depolarisation, at a Landauer cost of 2 k_BT ln 3.
  * The forgotten environment is the partner plus its purification, a 9-dimensional space: two trits.

**Reading.** The paper's arrow of time comes from forgetting a two-trit record. Here that record is supplied by
entanglement: a maximally entangling tick hands A's information to the partner, and discarding the partner discards it.
A perfect gate does this while also leaving A's past invisible in B's future alone.

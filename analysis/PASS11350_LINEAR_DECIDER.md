# Pass 11350 — one magic gate is F₃ linear algebra: an exact linear decider

Producer: `analysis/w33_pass11350_linear_decider.py`
Certificate: `data/w33_pass11350_linear_decider.json`
Regression: `tests/test_w33_pass11350_11354.py`

## Derivation

Pass 11331: U = C·T₁ (C = W(a)·V_M) is substrate-time-reversible iff a Clifford E = W(e)·V_Q solves
T₁ E T₁⁻¹ = λ·C⁻¹ E Cᵀ. Write both sides in the normal form W(frame)·V_(symplectic); both identities were verified
numerically.

* **Right side:** C⁻¹ E Cᵀ ∝ W(M⁻¹(e − a) − P J a)·V_P, with P = M⁻¹ Q J M⁻¹ J and J = diag(1, −1, …).
* **Left side:** T₁ W(e) V_Q T₁⁻¹ ∝ W(g(e) + s^{−k} r_Q)·V_{s^{−k}Q}.
  * k = e_{x₁}, and s is the unit shear z₁ ↦ z₁ + x₁.
  * g(e) = e + t_k (verified).
  * r_Q is the frame of T₁ V_Q T₁⁻¹ V_Q†.

Matching symplectic parts and frames (phases are free) gives two conditions:

* **(S)** M·s^{−k}·Q·(JMJ) = Q, Q z₁ = z₁, Q symplectic. This is linear in Q, plus a symplecticity filter.
* **(F)** (I − M⁻¹) e = −(M⁻¹ + s^{−k} Q J) a − t_k − s^{−k} r_Q, with e_{x₁} = k. For fixed (Q, k) this is affine
  over F₃ in (e, a).

**Consequence.** The reversible frames of a class are a **union of affine subspaces**, one for each solution of (S).
The affine law of Passes 11309/11330 says this union is always a single affine subspace at n = 2. Pass 11352 shows it
is **not** always at n = 3.

## Validation

| | |
|---|---|
| two-qutrit classes decided by (S)+(F) | **51 838 / 51 840 agree with Pass 11330's exhaustive frame counts; 0 mismatches** |
| the other 2 classes | M = ±I. (S) has too many solutions to enumerate, but both are trivial: for M = I, C is a Pauli and W(a)·T₁ violates iff a shifts X₁ (54 frames, the one-qutrit rule); for M = −I no frame violates. Both agree with Pass 11330. |
| frame-level spot checks against the Weyl criterion | 300 / 300 |
| three qutrits | agrees with the exhaustive 729-frame results of Pass 11333; 660 / 660 Weyl checks in Pass 11352 |

**Speed.** One Weil unitary per candidate solution of (S), built as Σ_p W(Mp)|j⟩⟨0|W(p)† (a single matrix product,
checked on generators). A three-qutrit class takes about a second, against hours for 729 Weyl decisions. This is what
makes Pass 11352's 32 000-class three-qutrit census possible.

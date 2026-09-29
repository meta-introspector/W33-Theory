# Pass 11143 — space and time are two quadrangles on the same 40 points

Producer: `analysis/w33_pass11143_space_time_two_quadrangles.py` (self-contained, exact over F₃)
Certificate: `data/w33_pass11143_space_time_two_quadrangles.json`
Regression: `tests/test_w33_pass11143_space_time_two_quadrangles.py`

## Background in the repo

The paper builds W(3,3) twice. Section 2 builds it from **two qutrits**: the commutation of P_a ⊗ P_b, with form
[a,a′] + [b,b′]. Theorem 3.1 builds it from **one qutrit acting on itself**: the left/right action A ↦ P_u A P_v†, with
form [u,u′] − [v,v′]. It reads the left–right exchange as a built-in time reversal. The two constructions were never
compared.

## What they are together

Same 40 projective points of F₃⁴, two symplectic forms:

| | spatial w_s = ω ⊕ ω | temporal w_t = ω ⊕ (−ω) |
|---|---|---|
| structure | GQ(3,3), SRG(40,12,2,4) | GQ(3,3), SRG(40,12,2,4) |
| lines | 40 | 40 |
| **shared** | **16 product lines ℓ₁ ⊕ ℓ₂** | **the same 16** |
| own lines | 24 graphs {(a, N a)}, det N = −1 | 24 graphs {(a, M a)}, M ∈ SL(2,3) |
| states / processes | 216 maximally entangled stabilizer states (Schmidt ⅓,⅓,⅓); 144 product | 24 Clifford classes = the qutrit Clifford group mod Paulis (216 unitaries mod phase) |

* **Temporal lines are Clifford channels.** The two-time pseudo-density operator of a qutrit that starts maximally mixed
  and passes through a Clifford gate U is R_U = J(U)/3, with J the Jamiołkowski matrix. It is supported exactly on the
  temporal line graph(−M_U). Its spectrum is (+⅓)⁶(−⅓)³, so its temporal negativity is **1**, the same as the spatial
  negativity of a maximally entangled pair.
* **Partial transposition exchanges the two own-line sets.** Operator for operator, ρ^{T_B} of each of the 216 maximally
  entangled two-qutrit stabilizer states **is** R_U for one of the 216 qutrit Clifford unitaries. This is a **bijection**.
* **Only local operations preserve both.** An exhaustive search over all 103 680 similitudes of the spatial form finds
  exactly **2304** that also preserve the temporal one. These are the block-diagonal maps with equal multipliers and the
  factor swap: local operations and exchange of the parties. Every entangling Clifford operation distinguishes space from
  time.

## Reading

* Space and time agree on everything classical (the 16 product lines) and disagree on exactly the 24 entangled lines.
* A maximally entangled pair of qutrits and **one qutrit persisting through a gate** are the same line of points, read
  with the two forms. The self-entanglement of Section 3 (one qutrit's operator space, "input" and "output" legs) is the
  temporal half of this picture.
* **Spooky action at a distance, precisely.** The correlation frame of a spatial Bell pair is the partial transpose of the
  frame of a single qutrit correlated with itself across time. The spatial version cannot be used to signal; the temporal
  one can, since the system itself carries the information. Nothing moves backwards in time. What is shared is the
  *geometry* of the correlations, not a mechanism.

Prior art: the general statement "partial transpose = space–time swap" is Fullwood–Li (arXiv:2508.12256). The
bipartite ↔ temporal correlation map is Marcovitch–Reznik (arXiv:1107.2186). The count 360 = 216 + 144 is in BT821. New
here and specific to W(3,3):
* the two-quadrangle structure on one point set;
* the 16/24 split and its meaning (product ↔ classical, own lines ↔ entangled / Clifford);
* the state ↔ unitary bijection;
* the common automorphism group of order 2304.
